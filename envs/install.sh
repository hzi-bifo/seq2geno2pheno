#!/usr/bin/env bash
# Install every tool Seq2Geno2Pheno runs, without containers: one conda environment per tool
# from the lock files next to this script, plus three checksummed upstream downloads, and
# wrapper scripts in <prefix>/bin that the workflow puts on PATH (-profile conda).
#
#   bash envs/install.sh [--prefix DIR] [--force]
#
# DIR defaults to $S2G2P_ENV_DIR, else envs/installed/ in this checkout (what the workflow's
# --env_dir defaults to). Re-running is cheap: an environment whose lock file has not changed
# is kept. Needs micromamba, mamba or conda on PATH, internet access (conda-forge, bioconda,
# github.com, ftp.ncbi.nlm.nih.gov, ftp.gnu.org, CPAN) and ~6 GB. Linux x86-64 only: the locks
# are linux-64 and BLAST+ 2.4.0 is NCBI's linux binary.
#
# What is installed, and why at these versions (docs/design.md, "The published method"):
#   prokka      Prokka 1.14.6 with the Bioconda tool versions of the Hu et al. 2024 reproduction
#   roary-perl  Perl 5.22 + BioPerl + mafft 7.310: Seq2Geno's Roary environment, exact package list
#   roary-deps  CD-HIT 4.6.6, MCL 14-137, bedtools 2.26, FastTree, PRANK (Roary's external programs)
#   roary       Seq2Geno's Roary fork, hzi-bifo/Roary @ 7402d1b (Roary 3.8.2 + a coverage filter)
#   blast-2.4.0 NCBI BLAST+ 2.4.0, the version the fork installs
#   parallel    GNU parallel 20160722, the version the fork installs
#   mafft       mafft 7.310 (Seq2Geno's indel alignment environment, exact package list)
#   msa2vcf     jvarkit msa2vcf 1.0 (Seq2Geno's indel calling environment, exact package list)
#   s2g2p       this package, editable, on Python 3.12
#   hu2024ml    the published model environment (Python 3.7, scikit-learn 0.24.1) + this package
#   folds       R 4.2 + ape + phangorn (phylogeny-aware folds) and kma_clustering, compiled here
#               from envs/vendor/kma_clustering.c as Hu et al. 2024 shipped it (homology-aware folds)
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(dirname "$here")
prefix=${S2G2P_ENV_DIR:-$here/installed}
force=0
while [ $# -gt 0 ]; do
    case $1 in
        --prefix) prefix=$2; shift 2 ;;
        --prefix=*) prefix=${1#*=}; shift ;;
        --force) force=1; shift ;;
        -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done
mkdir -p "$prefix"
prefix=$(cd "$prefix" && pwd)

die() { echo "install.sh: $*" >&2; exit 1; }
say() { echo "== $*"; }

[ "$(uname -s)/$(uname -m)" = "Linux/x86_64" ] || die "Linux x86-64 only (the locks are linux-64)"
solver=""
for s in micromamba mamba conda; do
    if command -v "$s" > /dev/null 2>&1; then solver=$s; break; fi
done
[ -n "$solver" ] || die "needs micromamba, mamba or conda on PATH"
say "solver: $solver; prefix: $prefix"

# Upstream tarballs (sha256 from dev/lock-envs.sbatch)
ROARY_URL=https://github.com/hzi-bifo/Roary/archive/7402d1bab868d73b9c1e7cf30e8af8c5b2578313.tar.gz
ROARY_SHA=caa28ece0c9aa4011aedd647dcc5eec5794ef7f2bb7c0336751955939f32bda9
BLAST_URL=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-x64-linux.tar.gz
BLAST_SHA=21aa7ca60954ce9c6d3f572e427fee804291fcad41447d554033a81d5af96a2b
PARALLEL_URL=https://ftp.gnu.org/gnu/parallel/parallel-20160722.tar.bz2
PARALLEL_SHA=e391ebd081e8ba13e870be68106d1beb5def2b001fa5881f46df0ab95304f521
# envs/vendor/kma_clustering.c = hzi-bifo/AMR_benchmarking@2850078
# AMR_software/AytanAktug/cluster_preparation/kma_clustering.c (KMA_index 1.0.3, Apache-2.0)
KMA_SRC_SHA=8e2960827f652ffe0fb08754b0b9376742781bb547e96e80172a9409a4aee48d

# -- conda environments --------------------------------------------------------------------
# An environment is (re)created when its lock file changed; the lock's sha256 is the stamp.
create_env() {
    local name=$1 lock=$here/$1.lock.txt
    local stamp=$prefix/$name/.s2g2p-lock
    local want; want=$(sha256sum "$lock" | cut -d' ' -f1)
    if [ $force = 0 ] && [ -f "$stamp" ] && [ "$(cat "$stamp")" = "$want" ]; then
        say "$name: up to date"; return 1
    fi
    say "$name: creating from $(basename "$lock")"
    rm -rf "${prefix:?}/$name"
    case $solver in
        micromamba) micromamba create -y -q -p "$prefix/$name" --file "$lock" ;;
        *) "$solver" create -y -q -p "$prefix/$name" --file "$lock" ;;
    esac
    echo "$want" > "$stamp.pending"
    return 0
}
done_env() { mv "$prefix/$1/.s2g2p-lock.pending" "$prefix/$1/.s2g2p-lock"; }

# Java without its memory-mapped /tmp/hsperfdata file: on a node whose /tmp is full the JVM
# dies at start with SIGBUS, and Prokka reads that as minced failing its version check.
# -XX:-UsePerfData only turns off JVM performance counters. --follow-symlinks: the launcher
# is a link, and replacing it breaks the script's own jar lookup.
no_perfdata() {
    local f=$1 pat=$2
    sed -i --follow-symlinks "s|$pat|& -XX:-UsePerfData|" "$f"
    grep -q -- '-XX:-UsePerfData' "$f" || die "could not patch $f"
}

if create_env prokka; then
    no_perfdata "$prefix/prokka/bin/minced" '^exec $java'
    PATH=$prefix/prokka/bin:$PATH prokka --setupdb > /dev/null
    done_env prokka
fi

if create_env roary-perl; then
    # Roary's Perl modules at the versions the reproduction loaded. Seq2Geno's
    # install_perl_mods.sh is not run: it pins newer versions, and in the reproduction every
    # one of its cpanm calls failed, so Roary ran with these. A module the environment does not
    # provide on Perl 5.22's path is installed with cpanm at the reproduction's version.
    # cpanm installs into an explicit library inside the environment (roary-perl/cpan). Left to
    # itself it decides whether site_perl is writable from the directory's owner -- wrong on
    # storage that reports another owner -- and then falls back to ~/perl5, where a module it
    # finds from some other project counts as "up to date". The caller's own Perl settings are
    # cleared for the same reason: Roary must load exactly what is checked here.
    perl=$prefix/roary-perl/bin/perl
    cpan_lib=$prefix/roary-perl/cpan
    clean_perl=(env -u PERL5LIB -u PERLLIB -u PERL_LOCAL_LIB_ROOT -u PERL_MB_OPT -u PERL_MM_OPT
                PERL5LIB="$cpan_lib/lib/perl5")
    while read -r m v; do
        if ! "${clean_perl[@]}" "$perl" -M"$m" -e1 2> /dev/null; then
            "${clean_perl[@]}" PATH="$prefix/roary-perl/bin:$PATH" cpanm --notest -l "$cpan_lib" "$m@$v"
        fi
        got=$("${clean_perl[@]}" "$perl" -M"$m" -e "print \$${m}::VERSION")
        [ "$got" = "$v" ] || die "Perl module $m is $got, the reproduction loaded $v"
    done <<'MODS'
Bio::Perl 1.006924
Moose 2.1804
Log::Log4perl 1.47
Array::Utils 0.5
File::Grep 0.02
File::Which 1.20
File::Slurper 0.008
Text::CSV 1.33
XML::Parser 2.44
Graph 0.9704
Graph::Writer::Dot 2.09
File::Find::Rule 0.34
Exception::Class 1.40
MODS
    done_env roary-perl
fi

create_env roary-deps && done_env roary-deps
create_env mafft && done_env mafft

if create_env msa2vcf; then
    no_perfdata "$prefix/msa2vcf/bin/msa2vcf" '^java'
    done_env msa2vcf
fi

if create_env s2g2p; then
    "$prefix/s2g2p/bin/pip" install -q --no-deps -e "$repo"
    done_env s2g2p
fi

if create_env hu2024ml; then
    "$prefix/hu2024ml/bin/pip" install -q --no-deps -r "$here/hu2024ml.pip.txt"
    "$prefix/hu2024ml/bin/pip" install -q --no-deps "pyyaml==6.0.1"
    # s2g2p by path: its build backend needs Python >= 3.8, this environment is 3.7
    site=$("$prefix/hu2024ml/bin/python" -c "import site; print(site.getsitepackages()[0])")
    echo "$repo/src" > "$site/s2g2p.pth"
    done_env hu2024ml
fi

create_env folds && done_env folds

# kma_clustering, compiled as AytanAktug_SSSA.sh did (gcc -O3 ... -lm -lz), with the
# environment's compiler and zlib. Rebuilt when the vendored source changes.
kma_src=$here/vendor/kma_clustering.c
kma_bin=$prefix/folds/bin/kma_clustering
[ "$(sha256sum < "$kma_src" | cut -d' ' -f1)" = "$KMA_SRC_SHA" ] \
    || die "envs/vendor/kma_clustering.c is not the published file (sha256 $KMA_SRC_SHA)"
if [ $force = 1 ] || [ ! -x "$kma_bin" ] || [ "$(cat "$kma_bin.sha256" 2> /dev/null)" != "$KMA_SRC_SHA" ]; then
    say "kma_clustering: compiling $kma_src"
    cc=$(ls "$prefix"/folds/bin/x86_64-conda*-linux-gnu-cc 2> /dev/null | head -1)
    [ -n "$cc" ] || die "folds: the environment has no C compiler"
    "$cc" -O3 -o "$kma_bin" "$kma_src" -I"$prefix/folds/include" -L"$prefix/folds/lib" \
        -Wl,-rpath,"$prefix/folds/lib" -lm -lz
    echo "$KMA_SRC_SHA" > "$kma_bin.sha256"
fi

# -- upstream tarballs ---------------------------------------------------------------------
fetch() {  # url sha256 destdir
    local url=$1 sha=$2 dest=$3
    if [ $force = 0 ] && [ -f "$dest/.s2g2p-sha256" ] && [ "$(cat "$dest/.s2g2p-sha256")" = "$sha" ]; then
        say "$(basename "$dest"): up to date"; return
    fi
    say "$(basename "$dest"): downloading $url"
    local tmp; tmp=$(mktemp -d "$prefix/.dl.XXXXXX")
    curl -fsSL --retry 3 -o "$tmp/archive" "$url"
    echo "$sha  $tmp/archive" | sha256sum -c --quiet - || die "checksum mismatch for $url"
    rm -rf "$dest"; mkdir -p "$dest"
    tar -xf "$tmp/archive" --no-same-owner --strip-components=1 -C "$dest"
    rm -rf "$tmp"
    echo "$sha" > "$dest/.s2g2p-sha256"
}
fetch "$ROARY_URL" "$ROARY_SHA" "$prefix/roary"
fetch "$BLAST_URL" "$BLAST_SHA" "$prefix/blast-2.4.0"
fetch "$PARALLEL_URL" "$PARALLEL_SHA" "$prefix/parallel"

missing=$(ldd "$prefix/blast-2.4.0/bin/blastp" | grep 'not found' || true)
[ -z "$missing" ] || die "BLAST+ 2.4.0 needs system libraries this host lacks:
$missing
(on RHEL/Alma/Rocky: dnf install libnsl; on Debian/Ubuntu: the libc's libnsl.so.1 package)"

# -- wrappers --------------------------------------------------------------------------------
# Each tool runs with the environment its Seq2Geno rule gave it. Roary sees only the programs
# it calls from roary-deps (that environment's own perl must not shadow Perl 5.22).
mkdir -p "$prefix/roary-tools" "$prefix/bin" "$prefix/bin-hu2024ml"
# mcxdeblast is MCL's BLAST-to-graph step, which Roary calls before mcl; `roary -a` does not list
# it, so its absence only shows when clustering starts (and Roary then exits 0 with no tables).
for t in cd-hit mcl mcxdeblast bedtools FastTree prank; do
    [ -x "$prefix/roary-deps/bin/$t" ] || die "roary-deps lacks $t"
    ln -sfn "$prefix/roary-deps/bin/$t" "$prefix/roary-tools/$t"
done

wrap() {  # name, body
    printf '#!/usr/bin/env bash\n# written by envs/install.sh\nexport SHELL=/bin/bash LC_ALL=C\n%s\n' "$2" > "$prefix/bin/$1"
    chmod +x "$prefix/bin/$1"
}
P=$prefix
wrap prokka "export PATH=$P/prokka/bin:\$PATH
exec $P/prokka/bin/prokka \"\$@\""
# Everything Roary runs comes from this PATH, not from the caller's: its own scripts, the programs
# above, BLAST+ 2.4.0, GNU parallel, mafft 7.310 (-e --mafft) and Perl 5.22.
ROARY_PATH=$P/roary/bin:$P/roary-tools:$P/blast-2.4.0/bin:$P/parallel/src:$P/mafft/bin:$P/roary-perl/bin
wrap roary "export PATH=$ROARY_PATH:\$PATH
unset PERLLIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
export PERL5LIB=$P/roary/lib:$P/roary-perl/cpan/lib/perl5
exec $P/roary-perl/bin/perl $P/roary/bin/roary \"\$@\""
wrap blastp "exec $P/blast-2.4.0/bin/blastp \"\$@\""
wrap makeblastdb "exec $P/blast-2.4.0/bin/makeblastdb \"\$@\""
wrap mafft "export PATH=$P/mafft/bin:\$PATH
exec $P/mafft/bin/mafft \"\$@\""
wrap msa2vcf "export PATH=$P/msa2vcf/bin:\$PATH
exec $P/msa2vcf/bin/msa2vcf \"\$@\""
wrap s2g2p "exec $P/s2g2p/bin/python -m s2g2p \"\$@\""
# R sees only the environment's library: no site or user library, profile or environ file.
wrap nj_tree "unset R_LIBS R_LIBS_SITE R_HOME
export R_LIBS_USER=$P/folds/lib/R/library
exec $P/folds/bin/Rscript --vanilla $repo/src/s2g2p/R/nj_tree.R \"\$@\""
wrap kma_clustering "exec $P/folds/bin/kma_clustering \"\$@\""
printf '#!/usr/bin/env bash\n# written by envs/install.sh\nexec %s/hu2024ml/bin/python -m s2g2p "$@"\n' "$P" > "$P/bin-hu2024ml/s2g2p"
chmod +x "$P/bin-hu2024ml/s2g2p"

# -- checks ------------------------------------------------------------------------------------
say "checking"
B=$P/bin
"$B/prokka" --version 2>&1 | grep -q "1.14.6" || die "prokka is not 1.14.6"
"$P/prokka/bin/minced" --version | grep -q "0.4.2" || die "minced does not start"
"$B/blastp" -version | grep -q "2.4.0" || die "blastp is not 2.4.0"
"$B/makeblastdb" -version | grep -q "2.4.0" || die "makeblastdb is not 2.4.0"
"$B/mafft" --version 2>&1 | grep -q "7.310" || die "mafft is not 7.310"
printf ">a_1\nACGT\n>b_1\nACGA\n" | "$B/msa2vcf" > /dev/null || die "msa2vcf does not run"
# The fork's patch: BLAST reports both lengths, for its coverage filter (upstream Roary does not)
grep -q 'slen qlen' "$P/roary/lib/Bio/Roary/External/Blastp.pm" || die "the Roary fork's coverage filter is missing"
# Roary's own dependency report must find every program at the pinned version
"$B/roary" -a > "$P/roary-dependencies.txt" 2>&1 || true
for want in "blastp version is 2.4.0" "makeblastdb version is 2.4.0" "mcl version is 14-137" \
            "cd-hit version is 4.6" "bedtools version is 2.26" "mafft version is 7.310" \
            "parallel version is 20160722"; do
    grep -q "$want" "$P/roary-dependencies.txt" || die "roary: expected '$want' (see $P/roary-dependencies.txt)"
done
PATH=$P/roary/bin:$PATH command -v extract_proteome_from_gff > /dev/null || die "Roary's helper scripts are missing"
# Every external program Roary invokes (Bio/Roary/External/*.pm) must resolve from the wrapper's
# own PATH, with nothing from the caller's environment to fall back on.
for t in blastp makeblastdb cd-hit mcl mcxdeblast FastTree prank mafft bedtools parallel; do
    env -i PATH="$ROARY_PATH:/usr/bin:/bin" bash -c "command -v $t" > /dev/null \
        || die "roary: '$t' is not on the wrapper's PATH"
done
"$B/s2g2p" --version
"$P/bin-hu2024ml/s2g2p" --version
"$P/hu2024ml/bin/python" -c "import sklearn; assert sklearn.__version__ == '0.24.1', sklearn.__version__"
"$B/kma_clustering" -v | grep -q "KMA_index-1.0.3" || die "kma_clustering is not KMA_index 1.0.3"
chk=$(mktemp -d "$P/.check.XXXXXX")
printf '>a\nACGTACGTACGTAAGG\n>b\nACGTACGTACGTAAGC\n>c\nTCGAACGTACGAAAGC\n>d\nTCGAACGTACGAATGC\n' > "$chk/aln.fa"
"$B/nj_tree" "$chk/aln.fa" "$chk/t.nwk" > /dev/null || die "nj_tree does not run"
for s in a b c d; do grep -q "$s" "$chk/t.nwk" || die "nj_tree wrote no leaf $s"; done
"$B/nj_tree" "$chk/aln.fa" "$chk/t.nwk" | grep -q "phangorn 2.10.0, ape 5.6.2" || die "nj_tree: phangorn/ape are not 2.10.0/5.6-2"
rm -rf "$chk"
command -v ps > /dev/null || echo "warning: no 'ps' on PATH; Nextflow needs it for task metrics"

{
    printf 'component\tsource\n'
    for n in prokka roary-perl roary-deps mafft msa2vcf s2g2p hu2024ml folds; do
        printf '%s\t%s sha256:%s\n' "$n" "envs/$n.lock.txt" "$(cat "$P/$n/.s2g2p-lock")"
    done
    printf 'kma_clustering\tenvs/vendor/kma_clustering.c sha256:%s\n' "$KMA_SRC_SHA"
    printf 'roary\t%s\nblast-2.4.0\t%s\nparallel\t%s\n' "$ROARY_URL" "$BLAST_URL" "$PARALLEL_URL"
} > "$P/manifest.tsv"
date -u +%FT%TZ > "$P/.installed"
say "done. Run the workflow with: -profile conda --env_dir $P"
