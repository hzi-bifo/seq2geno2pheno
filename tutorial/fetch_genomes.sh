#!/usr/bin/env bash
# Download the tutorial's 35 Enterococcus faecium assemblies from BV-BRC (formerly PATRIC) into
# tutorial/genomes/, and check them against the copies the Hu et al. 2024 benchmark used
# (genomes.md5). About 100 MB; already downloaded files are kept.
#
#   bash tutorial/fetch_genomes.sh
#
# The genomes come from BV-BRC's data API, which returns the contigs of the genome FASTA files
# PATRIC distributed in another order and layout. They are written back in that layout (contigs
# in accession order, `accn|` dropped, a missing description left empty rather than the API's
# `undefined`, 60 bases per line, a blank line after each contig), so the
# files are byte-identical to the benchmark's and Prokka numbers the genes the same way.
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$here/genomes"
cd "$here"
api="https://www.bv-brc.org/api/genome_sequence/"
sed 's|.*genomes/||; s|\.fna$||' genomes.md5 | while read -r id; do
    f=genomes/$id.fna
    [ -s "$f" ] && continue
    echo "fetching $id"
    curl -fsSL --retry 5 --retry-delay 10 --retry-connrefused -m 300 \
        "${api}?eq(genome_id,${id})&limit(25000)&http_accept=application/dna+fasta" |
    awk '
        /^>/ { if (h != "") print h "\t" s; h = substr($0, 2); sub(/^accn\|/, "", h)
               sub(/   undefined   \[/, "      [", h); s = ""; next }
        { gsub(/[ \t\r]/, ""); s = s tolower($0) }
        END { if (h != "") print h "\t" s }' |
    LC_ALL=C sort -t $'\t' -k1,1 |
    awk -F'\t' '{ print ">" $1; for (i = 1; i <= length($2); i += 60) print substr($2, i, 60); print "" }' \
        > "$f.part"
    [ -s "$f.part" ] || { echo "no sequence returned for $id" >&2; exit 1; }
    mv "$f.part" "$f"
done
if md5sum -c --quiet genomes.md5; then
    echo "all $(wc -l < genomes.md5) genomes downloaded and identical to the benchmark's copies"
else
    echo "warning: the files above differ from the copies the benchmark used (BV-BRC may have" >&2
    echo "updated them); the tutorial still runs, but numbers can differ from docs/tutorial.md" >&2
fi
