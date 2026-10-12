#!/usr/bin/env bash
# Download the tutorial's 35 Enterococcus faecium assemblies from BV-BRC (formerly PATRIC) into
# tutorial/genomes/, and check them against the copies the Hu et al. 2024 benchmark used
# (genomes.md5). About 100 MB; already downloaded files are kept.
#
#   bash tutorial/fetch_genomes.sh
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$here/genomes"
cd "$here"
sed 's|.*genomes/||; s|\.fna$||' genomes.md5 | while read -r id; do
    f=genomes/$id.fna
    [ -s "$f" ] && continue
    echo "fetching $id"
    curl -fsSL --retry 3 -o "$f.part" "https://ftp.bvbrc.org/genomes/$id/$id.fna"
    mv "$f.part" "$f"
done
if md5sum -c --quiet genomes.md5; then
    echo "all $(wc -l < genomes.md5) genomes downloaded and identical to the benchmark's copies"
else
    echo "warning: the files above differ from the copies the benchmark used (BV-BRC may have" >&2
    echo "updated them); the tutorial still runs, but numbers can differ from docs/tutorial.md" >&2
fi
