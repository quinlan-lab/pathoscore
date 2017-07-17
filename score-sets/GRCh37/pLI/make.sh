set -euo pipefail

name=forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz
if [[ ! -e "$name" ]]; then
    wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/manuscript_data/$name
fi
python make.py $name | bgzip -c > pLI.bed.gz
tabix pLI.bed.gz
