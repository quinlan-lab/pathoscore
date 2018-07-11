if [ ! -s mpc.txt.gz ]; then
    (echo -n '#' && wget -O - ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/regional_missense_constraint/fordist_constraint_official_mpc_values.txt.gz \
        | zcat - \
        | cut -f 1-4,19 ) \
        | bgzip -c > mpc.txt.gz
    tabix -b2 -e2 mpc.txt.gz
fi
