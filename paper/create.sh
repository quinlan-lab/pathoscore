#create score sets

cd ../score-sets/GRCh37/CADD
bash make.sh
cd -

cd ../score-sets/GRCh37/CCR
bash make.sh
cd -

cd ../score-sets/GRCh37/DANN
bash make.sh
cd -

cd ../score-sets/GRCh37/GERP
bash make.sh
cd -

cd ../score-sets/GRCh37/MCAP
bash make.sh
cd -

cd ../score-sets/GRCh37/MPC
bash make.sh
cd -

cd ../score-sets/GRCh37/MTR
bash make.sh
cd -

cd ../score-sets/GRCh37/REVEL
bash make.sh
cd -

cd ../score-sets/GRCh37/RVIS
bash make.sh
cd -

cd ../score-sets/GRCh37/VVP
bash make.sh
cd -

cd ../score-sets/GRCh37/aloft
bash make.sh
cd -

cd ../score-sets/GRCh37/fathmm
bash make.sh
cd -

cd ../score-sets/GRCh37/fitcons
bash make.sh
cd -

cd ../score-sets/GRCh37/metasvm
bash make.sh
cd -

cd ../score-sets/GRCh37/mis_Z
bash make.sh
cd -

cd ../score-sets/GRCh37/pLI
bash make.sh
cd -

cd ../score-sets/GRCh37/phastCons
bash make.sh
cd -

cd ../score-sets/GRCh37/phylop
bash make.sh
cd -

cd ../score-sets/GRCh37/polyphen2  
bash make.sh
cd -

cd ../score-sets/GRCh37/sift
bash make.sh
cd -

# make truth sets

cd ../truth-sets/GRCh37/clinvar
bash make.sh
cd -

cd ../truth-sets/GRCh37/homsy
bash make.sh
cd -

cd ../truth-sets/GRCh37/samocha
bash make.sh
cd -

cd ../truth-sets/GRCh37/wellderly
bash make.sh
cd -

cd ../truth-sets/GRCh37/wigler-autism
bash make.sh
cd -

# make gene sets

cd ../gene-sets/GRCh37/ad_genes/
bash make.sh
cd -

cd ../gene-sets/GRCh37/hi_genes/
bash make.sh
cd -

# make population sets

cd ../scripts/gnomad
bash makeexac.sh
bash makegnomad.sh
cd -
