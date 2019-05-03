mkdir -p $HOME/public_html/pathoscorepaper

#########
#CLINVAR#
#########

# metrics of choice: VVP, CCR, MPC, CADD, REVEL, RVIS, MCAP, GERP++, MetaSVM, aLoFT_Dominant, MTR, Grantham (chosen based on performance and historical value)

# to get AD variant data and ROCs
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvarad -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf adgene.pathogenic.vcf.gz adgene.benign.vcf.gz

# to get AR variant data and ROCs
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvarar -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf argene.pathogenic.vcf.gz argene.benign.vcf.gz

# regular clinvar
python ../pathoscore.py evaluate --prefix $HOME/public_html/pathoscorepaper/clinvarall -s CCR -s CADD -s DANN -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s aloft_het -i aloft_lof -i aloft_rec -s fathmm_non -s fathmm_coding -s fitCons -s MetaSVM -s missense_z -s pLI -s phastCons -s polyphen2_hvar -s SIFT -s mis_badness -i BLOSUM -s Grantham --suffix pdf --goi pathogenicgenes.txt pathogenic.vcf.gz benign.vcf.gz

python ../pathoscore.py evaluate --prefix $HOME/public_html/pathoscorepaper/clinvar -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf --goi pathogenicgenes.txt pathogenic.vcf.gz benign.vcf.gz

# pathogenic filtered clinvar on functional and gnomAD
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvarfilter -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf pathogenic.filter.vcf.gz benign.vcf.gz

# clinvar 100% gnomad filtered plus functional
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvarbenignfilter -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf pathogenic.filter.vcf.gz benign.filter.vcf.gz

#########
#SAMOCHA#
#########

# regular samocha; functional filter doesn't matter, since there are no non-protein-altering variants.
python ../pathoscore.py evaluate --prefix $HOME/public_html/pathoscorepaper/samocha -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf --goi neurodevgenes.txt neurodev.vcf.gz control.vcf.gz

# samocha benigns, clinvar pathogenics, filtered functionally
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvarsamocha -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf pathogenic.vcf.gz control.vcf.gz

# filtered samocha
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/samochafilter -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf neurodev.filter.vcf.gz control.vcf.gz

#########
#HOMSY###
#########

# regular homsy
python ../pathoscore.py evaluate --prefix $HOME/public_html/pathoscorepaper/homsy -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf hompathogenic.vcf.gz hombenign.vcf.gz

# filtered homsy
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/homsyfilter -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf hompathogenic.filter.vcf.gz hombenign.vcf.gz

# homsy benigns, clinvar pathogenics, filtered functionally
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvarhomsy -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf pathogenic.vcf.gz hombenign.vcf.gz

#########
#WIGLER##
#########

# regular wigler
python ../pathoscore.py evaluate --prefix $HOME/public_html/pathoscorepaper/wigler -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf wigpathogenic.vcf.gz wigbenign.vcf.gz

# filtered wigler
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/wiglerfilter -s CCR -s CADD -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s MetaSVM -s Grantham --suffix pdf wigpathogenic.filter.vcf.gz wigbenign.vcf.gz
