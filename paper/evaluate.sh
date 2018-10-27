#########
#CLINVAR#
#########

# to get AD variant data and ROCs
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvarad -s CCR -s CADD -s DANN -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s aloft_het -i aloft_lof -i aloft_rec -s fathmm_non -s fathmm_coding -s fitCons -s MetaSVM -s missense_z -s pLI -s phastCons -s polyphen2_hvar -s SIFT --suffix pdf adgene.pathogenic.vcf.gz adgene.benign.vcf.gz

# to get AR variant data and ROCs
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvarar -s CCR -s CADD -s DANN -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s aloft_het -i aloft_lof -i aloft_rec -s fathmm_non -s fathmm_coding -s fitCons -s MetaSVM -s missense_z -s pLI -s phastCons -s polyphen2_hvar -s SIFT --suffix pdf argene.pathogenic.vcf.gz argene.benign.vcf.gz

# regular clinvar
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvar -s CCR -s CADD -s DANN -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s aloft_het -i aloft_lof -i aloft_rec -s fathmm_non -s fathmm_coding -s fitCons -s MetaSVM -s missense_z -s pLI -s phastCons -s polyphen2_hvar -s SIFT --suffix pdf pathogenic.vcf.gz benign.vcf.gz

# filtered clinvar
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvarfilter -s CCR -s CADD -s DANN -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s aloft_het -i aloft_lof -i aloft_rec -s fathmm_non -s fathmm_coding -s fitCons -s MetaSVM -s missense_z -s pLI -s phastCons -s polyphen2_hvar -s SIFT --suffix pdf pathogenic.filter.vcf.gz benign.vcf.gz

#########
#SAMOCHA#
#########

# regular samocha
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/samocha -s CCR -s CADD -s DANN -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s aloft_het -i aloft_lof -i aloft_rec -s fathmm_non -s fathmm_coding -s fitCons -s MetaSVM -s missense_z -s pLI -s phastCons -s polyphen2_hvar -s SIFT --suffix pdf neurodev.vcf.gz control.vcf.gz

# samocha benigns, clinvar pathogenics
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/clinvarsamocha -s CCR -s CADD -s DANN -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s aloft_het -i aloft_lof -i aloft_rec -s fathmm_non -s fathmm_coding -s fitCons -s MetaSVM -s missense_z -s pLI -s phastCons -s polyphen2_hvar -s SIFT --suffix pdf pathogenic.vcf.gz control.vcf.gz

# filtered samocha
python ../pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscorepaper/samochafilter -s CCR -s CADD -s DANN -s GERP -s MCAP -s MPC -i MTR -s REVEL -i RVIS -s VVP -s aloft_het -i aloft_lof -i aloft_rec -s fathmm_non -s fathmm_coding -s fitCons -s MetaSVM -s missense_z -s pLI -s phastCons -s polyphen2_hvar -s SIFT --suffix pdf neurodev.filter.vcf.gz control.vcf.gz
