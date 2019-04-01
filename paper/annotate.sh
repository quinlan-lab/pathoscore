#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o %j-%N.out
#SBATCH -e %j-%N.err
#SBATCH --time=48:00:00

#########
#CLINVAR#
#########

python ../pathoscore.py annotate --scores ../score-sets/GRCh37/CCR/ccrs.autosomes.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CCR/ccrs.xchrom.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/CADD/InDels.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/DANN/DANN_whole_genome_SNVs.tsv.bgz:DANN:5:max --scores ../score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores ../score-sets/GRCh37/MCAP/mcap.txt.gz:MCAP:5:max --scores ../score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores ../score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores ../score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores ../score-sets/GRCh37/RVIS/rvis.bed.gz:RVIS:5:max --scores ../score-sets/GRCh37/VVP/VVP_scores_allChr.txt.gz:VVP:5:max --scores ../score-sets/GRCh37/aloft/aloft.txt.gz:aloft_het,aloft_lof,aloft_rec:5,6,7:max,max,max --scores ../score-sets/GRCh37/fathmm/fathmm/fathmm.txt.gz:fathmm_non,fathmm_coding:5,6:max,max --scores ../score-sets/GRCh37/fitcons/fitcons/fitcons.bed.gz:fitCons:4:max --scores ../score-sets/GRCh37/metasvm/metasvm.txt.gz:MetaSVM:5:max --scores ../score-sets/GRCh37/mis_Z/missensez.bed.gz:missense_z:5:max --scores ../score-sets/GRCh37/pLI/pLI.bed.gz:pLI:5:max --scores ../score-sets/GRCh37/phastCons/phastcons_vertebrate.txt.gz:phastCons:4:max --scores ../score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:polyphen2_hvar:6:max --scores ../score-sets/GRCh37/sift/sift.txt.gz:SIFT:12:max ../truth-sets/GRCh37/clinvar/clinvar-benign.20170905.vcf.gz --prefix benign

python ../pathoscore.py annotate --scores ../score-sets/GRCh37/CCR/ccrs.autosomes.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CCR/ccrs.xchrom.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/CADD/InDels.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/DANN/DANN_whole_genome_SNVs.tsv.bgz:DANN:5:max --scores ../score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores ../score-sets/GRCh37/MCAP/mcap.txt.gz:MCAP:5:max --scores ../score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores ../score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores ../score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores ../score-sets/GRCh37/RVIS/rvis.bed.gz:RVIS:5:max --scores ../score-sets/GRCh37/VVP/VVP_scores_allChr.txt.gz:VVP:5:max --scores ../score-sets/GRCh37/aloft/aloft.txt.gz:aloft_het,aloft_lof,aloft_rec:5,6,7:max,max,max --scores ../score-sets/GRCh37/fathmm/fathmm/fathmm.txt.gz:fathmm_non,fathmm_coding:5,6:max,max --scores ../score-sets/GRCh37/fitcons/fitcons/fitcons.bed.gz:fitCons:4:max --scores ../score-sets/GRCh37/metasvm/metasvm.txt.gz:MetaSVM:5:max --scores ../score-sets/GRCh37/mis_Z/missensez.bed.gz:missense_z:5:max --scores ../score-sets/GRCh37/pLI/pLI.bed.gz:pLI:5:max --scores ../score-sets/GRCh37/phastCons/phastcons_vertebrate.txt.gz:phastCons:4:max --scores ../score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:polyphen2_hvar:6:max --scores ../score-sets/GRCh37/sift/sift.txt.gz:SIFT:12:max  ../truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170905.vcf.gz --prefix pathogenic

# ad gene files
python ../pathoscore.py annotate pathogenic.vcf.gz --exclude ../gene-sets/GRCh37/ad_genes/ad_gene_complement.bed.gz  --prefix adgene.pathogenic
python ../pathoscore.py annotate benign.vcf.gz --exclude ../gene-sets/GRCh37/ad_genes/ad_gene_complement.bed.gz  --prefix adgene.benign

# ar gene files
python ../pathoscore.py annotate pathogenic.vcf.gz --exclude ../gene-sets/GRCh37/ar_genes/ar_gene_complement.bed.gz  --prefix argene.pathogenic
python ../pathoscore.py annotate benign.vcf.gz --exclude ../gene-sets/GRCh37/ar_genes/ar_gene_complement.bed.gz  --prefix argene.benign

# filtered pathogenics/benigns on ClinVar
python ../pathoscore.py annotate pathogenic.vcf.gz --exclude ../scripts/gnomad/ExAC.vcf.gz --exclude ../scripts/gnomad/gnomad.vcf.gz --prefix pathogenic.filter
python ../pathoscore.py annotate benign.vcf.gz --exclude ../scripts/gnomad/ExAC.vcf.gz --exclude ../scripts/gnomad/gnomad.vcf.gz --prefix benign.filter

#########
#SAMOCHA#
#########

python ../pathoscore.py annotate --scores ../score-sets/GRCh37/CCR/ccrs.autosomes.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CCR/ccrs.xchrom.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/CADD/InDels.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/DANN/DANN_whole_genome_SNVs.tsv.bgz:DANN:5:max --scores ../score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores ../score-sets/GRCh37/MCAP/mcap.txt.gz:MCAP:5:max --scores ../score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores ../score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores ../score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores ../score-sets/GRCh37/RVIS/rvis.bed.gz:RVIS:5:max --scores ../score-sets/GRCh37/VVP/VVP_scores_allChr.txt.gz:VVP:5:max --scores ../score-sets/GRCh37/aloft/aloft.txt.gz:aloft_het,aloft_lof,aloft_rec:5,6,7:max,max,max --scores ../score-sets/GRCh37/fathmm/fathmm/fathmm.txt.gz:fathmm_non,fathmm_coding:5,6:max,max --scores ../score-sets/GRCh37/fitcons/fitcons/fitcons.bed.gz:fitCons:4:max --scores ../score-sets/GRCh37/metasvm/metasvm.txt.gz:MetaSVM:5:max --scores ../score-sets/GRCh37/mis_Z/missensez.bed.gz:missense_z:5:max --scores ../score-sets/GRCh37/pLI/pLI.bed.gz:pLI:5:max --scores ../score-sets/GRCh37/phastCons/phastcons_vertebrate.txt.gz:phastCons:4:max --scores ../score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:polyphen2_hvar:6:max --scores ../score-sets/GRCh37/sift/sift.txt.gz:SIFT:12:max ../truth-sets/GRCh37/samocha/samocha.benign.vcf.gz --prefix control

python ../pathoscore.py annotate --scores ../score-sets/GRCh37/CCR/ccrs.autosomes.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CCR/ccrs.xchrom.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/CADD/InDels.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/DANN/DANN_whole_genome_SNVs.tsv.bgz:DANN:5:max --scores ../score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores ../score-sets/GRCh37/MCAP/mcap.txt.gz:MCAP:5:max --scores ../score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores ../score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores ../score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores ../score-sets/GRCh37/RVIS/rvis.bed.gz:RVIS:5:max --scores ../score-sets/GRCh37/VVP/VVP_scores_allChr.txt.gz:VVP:5:max --scores ../score-sets/GRCh37/aloft/aloft.txt.gz:aloft_het,aloft_lof,aloft_rec:5,6,7:max,max,max --scores ../score-sets/GRCh37/fathmm/fathmm/fathmm.txt.gz:fathmm_non,fathmm_coding:5,6:max,max --scores ../score-sets/GRCh37/fitcons/fitcons/fitcons.bed.gz:fitCons:4:max --scores ../score-sets/GRCh37/metasvm/metasvm.txt.gz:MetaSVM:5:max --scores ../score-sets/GRCh37/mis_Z/missensez.bed.gz:missense_z:5:max --scores ../score-sets/GRCh37/pLI/pLI.bed.gz:pLI:5:max --scores ../score-sets/GRCh37/phastCons/phastcons_vertebrate.txt.gz:phastCons:4:max --scores ../score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:polyphen2_hvar:6:max --scores ../score-sets/GRCh37/sift/sift.txt.gz:SIFT:12:max ../truth-sets/GRCh37/samocha/samocha.pathogenic.vcf.gz --prefix neurodev

# filtered pathogenics on samocha
python ../pathoscore.py annotate neurodev.vcf.gz --exclude ../scripts/gnomad/ExAC.vcf.gz --exclude ../scripts/gnomad/gnomad.vcf.gz --prefix neurodev.filter

#######
#HOMSY#
#######

python ../pathoscore.py annotate --scores ../score-sets/GRCh37/CCR/ccrs.autosomes.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CCR/ccrs.xchrom.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/CADD/InDels.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/DANN/DANN_whole_genome_SNVs.tsv.bgz:DANN:5:max --scores ../score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores ../score-sets/GRCh37/MCAP/mcap.txt.gz:MCAP:5:max --scores ../score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores ../score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores ../score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores ../score-sets/GRCh37/RVIS/rvis.bed.gz:RVIS:5:max --scores ../score-sets/GRCh37/VVP/VVP_scores_allChr.txt.gz:VVP:5:max --scores ../score-sets/GRCh37/aloft/aloft.txt.gz:aloft_het,aloft_lof,aloft_rec:5,6,7:max,max,max --scores ../score-sets/GRCh37/fathmm/fathmm/fathmm.txt.gz:fathmm_non,fathmm_coding:5,6:max,max --scores ../score-sets/GRCh37/fitcons/fitcons/fitcons.bed.gz:fitCons:4:max --scores ../score-sets/GRCh37/metasvm/metasvm.txt.gz:MetaSVM:5:max --scores ../score-sets/GRCh37/mis_Z/missensez.bed.gz:missense_z:5:max --scores ../score-sets/GRCh37/pLI/pLI.bed.gz:pLI:5:max --scores ../score-sets/GRCh37/phastCons/phastcons_vertebrate.txt.gz:phastCons:4:max --scores ../score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:polyphen2_hvar:6:max --scores ../score-sets/GRCh37/sift/sift.txt.gz:SIFT:12:max ../truth-sets/GRCh37/homsy/homsy.benign.vcf.gz --prefix hombenign

python ../pathoscore.py annotate --scores ../score-sets/GRCh37/CCR/ccrs.autosomes.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CCR/ccrs.xchrom.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/CADD/InDels.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/DANN/DANN_whole_genome_SNVs.tsv.bgz:DANN:5:max --scores ../score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores ../score-sets/GRCh37/MCAP/mcap.txt.gz:MCAP:5:max --scores ../score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores ../score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores ../score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores ../score-sets/GRCh37/RVIS/rvis.bed.gz:RVIS:5:max --scores ../score-sets/GRCh37/VVP/VVP_scores_allChr.txt.gz:VVP:5:max --scores ../score-sets/GRCh37/aloft/aloft.txt.gz:aloft_het,aloft_lof,aloft_rec:5,6,7:max,max,max --scores ../score-sets/GRCh37/fathmm/fathmm/fathmm.txt.gz:fathmm_non,fathmm_coding:5,6:max,max --scores ../score-sets/GRCh37/fitcons/fitcons/fitcons.bed.gz:fitCons:4:max --scores ../score-sets/GRCh37/metasvm/metasvm.txt.gz:MetaSVM:5:max --scores ../score-sets/GRCh37/mis_Z/missensez.bed.gz:missense_z:5:max --scores ../score-sets/GRCh37/pLI/pLI.bed.gz:pLI:5:max --scores ../score-sets/GRCh37/phastCons/phastcons_vertebrate.txt.gz:phastCons:4:max --scores ../score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:polyphen2_hvar:6:max --scores ../score-sets/GRCh37/sift/sift.txt.gz:SIFT:12:max ../truth-sets/GRCh37/homsy/homsy.pathogenic.vcf.gz --prefix hompathogenic

# filtered pathogenics on homsy
python ../pathoscore.py annotate hompathogenic.vcf.gz --exclude ../scripts/gnomad/ExAC.vcf.gz --exclude ../scripts/gnomad/gnomad.vcf.gz --prefix hompathogenic.filter

########
#WIGLER#
########

python ../pathoscore.py annotate --scores ../score-sets/GRCh37/CCR/ccrs.autosomes.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CCR/ccrs.xchrom.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/CADD/InDels.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/DANN/DANN_whole_genome_SNVs.tsv.bgz:DANN:5:max --scores ../score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores ../score-sets/GRCh37/MCAP/mcap.txt.gz:MCAP:5:max --scores ../score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores ../score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores ../score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores ../score-sets/GRCh37/RVIS/rvis.bed.gz:RVIS:5:max --scores ../score-sets/GRCh37/VVP/VVP_scores_allChr.txt.gz:VVP:5:max --scores ../score-sets/GRCh37/aloft/aloft.txt.gz:aloft_het,aloft_lof,aloft_rec:5,6,7:max,max,max --scores ../score-sets/GRCh37/fathmm/fathmm/fathmm.txt.gz:fathmm_non,fathmm_coding:5,6:max,max --scores ../score-sets/GRCh37/fitcons/fitcons/fitcons.bed.gz:fitCons:4:max --scores ../score-sets/GRCh37/metasvm/metasvm.txt.gz:MetaSVM:5:max --scores ../score-sets/GRCh37/mis_Z/missensez.bed.gz:missense_z:5:max --scores ../score-sets/GRCh37/pLI/pLI.bed.gz:pLI:5:max --scores ../score-sets/GRCh37/phastCons/phastcons_vertebrate.txt.gz:phastCons:4:max --scores ../score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:polyphen2_hvar:6:max --scores ../score-sets/GRCh37/sift/sift.txt.gz:SIFT:12:max ../truth-sets/GRCh37/wigler-autism/wigler-autism.benign.vcf.gz --prefix wigbenign

python ../pathoscore.py annotate --scores ../score-sets/GRCh37/CCR/ccrs.autosomes.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CCR/ccrs.xchrom.v2.20180420.bed.gz:CCR:4:max --scores ../score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/CADD/InDels.tsv.gz:CADD:6:max --scores ../score-sets/GRCh37/DANN/DANN_whole_genome_SNVs.tsv.bgz:DANN:5:max --scores ../score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores ../score-sets/GRCh37/MCAP/mcap.txt.gz:MCAP:5:max --scores ../score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores ../score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores ../score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores ../score-sets/GRCh37/RVIS/rvis.bed.gz:RVIS:5:max --scores ../score-sets/GRCh37/VVP/VVP_scores_allChr.txt.gz:VVP:5:max --scores ../score-sets/GRCh37/aloft/aloft.txt.gz:aloft_het,aloft_lof,aloft_rec:5,6,7:max,max,max --scores ../score-sets/GRCh37/fathmm/fathmm/fathmm.txt.gz:fathmm_non,fathmm_coding:5,6:max,max --scores ../score-sets/GRCh37/fitcons/fitcons/fitcons.bed.gz:fitCons:4:max --scores ../score-sets/GRCh37/metasvm/metasvm.txt.gz:MetaSVM:5:max --scores ../score-sets/GRCh37/mis_Z/missensez.bed.gz:missense_z:5:max --scores ../score-sets/GRCh37/pLI/pLI.bed.gz:pLI:5:max --scores ../score-sets/GRCh37/phastCons/phastcons_vertebrate.txt.gz:phastCons:4:max --scores ../score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:polyphen2_hvar:6:max --scores ../score-sets/GRCh37/sift/sift.txt.gz:SIFT:12:max ../truth-sets/GRCh37/wigler-autism/wigler-autism.pathogenic.vcf.gz --prefix wigpathogenic

# filtered pathogenics on wigler
python ../pathoscore.py annotate wigpathogenic.vcf.gz --exclude ../scripts/gnomad/ExAC.vcf.gz --exclude ../scripts/gnomad/gnomad.vcf.gz --prefix wigpathogenic.filter