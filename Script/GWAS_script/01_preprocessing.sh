#!/bin/bash -l
#SBATCH -A 
#SBATCH -o 
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 06:00:00
#SBATCH -J 
#SBATCH --mail-type=ALL
#SBATCH --mail-user matias.becker-burgos.1399@student.uu.se

## USAGE:

## 01_preprocessing.sh <INPUTFILE> <IDFILE> <OUTPUTLOCATION>
############################################################################
################################## SCRIPT ##################################
############################################################################
##Author: Matias Becker Burgos (matias.becker-burgos.1399@student.uu.se).
############################################################################

############################################################################
# Load modules:

module load bioinfo-tools
module load vcftools
module load plink

############################################################################
# Input:

INPUTFILE=$(readlink -f $1)
IDFILE=$(readlink -f $2)
OUTPUTLOCATION=$(readlink -f $3)

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
# Output:

OUTPUTFILEPREFIX=$(echo "${OUTPUTLOCATION}/Microbiome_variant_sites")

TMPFILEPREFIX=$(echo "tmp-file")
TMPFILEPREFIX=$(echo "${OUTPUTLOCATION}/${TMPFILEPREFIX}")

############################################################################
# Commands:

# Keep only microbiome data individuals
vcftools --vcf ${INPUTFILE} \
 --keep ${IDFILE}\
 --recode --out ${TMPFILEPREFIX}_269individuals.vcf.gz

# Filtering low frequency allèles and missigness
vcftools --vcf ${TMPFILEPREFIX}_269individuals.vcf.gz \
 --maf 0.05 --max-missing 0.95 \
 --recode --out ${TMPFILEPREFIX}_filt.vcf.gz

 # Converting the vcf to PLINK format
plink --vcf ${TMPFILEPREFIX}_filt.vcf.gz --make-bed --allow-extra-chr --out ${OUTPUTFILEPREFIX}
