#!/usr/bin/env bash

SNPSIFT_JAR_PATH=$1
STRING_FOR_SNPSIFT_FILTERING="$2"  ## double quote mandatory here
VCF=$3

if [[ ${VCF##*.} == "vcf" ]] ; then
    VCF_OUT=$(basename ${VCF} ".vcf").filt.vcf
elif [[ ${VCF##*.} == "gz" ]] ; then
    VCF_OUT=$(basename ${VCF} ".vcf.gz").filt.vcf
else
    echo -e "ERROR: Unexpected file extension. File Extension is neither .vcf or .vcf.gz; \
    Please Comply to the expected file extension of the VCF file; Aborting."
   exit 1
fi

set -eu

zcat -f ${VCF} | java -jar ${SNPSIFT_JAR_PATH} filter " ${STRING_FOR_SNPSIFT_FILTERING} " > ${VCF_OUT}
if [[ $? -ne 0 ]] ;
then
    echo -e "ERROR: Filtering VCF step FAILED;\n Command  << zcat -f ${VCF} | java -jar ${SNPSIFT_JAR_PATH} filter ' ${STRING_FOR_SNPSIFT_FILTERING} ' > ${VCF_OUT} >> FAILED"
    exit 1
fi