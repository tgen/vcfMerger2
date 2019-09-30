#!/usr/bin/env bash

### vcfMerger2
###
### MIT License
###
### Copyright (c) 2018 Translational Genomics Research Institute
###
### Permission is hereby granted, free of charge, to any person obtaining a copy
### of this software and associated documentation files (the "Software"), to deal
### in the Software without restriction, including without limitation the rights
### to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
### copies of the Software, and to permit persons to whom the Software is
### furnished to do so, subject to the following conditions:
###
### The above copyright notice and this permission notice shall be included in all
### copies or substantial portions of the Software.
###
### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
### IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
### FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
### AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
### LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
### OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
### SOFTWARE.
###
### Major Contributors: Christophe Legendre 
### Minor Contributors:


SNPSIFT_JAR_PATH=$1
CODE_EXT_FILT=$2  # can be either "filt" or "pass"
STRING_FOR_SNPSIFT_FILTERING="$3"  ## double quote mandatory here
DIR_TEMP=$4
VCF=$5

if [[ ${VCF##*.} == "vcf" ]] ; then
    VCF_OUT=$(basename ${VCF} ".vcf").${CODE_EXT_FILT}.vcf
elif [[ ${VCF##*.} == "gz" ]] ; then
    VCF_OUT=$(basename ${VCF} ".vcf.gz").${CODE_EXT_FILT}.vcf
else
    echo -e "ERROR: Unexpected file extension. File Extension is neither .vcf or .vcf.gz; \
    Please Comply to the expected file extension of the VCF file; Aborting."
   exit 1
fi

set -eu

zcat -f ${VCF} | java -jar ${SNPSIFT_JAR_PATH} filter " ${STRING_FOR_SNPSIFT_FILTERING} " > ${DIR_TEMP}/${VCF_OUT}
if [[ $? -ne 0 ]] ;
then
    echo -e "ERROR: Filtering VCF step FAILED;\n Command  << zcat -f ${VCF} | java -jar ${SNPSIFT_JAR_PATH} filter ' ${STRING_FOR_SNPSIFT_FILTERING} ' > ${DIR_TEMP}/${VCF_OUT} >> FAILED"
    exit 1
fi