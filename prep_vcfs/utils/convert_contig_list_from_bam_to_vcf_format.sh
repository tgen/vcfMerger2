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


function usage(){
	echo -e "This script creates the contig lines needed to be added to the Lancet's VCF using data from a BAM file"
	echo -e "\nUSAGE:"
	echo -e "$0 --bam \${BAM} --outfilename \${OUTFILENAME}"
	echo -e "OR"
	echo -e "$0 -i \${BAM} -o \${OUTFILENAME}"
	
}

function check_ev(){
	if [[ $1 -ne 0 ]] ; then echo -e "ERROR: ${2} FAILED ;  found ev=${1} ; Aborting! " ; exit -1 ; fi
}

function checkFile(){
	local F=$1
	if [[ ! -e ${F} ]] ; then echo -e "FILE NOT FOUND << ${F} >>; Aborting!" ; usage ; exit -1 ; fi ;
}

function init_some_vars(){
	LI="RECAP_INPUTS_USED:"
	BAM=""
	OUTFILENAME=""	
}

function getOptions(){ 


# options may be followed by one colon to indicate they have a required argument
if ! options=`getopt -o hi:o: -l help,bam:,outfilename: -- "$@" `
	then
	# something went wrong, getopt will put out an error message for us
		echo "ERROR in Arguments" ; usage
		exit -1
	fi
	eval set -- "$options"
	while [[ $# -gt 0 ]]
	do
		# for options with required arguments, an additional shift is required
		case $1 in
		-i|--bam) BAM=$2 ; LI="${LI}\nBAM==\"${BAM}\"";  shift ;; 
		-o|--outfilename) OUTFILENAME=$2 ; LI="${LI}\nOUTFILENAME==\"${OUTFILENAME}\"";  shift ;;
		-h|--help) usage ; exit ;;
		(--) shift ;;
		(-*) echo -e "$0: error - unrecognized option $1\n\n" 1>&2   ; usage;  exit -1  ;;
		(*) break ; echo "$0: error --- unrecognized option $1" 1>&2 ; usage;  exit -1  ;;
		esac
		shift
	done
	
	#input recap
	LI="${LI}\nCURR_DIR==\"${PWD}\""
	echo -e "\n\n+------------------------------------------------+\n${LI[@]}\n+------------------------------------------------+\n\n"

}

##@@@@@@@@##
## MAIN
##@@@@@@@@##

## init variables
init_some_vars
## get options
getOptions $@



for F in "${BAM}" "${OUTFILENAME}"
do
	if [[ ${F} == "" ]] ;
	then
		echo -e "ERROR: Missing filenames ; options --bam and/or --outfilename REQUIRED; Aborting.\nuse $0 --help for usage."
		exit -1
	fi
	
done

type samtools >/dev/null 2>&1 || { echo >&2 "Require \"\nsamtools\" executable but it's not in the PATH.  Aborting."; exit -1; }

checkFile ${BAM}

echo -e "running ...\nsamtools view -H ${BAM} | grep @SQ | sed 's/@SQ\t/##contig=<ID=/ ; s/SN:// ; s/\tLN:/,length=/ ; s/$/>/' > ${OUTFILENAME}" 1>&2
samtools view -H ${BAM} | grep @SQ | sed 's/@SQ\t/##contig=<ID=/ ; s/SN:// ; s/\tLN:/,length=/ ; s/$/>/' > ${OUTFILENAME}

check_ev $? "samtools command " 1>&2

echo -e "\nCreate Intermediate contigs file for vcfMerger2.0 pre-processing step completed succesfully\n"

