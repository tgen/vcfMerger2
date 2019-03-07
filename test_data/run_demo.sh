#!/bin/bash

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

##
## DO NOT MOVE this script from its original location [[ aka vcfMerger2/text_data ]]
## 


## capturing directory from where this script is run
DIR_WORK="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
echo DIR_WORK == ${DIR_WORK}
cd ${DIR_WORK}
echo "curr_dir: ${PWD}"


if [[ ! -e "raw_tool_vcfs/strelka2.raw.vcf" ]]
then
	echo "We are not in the Expected directory to run current script. Please go manually to the vcfMerger2 installation directory and 'cd' into 'test_data' directory. Thanks. Sorry for the inconvenience."
	exit 1
fi

type python >/dev/null 2>&1 || { echo >&2 "Require \"\python\" executable but it's not in the PATH.  Aborting."; exit 1; } || python -V 
if [[ ! `python3 -V 2>&1 | sed 's/Python //g' | cut -d"." -f1 ` == "3" ]] ; then echo -e "ERROR: Python 3 or up
Expected in PATH; Aborting " ; exit 1 ; fi

type vt >/dev/null 2>&1 || { echo >&2 "Require \"vt\" executable but it's not in the PATH.  Aborting."; exit 1; } || vt --version
type bcftools >/dev/null 2>&1 || { echo >&2 "Require \"bcftools\" executable but it's not in the PATH.  Aborting."; exit 1; }
if [[ $( echo "`bcftools --version-only  2>&1 | sed 's/+.*//g'` <  1.7 " | bc -l ) -eq 1  ]] ; then echo -e "ERROR: bcftools 1.7 or up Expected in PATH; Aborting " ; exit 1 ; fi 
type samtools >/dev/null 2>&1 || { echo >&2 "Require \"samtools\" executable but it's not in the PATH.  Aborting."; exit 1; }

## uncompressing reference genome
cd ref_genome
if [[  -e grch37.22.fa && -e grch37.22.fa.fai  ]] ; 
then 
	echo "reference grch37.22.fa and grch37.22.fa.fai FOUND";
else
	tar -xzf genome.tar.gz
	if [[ $? -ne 0 ]] ; then echo "untar genome FAILED. check inputs" ; exit 1 ; fi 
fi
if [[ ! -e grch37.22.fa || ! -e grch37.22.fa.fai  ]] ; then echo -e "ERROR: grch37.22.fa or grch37.22.fa.fai NOT FOUND ; check directories and files" ; exit 1 ; fi 
cd ..


mkdir demo_out
cd demo_out

echo "Directory 'demo_out' has been created and has become our working directory"
echo "All the files to be used should be relative to that directory now;"
echo "TIPS: symlink or copy all the needed files in here or use relative paths to the file of interest [we are going to use relative paths in our example]"
echo
echo "Checking if vcfMerger2 executable is in expected folder relative to our current directory"
BIN_VM2="../../bin/vcfMerger2.py"
if [[ ! -e "${BIN_VM2}" ]]
then
	echo "vcfMerger2.py file NOT FOUND. Please go manually to the vcfMerger2 installation directory and 'cd' into 'test_data' directory, and restart the script. Thanks. Sorry for the inconvenience."
	exit 1
else
	echo "vcfMerger2.py FOUND"
fi

echo "BASH `date` ---> Starting vcfMerger2 ..."

python ../../bin/vcfMerger2.py --toolnames "strelka2|mutect2|lancet|octopus" --vcfs "../raw_tool_vcfs/strelka2.raw.vcf|../raw_tool_vcfs/mutect2.raw.vcf|../raw_tool_vcfs/lancet.raw.vcf|../raw_tool_vcfs/octopus.raw.vcf" -a "SLK|MUT|LAN|OCT"  -g ../ref_genome/grch37.22.fa --prep-outfilenames "strelka2.prepped.vcf|mutect2.prepped.vcf|lancet.prepped.vcf|octopus.prepped.vcf" --normal-sname NORMAL  --tumor-sname TUMOR -o merged.vcf --contigs-file "||../contigs/contigs.txt|../contigs/contigs.txt"

if [[ $? -eq 0 ]] ; 
then
	echo "BASH `date` ---> ALL RESULTS can be found in:  ${PWD}"
	echo "BASH `date` ---> vcfMerger2 completed Successfully"
else
	echo "BASH `date` ---> see LOGS in:  ${PWD}"
	echo "BASH `date` ---> vcfMerger2 FAILED" ; 
	exit 1
fi

