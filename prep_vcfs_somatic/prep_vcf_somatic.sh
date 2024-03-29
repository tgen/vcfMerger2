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

## set -ue ## DO not uncomment this line; it you do so, the script does not work anymore. Why? unknown
set -uo pipefail

## trap to capture the exit value from a function and from within a function
trap "exit 1" TERM
export TOP_PID=$$


## CONSTANT VARIABLE (modified accordingly)
export DIR_PATH_TO_SCRIPTS="$( dirname `readlink -f $0` )"
echo -e "DIR_PATH_TO_SCRIPTS=${DIR_PATH_TO_SCRIPTS}"

## CONSTANT VARIABLE : add of modify toolnames accordingly
VALID_TOOLNAMES="lancet,  mutect2,  octopus,  strelka2, vardict, vardictjava, VarDictJava, or their corresponding abbreviations, LAN|LCT, MUT,
OCT, SLK, VDJ|VDT, respectively [case Insensitive]"  ## if tools are later added, we will update this variable along
# with the
# function run_tool(), where the case statement will need to be updated.

###@@@@@@@@@@@@@@
### START HERE
###@@@@@@@@@@@@@@

type python3 >/dev/null 2>&1 || { echo >&2 "Require \"python3\" executable but it's not in the PATH.  Aborting."; exit
1; } || python3 -V

for NEEDED_MODULE in cyvcf2 intervene collections argparse getopt gzip json logging PIL natsort shutil subprocess warnings
do
    python3 -c "import ${NEEDED_MODULE}"
    if [[ $? -ne 0 ]] ; then echo -e "ERROR: module << ${NEEDED_MODULE} NOT FOUND >> with this python version `python3 -V`; \nPlease Install module ${NEEDED_MODULE} or change
    the python version accordingly; ABORTING\n" ; fi
done

#python_main_version_number=`python3 -V 2>&1 | sed 's/Python //g' | cut -d"." -f1 `
#if [[ ! "${python_main_version_number}" == "3" ]] ; then echo -e "ERROR: Python 3 or up Expected in PATH; Aborting "
# ; exit 1 ; fi
type vt >/dev/null 2>&1 || { echo >&2 "Require \"vt\" executable but it's not in the PATH.  Aborting."; exit 1; } ||
vt --version
type bcftools >/dev/null 2>&1 || { echo >&2 "Require \"bcftools\" executable but it's not in the PATH.  Aborting."; exit 1; }
#if [[ $( echo "`bcftools --version-only  2>&1 | sed 's/+.*//g'` <  1.7 " | bc -l ) -eq 1  ]] ; then echo -e "ERROR: bcftools 1.7 or up Expected in PATH; Aborting " ; exit 1 ; fi

# sourcing functions
source ${DIR_PATH_TO_SCRIPTS}/prep_vcf_somatic_functions.sh
## init variables
init_some_vars
## get options given by user
getOptions $@
## recap inputs
recap_input
## now we have all our inputs set up, let's process data ...
main
check_ev $? "main function "
exit
