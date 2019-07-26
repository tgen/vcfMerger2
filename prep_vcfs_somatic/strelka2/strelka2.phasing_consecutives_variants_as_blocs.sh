#!/usr/bin/env bash

## NOTE: the VCF MUST HAVE BEEN PREPARED for VCFMERGER2 or having at least AR, DP and GT flags in the FORMAT columns
VCF=$1
TBAM=$2
SNAME_T=$3
CPUS=$4

set -euo pipefail
#module add python/3.6.0 samtools/1.9 R/3.4.1 BEDTools/2.26.0
MAX_CPUS_IN_CPUINFO=$(cat /proc/cpuinfo | grep processor | wc -l)
PHASER_EXE=/home/clegendre/tools/phASER/v1.1.1/phaser/phaser/phaser.py  ; ## phASER requires python 2.xx  (not 3)
echo ${PHASER_EXE}
DIR_PATH_TO_SCRIPTS="$( dirname `readlink -f $0` )"
SCRIPT_GET_CONSPOS=${DIR_PATH_TO_SCRIPTS}/get_consecutive_list_of_numbers_from_vcfFileInput.py
SCRIPT_PYTHON_RECOMPOSE_VARIANTS=${DIR_PATH_TO_SCRIPTS}/strelka2.recompose_phased_variants.py



## FUNCTIONS
function check_exe_in_path(){
	EXE=$1
	type $EXE >/dev/null 2>&1 || { echo >&2 "Require \"$EXE\" executable but it's not in the PATH.  Aborting."; exit -1; } || $EXE --help
}

function check_ev(){
	local ev=$1
	local msg=$(echo "$2" | sed 's/[ \+]/_/g')
	if [[ ${ev} -ne 0 ]] ; then touch FAILED_Recomposition_Strelka2_${msg}.flag ; echo -e "ERROR: ${msg} FAILED ;\nexit_value:${ev} ; Aborting! " ; exit -1 ; fi
}

function checkDir(){
	local D=$1
	if [[ ! -e ${D} ]] ; then echo -e "DIR NOT FOUND << ${D} >>; Aborting!" ; usage ; exit -1 ; fi ;
}

function checkFile(){
	local F=$1
	if [[ "" == "${VCF}" || ! -e ${F} ]] ; then echo -e "FILE NOT FOUND << ${F} >>; Aborting!" ; usage ; exit -1 ; fi ;
}

function usage(){
echo -e "USAGE: $0 \$VCF.gz \$TBAM \$SNAME_T \$CPUS \n1) compressed_VCF\n2) BAM file used to call variant\n3) Sample Name found in VCF\n4) Number of Threads used by phASEr tool (as much as you can unless I/O is slow on system ; limitation is the speed of reading BAM file ;   1 cpu will be 1 contig processed; so 10 cpus will process 10 contigs concurrently (aka 10 threads)\n
## NOTE: the VCF MUST HAVE BEEN PREPARE for VCFMERGER2 or having at least AR, DP and GT flags in the FORMAT columns"
}

## checkings section
if [[ $# -ne 4 ]] ; then usage ; exit -1 ; fi 
for F in ${SCRIPT_GET_CONSPOS} ${VCF} ${TBAM} ; do checkFile ${F} ; done
for EXE in samtools bcftools phaser.py ; do check_exe_in_path ${EXE} ; done
if [[ ${CPUS} -ge ${MAX_CPUS_IN_CPUINFO} ]] ; then CPUS=$((${MAX_CPUS_IN_CPUINFO}-1)) ; fi

VCF_ORIGINAL_INPUT=${VCF}

## WARNING: INPUT VCF MUST BE block-compressed vcf with extension .vcf.gz and associated with an index file
if [[ "23236669" == $( xxd ${VCF_ORIGINAL_INPUT} | head -n 1 | cut -d" " -f2-3 | sed 's/ \+//' ) ]]
then
    echo -e "we are dealing with a VCF file ..."
    if [[ ${VCF_ORIGINAL_INPUT##*.} == "vcf" ]]
    then
        bcftools view --threads 2 -O z -o ${VCF_ORIGINAL_INPUT}.gz ${VCF_ORIGINAL_INPUT}
        bcftools index --tbi ${VCF_ORIGINAL_INPUT}.gz
        VCF_ORIGINAL_INPUT=${VCF_ORIGINAL_INPUT}.gz
    else
        bcftools view --threads 2 -O z -o ${VCF_ORIGINAL_INPUT}.vcf.gz ${VCF_ORIGINAL_INPUT}
        bcftools index --tbi ${VCF_ORIGINAL_INPUT}.vcf.gz
        VCF_ORIGINAL_INPUT=${VCF_ORIGINAL_INPUT}.vcf.gz
    fi
elif [[ "1f8b0804" != $(xxd ${VCF_ORIGINAL_INPUT} | head -n 1 | cut -d" " -f2-3 | sed 's/ \+//' ) ]]
then
    echo -e "we deal with a .vcf.gz file type ...; checking extension "
    if [[ $(echo -e ${VCF_ORIGINAL_INPUT} | grep -o -m 1 -E ".vcf.gz$") != ".vcf.gz" ]]
    then
        mv ${VCF_ORIGINAL_INPUT} ${VCF_ORIGINAL_INPUT}.vcf.gz
        VCF_ORIGINAL_INPUT=${VCF_ORIGINAL_INPUT}.vcf.gz
    fi
else
    echo -e "ERROR: INPUT FILE is NOT a VALID VCF ; Aborting; "
    exit -1
fi

VCF=${VCF_ORIGINAL_INPUT}

if [[ 1 == 1 ]] ;then

echo -e "get consecutive positions ... as tabulated text file for bcftools ..."
python ${SCRIPT_GET_CONSPOS} ${VCF}
check_ev $? "$(basename ${SCRIPT_GET_CONSPOS})"

if [[ $(cat ${VCF}.consPos.txt | wc -l ) -lt 2 ]] ;
then
    echo -e "\nNo Consecutive Position found in VCF; ending Phasing section now"
    echo -e "renaming input file to match expected outfile from phasing section\n"
    mycmd="cp ${VCF_ORIGINAL_INPUT} ${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf.gz}"
    echo ${mycmd} ; eval ${mycmd}
    check_ev $? "cp command"
    echo "${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf.gz}"
    bcftools view -O v --threads 2 -o ${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf} ${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf.gz}
    exit 0
fi

echo -e "Subset VCF file with only the captured position form step above ..."
bcftools filter --threads 2 -O z -T ${VCF}.consPos.txt -o ${VCF/vcf.gz/subByConsPos.vcf.gz} ${VCF}
check_ev $? "bcftools filter #1"
bcftools index --threads 2 --tbi ${VCF/vcf.gz/subByConsPos.vcf.gz}
check_ev $? "bcftools index #1"

bcftools filter --threads 2 -O z -T ^${VCF}.consPos.txt -o ${VCF/vcf.gz/TempNoConsPos.vcf.gz} ${VCF}
check_ev $? "bcftools filter #2"
bcftools index --threads 2 --tbi ${VCF/vcf.gz/TempNoConsPos.vcf.gz}
check_ev $? "bcftools index #2"

fi

## https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
## https://github.com/secastel/phaser/tree/master/phaser
##@@@@@@@@@@@@@
### WARNING: Make sure the phaser.py file is in your path , is executable AND has the correct shebang (#!/usr/bin/env python2 )
##@@@@@@@@@@@@@
echo -e "Phasing the Consecutive positions using phASER ..."
echo -e "RUNNING phASER ... "
echo -e "Command:"
module load python/2.7.13
VCF_IN=${VCF/vcf.gz/subByConsPos.vcf.gz}
VCF_OUT_PHASER=phased_${VCF_IN}
TEMP_DIR=$(dirname ${VCF_IN} )/temp_phASER
mkdir -p ${TEMP_DIR}

mycmd="python ${PHASER_EXE} --bam ${TBAM} --sample ${SNAME_T}  --vcf ${VCF_IN}  --o ${VCF_OUT_PHASER}  --mapq 20 --baseq 20 --paired_end 1 --gw_phase_vcf 0  --remove_dups 1 --write_vcf 1  --gw_af_field AR --gw_phase_vcf 2 --temp_dir ${TEMP_DIR} --max_block_size 10 --threads ${CPUS}"  ## HARDCODED value for options here ; NOTE: 10 CPUS is a good trade-off between I/O and CPUS; Too amany cpus leads to too many samtools view which increases I/O and waiting time
echo -e "${mycmd}"
eval ${mycmd} | tee log_for_phASER_$(basename ${VCF_OUT_PHASER} ).log
check_ev $? "phASER using ${PHASER_EXE} with python-2.7.13"
module unload python/2.7.13



echo -e "Annotation of the subByConPos VCF with the vcf outputted by Phaser in order to keep a vcf with both NORMAL and TUMOR sample" ; ## INDEED, phASER only phases ONE sample at a time and we provided only the TUMOR sample
VCF_IN=${VCF_IN}
VCF_FOR_ANNO=${VCF_OUT_PHASER}.vcf.gz
VCF_OUT=${VCF_IN/.vcf.gz/.phased.vcf.gz}
mycmd="bcftools annotate --threads 2 -a ${VCF_FOR_ANNO} -c FORMAT/PG,FORMAT/PB,FORMAT/PW,FORMAT/PC,FORMAT/PM,FORMAT/PS,FORMAT/PI -m \"phASER\" -O z -o ${VCF_OUT} ${VCF_IN}"
echo ${mycmd}
eval ${mycmd}
check_ev $? "bcftools annotate "

bcftools index --threads 2 --tbi ${VCF_OUT}
check_ev $? "bcftools index"


echo -e "Merged as one Block the Consecutive newly-Phased Variants ..."
VCF_IN=${VCF_OUT} ## can be vcf or block-compressed vcf
VCF_OUT=${VCF_IN/vcf.gz/blocs.vcf}
python  ${SCRIPT_PYTHON_RECOMPOSE_VARIANTS} --tumor-col 11 --normal-col 10 -i ${VCF_IN} -o ${VCF_OUT}
check_ev $? "python strelka2.recompose_phased_variants.py"


echo -e "Block Size Distribution in new recomposed variants ..."
echo -e "block_size\tcount\tcount*blocksize\taccumulated_count"
for C in $(seq 1 10 ); do echo -e "${C}\t$(bcftools view -H ${VCF_OUT} | awk -v C=${C} 'length($4)==C && length($5)==C ' | wc -l )" ; done | awk ' {SUM+=$1*$2} { print $0,$1*$2,SUM } '
echo -e "indels_count\t$(bcftools view -H ${VCF_OUT} | awk 'length($4)<length($5) || length($4)>length($5) ' | wc -l )"


echo -e "sort and block-compressed vcf then index it ..."
VCF_IN=${VCF_OUT}
VCF_OUT=${VCF_IN}.gz
bcftools sort -Oz --max-mem 2G -o ${VCF_OUT} ${VCF_IN}
bcftools index --tbi --threads 2 ${VCF_OUT}
check_ev $? "bcftools index ${VCF_OUT}"



echo -e "Concatenating unphased and phased vcfs ... "
VCF_IN_UNPHASED=${VCF_ORIGINAL_INPUT/vcf.gz/TempNoConsPos.vcf.gz}
VCF_IN_PHASED=${VCF_OUT}
## WE OUTPUT An UNcompressed VCF;
VCF_OUT=${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf}
mycmd="bcftools concat -a -O v ${VCF_IN_UNPHASED} ${VCF_IN_PHASED} | bcftools sort -O v -T ${PWD} --max-mem 2G  - > ${VCF_OUT}"
echo ${mycmd}
eval ${mycmd}
check_ev $? "bcftools concat ${VCF_OUT}"

## WE OUTPUT An UNcompressed VCF; if we nee to go back to block-compressedwe will need to uncomment the two lines below
#bcftools index --tbi --threads 2 ${VCF_OUT}
#check_ev $? "bcftools index ${VCF_OUT}"

touch Completed_Recomposition_Strelka2_for_VCF_$(basename ${VCF_ORIGINAL_INPUT}).flag

echo ${VCF_OUT} 2>&1

exit 0
