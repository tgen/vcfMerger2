#!/usr/bin/env bash

## NOTE: the VCF MUST HAVE BEEN PREPARED for VCFMERGER2 or having at least AR, DP and GT flags in the FORMAT columns
VCF=$1
TBAM=$2
SNAME_T=$3
CPUS=$4

set -euo pipefail
#module add python/3.6.0 samtools/1.9 R/3.4.1 BEDTools/2.26.0
MAX_CPUS_IN_CPUINFO=$(cat /proc/cpuinfo | grep processor | wc -l)
PHASER_EXE=phaser.py  ; ## phASER requires python2.xx  (not 3)
echo -e "expected phASER executable: ${PHASER_EXE}"
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
## NOTE: the VCF MUST HAVE BEEN PREPARED for VCFMERGER2 or having at least AR, DP and GT flags in the FORMAT columns"
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

		echo -e "Keeping out the Indels as phASER exclude them from phasing anyway"   1>&2
		## We had to do this due to encountering an edge case mentioned in issue #27 in github
    bcftools filter -O z -i 'TYPE="indel"' -o ${VCF/vcf.gz/indels.vcf.gz} ${VCF}
    check_ev $? "bcftools filter out indels"
    bcftools index --tbi  ${VCF/vcf.gz/indels.vcf.gz}
    check_ev $? "bcftools index indels calls"
    VCF_INDELS_ONLY=${VCF/vcf.gz/indels.vcf.gz}

    ## now Subsetting the VCF_SBCP to continue without the indels
    bcftools filter -O z -e 'TYPE="indel"' -o ${VCF/vcf.gz/noindels.vcf.gz} ${VCF}
    check_ev $? "bcftools filter out indels"
    bcftools index --tbi  ${VCF/vcf.gz/noindels.vcf.gz}
    check_ev $? "bcftools index indels calls"
    cp ${VCF} ${VCF}.original_input.vcf
    VCF_NO_INDELS_FOR_PHASER=${VCF/vcf.gz/noindels.vcf.gz}
    #cp ${VCF/vcf.gz/noindels.vcf.gz} ${VCF}

    echo -e "get consecutive positions ... as tabulated text file for bcftools ..."
    python3 ${SCRIPT_GET_CONSPOS} ${VCF_NO_INDELS_FOR_PHASER}
    check_ev $? "$(basename ${SCRIPT_GET_CONSPOS})"

    if [[ $(cat ${VCF_NO_INDELS_FOR_PHASER}.consPos.txt | wc -l ) -lt 2 ]] ;
    then
        echo -e "\nNo Consecutive Position found in VCF; ending Phasing section now"
        echo -e "renaming input file to match expected outfile from phasing section\n"
        mycmd="cp ${VCF_ORIGINAL_INPUT} ${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf.gz}"
        echo ${mycmd} ; eval ${mycmd}
        check_ev $? "cp command"
        echo "${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf.gz}  ; Converting vcf.gz to vcf ..."
        bcftools view -O v --threads 2 -o ${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf} ${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf.gz}
        check_ev $? "bcftools view"
        exit 0
    fi

    echo -e "Subset VCF file with only the captured position from step above ..."
    bcftools filter --threads 2 -O z -T ${VCF_NO_INDELS_FOR_PHASER}.consPos.txt -o ${VCF_NO_INDELS_FOR_PHASER/vcf.gz/subByConsPos.vcf.gz} ${VCF_NO_INDELS_FOR_PHASER}
    check_ev $? "bcftools filter #1"
    bcftools index --threads 2 --tbi ${VCF_NO_INDELS_FOR_PHASER/vcf.gz/subByConsPos.vcf.gz}
    check_ev $? "bcftools index #1"

    bcftools filter --threads 2 -O z -T ^${VCF_NO_INDELS_FOR_PHASER}.consPos.txt -o ${VCF_NO_INDELS_FOR_PHASER/vcf.gz/TempNoConsPos.vcf.gz} ${VCF_NO_INDELS_FOR_PHASER}
    check_ev $? "bcftools filter #2"
    bcftools index --threads 2 --tbi ${VCF_NO_INDELS_FOR_PHASER/vcf.gz/TempNoConsPos.vcf.gz}
    check_ev $? "bcftools index #2"

    ## As phASER dose not deal with HOMOZYGOUS, we need to filter out the homs from subByConsPos.vcf.gz VCF_NO_INDELS_FOR_PHASER and recheck of number of ConsPOs lt 2
    ## due to the edge case encounter with MMRF_1073, we need to exclude manually ALL the Homozygous variant as phaser
    ## does not phase them and raises an error doing so if only homozygous variant are within the input vcf
    ## value for GT according to BCFTOOLS documentation
    # sample genotype: reference (haploid or diploid), alternate (hom or het, haploid or diploid), missing genotype, homozygous, heterozygous, haploid, ref-ref hom, alt-alt hom, ref-alt het, alt-alt het, haploid ref, haploid alt (case-insensitive)
    #GT="ref"
    #GT="alt"
    #GT="mis"
    #GT="hom"
    #GT="het"
    #GT="hap"
    #GT="RR"
    #GT="AA"
    #GT="RA" or GT="AR"
    #GT="Aa" or GT="aA"
    #GT="R"
    #GT="A"

    VCF_SBCP=${VCF_NO_INDELS_FOR_PHASER/vcf.gz/subByConsPos.vcf.gz}

    echo -e "removing ALT-ALT homozygous from subByConsPos VCF ...  using bcftools ..."
    bcftools filter -O z -e 'GT="AA"' -o ${VCF_SBCP/vcf.gz/hets.vcf.gz} ${VCF_SBCP}
    check_ev $? "bcftools filter out homz sites"

    VCF_HETS=${VCF_SBCP/vcf.gz/hets.vcf.gz}
    bcftools index --tbi ${VCF_HETS}
    check_ev $? "bcftools index"

    echo -e "keeping ALT-ALT homozygous from subByConsPos VCF ...  using bcftools ..."  1>&2
    bcftools filter -O z -i 'GT="AA"' -o ${VCF_SBCP/vcf.gz/homs.vcf.gz} ${VCF_SBCP}
    check_ev $? "bcftools filter out homz sites"

    VCF_HOMS=${VCF_SBCP/vcf.gz/homs.vcf.gz}
    bcftools index --tbi ${VCF_HOMS}
    check_ev $? "bcftools index"


    if [[ $( bcftools view -H ${VCF_HOMS} | wc -l ) -gt 0 ]] ;
    then
        echo -e "concat variants from homs.vcf.gz VCF with tempNoConsPos.vcf.gz VCF ... "  1>&2
        bcftools concat -a --threads 2 -O z -o ${VCF_NO_INDELS_FOR_PHASER/vcf.gz/TempNoConsPos_with_homs.vcf.gz} ${VCF_NO_INDELS_FOR_PHASER/vcf.gz/TempNoConsPos.vcf.gz}  ${VCF_HOMS}
        check_ev $? "bcftools concat #3"
        mv ${VCF_NO_INDELS_FOR_PHASER/vcf.gz/TempNoConsPos_with_homs.vcf.gz} ${VCF_NO_INDELS_FOR_PHASER/vcf.gz/TempNoConsPos.vcf.gz}
        check_ev $? "move #3"

        bcftools index --threads 2 -f --tbi ${VCF_NO_INDELS_FOR_PHASER/vcf.gz/TempNoConsPos.vcf.gz}
        check_ev $? "bcftools index #3"
    fi



    ## from now we assume that even if some variants were homozygous within the consecutives extracted varianats, not all of
    ## them will be and we run phASER from the file VCF named: VCF_HETS=${VCF_SBCP/vcf.gz/hets.vcf.gz}
    ## but we check if there is any variants left in the VCF (taht will take care of edge case found with MMRF_1073)

    if [[ $(bcftools view -H  ${VCF_HETS} | wc -l ) -eq 0 ]] ;
    then
        echo -e "\nNo Variant Left after removing Homozygous ; ending Phasing section now"
        echo -e "renaming input file to match expected outfile from phasing section\n"
        mycmd="cp ${VCF_ORIGINAL_INPUT} ${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf.gz}"
        echo ${mycmd} 1>&2 ;
        eval ${mycmd}
        check_ev $? "cp command"
        echo "${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf.gz}  ; Converting vcf.gz to vcf ..."  1>&2
        bcftools view -O v --threads 2 -o ${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf} ${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf.gz}
        check_ev $? "bcftools view convert vcf.gz too vcf" 1>&2
        echo "${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf.gz}"
        exit 0
    else
        VCF_FOR_PHASER=${VCF_HETS}
        echo -e "vcf with hets only and subsConsPOS is: ${VCF_FOR_PHASER}"
    fi

fi

# # Despite the check above for the absence of HETS, looks like it only check if there is no variant at all instead of just
# # checking if there is no HETS at all. So if we still have indels only it would keep going and phASER will FAIL since there
# # is no HETS to use to perform phasing; if this is the case, phASER returns an exit value of 1
echo "checking if we have any HETS in the VCF using bcftools ... " 1>&2
C=$(bcftools filter -i 'TYPE=="snp"' "${VCF_FOR_PHASER}" | bcftools view -H | wc -l)
if [[ ${C} -eq 0 ]]
then
  echo -e "No Hets in VCF, so No Phasing to perform; Skipping phASER" 1>&2
  echo "Making the expected VCF filename as if phASER had run:"
  echo "zcat \"${VCF_ORIGINAL_INPUT}\"  \"${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf}\" " 1>&2
  zcat "${VCF_ORIGINAL_INPUT}" > "${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf}"
  echo "${VCF_ORIGINAL_INPUT/.vcf.gz/.blocs.vcf}" 2>&1
  echo "Exiting $0 script without having phASER ran. ev = 0" 1>&2
  exit 0
fi
## https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
## https://github.com/secastel/phaser/tree/master/phaser
##@@@@@@@@@@@@@
### WARNING: Make sure the phaser.py file is in your path , is executable AND has the correct shebang (#!/usr/bin/env python2 )
##@@@@@@@@@@@@@
echo -e "Phasing the Consecutive positions using phASER ..."
echo -e "RUNNING phASER ... "
echo -e "Command:"
#module load python/2.7.13
VCF_IN=${VCF_FOR_PHASER}
VCF_OUT_PHASER=phased_${VCF_IN}
TEMP_DIR=$(dirname ${VCF_IN} )/temp_phASER
mkdir -p ${TEMP_DIR}

## HARDCODED value for options here ; NOTE: 8 CPUS is a good trade-off between I/O and CPUS; Too many cpus leads to too many samtools view which increases I/O and waiting time
mycmd="${PHASER_EXE} --bam ${TBAM} --sample ${SNAME_T}  --vcf ${VCF_IN}  --o ${VCF_OUT_PHASER}  --mapq 20 --baseq 20 --paired_end 1 --gw_phase_vcf 0  --remove_dups 1 --write_vcf 1  --gw_af_field AR --gw_phase_vcf 2 --temp_dir ${TEMP_DIR} --max_block_size 10 --threads ${CPUS} --id_separator '%' --unique_ids 1 "
echo -e "${mycmd}"
eval ${mycmd} | tee log_for_phASER_$(basename ${VCF_OUT_PHASER} ).log
check_ev $? "phASER using ${PHASER_EXE} with python-2.7.13"
echo -e "PHASER output file: ${VCF_OUT_PHASER}"
#module unload python/2.7.13



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
echo "python3  ${SCRIPT_PYTHON_RECOMPOSE_VARIANTS} --tumor-col 11 --normal-col 10 -i ${VCF_IN} -o ${VCF_OUT}"
python3  ${SCRIPT_PYTHON_RECOMPOSE_VARIANTS} --tumor-col 11 --normal-col 10 -i ${VCF_IN} -o ${VCF_OUT}
check_ev $? "python3 strelka2.recompose_phased_variants.py"


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
mycmd="bcftools concat -a -O v ${VCF_INDELS_ONLY} ${VCF_IN_UNPHASED} ${VCF_IN_PHASED} | bcftools sort -O v -T ${PWD} --max-mem 2G  - > ${VCF_OUT}"
echo ${mycmd}
eval ${mycmd}
check_ev $? "bcftools concat ${VCF_OUT}"

## WE OUTPUT An UNcompressed VCF; if we nee to go back to block-compressedwe will need to uncomment the two lines below
#bcftools index --tbi --threads 2 ${VCF_OUT}
#check_ev $? "bcftools index ${VCF_OUT}"

touch Completed_Recomposition_Strelka2_for_VCF_$(basename ${VCF_ORIGINAL_INPUT}).flag

echo "${VCF_OUT}" 2>&1

exit 0