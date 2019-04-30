<p align="center">
<img src="vcfMerger2.logo.png"/>
</p>
DEVEL Branch

# vcfMerger2 
vcfMerger2 is a tool merger for 2 to N somatic variants vcf files. 


#### What is vcfMerger2?
This tool allows merging from 2 to N somatic vcfs obtained after running different variant callers using the same bam file(s). 
**vcfMerger2** will generates one vcf.

#### How vcfMerger2 operates?
vcfMerger2 inputs are raw vcfs outputted from variant callers.

Before merging the vcfs' data, each vcf file MUST be brought up to vcfMerger2 specifications. These specs follow vcf specifications but has one specific requirement.
The `FORMAT` fields must contains expected flags.

vcfMerger2 team has provided scripts to bring vcfs up to vcfMerger2 specs. So far, there are only few supported variant caller tools (see section: **_what vcfs and tools are supported?_** below)   

vcfMerger2-upto-spces vcfs are merged. The merging is lossless by default (bigger vcfs). Lossless means the merged vcf 
will contain ALL the information from **ALL** the inputs vcfs


# How to Install 
You can install vcfMerger2 in three ways:
1. download the git repository as a compressed archived, uncompressed and untar the archive, then optionally update your `PATH` with the command  
`export PATH=/install_directory/vcfMerger2/bin/:${PATH}` 

2. clone the repository and add the bin directory to your `PATH` as described above. 
3. clone a specific release version of vcfMerger2, such as for instance for the versio 0.5.0-beta, the following command:
`git clone git@github.com:tgen/vcfMerger2.git --release  0.5.0-beta` 
or
`git clone https://github.com/tgen/vcfMerger2.git --release  0.5.0-beta`


After installing, you may delete the test_data directory unless you want to run the demo script. 


# How to run vcfMerger2?
#### The All-in-One way [ prep + merge ]

The easiest way to use vcfMerger2 is to call on the `vcfMerge2.py` script in `bin` directory (normally in your `PATH`) 

vcfMerger2 can merge from 2 to N vcfMerger2-upto-specs somatic variants VCFs files. 
 
###### Example_1:  
``` 
python vcfMerger2.py    
-g hs37d5.fa  
--toolnames "strelka2|mutect2|lancet|octopus" 
--vcfs "./raw_tool_vcfs/strelka2.raw.vcf|./raw_tool_vcfs/mutect2.raw.vcf|./raw_tool_vcfs/lancet.raw.vcf|./raw_tool_vcfs/octopus.raw.vcf" 
--prep-outfilenames "strelka2.prepped.vcf|mutect2.prepped.vcf|lancet.prepped.vcf|octopus.prepped.vcf" 
--normal-sname NORMAL  
--tumor-sname TUMOR 
-o merged.vcf 
--contigs-file "||contigs/contigs.txt|contigs/contigs.txt"
-a "SLK|MUT|LAN|OCT"
--precedence "MUT|LAN|SLK|OCT"    
```
 
##  
_**Note_1**: empty string information for the value of `--contigs-file` is necessary to match the number of values given to the `--vcfs` option even though the information for the 1st and 2nd tools is absent in this case._  
_**Note_2**: run command above from top folder vcfMerger2._  
_**Note_3**: hs37d5.fa filename for the reference genome should be replaced with your own version of the GRCh37._  
_**Note_4**: The given command line example uses 4 input vcfs. [ you may also test it using only 2, or 3 input vcfs ; for instance, run the same command by removing lancet and octopus tool from the lists ]_  

[ run `python vcfMerger2.py --help` for more options ]

-----------

## Options to vcfMerger2.py 

##### required arguments:  
  `--vcfs VCFS`           List of vcfs file delimited by DELIM character;
                        default DELIM is "|" (piped character) ; if needed re-
                        assign DELIM using --delim option
                         
  `--toolnames TOOLNAMES` 
                        List of vcfs file delimited by DELIM character;
                        default DELIM is pipe unless --delim option is used
                        using --delim option
                        
  `-g REFGENOME, --refgenome REFGENOME` 
                        reference genome `.fa` and its `.fai` index are used with bcftools norm ; must match
                        reference used for alignment

                     
  `--prep-outfilenames PREP_OUTFILENAMES` 
                        delim-separated names to the tool-specific prepared
                        vcf files
                        
  `-o MERGED_VCF_OUTFILENAME, --merged-vcf-outfilename MERGED_VCF_OUTFILENAME` 
                        outfilename for the merge vcf (can be relative or full
                        path)

##### additional required arguments if merging SOMATIC Vcfs:                         
                      
  `--tumor-sname TUMOR_SNAME` 
                        expected name of tumor sample in vcf file (already present in VCF; mostly used to make sure of the right column)
  
  `--normal-sname NORMAL_SNAME`
                        expected name of normal sample in vcf file (already present in VCF; mostly used to make sure of the right column)
   

##### optional arguments: 

  `--germline` 
                        consider the inputs VCFs as GERMLINE only VCFs (and not somatic)
                        This implies using the right toolnames associated to germline calls. as of 2019, we implemented
                        preprocessing for three (3) germline variants callers: HaplotypeCaller, Freebayes and and Deepvariant.
                        (see exact name using --list-tools )

  `--germline-snames`
                        ordered list of sample name(s) that should be present in the VCF beyond column 9.
                        As of 2019-04-04, ONLY one sample per VCF is taken into account; So this option should be a list of ONE name only
                     


  `-c PRECEDENCE, --precedence PRECEDENCE`
                        sorted delim-separated list of the toolnames as listed
                        in --toolnames; This list stipulates an order of precedence for the tools different from the 
                        default order given by the `--toolnames` list 
   
   ---                      

**NOTE**: What does `PRECEDENCE` mean?  A vcf contains information in INFO and FORMAT columns. Unfortunately, redundant information exist from one tools to another in a vcf.
For instance, the AR field may exist in ALL the given vcf in the FORMAT columns, but the values may vary from one tool to another. Unfortunately, only one value can be kept in the AR field within the merged vcf. 
So which one would the user preferably keep in the merged vcf? What tool does have our "liking" the most. This is where the precedence is used. It gives an order of preference for the tool when the vairant is called 
by more than one tool. This Precedence is subjective to the user.   

   ---                       

  `--bams BAMS`          List of bams necessary for capturing contigs if not
                        present in input vcf; otherwise put empty_string as
                        value for each tool [maybe deprecated in future as ALL the vcfs MUST contain the exact same contigs in their headers (if not check you input vcfs and update them)]
                        
  `--contigs-file-for-vcf-header CONTIGS_FILE_FOR_VCF_HEADER`
                        List of contigs necessary for capturing adding them to
                        tool vcf header if needed; otherwise put empty_string
                        as value for each tool ;do not provide if bam file is
                        given instead [maybe deprecated in future as ALL the vcfs MUST contain the exact same contigs in their headers (if not check you input vcfs and update them)]
                        
  `-a TOOLACRONYMS, --toolacronyms TOOLACRONYMS`
                        List of Acronyms for toolnames to be used as PREFIXES
                        in INFO field ; same DELIM as --vcfs
                        
  `--delim DELIM`         delimiter which will be use to create the arguments
                        value for the vcfMerger2.0 tool ; default is "|"
                        (a.k.a pipe character)
                        
  `--lossy`               This will create a lossy merged vcf by only keeping
                        the infromation from the tool with precedence
                        
  `--skip-prep-vcfs`      skip the step for preparing vcfs up to specs and only
                        run the merge step; implies all << given >> vcfs are 
                        already up-to-specs (uts); whether or not this option is used, user has to provide the same 
                        options and inputs required to run the prep step, so as if the prep step is run.
                        
  `--skip-merge`          enabling this flag prevents doing the merging step
                        [useful if only the prep step needs to be done ]
                        
  `--threshold-AR`        Threshold Allele Ratio (AR) to assign genotype GT value with 0/1 or 1/1 ; 
                          GT=0/1 if below threshold, GT=1/1 if equal or above threshold [ default is 0.90 ; range ]0,1] ] 
                           
  `-n, --dry-run`         perform a dry run to just see the commands lines that
                        will be run behind the scenes
                         
`--filter-by-pass`      filter the variants by PASS. This implies that the keyword `PASS` must be present in column 7 of the VCF, and not a dot.
                        This filtering step is specifically performed before preparing the vcfs to specs steps ; if you want to filter PASS after prep_step and before merging, use `--filter FILTER`; 

`--filter FILTER`   filter variants using snpSift in the backend; This filtering process is performed ONLY on up-to-specs vcf (so after the stage prep-vcf [if run] );
                    A string must be provided as if you were using snpSift (see snpSift user manual); For each tool, a string must be provided; If you have 3 tools, three string must be provided;
                    Hint: if you want to apply the exact same filtering to all the tools because they do possess the same Flags, you may provide only one string. 

`--path-jar-snpsift`    FULL PATH to the SnpSift JAR file MUST be provided if any of the filter option is used. 

 `--do-venn`    Will make a Venn diagram using the vcf files provided to the merging step. This is a simple vennDiagram.
 
 `--beds `      If the user already have the data in bed format to make the Venn Diagram, the user can provide these bed files here; The number of bed files MSUT be the same as the number of toolnames or VCFs files. 
                If bed files are not provided and `--do-venn` is enabled, the bed files are created on the fly from the vcf files;
                If `--beds` is used, `--do-venn` must also be used otherwise the venn won't be created.   


--- 


## vcfMerger2 requirements

All the following tools **must** be in your `PATH` before running `vcfMerger2` scripts 

- linux system 
- grep, awk, sed, echo, find, ...
- python 3.6.0 or up
    - cyvcf2 (tested with versions: 0.8.2 or 0.8.4 ; other versions not tested and not recommended)
    - intervene
    - collections
    - argparse
    - getopt
    - gzip
    - json
    - logging
    - Pillow
    - natsort
    - shutil
    - subprocess
    - sys, os, time, re 
    - warnings
    
- samtools 1.7 or up
- bcftools 1.7 or up
- bedtools 2.26.0 or up
- vt v0.57721 or up (tests done with v0.57721)
- snpSift 4.3t or up (snpSift.jar)


## FILTERING VCF
### PASS variant and FILTER column
vcf files must contain the keyword "PASS" in the FILTER column if you want to filter (i.e., keep) the variants the tool considers as of "PASS" or 'OK' 

### Filtering options [[ --filter_by_pass and --filter options ]]
VCF Filtering is performed here using a third party tool called SnpSift from the snpEff tool. 
The full path to "SnpSift.jar" file MUST be provided to vcfMerger2 if any of the two aforementioned options is given.
`--path_jar_snpsift` is used to provide the full path to the jar file
Two options are available to filter the variant before merging the vcf files. 
`--filter-by-pass` and `--filter` 
1. `filter-by-pass` filter the raw vcf outputted by the variant caller using the `FILTER` column and looking for the keyword `PASS` in that column; If a dot is present instead of `PASS` to signify a PASS variant, all the good variant will be skipped and no variant are going to be found in th vcf.  


## Running demo
Once all the requiremetns are available, it is advised to run the **demo** script.  
You must run the script from within the **test_data** directory [in theory, it could be run from anywhere]. The results will be outputted into a **demo_out** folder created by the script itself into the **test_data** directory.
Any previous **demo_out** directory from a previous run will be deleted.

----
----
# For advanced users
vcfMerger2 is composed of two main steps:
1. Preparing vcf to specs
2. Merging vcfs
3. Filtering any vcf

Users can run these two steps independently if needed, even though this can be achieved by using the `--skip-prep-vcfs` or `--skip-merge` optional arguments from the `all-in-one` way. 

#### How to _**only**_ prepare the vcf files using vcfMerger2?
Users can prepare tool-specific vcf by only running the sub-step of vcfMerger2 called `prep_vcf` independently of the `all-in-one` way described above. 
The script called `prep_vcf.sh` allows you to specifically bring up to vcfMerger2-specs vcfs for the supported tools *(see list of supported tools below)*.
This script can be called directly and run independently for each vcf and tool. 

###### Example_2:  
`prep_vcf.sh --toolname lancet --vcf ./test_data/raw_tool_vcfs/lancet.raw.vcf -d ./test_data -g  ./ref_genome/grch37.22.fa` 

###### Example_3:  
`prep_vcf.sh --toolname lancet --normal-sname NORMAL --tumor-sname TUMOR --vcf raw_tool_vcfs/lancet.raw.vcf -g ref_genome/grch37.22.fa -o lancet.prepped.vcf --contigs-file ./contigs/contigs/txt`
 

_**HINT**_: if you already have the command line from `all-in-one` way style, you just add the option `--skip-merge` to that command and only the `prep` step will be performed for all the given vcfs instead of running one by one the `prep_vcf.sh` script


#### How to **_only_** merge vcfs? 

The script called `vcfMerger.py` allows you to merge the data of several already ready-prepped-vcf files .
This script can be called and run independently.  
vcfMerger2 can merge from 2 to N vcfMerger2-upto-specs somatic vcfs;  
Bringing the vcfs up to vcfMerger2 specs is **mandatory** before running `vcfMerger.py` script ; this script can be found in the vcfMerger sub-directory 

###### Example_4 (simplest way of running vcfMerger step only):  
`vcfMerger.py --toolnames "strelka2|mutect2|lancet|octopus" --vcfs "strelka2.prepped.vcf|mutect2.prepped.vcf|lancet.prepped.vcf|octopus.prepped.vcf"`  
###### Example_5 (adding the name of the output file; otherwise default is `input_vcf_filename.merge.vcf`):  
`vcfMerger.py --toolnames "strelka2|mutect2|lancet|octopus" --vcfs "strelka2.prepped.vcf|mutect2.prepped.vcf|lancet.prepped.vcf|octopus.prepped.vcf" -o merged.vcf `
###### Example_6 (adding acronyms to reduce file size):  
`vcfMerger.py --toolnames "strelka2|mutect2|lancet|octopus" --vcfs "strelka2.prepped.vcf|mutect2.prepped.vcf|lancet.prepped.vcf|octopus.prepped.vcf" -o merged.vcf -a "SLK|MUT|LAN|OCT" `   

_**HINT**_: You also can run the same command listed in `all-in-one` way by just adding the option `--skip-prep-vcfs`, and only vcfMerger step will be performed.

###### Example_7 (Germline calls instead of somatic ; adding acronyms to reduce file size):  
`python3 /home/clegendre/qscripts/gits/vcfMerger2_devel_branch/bin/vcfMerger2.py -g ${REF_GENOME}  --toolnames "haplotypecaller|freebayes|samtools" --vcfs "testFile_HC.100000lines.vcf|testFile_FB.100000lines.vcf|testFile_ST.100000lines.vcf" --prep-outfilenames "HC_prep.vcf|FB_prep.vcf|ST_prep.vcf" --germline  --germline-snames "HAPI_0001_000001_OV_Whole_T1_TSWGS_A28333" -o "merged_germline_calls_3tools.vcf" -a "HC|FB|ST"`

###### Example_8 (using filtering options [here both filtering options are being used]):
```
python3  ${VCFMERGER2_INSTALL_DIR}/bin/vcfMerger2.py 
--toolnames "strelka2|mutect2|lancet|octopus" 
--vcfs "strelka2.somatic.snvs_indels.vcf|mutect2.somatic.snvs_indels.FiltMutCallsTool.vcf|lancet.commpressed.somatic.snvs_indels.vcf.gz|octopus.legacy.vcf" 
--normal-sname "COLO829_C2" 
--tumor-sname "COLO829_T1"  
-g ${REF_GENOME_FASTA_FILE} 
-o colo829.merged.vcf 
--path-jar-snpsift /tools/snpEff/4.3/snpEff/SnpSift.jar 
--filter-by-pass 
--prep-outfilenames "SLK.prep.vcf|MUT.prepped.vcf|LAN.prepped.vcf|OCT_prepped.vcf"  
--filter "( GEN[0].DP>=10 && GEN[0].DP>=10 ) && ( GEN[0].AR<=0.02 && GEN[1].AR >= 0.05 )" 
--do-venn`
```

#### what tool-specific raw vcfs are currently supported?
To bring each vcf up to vcfMerger2 specs, a tool-specific script that modifies and update tool-specific vcf is provided. (see `prep_vcfs` directory) 
To date, 4 somatic tools are supported.
- strelka2
- mutect2
- lancet
- octopus
- vardict

To date, 3 germline callers are supported:
- haplotype caller
- freebayes
- samtools mpileup
###### Coming soon:
- strelka2
- deepvariant


If the variant caller is not listed above, users can make their own `prep` script for their tool then create a pull-request on github, or open an issue with keyword `improvement` as first line in bloc text. 

### Licence
vcfMerger2 is under MIT licence.

## vcfMerger2 Core tool 
###### This is the **CORE** functionnality of the vcfMerger2 tool
###### when user already has vcfMerger2-merge-ready VCFs, running vcfMerger2 executable is as easy as pie
###### The example [Example_1](#Example_1:) above show the minimum inputs the user has to give.
![flowchart](/vcfMerger2.Core.Functionality.png)


## vcfMerger2 and wrapper utilities
###### for know variant caller, we provide scripts allowing preparing the tool-specific vcf to specifications for the merger tool
###### We provide a prep_vcf.sh script that allows to prep vcf independently from the main executable that is vcfMerger2.py located in the bin directory
###### Furthermore, we also provide a simplified way of filtering input vcf by PASS and to filter prepped vcfs the same way you would using snpSift directly. 
![flowchart](/vcfMerger2.flowchart.png)