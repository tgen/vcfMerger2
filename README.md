# vcfMerger2 
Description : Dynamic vcfMerger for 2 to N somatic variants vcf files. 


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
You can install vcfMerger2 in two ways:
1. download the git repository as a compressed archived, uncompressed and untar the archive, and update your `PATH` with the command  
`export PATH=/install_directory/vcfMerger2/bin/:${PATH}` 

2. clone the repository and add the bin directory to your `PATH` as described above.

After installing, you may delete the test_data directory unless you want to run the demo script. 


# How to run vcfMerger2?
#### The All-in-One way [ prep + merge ]

The easiest way to use vcfMerger2 is to call on the `vcfMerge2.py` script in `bin` directory (normally in your `PATH`) 

vcfMerger2 can merge from 2 to N vcfMerger2-upto-specs somatic variants vcfs. 
 
Example:  
``` 
python vcfMerger2.py \   
--toolnames "strelka2|mutect2|lancet|octopus" 
--vcfs "./raw_tool_vcfs/strelka2.raw.vcf|./raw_tool_vcfs/mutect2.raw.vcf|./raw_tool_vcfs/lancet.raw.vcf|./raw_tool_vcfs/octopus.raw.vcf" 
-a "SLK|MUT|LAN|OCT"  
-g hs37d5.fa \ 
--prep-outfilenames "strelka2.prepped.vcf|mutect2.prepped.vcf|lancet.prepped.vcf|octopus.prepped.vcf" 
--normal-sname NORMAL  
--tumor-sname TUMOR -o merged.vcf 
--contigs-file "||contigs/contigs.txt|contigs/contigs.txt"  
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
                        
  `--normal-sname NORMAL_SNAME` 
                        expected name of normal sample in vcf file
                        
  `--tumor-sname TUMOR_SNAME` 
                        expected name of tumor sample in vcf file
                        
  `--prep-outfilenames PREP_OUTFILENAMES` 
                        delim-separated names to the tool-specific prepared
                        vcf files
                        
  `-o MERGED_VCF_OUTFILENAME, --merged-vcf-outfilename MERGED_VCF_OUTFILENAME` 
                        outfilename for the merge vcf (can be relative or full
                        path)

##### optional arguments: 

  `-c PRECEDENCE, --precedence PRECEDENCE`
                        sorted delim-separated list of the toolnames as listed
                        in --toolnames; This list stipulates an order of precedence for the tools different from the 
                        default order given by the `--toolnames` list 
   
   ---                      

**NOTE**: What does `PRECEDENCE` mean here?  A vcf contains information in INFO and FORMAT columns. Unfortunately, redundant information exist from one tools to another in a vcf.
For instance, the AR field may exist in ALL the given vcf in the FORMAT columns, but the values may vary from one tool to another. Unfortunately, only one value can be kept in the AR field within the merged vcf. 
So which one would the user preferably keep in the merged vcf? What tool does have our "liking" the most. This is where the precedence is used. It gives an order of preference for the tool when the vairant is called 
by more than one tool. This Precedence is subjective to the user.   

   ---                       

  `--bams BAMS`          List of bams necessary for capturing contigs if not
                        present in input vcf; otherwise put empty_string as
                        value for each tool
                        
  `--contigs-file-for-vcf-header CONTIGS_FILE_FOR_VCF_HEADER`
                        List of contigs necessary for capturing adding them to
                        tool vcf header if needed; otherwise put empty_string
                        as value for each tool ;do not provide if bam file is
                        given instead
                        
  `-a TOOLACRONYMS, --toolacronyms TOOLACRONYMS`
                        List of Acronyms for toolnames to be used as PREFIXES
                        in INFO field ; same DELIM as --vcfs
                        
  `--delim DELIM`         delimiter which will be use to create the arguments
                        value for the vcfMerger2.0 tool ; default is "|"
                        (a.k.a pipe character)
                        
  `--lossy`               This will create a lossy merged vcf by only keeping
                        the infromation from the tool with precedence
                        
  `--skip-prep-vcfs`      skip the step for preparing vcfs up to specs and only
                        run the merge step; implies all prep-vcfs are ready
                        already ; same options and inputs required as if prep
                        step was run
                        
  `--skip-merge`          enabling this flag prevents doing the merging step
                        [useful if only the prep step needs to be done ]
                        
  `-n, --dry-run`         perform a dry run to just see the commands lines that
                        will be run behind the scenes
                         

--- 


## vcfMerger2 requirements

All the following tools **must** be in your `PATH` before running `vcfMerger2` scripts 

- linux system 
- grep, awk 
- python 3.6.0 or up
    - cyvcf2
    - collections
    - argparse
    - logging
    - natsort
    - sys, os, time, re 
    
- samtools 1.7 or up
- bcftools 1.7 or up
- vt v0.57721 or up (tests done with v0.57721)


## Running demo
Once all the requiremetns are available, it is advised to run the **demo** script.  
You can run the script from anywwhere, and the results will be outputted into a **demo_out** folder created by the script itself into the **test_data** directory.

----
----
# For advanced users
vcfMerger2 is composed of two main steps:
1. Preparing vcf to specs
2. Merging vcfs

Users can run these two steps independently if needed, even though this can be achieved by using the `--skip-prep-vcfs` or `--skip-merge` optional arguments from the `all-in-one` way. 

#### How to _**only**_ prepare the vcf files using vcfMerger2?
Users can prepare tool-specific vcf by only running the sub-step of vcfMerger2 called `prep_vcf` independently of the `all-in-one` way described above. 
The script called `prep_vcf.sh` allows you to specifically bring up to vcfMerger2 specs vcfs for the supported tools *(see list of supported tools below)*.
This script can be called directly and run independently for each vcf and tool. 

Example_1:  
`prep_vcf.sh --toolname lancet --vcf ./test_data/raw_tool_vcfs/lancet.raw.vcf -d ./test_data -g  ./ref_genome/grch37.22.fa` 

Example_2:  
`prep_vcf.sh --toolname lancet --normal-sname NORMAL --tumor-sname TUMOR --vcf raw_tool_vcfs/lancet.raw.vcf -g ref_genome/grch37.22.fa -o lancet.prepped.vcf --contigs-file ./contigs/contigs/txt`
 

_**HINT**_: if you already have the command line from `all-in-one` way style, you can just add the option `--skip-merge` to that command and only the `prep` step will be performed for all the given vcfs instead of running one by one the `prep_vcf.sh` script


#### How to **_only_** merge vcfs? 

The script called `vcfMerger.py` allows you to merge the data of several already ready-prepped-vcf files .
This script can be called and run independently.  
vcfMerger2 can merge from 2 to N vcfMerger2-upto-specs somatic vcfs;  
Bringing the vcfs up to vcfMerger2 specs is **mandatory** before running `vcfMerger.py` script ; this script can be found in the vcfMerger sub-directory 

Example_3 (simplest way of running vcfMerger step only):  
`vcfMerger.py --toolnames "strelka2|mutect2|lancet|octopus" --vcfs "strelka2.prepped.vcf|mutect2.prepped.vcf|lancet.prepped.vcf|octopus.prepped.vcf"`  
Example_4 (adding the name of the output file; otherwise default is `input_vcf_filename.merge.vcf`):  
`vcfMerger.py --toolnames "strelka2|mutect2|lancet|octopus" --vcfs "strelka2.prepped.vcf|mutect2.prepped.vcf|lancet.prepped.vcf|octopus.prepped.vcf" -o merged.vcf `
Example_3 (adding acronyms to reduce file size):  
`vcfMerger.py --toolnames "strelka2|mutect2|lancet|octopus" --vcfs "strelka2.prepped.vcf|mutect2.prepped.vcf|lancet.prepped.vcf|octopus.prepped.vcf" -o merged.vcf -a "SLK|MUT|LAN|OCT" `   

_**HINT**_: You also can run the same command listed in `all-in-one` way by just adding the option `--skip-prep-vcfs`, and only vcfMerger step will be performed.

#### what tool-specific raw vcfs are currently supported?
To bring each vcf up to vcfMerger2 specs, a tool-specific script that modifies and update tool-specific vcf is provided. (see `prep_vcfs` directory) 
To date, 4 tools are supported.
- strelka2
- mutect2
- lancet
- octopus

If the variant caller users use is not listed above, users can make their own `prep` script for their tool then create a pull-request on github, or open an issue with keyword `improvement` as first line in bloc text. 

### Licence
vcfMerger2 is under MIT licence.