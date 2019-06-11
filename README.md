<p align="center">
<img src="images/vcfMerger2.logo.png"/>
</p>

# vcfMerger2 
vcfMerger2 is a tool merger for 2 to N somatic variants vcf files. 


#### What is vcfMerger2?
This tool allows merging from 2 to N somatic vcfs obtained after running different variant callers using the same bam file(s). 
**vcfMerger2** will generate one vcf. 

#### How vcfMerger2 operates?
vcfMerger2 can use two types of inputs vcfs;  
Either the input vcfs are `RAW vcfs` (case_1) - OR -  the vcfs are `vcfMerger2-prepped-ready` vcfs (case_2)

##### Case_1
See [Example_1](https://github.com/tgen/vcfMerger2/wiki/Examples#Example_1)
vcfMerger2 inputs are raw vcfs from different variant callers. (see list of compatible variant caller vcfs in wiki page "Variant_Callers_compatible"). 
In this case, before merging, the vcfs must be prepared to get them to vcfMerger2 specifications;
this part is transparent to the user but can be achieved only with certain variant callers' vcfs 

##### Case_2
See [Example_2](https://github.com/tgen/vcfMerger2/wiki/Examples#Example_2)
vcfMerger2 inputs are ready-to-be-merged vcfs from different variant callers. Before running vcfMerger2, these vcfs have been prepared 
to be vcfMerger2-ready for merging. User modifies the vcfs to make them compatible or uses the vcfMerger2 utilities 
script provided for some variant caller (see "Variant_Callers_compatible" or "Examples" wiki pages) 

#### Somatic or Germline? 
VcfMerger2 has been developped to merge SOMATIC vcfs into one vcf. The current release can merge from 2 to N vcfs.  
To see the list  variant callers compatible with vcfMerger2 go to "Variant_Callers_compatible" wiki page

Merging GERMLINE vcfs will be available in the next version; some features are already there but are still under development and should be consider as in Alpha version. 


# How to Install 
See wiki [Installation](https://github.com/tgen/vcfMerger2/wiki/Installation#Installation) page

# How to run vcfMerger2?
The core utility and main purpose of that tool is to MERGE VCFs information into one VCF without losing any data from any of the input vcfs files
See [QuickStart](https://github.com/tgen/vcfMerger2/wiki/QuickStart#QuickStart) wiki page.



### License
vcfMerger2 is under MIT licence.
