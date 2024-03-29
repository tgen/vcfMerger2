#!/usr/bin/env python3

# ## vcfMerger2
# ##
# ## MIT License
# ##
# ## Copyright (c) 2018 Translational Genomics Research Institute
# ##
# ## Permission is hereby granted, free of charge, to any person obtaining a copy
# ## of this software and associated documentation files (the "Software"), to deal
# ## in the Software without restriction, including without limitation the rights
# ## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# ## copies of the Software, and to permit persons to whom the Software is
# ## furnished to do so, subject to the following conditions:
# ##
# ## The above copyright notice and this permission notice shall be included in all
# ## copies or substantial portions of the Software.
# ##
# ## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# ## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# ## FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
# ## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# ## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# ## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# ## SOFTWARE.
# ##
# ## Major Contributors: Christophe Legendre
# ## Minor Contributors:

from sys import exit
from sys import argv
from os import path
import getopt
from cyvcf2 import VCF
import logging as log
import warnings
from collections import defaultdict


class Genotype(object):
    __slots__ = ('alleles', 'phased')
    
    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]
    
    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123456789."[a] for a in self.alleles)
    
    __repr__ = __str__


def usage(scriptname):
    print("\nUSAGE: \npython3 " + scriptname + '  --help')
    print("python3 " + scriptname + " -i octopus.somatic.snvs.pass.vcf  --normal_column 10 --tumor_column 11 -o "
                                    "updated_vcf.vcf\n")
    print("#" * 40 + "\nWARNING WARNING: This script is to be used only with somatic snvs vcf having NORMAL sample in column 10 and TUMOR sample in column 11; "
                     "if not like this, update your vcf file to these specs; and the vcf has to be decomposed as well.\n" + "#" * 40)


def parseArgs(scriptname, argv):
    new_vcf_name, column_tumor, column_normal = None, None, None
    
    warnings.simplefilter(action='ignore', category=FutureWarning)
    FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
    log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)
    opts = None
    column_tumor, column_normal, new_vcf_name = None, None, None
    
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["vcf=", "tumor_column=", "normal_column=", "debug", "outfilename=", "help"])
        log.info(opts)
        print("%" * 20 + "\n" + str(opts) + "\n" + "%" * 20)
    except getopt.GetoptError:
        usage(scriptname)
        exit(2)
    for opt, arg in opts:
        if opt == '-h' or opt == "--help":
            usage(scriptname)
            exit()
        elif opt in ("", "--tumor_column"):
            column_tumor = arg
        elif opt in ("", "--normal_column"):
            column_normal = arg
        elif opt in ("", "--debug"):
            log.info("DEBUG MODE Enabled")
            log.getLogger().setLevel(log.DEBUG)
        elif opt in ("-o", "--outfilename"):
            new_vcf_name = arg.strip()
        elif opt in ("-i", "--vcf"):
            fvcf = arg
            if not path.exists(fvcf):
                exit("ERROR: FNF --> " + fvcf)
        else:
            exit("Unknown Option --> " + opt)
    
    if column_tumor is None or column_normal is None:
        exit(
            "Please Provide column number for tumor and normal Samples; should be 10 and 11  - or -  11 and 10 respectively;\n"
            "options are: --tumor_column and --normal_column;\nAborting. ")
    if column_normal == column_tumor:
        exit("ERROR: number for the columns Tumor and Normal MUST be different")
    
    return fvcf, column_tumor, column_normal, new_vcf_name


def collapse_variants(LOV):
    """
    After going through several tests, it looks like the variants in the LOV object can be considered as
    a block-substitution variant, and therefore need to be collapsed into one line;
    To keep track of which variants record were collapsed we add the INFO field  'BLOCSUBSFROM'  which gathers
    a comma-separated list of CHROM_POS  variant record information found in the original octopus output file.
    The First record is taken as a template and the ALT is going to be replaced with the concatenated string
    of all the ALTs in the LOV records.
    :param LOV:  List Variant Object (a.k.a. variant record in vcf file
    :return: variant as a list of String ready to be written into the output vcf
    """
    newREF = ""
    newALT = ""
    BLOCSUBSFROM = ""
    delim = ","
    for rec in LOV:
        newREF = newREF + rec.REF
        newALT = newALT + rec.ALT[0]
        locus = str(rec.CHROM) + "_" + str(rec.POS)
        BLOCSUBSFROM = delim.join([BLOCSUBSFROM, locus]) if BLOCSUBSFROM != "" else locus
    newV = str(LOV[0]).split("\t")
    newV[4 - 1] = newREF
    newV[5 - 1] = newALT
    newV[8 - 1] = ';'.join([newV[8 - 1], "BLOCSUBSFROM=" + BLOCSUBSFROM])
    return str('\t'.join(newV))


def processing_variants_as_block_substitution(LOV, wo):
    """
    We consider that ALL the variants within the List_Of_Variant Object in argument belong to the same BLOCK
    We need to "merge" the records to create only one
    cyvcf2 does NOT allow to modify the ALT (attribute 'ALT' of 'cyvcf2.cyvcf2.Variant' objects is not writable. ## At least at the time of this implementation)
    The only thing we could think of so far with our basic knowledge of Python, is to convert the Variant object
    into a list of Strings (variant record is tab-separated) and modify the column 5 (index 4 in python) with
    the value of the ALT alleles from the variant in the list_Of_Variant_Object argument.
    The LAST test we need to perform before considering the variants as a BLOCK is to look at the AR
    If their AR is too FAR apart, should we consider it as a block?.
    :param LOV: List Of Variant object (a.k.a. variant record in vcf file)
    :param wo: writer object to write processed record into output vcf file
    """
    # we will need to check if:
    # 1) we only deal with SNV in all ALTs or DEL in all ALTs or INS in all ALTs
    # 2) need to check if AR within x% from each other calls otherwise might mean it is two UN-related events
    # 3) should we check for a minimum DEPTH to consider worth combining records or accurate and specs?
    # 4) create a new updated record using the first variant in the block
    #      -- update ALT field
    #      -- update DP field or ADD a field in info field that might help keeping track of data
    
    tests = []
    for rec in LOV:  # HARDCODED Here Below for TESTS
        if 'MAP_HF' in rec.FORMAT and 'DP' in rec.FORMAT:
            log.debug("v.format('MAP_HF')[1] == {}".format(rec.format('MAP_HF')[1]))
            # we hardcoded index 1 in max(rec.format('MAP_HF')[1]), because the 1 is supposed to represent the tumor sample (column 11 in vcf)
            if int(rec.format('DP')[1]) >= 10 \
                    and float(max(rec.format('MAP_HF')[1])) >= 0.20 \
                    and not rec.is_indel:
                tests.append(True)
            else:
                tests.append(False)
        else:
            log.debug("MAF_HF or DP filed does not exist in current record ...")
            tests.append(False)
    
    if all(tests):
        log.info('\n'.join(["Potential BlockSub Found:"] + [str(rec).strip() for rec in LOV]))
        wo.write(collapse_variants(LOV))
    else:
        log.debug("BLOCK FAILED DP_TUMOR>=10 and/or MAP_HF_TUMOR >= 0.05 ")
        for rec in LOV:
            wo.write(str(rec))


def check_for_block_substitution(vcf, column_tumor, wo):
    log.info("looping over records to capture and concatenate Block Substitutions Variants ...")
    
    idxT = 0 if int(column_tumor) == 10 else 1
    
    dico_PS = defaultdict(list)
    k = -1  # init the key value to -1 for the first next if statement
    
    for v in vcf:  # v for variant which represents one "variant record" ; !!! WARNING: We can consume the iterator only once. Once that loop is done, the VCF object will be
        # empty and we cannot loop over vcf object anymore.
        # we gather the variant with the same PhaseSet Value
        log.debug("v.format('PS') : {}".format(v.format('PS')))
        log.debug("v.format('PS')[idxT] : {}".format(v.format('PS')[idxT]))
        log.debug("len(v.format('PS')[idxT]) : {}".format(len([v.format('PS')[idxT]])))
        if len([v.format('PS')[idxT]]) != 1:
            log.error("Found more than 1 PS value for the Tumor; The script does not handle that case; Check your input or/and Submit an issue to the vcfMerger2 github; Aborting.")
            log.error("PS value Found: {}".format(v.format('PS')[idxT]))
            exit(2)
        try:
            ps_value = int(v.format('PS')[idxT])
        except TypeError as te:
            log.error("Could not convert the PS value to an intege. PS value found is: {}".format(str(ps_value)))
            log.error(te)
            exit(2)
        new_k = '_'.join([v.CHROM, str(ps_value)])

        if k != -1 and k != new_k and len(dico_PS[k]) == 1:
            # dico has only one item, so we write it to output vcf right away and re-init dico for olk key
            for rec in dico_PS[k]:
                wo.write(str(rec))
            k = new_k
            dico_PS[k] = []
        elif k != -1 and k != new_k and len(dico_PS[k]) > 1:
            # dico has only one item, so we write it to output vcf right away and re-init dico for olk key
            # this mean we already processed value for the previous PS value, and we reached a new PS value, so the key k is different
            # we also check if the dico has more than one records if not we print record
            processing_variants_as_block_substitution(dico_PS[k], wo)
            k = new_k
            dico_PS[k] = []  # we reinit the dico for the current k value
        
        if len(dico_PS[k]) == 0:  # if dico[k] empty we add the current variant
            dico_PS[k].append(v)
        
        elif len(dico_PS[k]) == 1:  # this mean we already had one variant with same PS
            # we compare the variant POSITION(implies that the VCF was sorted by variant records)
            # if position is shifted by one, this means the variant calls are right next to each other, and we keep both in the dictionary
            # if not, we ony keep the current variant in the dico for that PS and continue to next variant in the loop
            
            # let's get 'GT' for the current variant and the variant in position [-1] in dico
            # the variant of the Tumor should only be considered here
            # same genotype means similar to same PHASE
            lgenotype_current = [str(Genotype(li)) for li in v.genotypes][1]
            lgenotype_previous = [str(Genotype(li)) for li in dico_PS[k][-1].genotypes][1]
            # we then compare the current locus position (v.POS) with the last locus added to the dictionary
            # if consecutive loci and genotype is identical, good chances we deal with variants that can
            # constitute a block-substitution
            if int(v.POS) == int(dico_PS[k][-1].POS) + 1 and lgenotype_current == lgenotype_previous:
                dico_PS[k].append(v)  # we gather the variant with the same PhaseSet Value, Same genotype in Tumor and consecutive POS
            else:
                # the variants in the dico either do not have the same genotype or were not consecutive
                # so, we have to process what we have, and as it is only one variant in the dico,
                # we write it here, directly to the output file;
                for rec in dico_PS[k]:
                    wo.write(str(rec))
                dico_PS[k] = [v]  # we re-init the value to the current variant as the previous variant is definitely not right next to the current one even though in the same
        # PhaseSet
        
        elif len(dico_PS[k]) > 1:
            # this mean we already have at least two variants right next to each other in the dico
            # Either the 3rd or nth variant is right after the position of the latest variant in the list
            # or if not this implies two statements:
            # 1) the current variant does not belong to the same block-substitution as the last variant already in the dictionary
            # 2) we have to process the variants currently in the dictionary has "potential" block-substitution and the
            # current variant will be added to the dico_PS[k] replacing the old variants
            # 2a) if Genotype of each consecutive position is the same; it should be considered as a block
            # 2b) if Genotype of each consecutive position is NOT the same; we have to Process the variants that
            #      are already in the dictionary has belonging to the same Block-Substitution Variant
            lgenotype_current = [str(Genotype(li)) for li in v.genotypes][1]
            lgenotype_previous = [str(Genotype(li)) for li in dico_PS[k][-1].genotypes][1]
            if int(v.POS) == int(dico_PS[k][-1].POS) + 1 and lgenotype_current == lgenotype_previous:
                dico_PS[k].append(v)  # we gather the variant with the same PhaseSet Value
            else:
                processing_variants_as_block_substitution(dico_PS[k], wo)  # dicoPS[k] is a list of Variants
                dico_PS[k] = [v]
    
    # we need to write the data from the last phase set captured
    if len(dico_PS[k]) == 1:
        for rec in dico_PS[k]:
            wo.write(str(rec))
    else:
        processing_variants_as_block_substitution(dico_PS[k], wo)


def check_if_PS_in_FORMAT_field(vcf_cyobj, input_vcf_path, p_new_vcf_name):
    iterVCF = iter(vcf_cyobj)
    v1 = next(iterVCF)
    log.info("Checking PS flag presence in FORMAT ...")
    if 'PS' not in v1.FORMAT:
        log.warning(
            "PS tag s not present in the FORMAT field of OCTOPUS; We assume that the VCF has already been modified from its original copy.")
        if ':'.join(v1.FORMAT) == "GT:DP:AR:AD":
            log.warning("FORMAT field is equivalent to 'GT:DP:AR:AD'")
            log.warning(
                "We assume the vcf has already been prepared for vcfMerger2 and therefore just copy the vcf by assigning the decomposed expected filename output")
            from shutil import copyfile
            copyfile(input_vcf_path, p_new_vcf_name)
            exit()
        else:
            log.error("PS flag NOT found in OFRMAT field; Aborting VCF preparation.")
            exit("PS flag Absent")
    else:
        log.warning("PS flag Found")


# @#########
# @ MAIN  ##
# @#########
if __name__ == "__main__":
    
    vcf_path, column_tumor, column_normal, new_vcf_name = parseArgs(argv[0], argv[1:])
    
    vcf = VCF(vcf_path)
    tot_number_samples = len(vcf.samples)
    if tot_number_samples > 2 or tot_number_samples < 2:
        exit("ERROR: Number of Samples is different from 2; Expected 2 samples only NORMAL and TUMOR (in that order in the VCF) ; Aborting Here.")
    
    if new_vcf_name is None:
        new_vcf = '.'.join([str(vcf_path), "blocsubs.vcf"])
    else:
        new_vcf = new_vcf_name
    
    # checking if PS flag is still present in the VCF genotype fields
    check_if_PS_in_FORMAT_field(vcf, vcf_path, new_vcf_name)
    
    vcf = VCF(vcf_path)
    # Adding Fields to INFO field
    vcf.add_info_to_header({'ID': 'BLOCSUBSFROM', 'Description': 'List of Original Position in OctopusVcf That were used to make current block substitution variant',
                            'Type': 'String', 'Number': '.'})
    
    # create a new vcf Writer using the input vcf as a template.
    w = open(new_vcf, 'w')
    # writing the updated header immediately to the new output vcf file
    w.write(str(vcf.raw_header))
    
    log.info("looping over records ...")
    check_for_block_substitution(vcf, column_tumor, w)
    
    w.close()
    vcf.close()
    
    exit()
