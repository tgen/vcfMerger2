#!/usr/bin/env python3

# ##vcfMerger2
# ##
# ##MIT License
# ##
# ##Copyright (c) 2018 Translational Genomics Research Institute
# ##
# ##Permission is hereby granted, free of charge, to any person obtaining a copy
# ##of this software and associated documentation files (the "Software"), to deal
# ##in the Software without restriction, including without limitation the rights
# ##to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# ##copies of the Software, and to permit persons to whom the Software is
# ##furnished to do so, subject to the following conditions:
# ##
# ##The above copyright notice and this permission notice shall be included in all
# ##copies or substantial portions of the Software.
# ##
# ##THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# ##IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# ##FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
# ##AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# ##LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# ##OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# ##SOFTWARE.
# ##
# ##Major Contributors: Christophe Legendre'a0
# ##Minor Contributors:

# WHAT DOES THIS SCRIPT DO?
# In order to <<harmonize>> the GENOTYPE flags, we add, modify, remove or replace the Genotype Octopus' Flags
# We want some common following flags in the FORMAT fields for each of our vcf that need to be merged:
# GT:DP:AR:AD: they are the flags we want to be common in the FORMAT columns in the merged VCF ;
# We are not removing any flags here; we might modify some of them or adding some if missing such as AR here


from sys import exit
from sys import argv
from os import path
import getopt
from cyvcf2 import VCF, Writer
import numpy as np
import logging as log
import warnings
from math import isnan

global AR_threshold_for_GT
AR_threshold_for_GT = 0.90  # value HARDCODED but can be dynamically modify with --threshold_AR option


class Genotype(object):
    """
    # genotypes = [Genotype(li) for li in variant.genotypes ]
    # which shows: [./., ./., ./., 1/1, 0/0, 0/0, 0/0]
    """
    __slots__ = ('alleles', 'phased')
    
    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]
    
    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123456789."[a] for a in self.alleles)
    
    __repr__ = __str__


class GenotypeInv(object):
    
    def __init__(self, li):
        self.li = li
        try:
            # we added the if statements because Octopus may put only one letter in GT when dealing with chrY or chrM or any haploid calls
            # so, we decided to re-built the Genotype to be consistent here
            log.debug("LI==" + str(li))
            log.debug("LI_length==" + str(len(li)))
            
            if len(li) < 3:
                log.debug("LI[0]==" + str(li[0]))
                log.debug("length LI less than 3 ")
                if self.li[0] == "." or str(self.li[0]) == "0":
                    self.li = ['0', '|', '0']
                elif str(self.li[0]) == "1":
                    self.li = ['1', '|', '1']
                else:
                    self.li = [0, "|", li[1:]]
                self.allele1 = self.li[0]
                self.phased = bool(0) if self.li[1] == "/" else bool(1)
                self.allele2 = self.li[2]
            else:
                self.allele1 = li[0]
                self.phased = bool(0) if li[1] == "/" else bool(1)
                self.allele2 = li[2]
                if len(li) >= 5:
                    self.allele3 = li[4]
                if len(li) >= 7:
                    self.allele4 = li[6]
                self.GT = []
        except Exception as e:
            log.error(e)
            log.error("LI==" + str(li))
            exit(1)
    
    def get_gt_numpy_compatible(self):
        self.GT = []  # we need to reinit the GT list here otherwise shared by all instances.
        # Weird because we re-initiated it already in the _init_ ;
        # I am probably missing knowledge in some python features behaviour for classes.
        if self.phased:
            if len(self.li) == 1:  # this means we deal with Haploid call such as with chrY, chrM
                if self.allele1 == 0:
                    self.GT.append(0)
                    self.GT.append(0)  # we need two values for the GT to be compatible with the cyvcf2's numpy_array
                elif self.allele1 == 1:
                    self.GT.append(1)
                    self.GT.append(1)  # we need two values for the GT to be compatible with the cyvcf2's numpy_array
            else:
                if self.allele1 != ".":
                    self.GT.append((2 * int(self.allele1)) + 3)
                else:
                    self.GT.append(1)
                if self.allele2 != ".":
                    self.GT.append((2 * int(self.allele2)) + 3)
                elif len(self.li) >= 3:
                    self.GT.append(1)
                if len(self.li) >= 5 and self.allele3 != ".":
                    self.GT.append((2 * int(self.allele3)) + 3)
                else:
                    if len(self.li) >= 6:
                        self.GT.append(1)
                if len(self.li) >= 7 and self.allele4 != ".":
                    self.GT.append((2 * int(self.allele4)) + 3)
        else:  # Octopus is ALL Phased starting from version 0.7.4, so no need for not phased information
            if self.allele1 != ".":
                self.GT.append((2 * int(self.allele1)) + 2)
            else:
                self.GT.append(0)
            if self.allele2 != ".":
                self.GT.append((2 * int(self.allele2)) + 2)
            else:
                self.GT.append(0)
        log.debug("returning GT? : " + str(self.GT))
        return self.GT


def usage(scriptname, opts):
    print("\nUSAGE: \npython3 " + scriptname + '  --help')
    print("python3 " + scriptname + " -i octopus.somatic.snvs.pass.vcf  --tumor_column 11 --normal_column 10  -o updated_vcf.vcf\n")
    
    print("")
    print("options available:")
    print(" -i|--fvcf  [ Mandatory, no default value, String Filename full or relative path expected ]\n",
          "--tumor_column  [ Mandatory, no default value, Integer Expected ]\n",
          "--normal_column [ Mandatory, no default value, Integer Expected ]\n",
          "-o|--outfilename  [ Optional, no default value, String Expected ]\n",
          "--threshold_AR [ Optional; default value:0.9 ; float expected ]\n")
    print("")
    
    print("#" * 40 + "\nWARNING WARNING\n" + "#" * 40)
    print("1) This script is to be used only with somatic snvs vcf having NORMAL sample in "
          "column 10 and TUMOR sample in column 11; if not like this, update your vcf file to these "
          "specs;\n2) and the vcf has to be decomposed as well.\n")


def parseArgs(scriptname, argv):
    new_vcf_name = None
    column_tumor, column_normal = None, None
    
    warnings.simplefilter(action='ignore', category=FutureWarning)
    FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
    log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)
    # log.basicConfig(format=FORMAT_LOGGING, level=log.DEBUG)
    
    try:
        opts, args = getopt.getopt(argv, "hi:t:o:", ["vcf=", "tumor_column=", "normal_column=", "debug", "outfilename=", "threshold_AR=", "help"])
        log.info(opts)
    
    except getopt.GetoptError:
        usage(scriptname, opts)
        exit(2)
    for opt, arg in opts:
        if opt == '-h' or opt == "--help":
            usage(scriptname, opts)
            exit()
        elif opt in ("", "--threshold_AR"):
            try:
                global AR_threshold_for_GT
                AR_threshold_for_GT = round(float(arg), 4)
            except TypeError:
                log.info("ERROR: threshold values MUST be integer or float")
                exit(2)
        elif opt in ("", "--tumor_column"):
            column_tumor = int(arg)
        elif opt in ("", "--normal_column"):
            column_normal = int(arg)
        elif opt in ("-o", "--outfilename"):
            new_vcf_name = arg.strip()
        elif opt in ("-i", "--vcf"):
            fvcf = arg
            if not path.exists(fvcf):
                exit("ERROR: FNF --> " + fvcf)
        else:
            exit("Unknown Option --> " + opt)
    
    if column_tumor is None or column_normal is None:
        usage(scriptname, opts)
        exit(
            "Please Provide column number for tumor and normal Samples; should be 10 and 11  - or -  11 and 10;\n"
            "options are: --tumor_column and --normal_column;\nAborting. ")
    if column_normal == column_tumor:
        exit("ERROR: number for the columns Tumor and Normal MUST be different")
    
    return fvcf, column_tumor, column_normal, new_vcf_name


def update_header(vcf):
    """
    We modify the current header in the vcf object with the new fields or modifying old fields
    :param vcf: cyvcf2 VCF object
    :return vcf: cyvcf2 VCF object
    """
    # if Adding Fields to INFO field
    vcf.add_info_to_header(
        {'ID': 'OGT', 'Description': ''.join([
            'Original Octopus GT fields for each sample before reassigning the GT value based on AR threshold (GT_0/1 < AR_',
            str(AR_threshold_for_GT), ' and GT_1/1 >= AR_', str(AR_threshold_for_GT), ' )']), 'Type': 'String', 'Number': '.'})
    
    # adding field to the FORMAT columns
    # NOTE: if the field already exist in the Header, it will not be replaced or update; You must rename the Field already present in the VCF to add specifically the
    # following fields to the vcf HEADER
    # Adding AR
    vcf.add_format_to_header(
        {'ID': 'AR', 'Description': 'Alt Allelic Ratios for each sample in same order as list of samples found in VCF beyond column FORMAT', 'Type': 'Float', 'Number': '1'})
    return vcf


def is_obj_nan(obj):
    log.debug("obj test for nan: {}".format(str(obj)))
    if isnan(obj):
        return True
    return False


def get_GT_value_from_AR_Octopus_version_0_7_4_or_up(AR_value):
    """
        Because Octopus Version is ALL PHASED we do not need to consider whether it is or not phased for the GT
        return the GT value according to AR threshold values
        This is based off TGen's current thresholds of assigning the genotype
        1/1 if AR>=0.90 and 0/1 if AR<0.90
    
        Genotype representation for cyvcf2
        0 --> Unknown ; 1 --> Unknown_phased
        2 --> 0     ; 3 --> 0_PHASED_if_secondValue
        4 --> 1     ; 5 --> 1_PHASED_if_secondValue
        6 --> 2     ; 7 --> 2_PHASED_if_secondValue
        [0,0] == ./. ; [1,1] == .|.
        [1,0] == ./. ; [0,1] == .|.
        [2,2] == 0/0 ; [2,2] == 0/0
        [2,3] == 0|0 ; [3,2] == 0/0
        [2,4] == 0/1 ; [4,2] == 1/0
        [2,5] == 0|1 ; [5,2] == 1/0
        [2,6] == 0/2 ; [6,2] == 2/0
        [3,6] == 0/2 ; [6,3] == 2|0
        [4,6] == 1/2 ; [6,4] == 2/1
        [5,6] == 1/2 ; [6,5] == 2|1
        [9,8] == 3/3 ; [8,9] == 3|3
    
        [2,2] == 0/0 ; [4,4] == 1/1
        [3,3] == 0|0 ; [5,5] == 1|1
        [6,6] == 2/2 ; [7,7] == 2|2
        [8,8] == 3/3 ; [9,9] == 3|3
    
        return [int(2),int(4)] ; --> 0/1
        return [int(4),int(4)] ; --> 1/1
    
        """
    log.debug("0.7.4  AR =" + str(AR_value) + " ---  AR_threshold = " + str(AR_threshold_for_GT))
    
    try:
        if AR_value < AR_threshold_for_GT:
            # case we consider as HET or HOM_REF PHASED
            return [2, 5]
        else:
            return [5, 5]
    except ValueError:
        log.error("ERROR: AR value not a number")
    except TypeError:
        log.error("ERROR: AR value not of type Float")


def get_GT_value_from_GT_value(GT_value):
    """
		return the numpy array compatible GT value according to string GT value
		See also Genotype representation for cyvcf2 in Class Genotype
	"""
    
    dico_mapping_GT = {
        "./.": [0, 0],
        "0/0": [2, 2],
        "0/1": [2, 4], "1/0": [4, 2], "1/1": [4, 4],
        "0/2": [2, 6], "2/0": [6, 2], "2/2": [6, 6],
        "0/3": [2, 8], "3/0": [8, 2], "3/3": [8, 8],
        ".|.": [1, 1],
        "0|0": [3, 3],
        "0|1": [2, 5], "1|0": [4, 3],
        "0|2": [2, 7], "2|0": [6, 3],
        "0|3": [2, 9], "3|0": [8, 3], "3|3": [8, 9],
        "1|1": [4, 5], "1|2": [4, 7], "1|3": [4, 9], "1|4": [4, 11],
        "2|2": [6, 7], "2|3": [6, 9], "2|4": [6, 11], "2|5": [6, 13],
        "0|0|1": [3, 3, 5], "0|1|0": [3, 5, 3], "1|0|0": [5, 3, 3],
        "0|0|2": [3, 3, 7], "0|2|0": [3, 7, 3], "2|0|0": [7, 3, 3],
        "0|0|3": [3, 3, 9], "0|3|0": [3, 9, 3], "3|0|0": [9, 3, 3],
        "0|1|2": [3, 5, 7], "0|2|1": [3, 7, 5], "1|2|0": [5, 7, 3],
        "2|1|0": [7, 5, 3], "2|0|1": [7, 3, 5], "1|0|2": [5, 3, 7]
    }  # unused value ; here kept only for the mapping information
    
    x = GenotypeInv(list(GT_value))
    try:
        return x.get_gt_numpy_compatible()
    except ValueError:
        log.error("ERROR: GT value ")
    except TypeError:
        log.error("ERROR: GT value not of right type ")
    except Exception as e:
        log.error("ERROR: Unknown Error ; Check with the Author :-( ; " + str(e))


def process_GTs(tot_number_samples_p, v, col_tumor, col_normal):
    """
    Reassign GT value based on ala TGen threshold for AR value using _th_AR_for_GT_ CONSTANT
    :param tot_number_samples_p: (p is for parameter)
    :param v: variant record
    :param col_tumor:
    :param col_normal:
    :return: updated variant record
    """

    if tot_number_samples_p != 2:
        msg = "Expected 2 Samples in VCF found {}. We are suppose to process VCF file as a SOMATIC vcf and expect two samples;  Aborting.".format(tot_number_samples_p)
        raise Exception(msg)
    # capturing original GTs and adding them to INFO field
    # v.genotypes is a list of list example: v.genotypes: [[0, 0, 1, True], [0, 0, True]]  and it stype : <class 'list'>
    log.debug("v.genotypes: {}  and it stype : {}".format(str(v.genotypes), str(type(v.genotypes))))
    v.INFO["OGT"] = ','.join([str(Genotype(li)) for li in v.genotypes])
    # Re-Assigning GT with value based on AR thresholds comparison to CONSTANT threshold value
    GTs = [[0], [0]]  # need to init  list as we used index later for the list to replace values
    GTOs = [str(Genotype(li)) for li in v.genotypes]
    ARs = v.format('AR')
    idxN = 0 if col_normal == 10 else 1
    idxT = 1 if col_tumor == 11 else 0
    
    log.debug("xx" * 50 + "\nGTOs[idxN]: {}  and ARs[idxT]: {} and finally GTOs[idxT]:  {} ".format(GTOs[idxN], ARs[idxT], GTOs[idxT]))
    # we need to keep the order of the information based on the index; so the list GTs MUST be ordered;
    GTs[idxN] = get_GT_value_from_GT_value(GTOs[idxN])  # we do not modify the GT field for the Normal sample
    # GTs  # we do modify the GT field for the Tumor Sample based on defined threshold (version Octopus older than 0.7.4)
    GTs[idxT] = get_GT_value_from_AR_Octopus_version_0_7_4_or_up(ARs[idxT])
    
    v.set_format('GT', np.array(GTs))
    log.debug("v after reassigning GT: " + str(v))
    return v


def check_if_PS_in_FORMAT_field(vcf_cyobj, input_vcf_path, new_vcf_name, list_of_fields_to_check):
    try:
        # v1 = next(iter(vcf_cyobj))
        v1 = vcf_cyobj
        
        # Minimum_Expected_Fields_in_FORMAT_not_MANAGE_by_the_CODE
        ExpectedFlags = "GT:DP"
        for FIELD in list_of_fields_to_check:
            if FIELD not in v1.FORMAT:
                log.warning(
                    FIELD + " tag is ABSENT from the FORMAT field of OCTOPUS\nPlease Check you have run Octopus 0.7.4 or up "
                            "with appropriate options (i.e, with random forest option and no other filtering option; "
                            "random forest filtering add ALL the adequate fields to FORMAT columns )")
                if ':'.join(v1.FORMAT) == ExpectedFlags:
                    log.info("FORMAT field is equivalent to {}, but we request to have these flags at least: {}; "
                             "And we try to manage the rest, such as AD and AR; Would be better if you could run "
                             "random_Forest_Filtering; ".format(v1.FORMAT, ExpectedFlags))
                    log.warning(
                        "We assume the vcf has already been prepared for vcfMerger2 and therefore just copy the vcf by assigning the decomposed expected filename output")
                    from shutil import copyfile
                    copyfile(input_vcf_path, new_vcf_name)
                    exit()
                else:
                    log.error(FIELD + " flag NOT found in FORMAT field; Aborting VCF preparation.")
                    exit(FIELD + " flag Absent")
            else:
                log.info(FIELD + " flag Found")
    except StopIteration as sti:
        log.info(sti)
        log.warning("no records")
    log.info("Checking PS flag presence in FORMAT ...")


def if_dot_assign_value_zero(obj, idx):
    try:
        if obj == "." or is_obj_nan(float(obj)):
            return 0
        elif str(obj) != "." and str(obj) != "-2147483648" and str(obj) != "./.":
            return obj
        else:
            # Negative value was causing issue and adding more bytes to a file; So we switch to Zero number
            return 0
    except TypeError:
        return 0


def add_new_flags_v0_7_4(v, column_tumor, column_normal, tot_number_samples):
    """
    Calculate the AR for each sample in the variant record v
    The Total number of Sample in the VCF file is given by the variable tot_number_samples
    We assume that the number of sample does not vary form one record to another as recommended in VCF specs
    """

    idxT = 0 if int(column_tumor) == 10 else 1
    idxN = 1 if int(column_normal) == 11 else 0
    log.debug("___".join(str(x) for x in [idxT, idxN]))
    
    # in that version 0.7.4 of Octopus, the fields AD and AF and ADP are already present and AD, AF and ADP now follow the specs.
    # For each allele, there is a number now except for ADP with only one value which is supposed to represent Depth for an assigned allele but actually represent the depth for
    # the current sample in the FORMAT field whereas DP in INFO represent the total depth across all the samples (normally)
    # ##FORMAT=<ID=ADP,Number=1,Type=Integer,Description="Number of reads overlapping the position that could be assigned to an allele">
    # ##FORMAT=<ID=AF,Number=R,Type=Float,Description="Empirical allele frequency (AD / ADP)">
    # ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Empirical allele depth">
    # So we do not need to check if AD, AF or ADP are present, we assume they are there in our use case from the TGen pipeline
    
    if 'AD' in v.FORMAT:
        log.debug(str(v))
        log.debug("v.format('AF') : {}".format(v.format('AF')))
        log.debug("v.format('AF')[idxT] : {}".format(v.format('AF')[idxN]))
        log.debug("v.format('AF')[idxN] : {}".format(v.format('AF')[idxT]))
        
        AF_tumor = if_dot_assign_value_zero(v.format('AF')[idxT], 0)  # idxT to capture the list belonging to the tumor and 1 to capture the value for the ALT in case it is a
        # dot ## need to make it generic for both allele, i.e. testing both in case the ref is also a dot
        AF_normal = if_dot_assign_value_zero(v.format('AF')[idxN], 0)
        AR_tumor = float(AF_tumor)
        AR_normal = float(AF_normal)
        
        if idxT == 0:
            ARs = [AR_tumor, AR_normal]
        else:
            ARs = [AR_normal, AR_tumor]
    else:  # the else is entered only if there is no 'AD' flag in the octopus vcf
        log.error("AD not Found in the VCF; Are you sure you are using a vcf generated by a versino of Octopus 0.7.4 or up?; Aborting.")
        exit(2)
    
    # checking the values after processing and before adding them to the variant object v
    log.debug("ARs  are  : " + str(ARs))
    v.set_format('AR', np.array(ARs))
    log.debug("updated v object with new ARs: " + str(v))
    
    return process_GTs(tot_number_samples, v, column_tumor, column_normal)


if __name__ == "__main__":
    
    vcf_path, column_tumor, column_normal, new_vcf_name = parseArgs(argv[0], argv[1:])  # tth means tuple of thresholds
    vcf = VCF(vcf_path)
    vtest = None
    counter_rec = 0
    
    log.info(' '.join(["Constant AR threshold is: ", str(AR_threshold_for_GT)]))
    
    if new_vcf_name is None:
        new_vcf = '.'.join([str(vcf_path), "uts.vcf"])
    else:
        new_vcf = new_vcf_name
    
    # test if no variants in vcf:
    try:
        vtest = next(iter(vcf))  # we just try to check if we have no variant at all in VCF
    except StopIteration as si:
        # we bring the header up to specs and write down the output file
        log.warning("No Variants found in VCF; Creating Final Empty VCF now ..." + str(si))
        vcf = update_header(vcf)
        w = Writer(new_vcf, vcf)
        w.close()
        exit()
    
    # checking if PS flag is still present in the VCF genotype fields
    # vcf = VCF(vcf_path)  # as we have already consumed once the generator; if in case there is only one variant (edge case encountered already),
    # we need to re-init vcf here; need to think about another way of doing these tests without consuming the vcf object
    # vcf = update_header(vcf)
    check_if_PS_in_FORMAT_field(vtest, vcf_path, new_vcf_name, ["PS"])
    
    vcf = VCF(vcf_path)  # as we have already consumed twice the generator; we do not want to lose any variant, so we need to scan the VCF file again
    vcf = update_header(vcf)
    
    # create a new vcf Writer using the input vcf as a template.
    w = Writer(new_vcf, vcf)
    
    tot_number_samples = len(vcf.samples)
    if tot_number_samples != 2:
        exit("ERROR: Number of Sample greater than 2; Expected 2 samples only TUMOR and NORMAL")
    
    log.info("looping over records ...")
    for v in vcf:  # v for variant which represents one "variant record"
        log.debug(str(v))
        counter_rec = counter_rec + 1
        if counter_rec % 1000 == 0:
            log.info("processed rec: {}".format(str(counter_rec)))
        log.debug("v object: ORIGINAL RECORD: ...")
        log.debug(str(v))
        
        v = add_new_flags_v0_7_4(v, column_tumor, column_normal, tot_number_samples)
        
        log.debug("v object right before writing it to the VCF ...")
        log.debug(str(v))
        log.debug("---")
        if v is not None:
            w.write_record(v)
    
    w.close()
    vcf.close()
    log.info("work completed")
    log.info('new vcf is << {} >>'.format(new_vcf))
    exit()
