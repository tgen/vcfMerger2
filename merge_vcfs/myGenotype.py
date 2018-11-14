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
### Major Contributors: Christophe Legendre'a0
### Minor Contributors:


class Genotype(object):
    """
    Override class implementation borrowed from brentp @ https://github.com/brentp/cyvcf2/issues/58
    Override the string function of the Genotype Class
    """
    __slots__ = ('alleles', 'phased')

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)
    __repr__ = __str__

#gts = variant.genotypes
#genotypes = [Genotype(li) for li in gts]
#print genotypes
## which shows:
#[./., ./., ./., 1/1, 0/0, 0/0, 0/0]