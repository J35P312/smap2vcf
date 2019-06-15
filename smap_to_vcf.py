
description='''Stand-alone script to convert BioNano smap file format to vcf.'''

vcfheader_1='''##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=IRYSTYPE,Number=1,Type=String,Description="Type of structural variant"> 
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">   
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=QUERY,Number=1,Type=Integer,Description="Xmap query ID">
##ALT=<ID=DEL,Description="Deletion"> 
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Break-end">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'''


import os
import argparse
import math


#create .vcf at smappath (replacing suffix '.smap' with '.vcf')
#currently only processing insertions and deletions
def smap_to_vcf(smappath, sample, vcfh) :
    doconf = True #put Confidence in QUAL field
    defaultconf = "." #missing value in vcf format: used if !doconf
    maxconf = "20" #if doconf, this is maximum reported: equivalent to ppv of 0.99
    
    colhead = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t" #this is complete if no genotype
    colform = "FORMAT\t%s" % sample #if genotype, add this

    print vcfh
    print colhead+colform

    f1 = open(smappath)
    nent = 0
    conf = defaultconf #if !doconf, this is used
    for line in f1 :
        if line[0] == "#" :
            continue
        tokens = line.split()
        if len(tokens) < 17 :
            print "ERROR: line incomplete, terminating:\n%s\n" % line
            break
        svtype = tokens[9] 
        
        vcftype = "BND"
        if svtype == "deletion":
            vcftype="DEL"
        elif svtype == "insertion":
            vcftype = "INS"
        elif "inversion" in svtype:
            vcftype = "INV"
        
        # 0 : SmapEntryID
        qry = tokens[1] #QryContigID
        ref = int(tokens[2]) #RefcontigID1 -- must be int bc refcmap keys are ints
        # 3 : RefcontigID2 
        refstart = float(tokens[6]) #RefStartPos
        refstop  = float(tokens[7]) #RefEndPos

        conf = tokens[8] #Confidence
        
            #tokens[18] = Genotype
        if int(tokens[19]) > 0 : #GenotypeGroup
            if tokens[17] == "homozygous" : #zygosity
                gt = "1/1"
            else : #note: treating 'unknown' same as 'heterozygous'
                gt = "1/2"
        else : #note that GenotypeGroup should always be -1 or > 0
            gt = "0/1"

        refB=ref
        if not tokens[2] == tokens[3]:
            refB= int(tokens[3])
        if refB == -1:
           continue        

        pos=int(float(tokens[6]))
        end=int(float(tokens[7]))
        #note: two inner abs should not be necessary
        if refB == ref:
            svlen = int(abs(refstop-refstart))

        if ref== 23:
            ref="X"
        elif ref == 24:
            ref= "Y"

        nent += 1
        if not vcftype == "BND":
            pos=int(float(tokens[6]))
            end=int(float(tokens[7]))
            if end < pos:
               tmp = end
               end = pos
               pos = tmp
            INFO="END={};SVLEN={};SVTYPE={};IRYSTYPE={};QUERY={}".format(end,svlen,vcftype,svtype,qry)
            vcf_line=[str(ref),str(pos),str(nent),"N","<" + vcftype + ">",str(conf),"PASS",INFO,"GT",gt]
            print "\t".join(vcf_line)
        else:

            chrB= int(tokens[3])
            if chrB== 23:
                chrB="X"
            elif chrB == 24:
                chrB= "Y"
            posB=end
            INFO="SVTYPE={};IRYSTYPE={};QUERY={}".format(vcftype,svtype,qry)
            vcf_line=[str(ref),str(pos),str(nent),"N","N[" + str(chrB) + ":" + str(posB) + "[",conf,"PASS",INFO,"GT",gt]
            print "\t".join(vcf_line)
    #end loop on input smap
    f1.close()
#end smap_to_vcf

#copied from utilities.py
def checkFile(filepath, filesuff="", checkReadable=True) :
    try :
        valid = os.path.isfile(filepath)
        if filesuff != "" :
            valid = (valid and filepath.endswith(filesuff))
        if checkReadable :
            valid = (valid and os.access(filepath, os.R_OK))
    except :
        valid = False
    #print "checkFile:", filepath, valid #debug
    return valid



def getArgs() :    
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-s', dest='smappath'   , help='Path to smap file to convert (required)', type=str) 
    defsamp = "Sample"
    parser.add_argument('-n', dest='sample'     , help='Sample ID name for genotype data (optional, default "%s")'%defsamp, type=str, default=defsamp)
    result = parser.parse_args()

    smappath    = result.smappath
    sample      = result.sample

    if not checkFile(smappath,".smap") :
        print "ERROR: smap does not exist, is not readable, or does not end with '.smap':", smappath
        return None

    return smappath, sample


def run_smap_to_vcf():
    getargs = getArgs()
    if getargs != None :
        #smappath, refcmap, sample
        #smap_to_vcf(smappath, refcmap, sample, vcfh=vcfheader_1)
        smap_to_vcf(getargs[0], getargs[1], vcfh=vcfheader_1)


if __name__ == "__main__" :
    run_smap_to_vcf()
