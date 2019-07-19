import os, sys, glob

'''
    AUTHOR: SHWETHA HARA SRIDHAR
    TECHDEV TEAM, DEPARTMENT OF GENETICS AND GENOMICS
    ICAHN SCHOOL OF MEDICINE AT MOUNT SINAI, NY
    VERSION: 19.7.18

    INPUT: .CSV WITH COMBINED COUNTS FOR ALL THE GENES (OUTPUT FROM STEP1: COMBINE REPLICATES COUNTS)
    METHOD: USES GTF ANNOTATION TO RENAME GENE IDS WITH GENE NAMES AND COMBINED THE COUNTS FROM TWO CONDITIONS INTO ONE TABLE
    OUTPUT: _INPUT_FOR_DESEQ.CSV WITH ALL THE CONDITION COUNTS COMBINED (USE FOR DESEQ2 RSCRIPT INPUT)
    
'''


def make_gtf_dict(gtf_file):

    gtf_dictnry = {}

    with open(gtf_file, "r") as gtf:
        for _ in range(5):
            next(gtf)
        for line in gtf:
            line = line.strip().split("\t")[-1].split(";")
            geneid = line[1].split(" ")[2].lstrip("\"").rstrip("\"")
            genename = line[2].split(" ")[2].lstrip("\"").rstrip("\"")
            if geneid in gtf_dictnry:
                gtf_dictnry[geneid] = genename
            else:
                gtf_dictnry[geneid] = genename

    return gtf_dictnry


def make_finalcounts(conditionA, conditionB, outpath, gtf_dict):

    outfile = open(outpath + "_Input_for_DEseq.csv", "w")
    sys.stdout = outfile

    for lineA, lineB in zip(open(conditionA, "r"), open(conditionB, "r")):
        lineA = lineA.strip().split(",")
        lineB = lineB.strip().split(",")
        if lineA[0] == lineB[0]:
            ensemble_id = lineA[0]
            if gtf_dict and ensemble_id in gtf_dict:
                print(gtf_dict[ensemble_id] + ", " + ", ".join(lineA[1:len(lineA) - 1]) + ", " + ", ".join(lineB[1:len(lineB) - 1]))
            else:
                print(ensemble_id + ", " + ", ".join(lineA[1:len(lineA)-1]) + ", " + ", ".join(lineB[1:len(lineB)-1]))
    outfile.close()

    return



## MAIN ##


path = "/Candida/Candida_featurecounts"
files = os.listdir(path)
print(files)
for file in files:
    combine_counts(path+"/"+file)

# input files
conditionA = '/Candida/Candida_featurecounts/Resistant/Resistant_CombinedCounts.csv'
conditionB = '/Candida/Candida_featurecounts/Sensitive/Sensitive_CombinedCounts.csv'
gtf_file = "/Candida/Ref_GTF/GCA_002759435.2_Cand_auris_B8441_V2_genomic.gtf.txt"

gtf_dict = make_gtf_dict(gtf_file)

out_path = ("/".join(conditionA.split("/")[0:len(conditionA.split("/"))-2]))\
           +"/"+(conditionA.split("/")[-2]+'_vs_'+conditionB.split("/")[-2])

make_finalcounts(conditionA, conditionB, out_path, gtf_dict)
