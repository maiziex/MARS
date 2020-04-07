import pdb
#pdb.set_trace()
import glob
import os
from subprocess import Popen, PIPE

def multiseq_alignment(in_dir):
    files_all = glob.glob(in_dir +  "*.fasta")
    for one_file in files_all:
        print(one_file)
        output_file = one_file.split(".fasta")[0] + ".html"
        use_cmd = "/oak/stanford/groups/arend/Xin/Software/muscle3.8.31_i86linux64 -in " + one_file + " -out " + output_file + " -html"
        Popen(use_cmd,shell=True).wait()


        output_file = one_file.split(".fasta")[0] + ".afa"
        use_cmd = "/oak/stanford/groups/arend/Xin/Software/muscle3.8.31_i86linux64 -in " + one_file + " -out " + output_file 
        Popen(use_cmd,shell=True).wait()
    print("done~")


#multiseq_alignment("/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_Aquila/MSA_del_results_chr21/")
multiseq_alignment("/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_Aquila/MSA_SV_results_chr21/")
