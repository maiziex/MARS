#!/usr/bin/env python
import pdb
#pdb.set_trace()
import os
import sys
import numpy as np
from argparse import ArgumentParser
import pickle
from subprocess import Popen
from multiprocessing import Pool,cpu_count,active_children,Manager
import time
__author__ = "Maizie X.Zhou@Stanford"
parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--assembly_dir','-i_dir', help="The input folder where you store the diploid assembled contig files",required=True)
parser.add_argument('--out_dir','-o_dir', help="The folder name you can define to store the final results",default="./Results_ancestral_Aquila")
parser.add_argument('--ref_file','-r', help="The human reference fasta file",required=True)
parser.add_argument('--SV_len','-l',type=int,help="The SV size you can define", default=20)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--num_threads','-nt', help=" The number of threads you can define, which corresponds to number of samples",default=2)
parser.add_argument('--sample_list','-sl', help='The sample names corresponding to your contig files, which is the prefix of the contig files', type=str,required=True)
args = parser.parse_args()
if len(sys.argv) == 1:
    Popen("python Run_all_samples_Aquila.py -h",shell=True).wait()
else:
    sample_list = [item for item in args.sample_list.split(',')]
    script_path = os.path.dirname(os.path.abspath( __file__ ))
    code_path = script_path + "/" 
    other_tools_path = (os.path.abspath(os.path.join(code_path, '..'))) + "/"


def Call_one_sample(assembly_dir,out_dir,ref_file,sample_name,SV_len,num_of_threads,all_regions_flag,chr_start,chr_end,xin):
    use_cmd = code_path + "Aquila_assembly_based_variants_call.py " + " --assembly_dir " + assembly_dir + sample_name + "/"  + " --ref_file " + ref_file + " --num_of_threads " + str(num_of_threads) + " --all_regions_flag " + str(all_regions_flag) + " --out_dir " + out_dir + sample_name + "/" + " --var_size " + str(SV_len) + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end)
    Popen(use_cmd,shell=True).wait()


def Run_all(sample_list,out_dir,in_dir,ref_file,SV_len,num_of_threads,chr_start,chr_end,all_regions_flag):
    all_cmd = ""
    count = 1
    total_num = len(sample_list)
    pool = Pool(processes=num_of_threads)
    for sample_name in sample_list:
        pool.apply_async(Call_one_sample,(in_dir,out_dir,ref_file,sample_name,SV_len,1,all_regions_flag,chr_start,chr_end,"xin"))
        #Call_one_sample(in_dir,out_dir,ref_file,sample_name,SV_len,1,all_regions_flag,chr_start,chr_end,"xin")
        count += 1
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()
            if (count - 1) == total_num:
                print("finished all samples")
            else:
                pool = Pool(processes=num_of_threads)
    print("all done~")


def main():
    if len(sys.argv) == 1:
        Popen("python3 " + "MARS_step1.py -h",shell=True).wait()
    else:
        out_dir = args.out_dir + "/"
        if os.path.exists(out_dir):
            print("using existing output folder: " + out_dir)
        else:
            os.makedirs(out_dir)
        in_dir = args.assembly_dir + "/"
        ref_file = args.ref_file
        chr_start = args.chr_start
        chr_end = args.chr_end
        all_regions_flag = 1 
        SV_len = args.SV_len
        num_of_threads = int(args.num_threads)
        Run_all(sample_list,out_dir,in_dir,ref_file,SV_len,num_of_threads,chr_start,chr_end,all_regions_flag)




if __name__ == "__main__":
    main()

