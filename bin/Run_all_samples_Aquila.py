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
parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--in_dir','-i_dir', help="The input folder where you store the diploid assembled contig files",required=True)
parser.add_argument('--out_dir','-o_dir', help="The folder name you can define to store the final results",default="./Results_ancestral_Aquila")
parser.add_argument('--ref_file','-r', help="The human reference fasta file which can be download by running \"./install.sh\". ",required=True)
parser.add_argument('--SV_len','-l',type=int,help="The SV size you can define", default=20)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--all_regions_flag','-f',type=int,help="call variant in diploid regions", default=0)
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

def read_ref(fasta_file,chr_num,out_dir):
    f = open(fasta_file,"r")
    count = 0
    ref_seq = ""
    for line in f:
        if count > 0:
            data = line.rsplit()
            ref_seq += data[0]
        count += 1
    print("total_len for chr" + str(chr_num))
    pickle.dump(ref_seq, open(out_dir + "ref_seq_chr" + str(chr_num) +  ".p","wb"))


def extract_ref_chr(ref_file,chr_num,out_dir):
    fw = open(out_dir + "genome_ref_chr" + str(chr_num) + ".fasta","w")
    f = open(ref_file,"r")
    flag = 0
    total_len = 0
    for line in f:
        data = line.rsplit()
        if data[0] == ">chr" + str(chr_num):
            fw.writelines(">" + str(chr_num) + "\n")
            flag = 1
        elif data[0][0] == ">" and flag == 1:
            break
        else:
            if flag == 1:
                total_len += len(data[0])
                fw.writelines(data[0] + "\n")
    print("chr" + str(chr_num) + ":")
    print(total_len)


def Call_one_sample(assembly_dir,out_dir,ref_file,sample_name,SV_len,num_of_threads,all_regions_flag,chr_start,chr_end,xin):
    #use_cmd = code_path + "Assembly_based_variants_call.py " + " --contig_hp1_fasta " + hp1_file + " --contig_hp2_fasta " + hp2_file + " --out_dir " + out_dir + " --ref_file " + ref_file + " --sample_name " + sample_name  + " --SV_len " + str(SV_len)  + " --species_ref_list " + ref_list + " --species_name_list " + species_list 
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


def Merge_vcf_for_all_samples_del(out_dir):
    merge_cmd = "bcftools merge " 
    for sample_name in sample_list:
        origin_file = out_dir + sample_name + "_final_del.vcf"  
        sorted_file = out_dir + sample_name + "_final_del_sorted.vcf"  
        zip_file = out_dir + sample_name + "_final_del_sorted.vcf.gz"  
        sort_cmd = "bcftools sort " + origin_file + " > " + sorted_file
        zip_cmd = "bgzip -c " + sorted_file + " > " + zip_file
        idx_cmd = "tabix -p vcf " + zip_file
        Popen(sort_cmd,shell=True).wait()
        Popen(zip_cmd,shell=True).wait()
        Popen(idx_cmd,shell=True).wait()
        merge_cmd +=  zip_file + " " 
    merge_cmd += " -o " + out_dir + "Merged_del.vcf"
    Popen(merge_cmd,shell=True).wait()
    print(merge_cmd)
    print("all done")   


def Merge_vcf_for_all_samples_ins(out_dir):
    merge_cmd = "bcftools merge " 
    for sample_name in sample_list:
        origin_file = out_dir + sample_name + "_final_ins.vcf"  
        sorted_file = out_dir + sample_name + "_final_ins_sorted.vcf"  
        zip_file = out_dir + sample_name + "_final_ins_sorted.vcf.gz"  
        sort_cmd = "bcftools sort " + origin_file + " > " + sorted_file
        zip_cmd = "bgzip -c " + sorted_file + " > " + zip_file
        idx_cmd = "tabix -p vcf " + zip_file
        Popen(sort_cmd,shell=True).wait()
        Popen(zip_cmd,shell=True).wait()
        Popen(idx_cmd,shell=True).wait()
        merge_cmd +=  zip_file + " " 
    merge_cmd += " -o " + out_dir + "Merged_ins.vcf"
    Popen(merge_cmd,shell=True).wait()
    print(merge_cmd)
    print("all done")   


if __name__ == "__main__":
    if len(sys.argv) > 1:
        out_dir = args.out_dir + "/"
        if os.path.exists(out_dir):
            print("using existing output folder: " + out_dir)
        else:
            os.makedirs(out_dir)
        in_dir = args.in_dir + "/"
        ref_file = args.ref_file
        chr_start = args.chr_start
        chr_end = args.chr_end
        all_regions_flag = args.all_regions_flag
        SV_len = args.SV_len
        num_of_threads = int(args.num_threads)
        """
        for chr_num in range(1,23):
            extract_ref_chr(ref_file,chr_num,out_dir)
        extract_ref_chr(ref_file,"X",out_dir)
        for chr_num in range(1,23):
            read_ref(out_dir + "genome_ref_chr" + str(chr_num) + ".fasta",chr_num,out_dir)
        read_ref(out_dir + "genome_ref_chrX.fasta",23,out_dir)
        """
        Run_all(sample_list,out_dir,in_dir,ref_file,SV_len,num_of_threads,chr_start,chr_end,all_regions_flag)

        #Merge_vcf_for_all_samples_del(out_dir)
        #Merge_vcf_for_all_samples_ins(out_dir)

