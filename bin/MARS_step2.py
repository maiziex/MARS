#!/usr/bin/env python
import pdb
#pdb.set_trace()
from collections import defaultdict
import subprocess
from argparse import ArgumentParser
from subprocess import Popen, PIPE
from multiprocessing import Pool,cpu_count,active_children,Manager
import time
import pickle
import pysam
import os
import sys
import glob
from .Define_SV_and_Ref_within_break_points_3 import *
from .extract_TandemRepeats_dict_use_concensus import * 
from .Generate_SV_type_table_from_MSA_new import *
from .Extract_validated_SV_for_flankingseq import *
from .Evaluate_derived_sv_by_flankingseq_global_align import *
from .Extract_Apes_seqs import *
from .Print_linked_SNPs_for_SV import *
from .Merge_more_fields_to_final_table import *
from .generate_MSA_geno import *
parser = ArgumentParser(description="Author: maizie.zhou@vanderbilt.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--in_dir','-i_dir', help="The input folder where you store vcf files",required=True)
parser.add_argument('--out_dir','-o_dir', help="The folder name you can define to store the final results",default="./MARS_step2_results")
parser.add_argument('--assembly_dir','-a_dir', help="The folder to store assembly contigs",required=True)
parser.add_argument('--ref_file','-r', help="The human reference fasta file which can be download by running \"./install.sh\". ",required=True)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
#parser.add_argument('--num_threads_msa','-nt_msa', help=" The number of threads you can define to run MSA by muscle",default=20)
parser.add_argument('--num_threads_bychr','-nt_chr', help=" The number of threads you can define to run Step2 by chromosome",default=1)
parser.add_argument('--sample_list','-sl', help='The sample names corresponding to your contig files, which is the prefix of the contig files', type=str,required=True)
parser.add_argument('--Ape_ref_list','-lr', help='If --HARP_flag set to 1, the users need to iniatize the Ape reference genomes they want to use. "Gorilla_gorilla_ref.fasta", "pan_troglodytes_ref.fasta", "pongo_abelii_ref.fasta", and "macaca_mulatta_ref.fasta" are the reference fasta files for each Ape. Each reference file is seperately by comma (",") ', type=str,default="")
parser.add_argument('--HARP_flag', '-h_flag',help= 'If flag set to 1, it will output SV with ancestral state: derived deletions and derived insertions.',default= 0)
parser.add_argument('--gnomad_flag_for_linked_SNP', '-gnomad_flag',default= 0)
parser.add_argument('--gnomad_dir','-g_dir', help="If -gnomad_flag set to 1, the users need to download gnomad VCF files and give a path to the folder which stores these gnomAD VCF files",default="")

args = parser.parse_args()
if len(sys.argv) == 1:
    Popen("python MARS_step2.py -h",shell=True).wait()
else:
    sample_list = [item for item in args.sample_list.split(',')]
    num_of_samples = len(sample_list)
    HARP_flag = int(args.HARP_flag)
    gnomad_flag_for_linked_SNP = int(args.gnomad_flag_for_linked_SNP)
    if HARP_flag == 1:
        if args.Ape_ref_list == "":
            print("Error: Using HARP_flag = 1, please add --Ape_ref_list!!!")
            sys.exit()
        else:
            Ape_ref_list = [item for item in args.Ape_ref_list.split(',')]
    if gnomad_flag_for_linked_SNP == 1:
        if args.gnomad_dir == "":
            print("Error: Using gnomad_flag_for_linked_SNP = 1, please define gnomad folder for --gnomad_dir!")
            sys.exit()
        else:
            gnomad_dir = args.gnomad_dir + "/" 


    script_path = os.path.dirname(os.path.abspath( __file__ ))
    code_path = script_path + "/" 
    other_tools_path = (os.path.abspath(os.path.join(code_path, '..'))) + "/"

    print("Processing " + str(num_of_samples) + " samples...")

basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
def translate(seq):
    seq_r = []
    for character in seq:
        seq_r.append(basecomplement[character])
    seq_r_2 = "".join(seq_r)
    rseqn = seq_r_2[::-1]
    return rseqn


def load_and_save_Apes_ref(Ape_ref_dir,ref_file,Ape_name):
    count = 0
    ref_string = ""
    use_chr_num = ""
    ref_dict = defaultdict(str)
    with open(ref_file,"r") as f:
        for line in f:
            data = line.rsplit()
            if data[0][0] == ">":
                chr_num = data[0].split(">")[-1]
                if chr_num != use_chr_num and count > 0:
                    ref_dict[use_chr_num] = ref_string
                    pickle.dump(ref_string,open(Ape_ref_dir + Ape_name + "_chr" + use_chr_num + ".p","wb"))
                    ref_string = ""
                    ref_dict = defaultdict(str)
            else:
                use_chr_num = chr_num
                ref_string += data[0]
            count += 1
        ref_dict[use_chr_num] = ref_string
        pickle.dump(ref_string,open(Ape_ref_dir + Ape_name + "_chr" + use_chr_num + ".p","wb"))

    print("Extract and Save " + Ape_name +  "reference dict files done!")


def Extract_DEL_from_all_samples(sample_list,in_dir,out_dir,chr_start,chr_end):
    for chr_num_use in range(chr_start,chr_end + 1):
        print(chr_num_use)
        output_file = out_dir + "chr" + str(chr_num_use) + "_del.txt"
        output_file_sorted = out_dir + "chr" + str(chr_num_use) + "_del_sorted.txt"
        if chr_num_use == 23:
            chr_num_use = 23
        fw = open(output_file,"w")
        for sample_name in sample_list:
            vcf_file = in_dir  + sample_name + "/" + "Aquila_DEL_chr" + str(chr_num_use) + ".vcf"
            with open(vcf_file,"rb") as f:
                for line in f:
                    data = line.decode().rsplit()
                    if data[0][0] != "#":
                        chr_num = data[0]
                        if "chr" + str(chr_num_use) == chr_num:
                            start_ = int(data[1])
                            end_ = int(data[1]) + len(data[3])  # for deletion
                            SV_len = len(data[3])
                            event_name = data[2]
                            GT = data[-1].split(":")[0]
                            contig_info = data[-1].split(":")[1]
                            fw.writelines(sample_name + "_del" + "\t" + chr_num + "\t" + str(start_) + "\t" +  str(end_) + "\t" + str(SV_len) + "\t" +   event_name + "\t" + GT + "\t" + contig_info + "\n")
        fw.close()

        sort_cmd = "cat " + output_file + " | sort -k3n > " + output_file_sorted 
        Popen(sort_cmd,shell=True).wait()
    print("done")


def Extract_INS_from_all_samples(sample_list,in_dir,out_dir,chr_start,chr_end):
    for chr_num_use in range(chr_start,chr_end + 1):
        print(chr_num_use)
        output_file = out_dir + "chr" + str(chr_num_use) + "_ins.txt"
        output_file_sorted = out_dir + "chr" + str(chr_num_use) + "_ins_sorted.txt"
        if chr_num_use == 23:
            chr_num_use = 23
        fw = open(output_file,"w")
        for sample_name in sample_list:
            vcf_file = in_dir  + sample_name + "/" + "Aquila_INS_chr" + str(chr_num_use) + ".vcf"
            with open(vcf_file,"rb") as f:
                for line in f:
                    data = line.decode().rsplit()
                    if data[0][0] != "#":
                        chr_num = data[0]
                        if "chr" + str(chr_num_use) == chr_num:
                            start_ = int(data[1])
                            end_ = int(data[1])  # for insertions
                            SV_len = len(data[4])
                            event_name = data[2]
                            GT = data[-1].split(":")[0]
                            contig_info = data[-1].split(":")[1]
                            fw.writelines(sample_name + "_ins" + "\t" + chr_num + "\t" + str(start_) + "\t" +  str(end_) + "\t" + str(SV_len) + "\t" +   event_name + "\t" + GT + "\t" + contig_info + "\n")
        fw.close()

        sort_cmd = "cat " + output_file + " | sort -k3n > " + output_file_sorted 
        Popen(sort_cmd,shell=True).wait()
    print("done")


def Write_to_single_SV(save_data,fw,min_start,max_end,merge_DEL_dict,merge_INS_dict):
    new_save_data = []
    # if > 1: merge del and ins
    if len(save_data) > 1:
        sample_num_list = []
        event_name_list = []
        contig_info_list = []
        SV_size_list = []
        GT_list = []
        count = 1
        for data in save_data:
            # redefine the save_data after merge del and ins
            if data[1] == "del":
                new_save_data += merge_DEL_dict[(int(data[2]),int(data[3]))]
            elif data[1] == "ins":
                new_save_data += merge_INS_dict[(int(data[2]),int(data[3]))]

            chr_num = data[0]
            if count == 1:
                sv_type = data[1]
            else:
                sv_type += "_" + data[1]
            sample_num_list += data[4].split(":")
            event_name_list += data[5].split(":")
            contig_info_list += data[-2].split(":")
            GT_list += data[-1].split(":")
            SV_size_list += data[6].split(":")
            count += 1
        fw.writelines(chr_num + "\t" + sv_type + "\t" + str(min_start) + "\t" + str(max_end) + "\t" + ":".join(sample_num_list) + "\t" + ":".join(event_name_list) +"\t" + ":".join(SV_size_list) + "\t" +  ":".join(contig_info_list) + "\t" + ":".join(GT_list) + "\n")
    else:
        fw.writelines("\t".join(save_data[0]) + "\n")
        sv_type = save_data[0][1]
        for data in save_data:
            if data[1] == "del":
                new_save_data += merge_DEL_dict[(int(data[2]),int(data[3]))]
            elif data[1] == "ins":
                new_save_data += merge_INS_dict[(int(data[2]),int(data[3]))]


    return (fw, new_save_data,sv_type)


def Write_to_single_DEL(save_data,fw,min_start,max_end):
    sample_num_list = []
    event_name_list = []
    contig_info_list = []
    SV_size_list = []
    GT_list = []
    for data in save_data:
        chr_num = data[1]
        sample_num_list.append(data[0])
        event_name_list.append(data[5])
        contig_info_list.append(data[-1])
        GT_list.append(data[6])
        SV_size_list.append(data[4])
    fw.writelines(chr_num + "\t" + "del" + "\t" + str(min_start) + "\t" + str(max_end) + "\t" + ":".join(sample_num_list) + "\t" + ":".join(event_name_list) +"\t" + ":".join(SV_size_list) + "\t" +  ":".join(contig_info_list) + "\t" + ":".join(GT_list) + "\n")
    return fw


def Write_to_single_INS(save_data,fw,min_start,max_end):
    sample_num_list = []
    event_name_list = []
    contig_info_list = []
    SV_size_list = []
    GT_list = []
    for data in save_data:
        chr_num = data[1]
        sample_num_list.append(data[0])
        event_name_list.append(data[5])
        contig_info_list.append(data[-1])
        GT_list.append(data[6])
        SV_size_list.append(data[4])
    fw.writelines(chr_num + "\t" + "ins" + "\t" + str(min_start) + "\t" + str(max_end) + "\t" + ":".join(sample_num_list) + "\t" + ":".join(event_name_list) +"\t" + ":".join(SV_size_list) + "\t" +  ":".join(contig_info_list) + "\t" + ":".join(GT_list) + "\n")
    return fw


def Merge_DEL(out_dir,chr_start,chr_end):
    for chr_num_use in range(chr_start,chr_end + 1):
        merge_SV_dict = defaultdict(list)
        input_file = out_dir + "chr" + str(chr_num_use) + "_del_sorted.txt"
        output_file = out_dir + "chr" + str(chr_num_use) + "_merge_del.txt"
        fw = open(output_file,"w")
        count = 0
        save_data = []
        with open(input_file,"r") as f:
            for line in f:
                data = line.rsplit()
                _start = int(data[2])
                _end = int(data[3])
                if count == 0:
                    save_data.append(data)
                    min_start = _start
                    max_end = _end
                if _start <= max_end and count > 0:
                    save_data.append(data)
                    if _end > max_end:
                        max_end = _end
                elif _start > max_end and count > 0:
                    fw = Write_to_single_DEL(save_data,fw,min_start,max_end)
                    merge_SV_dict[(min_start,max_end)] = save_data
                    min_start = _start
                    max_end = _end
                    save_data = []
                    save_data.append(data)
                count += 1      
        fw = Write_to_single_DEL(save_data,fw,min_start,max_end)
        merge_SV_dict[(min_start,max_end)] = save_data
        pickle.dump(merge_SV_dict, open(out_dir + "merge_DEL_dict_chr" + str(chr_num_use)  + ".p", "wb"))
    return merge_SV_dict


def Merge_INS(out_dir,chr_start,chr_end):
    ins_threshold = 100
    for chr_num_use in range(chr_start,chr_end + 1):
        merge_SV_dict = defaultdict(list)
        input_file = out_dir + "chr" + str(chr_num_use) + "_ins_sorted.txt"
        output_file = out_dir + "chr" + str(chr_num_use) + "_merge_ins.txt"
        fw = open(output_file,"w")
        save_data = []
        save_start_list = []
        count = 0
        with open(input_file,"r") as f:
            for line in f:
                data = line.rsplit()
                _start = int(data[2])
                save_start_list.append(_start)
                if count == 0:
                    min_start = _start
                    save_data.append(data)
                if (_start - min_start) < ins_threshold and count > 0:
                    save_data.append(data)
                elif (_start - min_start) >= ins_threshold and count > 0:
                    max_end = save_start_list[-2]
                    fw = Write_to_single_INS(save_data,fw,min_start,max_end)
                    merge_SV_dict[(min_start,max_end)] = save_data
                    save_start_list = []
                    save_data = []
                    save_start_list.append(_start)
                    save_data.append(data)
                    min_start = _start
                count += 1
        max_end = save_start_list[-1]
        fw = Write_to_single_INS(save_data,fw,min_start,max_end)
        merge_SV_dict[(min_start,max_end)] = save_data
        pickle.dump(merge_SV_dict, open(out_dir + "merge_INS_dict_chr" + str(chr_num_use)  + ".p", "wb"))
    return merge_SV_dict


def Merge_DEL_and_INS(out_dir,chr_start,chr_end):
    for chr_num_use in range(chr_start,chr_end + 1):
        merge_DEL_dict = pickle.load(open(out_dir + "merge_DEL_dict_chr" + str(chr_num_use)  + ".p","rb"))
        merge_INS_dict = pickle.load(open(out_dir + "merge_INS_dict_chr" + str(chr_num_use)  + ".p","rb"))
        merge_SV_dict = defaultdict(list)
        del_file = out_dir + "chr" + str(chr_num_use) + "_merge_del.txt"
        ins_file = out_dir + "chr" + str(chr_num_use) + "_merge_ins.txt"
        sv_file = out_dir + "chr" + str(chr_num_use) + "_sv.txt"
        sv_file_sorted = out_dir + "chr" + str(chr_num_use) + "_sv_sorted.txt"
        output_file = out_dir + "chr" + str(chr_num_use) + "_merged_sv.txt"
        cat_cmd = "cat " + del_file + " " + ins_file + " > " + sv_file
        Popen(cat_cmd,shell=True).wait()
        sort_cmd = "cat " + sv_file + " | sort -k3n > " + sv_file_sorted 
        Popen(sort_cmd,shell=True).wait()
        fw = open(output_file,"w")
        count = 0
        save_data = []
        with open(sv_file_sorted,"r") as f:
            for line in f:
                data = line.rsplit()
                _start = int(data[2])
                _end = int(data[3])
                if count == 0:
                    save_data.append(data)
                    min_start = _start
                    max_end = _end
                if _start <= max_end and count > 0:
                    save_data.append(data)
                    if _end > max_end:
                        max_end = _end
                elif _start > max_end and count > 0:
                    fw,new_save_data,sv_type = Write_to_single_SV(save_data,fw,min_start,max_end,merge_DEL_dict,merge_INS_dict)
                    merge_SV_dict[(min_start,max_end,sv_type)] = new_save_data
                    min_start = _start
                    max_end = _end
                    save_data = []
                    save_data.append(data)
                count += 1      
        fw,new_save_data,sv_type = Write_to_single_SV(save_data,fw,min_start,max_end,merge_DEL_dict,merge_INS_dict)
        merge_SV_dict[(min_start,max_end,sv_type)] = new_save_data
        pickle.dump(merge_SV_dict, open(out_dir + "merge_SV_dict_chr" + str(chr_num_use)  + ".p", "wb"))
    return merge_SV_dict


def get_sam(fasta_file,sam_file,bam_file,chr_num,ref_dir,xin):
    ref_chr_mmi = ref_dir + "ref_chr" + str(chr_num) + ".mmi"
    ref_chr_fasta = ref_dir + "genome_ref_chr" + str(chr_num) + ".fasta"
    try:
        get_mmi_cmd = "minimap2 -d " + ref_chr_mmi + " " +  ref_chr_fasta
        align_cmd = "minimap2 -a " + ref_chr_mmi + " " + fasta_file + " > " + sam_file
        sort_cmd = "samtools view -Sb " + sam_file + " | samtools sort  > " + bam_file
        idx_cmd = "samtools index " + bam_file
        Popen(get_mmi_cmd,shell=True).wait()
        Popen(align_cmd,shell=True).wait()
        Popen(sort_cmd,shell=True).wait()
        Popen(idx_cmd,shell=True).wait()
    except:
        get_mmi_cmd = code_path + "minimap2 -d " + ref_chr_mmi + " " +  ref_chr_fasta
        align_cmd = code_path + "minimap2 -a " + ref_chr_mmi + " " + fasta_file + " > " + sam_file
        sort_cmd = code_path + "samtools view -Sb " + sam_file + " | " + code_path + "samtools sort  > " + bam_file
        idx_cmd = code_path + "samtools index " + bam_file
        Popen(get_mmi_cmd,shell=True).wait()
        Popen(align_cmd,shell=True).wait()
        Popen(sort_cmd,shell=True).wait()
        Popen(idx_cmd,shell=True).wait()


def Align_hp_fasta_to_bam(in_dir,out_dir,chr_num):
    for sample_name in sample_list:
        out_dir_sample = out_dir + sample_name + "/"
        if os.path.exists(out_dir_sample):
            print("using existing output folder: " + out_dir)
        else:
            os.makedirs(out_dir_sample)
        hp1_file = in_dir + sample_name + "/Assembly_Contigs_files/" +  "Aquila_Contig_chr" + str(chr_num) + "_hp1.fasta "
        sam1_file = out_dir_sample + sample_name + "_chr" + str(chr_num) + "_contig.1.sam"
        bam1_file = out_dir_sample + sample_name + "_chr" + str(chr_num) + "_contig.1.bam"
        hp2_file = in_dir + sample_name + "/Assembly_Contigs_files/" +  "Aquila_Contig_chr" + str(chr_num) + "_hp2.fasta "

        sam2_file = out_dir_sample + sample_name + "_chr" + str(chr_num) + "_contig.2.sam"
        bam2_file = out_dir_sample + sample_name + "_chr" + str(chr_num) + "_contig.2.bam"

        #pool = Pool(processes=2)
        #pool.apply_async(get_sam,(hp1_file,sam1_file,bam1_file,chr_num,ref_dir,"xin"))
        get_sam(hp1_file,sam1_file,bam1_file,chr_num,ref_dir,"xin")
        #pool.apply_async(get_sam,(hp2_file,sam2_file,bam2_file,chr_num,ref_dir,"xin"))
        get_sam(hp2_file,sam2_file,bam2_file,chr_num,ref_dir,"xin")
        """
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()
        """

def Delete_temp_files(out_dir_chr):
    _cmd = "rm " + out_dir_chr + "target.fa " + out_dir_chr + "target.mmi " + out_dir_chr + "query.fa " + out_dir_chr + "output.sam " + out_dir_chr + "output.bam*"
    Popen(_cmd,shell=True).wait()


## target_seq: contig ref_seq
def Align_refseq_to_contig(ref_seq,target_seq,out_dir_chr):
    target_fa = open(out_dir_chr + "target.fa","w")
    target_fa.writelines(">contig\n")
    target_fa.writelines(target_seq + "\n")
    target_fa.close()
    query_fa = open(out_dir_chr + "query.fa","w")
    query_fa.writelines(">refseq\n")
    query_fa.writelines(ref_seq + "\n")
    query_fa.close()

    _cmd = "minimap2 -d " + out_dir_chr + "target.mmi " + " " + out_dir_chr + "target.fa "
    Popen(_cmd,shell=True).wait()
    _cmd = "minimap2 -a " + out_dir_chr + "target.mmi " + " " + out_dir_chr + "query.fa " + " > " + out_dir_chr + "output.sam"
    Popen(_cmd,shell=True).wait()
    sam_to_bam_cmd = "samtools view -Sb " + out_dir_chr + "output.sam " + " | samtools sort  > " + out_dir_chr + "output.bam"
    Popen(sam_to_bam_cmd,shell=True).wait()
    index_cmd = "samtools index "  + out_dir_chr + "output.bam"
    Popen(index_cmd,shell=True).wait()


def Cal_extract_steps(cigar_temp):
    total_use = 0
    match_num = 0
    del_num = 0
    for num in range(len(cigar_temp)):
        one_cigar = cigar_temp[num]
        if one_cigar[0] == 0:
            match_num = one_cigar[1]
            total_use += match_num
        elif one_cigar[0] == 2:
            del_num = one_cigar[1]
            total_use += del_num
    return total_use


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


def extract_contig_fasta(fasta_file,lib_prefix,out_dir):
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    f = open(fasta_file,"r")
    contig_dict = defaultdict()
    count = 0
    for line in f:
        data = line.rsplit()
        if count%2 == 0:
            contig_num = data[0].split(">")[-1]
            contig_num_new = contig_num.replace(":","_")
        if count%2 == 1:
            try:
                contig_str = data[0]
            except:
                pass
            contig_dict[contig_num_new] =  contig_str
        count += 1
    pickle.dump(contig_dict, open(out_dir + "contig_dict_" + lib_prefix + ".p","wb"))
    return contig_dict


def Extract_fasta_for_SV(merge_SV_dict,out_dir,out_dir_chr,ref_dir,in_chr_num):
    ref_dict = pickle.load(open(ref_dir + "ref_seq_chr" + str(in_chr_num) + ".p","rb"))
    count_contig_0 = 0
    _count_miss = 0
    use_chr_num = "chr" + str(in_chr_num)
    for key,value in merge_SV_dict.items():
        if len(value) > 1:
            _start = key[0]
            _end = key[1]
            sv_type = key[2]
            fw = open(out_dir_chr + str(use_chr_num) + "_" + str(_start) + "_" + str(_end) + "_" + sv_type +  ".fasta","w")
            ref_seq = ref_dict[_start-500:_end+500]
            fw.writelines(">human_ref" + "\n")
            fw.writelines(ref_seq + "\n")
            for one_value in value:
                contig_info = one_value[-1]
                GT = one_value[6]
                if GT == "0/1":
                    if contig_info[0] == "0":
                        contig_1_num = 0
                        contig_2 = contig_info.split("&")[1]
                        contig_2_info = "_".join(contig_2.split("_")[:-3])
                        contig_2_num = int(contig_2_info.split("_")[0])
                        contig_2_start = int(contig_2.split("_")[-3])
                        contig_2_end = int(contig_2.split("_")[-2])
                        contig_2_strand = contig_2.split("_")[-1]
                    else:
                        contig_2_num = 0
                        contig_1 = contig_info.split("&")[0]
                        contig_1_info = "_".join(contig_1.split("_")[:-3])
                        contig_1_num = int(contig_1_info.split("_")[0])
                        contig_1_start = int(contig_1.split("_")[-3])
                        contig_1_end = int(contig_1.split("_")[-2])
                        contig_1_strand = contig_1.split("_")[-1]
                elif GT == "1/1":
                    contig_1 = contig_info.split("&")[0]
                    contig_1_info = "_".join(contig_1.split("_")[:-3])
                    contig_1_num = int(contig_1_info.split("_")[0])
                    contig_1_start = int(contig_1.split("_")[-3])
                    contig_1_end = int(contig_1.split("_")[-2])
                    contig_1_strand = contig_1.split("_")[-1]

                    contig_2 = contig_info.split("&")[1]
                    contig_2_info = "_".join(contig_2.split("_")[:-3])
                    contig_2_num = int(contig_2_info.split("_")[0])
                    contig_2_start = int(contig_2.split("_")[-3])
                    contig_2_end = int(contig_2.split("_")[-2])
                    contig_2_strand = contig_2.split("_")[-1]

                sample_name = one_value[0]
                sample_name_2 = sample_name.split("_")[0]
                event_name = one_value[5]
                chr_num = one_value[1]
                SV_size = one_value[4]
                ref_start = int(one_value[2])
                ref_end = int(one_value[3])
                _contig_1_dict = pickle.load(open(out_dir + sample_name_2 + "/" + "contig_dict_" + sample_name_2 + "_1_" + str(chr_num) + ".p" ,"rb"))
                _contig_2_dict = pickle.load(open(out_dir + sample_name_2 + "/" + "contig_dict_" + sample_name_2 + "_2_" + str(chr_num) + ".p","rb"))
                if contig_1_num != 0:
                    if contig_1_start >= 500 and (contig_1_end + 500) <= len(_contig_1_dict[contig_1_info]):
                        contig_1_seq = _contig_1_dict[contig_1_info][contig_1_start-500:contig_1_end+500]
                        fw.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c1_" + GT + "_1_" + SV_size + "\n")
                        if contig_1_strand == "+":
                            fw.writelines(contig_1_seq + "\n")
                        elif contig_1_strand == "-":
                            fw.writelines(translate(contig_1_seq) + "\n")
                elif contig_1_num == 0:
                    count_contig_0 += 1
                    print(contig_1_num)
                    ref_seq_1 = ref_dict[ref_start-500:ref_end+500]
                    bam_file = out_dir + sample_name_2 + "/" +  sample_name_2 + "_" + chr_num + "_contig.1.bam"
                    sam_file_temp = out_dir_chr + "temp.sam"
                    bam_file_temp = out_dir_chr + "temp.bam"
                    regions = "\'" + chr_num[3:] + ":" + str(ref_start) + "-" + str(ref_end) + "\'"
                    extract_cmd = "samtools view -h " + bam_file + " " + regions + " >  " + sam_file_temp
                    Popen(extract_cmd,shell=True).wait()
                    sam_to_bam_cmd = "samtools view -Sb " + sam_file_temp + " | samtools sort  > " + bam_file_temp
                    Popen(sam_to_bam_cmd,shell=True).wait()
                    index_cmd = "samtools index "  + bam_file_temp
                    Popen(index_cmd,shell=True).wait()
                    sam_file_curr = pysam.AlignmentFile(bam_file_temp, "rb")
                    for read in sam_file_curr.fetch(chr_num[3:]):
                        mapq = int(read.mapping_quality)
                        if mapq >= 20:
                            read_seq = read.seq
                            if read_seq == "*":
                                read_seq = _contig_1_dict[contig_1_info]
                            align_start = int(read.pos)
                            cigar = read.cigarstring
                            if read.is_reverse:
                                target_seq = translate(read_seq)
                                if len(target_seq) > len(ref_seq_1)*1.5:
                                    Align_refseq_to_contig(ref_seq_1,target_seq,out_dir_chr)
                                    sam_out = pysam.AlignmentFile(out_dir_chr + "output.bam", "rb")
                                    for read_temp in sam_out.fetch():
                                        mapq_temp = int(read_temp.mapping_quality)
                                        if mapq_temp >= 20:
                                            align_start_temp = int(read_temp.pos)
                                            cigar_temp = read_temp.cigar
                                            total_use = Cal_extract_steps(cigar_temp)
                                            #### contig sequence extract
                                            if total_use >= 0.9*1000:
                                                target_seq_extract = target_seq[align_start_temp:align_start_temp+total_use]
                                                fw.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c1_" + GT + "_0_"  + SV_size + "\n")
                                                fw.writelines(translate(target_seq_extract) + "\n")
                                                Delete_temp_files(out_dir_chr)
                                                break
                                        else:
                                            _count_miss += 1
                                    sam_out.close()
                            else:
                                target_seq = read_seq
                                if len(target_seq) > len(ref_seq_1)*1.5:
                                    Align_refseq_to_contig(ref_seq_1,target_seq,out_dir_chr)
                                    sam_out = pysam.AlignmentFile(out_dir_chr + "output.bam", "rb")
                                    for read_temp in sam_out.fetch():
                                        mapq_temp = int(read_temp.mapping_quality)
                                        if mapq_temp >= 20:
                                            align_start_temp = int(read_temp.pos)
                                            cigar_temp = read_temp.cigar
                                            total_use = Cal_extract_steps(cigar_temp)
                                            if total_use >= 0.9*1000:
                                                target_seq_extract = target_seq[align_start_temp:align_start_temp+total_use]
                                                fw.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c1_" + GT + "_0_"  + SV_size + "\n")
                                                fw.writelines(target_seq_extract + "\n")
                                                Delete_temp_files(out_dir_chr)
                                                break
                                        else:
                                            _count_miss += 1
                                    sam_out.close()
                        else:
                            _count_miss += 1

                    sam_file_curr.close()                                

                if contig_2_num != 0:
                    if contig_2_start >= 500 and (contig_2_end + 500) <= len(_contig_2_dict[contig_2_info]):
                        contig_2_seq = _contig_2_dict[contig_2_info][contig_2_start-500:contig_2_end+500]
                        fw.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c2_" + GT + "_1_" + SV_size + "\n")
                        if contig_2_strand == "+":
                            fw.writelines(contig_2_seq + "\n")
                        elif contig_2_strand == "-":
                            fw.writelines(translate(contig_2_seq) + "\n")
                elif contig_2_num == 0:
                    count_contig_0 += 1
                    ref_seq_2 = ref_dict[ref_start-500:ref_end+500]
                    print(contig_2_num)
                    bam_file = out_dir + sample_name_2 + "/" + sample_name_2 + "_" + chr_num + "_contig.2.bam"
                    sam_file_temp = out_dir_chr +  "temp.sam"
                    bam_file_temp = out_dir_chr + "temp.bam"
                    regions = "\'" + chr_num[3:] + ":" + str(ref_start) + "-" + str(ref_end) + "\'"
                    extract_cmd = "samtools view -h " + bam_file + " " + regions + " >  " + sam_file_temp
                    Popen(extract_cmd,shell=True).wait()
                    sam_to_bam_cmd = "samtools view -Sb " + sam_file_temp + " | samtools sort  > " + bam_file_temp
                    Popen(sam_to_bam_cmd,shell=True).wait()
                    index_cmd = "samtools index "  + bam_file_temp
                    Popen(index_cmd,shell=True).wait()
                    sam_file_curr = pysam.AlignmentFile(bam_file_temp, "rb")
                    for read in sam_file_curr.fetch(chr_num[3:]):
                        mapq = int(read.mapping_quality)
                        if mapq >= 20:
                            read_seq = read.seq
                            if read_seq == "*":
                                read_seq = _contig_2_dict[contig_2_num]
                            align_start = int(read.pos)
                            cigar = read.cigarstring
                            if read.is_reverse:
                                target_seq = translate(read_seq)
                                if len(target_seq) > len(ref_seq_2)*1.5:
                                    Align_refseq_to_contig(ref_seq_2,target_seq,out_dir_chr)
                                    sam_out = pysam.AlignmentFile(out_dir_chr + "output.bam", "rb")
                                    for read_temp in sam_out.fetch():
                                        mapq_temp = int(read_temp.mapping_quality)
                                        if mapq_temp >= 20:
                                            align_start_temp = int(read_temp.pos)
                                            cigar_temp = read_temp.cigar
                                            total_use = Cal_extract_steps(cigar_temp)
                                            if total_use >= 0.9*1000:
                                                target_seq_extract = target_seq[align_start_temp:align_start_temp+total_use]
                                                fw.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c2_" + GT + "_0_"  + SV_size + "\n")
                                                fw.writelines(translate(target_seq_extract) + "\n")
                                                Delete_temp_files(out_dir_chr)
                                                break
                                        else:
                                            _count_miss += 1
                                    sam_out.close()
                            else:
                                target_seq = read_seq
                                if len(target_seq) > len(ref_seq_2)*1.5:
                                    Align_refseq_to_contig(ref_seq_2,target_seq,out_dir_chr)
                                    sam_out = pysam.AlignmentFile(out_dir_chr + "output.bam", "rb")
                                    for read_temp in sam_out.fetch():
                                        mapq_temp = int(read_temp.mapping_quality)
                                        if mapq_temp >= 20:
                                            align_start_temp = int(read_temp.pos)
                                            cigar_temp = read_temp.cigar
                                            total_use = Cal_extract_steps(cigar_temp)
                                            if total_use >= 0.9*1000:
                                                target_seq_extract = target_seq[align_start_temp:align_start_temp+total_use]
                                                fw.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c2_" + GT + "_0_"  + SV_size + "\n")
                                                fw.writelines(target_seq_extract + "\n")
                                                Delete_temp_files(out_dir_chr)
                                                break
                                        else:
                                            _count_miss += 1
                                    sam_out.close()
                        else:
                            _count_miss += 1

                    sam_file_curr.close()                

            fw.close()
        
    print("count_contig_0: " + str(count_contig_0))
    print("_count_miss: " + str(_count_miss))


def Multiseq_alignment_by_muscle(in_dir):
    files_all = glob.glob(in_dir + "/" + "*.fasta")
    for one_file in files_all:
        output_file = one_file.split(".fasta")[0] + ".html"
        use_cmd = code_path + "muscle3.8.31_i86linux64 -in " + one_file + " -out " + output_file + " -html" + " -maxiters 1 -diags -sv -distance1 kbit20_3"
        Popen(use_cmd,shell=True).wait()

        output_file = one_file.split(".fasta")[0] + ".afa"
        use_cmd = code_path + "muscle3.8.31_i86linux64 -in " + one_file + " -out " + output_file + " -maxiters 1 -diags -sv -distance1 kbit20_3"
        Popen(use_cmd,shell=True).wait()
        
    print("msa by muscle done~")


def Multiseq_profile_alignment_by_muscle(in_dir,in_dir_2):
    files_all = glob.glob(in_dir + "/" + "*.fasta")
    for one_file in files_all:
        output_file = one_file.split(".fasta")[0] + ".afa"
        input_file_2 = in_dir_2 + "/" + output_file.split("/")[-1]
        use_cmd = code_path + "muscle3.8.31_i86linux64 -in " + one_file + " -out " + output_file  + " -maxiters 1 -diags -sv -distance1 kbit20_3"
        Popen(use_cmd,shell=True).wait()

        output_file = one_file.split(".fasta")[0] + ".html"
        use_cmd = code_path + "muscle3.8.31_i86linux64 -profile -in1 " + output_file + "  -in2 " + input_file_2 +  " -out " + output_file + " -html " + " -maxiters 1 -diags -sv -distance1 kbit20_3"
        Popen(use_cmd,shell=True).wait()
        
    print("msa by muscle done~")


def run_one_msa_by_muscle(one_file,output_html,output_afa,xin):
    use_cmd = code_path + "muscle3.8.31_i86linux64 -in " + one_file + " -out " + output_html + " -html" + " -maxiters 1 -diags -sv -distance1 kbit20_3"
    Popen(use_cmd,shell=True).wait()
    use_cmd = code_path  + "muscle3.8.31_i86linux64 -in " + one_file + " -out " + output_afa + " -maxiters 1 -diags -sv -distance1 kbit20_3"
    Popen(use_cmd,shell=True).wait()


def run_one_profile_msa_by_muscle(one_file,output_html,output_afa_1,input_afa_2,output_afa,xin):
    use_cmd = code_path + "muscle3.8.31_i86linux64 -in " + one_file + " -out " + output_afa_1 + " -maxiters 1 -diags -sv -distance1 kbit20_3"
    Popen(use_cmd,shell=True).wait()
    use_cmd = code_path + "muscle3.8.31_i86linux64 -profile -in1 " + input_afa_2 + " -in2 " + output_afa_1  + " -out " + output_html + " -html " + " -maxiters 1 -diags -sv -distance1 kbit20_3"
    Popen(use_cmd,shell=True).wait()
    use_cmd = code_path + "muscle3.8.31_i86linux64 -profile -in1 " + input_afa_2 + " -in2 " + output_afa_1  + " -out " + output_afa  + " -maxiters 1 -diags -sv -distance1 kbit20_3"
    Popen(use_cmd,shell=True).wait()


def Multiseq_alignment_by_muscle_Multithreads(in_dir,num_of_threads):
    files_all = sorted(glob.glob(in_dir + "/" +  "*.fasta"),key=os.path.getsize,reverse=True)
    count = 1
    total_num = len(files_all)
    #pool = Pool(processes=num_of_threads)
    for one_file in files_all:
        output_html = one_file.split(".fasta")[0] + ".html"
        output_afa = one_file.split(".fasta")[0] + ".afa"
        #pool.apply_async(run_one_msa_by_muscle,(one_file,output_html,output_afa,"xin"))
        run_one_msa_by_muscle(one_file,output_html,output_afa,"xin")
        """
        count += 1
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()
            if (count - 1) == total_num:
                print("finished")
            else:
                pool = Pool(processes=num_of_threads)
        """
    print("MSA to generate *html/*afa all done~")


def Multiseq_profile_alignment_by_muscle_Multithreads(in_dir,in_dir_2,num_of_threads):
    files_all = sorted(glob.glob(in_dir + "/" +  "*.fasta"),key=os.path.getsize,reverse=True)
    count = 1
    total_num = len(files_all)
    #pool = Pool(processes=num_of_threads)
    for one_file in files_all:
        output_html = one_file.split(".fasta")[0] + ".html"
        output_afa_1 = one_file.split(".fasta")[0] + "_1.afa"
        output_afa = one_file.split(".fasta")[0] + ".afa"
        input_afa_2 = in_dir_2 + "/" + output_afa.split("/")[-1]
        #pool.apply_async(run_one_profile_msa_by_muscle,(one_file,output_html,output_afa_1,input_afa_2,output_afa,"xin"))
        run_one_profile_msa_by_muscle(one_file,output_html,output_afa_1,input_afa_2,output_afa,"xin")
        """
        count += 1
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()
            if (count - 1) == total_num:
                print("finished")
            else:
                pool = Pool(processes=num_of_threads)
        """
    print("MSA to generate *html/*afa all done~")
        

def check_folder_existence(out_dir_derived_sv):
    if os.path.exists(out_dir_derived_sv):
        print("Delete existing output folder: " + out_dir_derived_sv)
        Popen("rm -rf " + out_dir_derived_sv,shell=True).wait()
        os.makedirs(out_dir_derived_sv)
    else:
        os.makedirs(out_dir_derived_sv)


def gnomad_snp_extract(gnomad_dir,vcf_bgz_file,use_chr_num):
    gnomad_snp_dict = defaultdict(list)
    with gzip.open(vcf_bgz_file,"r") as f:
        for line in f:
            data = line.decode().rsplit()
            if data[0][0] != "#":
                chr_num = data[0]
                _pos = int(data[1])
                rsid = data[2]
                ref = data[3]
                alt = data[4]
                gnomad_snp_dict[(chr_num,_pos)] = [rsid,ref,alt]
    pickle.dump(gnomad_snp_dict,open(gnomad_dir + "gnomad_snp_dict_chr" + str(use_chr_num) + ".p","wb"))
    return gnomad_snp_dict




def Run_step2_by_chr(out_dir,out_dir_all,out_dir_table,in_dir,assembly_dir,ref_dir,ref_file,chr_start,chr_end,gnomad_flag,HARP_flag,num_of_threads,xin):
    ### step 1: Merge SV for population
    Extract_DEL_from_all_samples(sample_list,in_dir,out_dir_all,chr_start,chr_end)
    Extract_INS_from_all_samples(sample_list,in_dir,out_dir_all,chr_start,chr_end)
    merge_DEL_dict = Merge_DEL(out_dir_all,chr_start,chr_end)
    merge_INS_dict = Merge_INS(out_dir_all,chr_start,chr_end)
    merge_SV_dict = Merge_DEL_and_INS(out_dir_all,chr_start,chr_end)
    for chr_num in range(chr_start,chr_end + 1):
        if gnomad_flag == 1:
            gnomad_dir = args.gnomad_dir + "/"
            gnomad_vcf_bgz_file = gnomad_dir  + "gnomad.genomes.r2.1.1.sites." + str(chr_num) + ".liftover_grch38.vcf.bgz"
            gnomad_snp_extract(gnomad_dir,gnomad_vcf_bgz_file,chr_num)
        out_dir_chr = out_dir + "/" + "MSA_SV_results_chr" + str(chr_num) + "/MSA_SV_files/"
        out_dir_derived_del = out_dir + "/" + "MSA_SV_results_chr" + str(chr_num) + "/HARP_derived_del/"
        out_dir_derived_ins = out_dir + "/" + "MSA_SV_results_chr" + str(chr_num) + "/HARP_derived_ins/"
        out_dir_derived_complex = out_dir + "/" + "MSA_SV_results_chr" + str(chr_num) + "/HARP_derived_complex/"
        check_folder_existence(out_dir_derived_del)
        check_folder_existence(out_dir_derived_ins)
        check_folder_existence(out_dir_derived_complex)
        if os.path.exists(out_dir_chr):
            print("using existing output folder: " + out_dir_chr)
        else:
            os.makedirs(out_dir_chr)
        #############################################################################
        ###############       load human ref for each chromosome ####################
        ############### load contig sequence files for each chromosome ##############
        #############################################################################
        if chr_num == 23:
            extract_ref_chr(ref_file,"X",ref_dir)
            read_ref(ref_dir + "genome_ref_chrX.fasta",chr_num,ref_dir)
        else:
            extract_ref_chr(ref_file,chr_num,ref_dir)
            read_ref(ref_dir + "genome_ref_chr" + str(chr_num) + ".fasta",chr_num,ref_dir)
        Align_hp_fasta_to_bam(assembly_dir,out_dir,chr_num)
        for sample_name in sample_list:
            hp1_file = assembly_dir + sample_name + "/" + "Assembly_Contigs_files/" + "Aquila_Contig_chr" + str(chr_num) + "_hp1.fasta"  
            hp2_file = assembly_dir + sample_name + "/" + "Assembly_Contigs_files/" + "Aquila_Contig_chr" + str(chr_num) + "_hp2.fasta"  
            contig_dict_1 = extract_contig_fasta(hp1_file,sample_name + "_1_chr" + str(chr_num),out_dir + sample_name + "/")
            contig_dict_2 = extract_contig_fasta(hp2_file,sample_name + "_2_chr" + str(chr_num),out_dir + sample_name + "/")
        ############################################################################
        merge_SV_dict = pickle.load(open(out_dir_all + "merge_SV_dict_chr" + str(chr_num) + ".p", "rb"))
        ### step 2: Extract SV sequences from contigs
        Extract_fasta_for_SV(merge_SV_dict,out_dir,out_dir_chr,ref_dir,chr_num)
        
        ### step 3: MSA by muscle
        Multiseq_alignment_by_muscle_Multithreads(out_dir_chr,num_of_threads)

        ### step 4: Define Break Point and only leave flanking region 10bp from left and right sides
        Define_SV_break_point(out_dir_chr)
        ### step 5: Calculate SV repeat percentage from trf, annotate Alu, LINEs
        SV_len_dict = Cal_SV_repeat_percentage_by_trf(out_dir_chr,out_dir_table,chr_num)
        ### step 6: Generate .txt/.excel table for SV
        output_file = out_dir_table + "SV_msa_table_chr" + str(chr_num) + ".txt"
        SV_len_dict = pickle.load(open(out_dir_chr + "SV_len_dict_chr" + str(chr_num) + ".p","rb"))
        SV_bk_dict = pickle.load(open(out_dir_chr + "SVs_breakpoints_dict.p","rb"))
        Generate_table(out_dir_chr,output_file,SV_len_dict,SV_bk_dict,1)
        if gnomad_flag == 1:
            snp_dict_file = out_dir_chr  + "/SNP_dict.p"
            snp_dict = pickle.load(open(snp_dict_file,"rb"))
            gnomad_dict_file = gnomad_dir  + "/gnomad_snp_dict_chr" + str(chr_num) + ".p"
            gnomad_dict = pickle.load(open(gnomad_dict_file,"rb"))
            Print_linked_SNPs_info(out_dir_chr,out_dir_table,gnomad_dict,chr_num,num_of_samples)
            
        #####################################################################################################
        #########################Hominid Ancestral Population analysis (HARP)################################
        #####################################################################################################
        if HARP_flag == 1:
            ### step 7: Extract flanking region for SVs to align to Apes
            Extract_flankingregion_for_SV_to_align_to_Apes(out_dir_all,ref_dir,Ape_ref_list,chr_num)

            ### step 8: Evaluate derived SVs by global alignment with Apes
            Evalute_derived_SV_by_flankingregion_global_align(out_dir_all,Ape_ref_list,chr_num)

            ### step 9: Extract derived SVs with Apes sequences
            Extract_derived_SV_with_Apes_seq(Ape_ref_list,out_dir,out_dir_all,chr_num)
           
            ### step 10: MSA by muscle
            ### step 11: Define Break Point and only leave flanking region 10bp from left and right sides
            ### step 12: Generate .txt/.excel table for derived del
            Multiseq_profile_alignment_by_muscle_Multithreads(out_dir_derived_del,out_dir_chr,num_of_threads)
            Define_SV_break_point(out_dir_derived_del)
            
            output_file = out_dir_table + "derived_del_msa_table_chr" + str(chr_num) + ".txt"
            Generate_table(out_dir_derived_del,output_file,SV_len_dict,SV_bk_dict,0)
            
            ### step 10: MSA by muscle
            ### step 11: Define Break Point and only leave flanking region 10bp from left and right sides
            ### step 12: Generate .txt/.excel table for derived ins
            Multiseq_profile_alignment_by_muscle_Multithreads(out_dir_derived_ins,out_dir_chr,num_of_threads)
            Define_SV_break_point(out_dir_derived_ins)

            output_file = out_dir_table + "derived_ins_msa_table_chr" + str(chr_num) + ".txt"
            Generate_table(out_dir_derived_ins,output_file,SV_len_dict,SV_bk_dict,0)
            
            ### step 10: MSA by muscle
            ### step 11: Define Break Point and only leave flanking region 10bp from left and right sides
            ### step 12: Generate .txt/.excel table for derived complex sv
            Multiseq_profile_alignment_by_muscle_Multithreads(out_dir_derived_complex,out_dir_chr,num_of_threads)
            Define_SV_break_point(out_dir_derived_complex)

            output_file = out_dir_table + "derived_complex_msa_table_chr" + str(chr_num) + ".txt"
            Generate_table(out_dir_derived_complex,output_file,SV_len_dict,SV_bk_dict,0)


def Run_MARS_step2_all(out_dir,out_dir_all,out_dir_table,in_dir,assembly_dir,ref_dir,ref_file,chr_start,chr_end,gnomad_flag,HARP_flag,num_of_threads,num_of_threads_bychr):
    count = 1
    total_num = chr_end - chr_start + 1
    pool = Pool(num_of_threads_bychr)
    for chr_num in range(chr_start,chr_end+1):
        count += 1
        pool.apply_async(Run_step2_by_chr,(out_dir,out_dir_all,out_dir_table,in_dir,assembly_dir,ref_dir,ref_file,chr_num,chr_num,gnomad_flag,HARP_flag,num_of_threads,"xin"))
        Run_step2_by_chr(out_dir,out_dir_all,out_dir_table,in_dir,assembly_dir,ref_dir,ref_file,chr_num,chr_num,gnomad_flag,HARP_flag,num_of_threads,"xin")
        if (count - 1)%num_of_threads_bychr == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()

            if (count - 1) == total_num:
                print("finished all!")
            else:
                pool = Pool(num_of_threads_bychr)
    print("all done~")




if __name__ == "__main__":
    out_dir = args.out_dir + "/"
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    in_dir = args.in_dir + "/" 
    assembly_dir = args.assembly_dir + "/" 
    ref_file = args.ref_file 
    chr_start = args.chr_start
    chr_end = args.chr_end
    gnomad_flag = int(args.gnomad_flag_for_linked_SNP)
    HARP_flag = int(args.HARP_flag)
    ref_dir = out_dir + "ref_dir" + "/"
    #num_of_threads = int(args.num_threads_msa)
    num_of_threads = 1
    num_of_threads_bychr = int(args.num_threads_bychr)
    if os.path.exists(ref_dir):
        print("using existing output folder: " + ref_dir)
    else:
        os.makedirs(ref_dir)
    out_dir_all = out_dir + "/" + "SV_files_all/"
    out_dir_table = out_dir + "/" + "Final_tables/"
    if os.path.exists(out_dir_all):
        print("using existing output folder: " + out_dir_all)
    else:
        os.makedirs(out_dir_all)
    if os.path.exists(out_dir_table):
        print("using existing output folder: " + out_dir_table)
    else:
        os.makedirs(out_dir_table)
    Ape_ref_dir = out_dir + "Apes_ref/" 
    if os.path.exists(Ape_ref_dir):
        print("using existing output folder: " + Ape_ref_dir)
    else:
        os.makedirs(Ape_ref_dir)
    if HARP_flag == 1:
        for Ape_ref_file in Ape_ref_list:
            Ape_name = Ape_ref_file.split("/")[-1].split(".")[0]
            load_and_save_Apes_ref(Ape_ref_dir,Ape_ref_file,Ape_name)

    Run_MARS_step2_all(out_dir,out_dir_all,out_dir_table,in_dir,assembly_dir,ref_dir,ref_file,chr_start,chr_end,gnomad_flag,HARP_flag,num_of_threads,num_of_threads_bychr)
            
            
    for chr_num in range(chr_start,chr_end + 1):        
        print("processing chr" + str(chr_num))
        input_table = out_dir_table + "SV_msa_table_chr" + str(chr_num) + ".txt"
        out_table = out_dir_table + "SV_msa_table_chr" + str(chr_num) + "_2.txt"
        out_dir_chr = out_dir + "/" + "MSA_SV_results_chr" + str(chr_num) + "/MSA_SV_files/"
        updatefile(input_table,out_table,out_dir_chr)

    if HARP_flag == 1:
        for chr_num in range(chr_start,chr_end + 1):
            print("processing chr" + str(chr_num))
            input_table = out_dir_table + "derived_del_msa_table_chr" + str(chr_num) + ".txt"
            derived_del_dict = read_table(input_table)
            input_table = out_dir_table + "derived_ins_msa_table_chr" + str(chr_num) + ".txt"
            derived_ins_dict = read_table(input_table)
            input_table = out_dir_table + "derived_complex_msa_table_chr" + str(chr_num) + ".txt"
            derived_complex_dict = read_table(input_table)

            input_table = out_dir_table + "SV_msa_table_chr" + str(chr_num) + "_2.txt"
            output_table = out_dir_table + "SV_msa_table_chr" + str(chr_num) + "_add_AncestralState.txt"
            Merge_table(input_table,output_table,derived_del_dict,derived_ins_dict,derived_complex_dict)
