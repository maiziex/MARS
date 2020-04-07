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
parser = ArgumentParser(description="Run MSA:")
parser.add_argument('--in_dir','-i_dir', help="The input folder where you store vcf files",required=True)
parser.add_argument('--out_dir','-o_dir', help="The folder name you can define to store the final results",default="./Results_MSA_Aquila")
parser.add_argument('--assembly_dir','-a_dir', help="The folder to store assembly contigs",required=True)
parser.add_argument('--ref_file','-r', help="The human reference fasta file which can be download by running \"./install.sh\". ",required=True)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
#parser.add_argument('--num_threads','-nt', help=" The number of threads you can define, which corresponds to number of samples",default=2)
parser.add_argument('--sample_list','-sl', help='The sample names corresponding to your contig files, which is the prefix of the contig files', type=str,required=True)
args = parser.parse_args()
if len(sys.argv) == 1:
    Popen("python Run_all_samples_Aquila.py -h",shell=True).wait()
else:
    sample_list = [item for item in args.sample_list.split(',')]


#sample_list = ['HG00250','HG00353','HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','HG00851','HG01971','HG02623','HG03115','HG03838','NA12878','NA18552','NA19068','NA19238','NA19239','NA19240','NA19440','NA19789','NA20587','NA24143','NA24149','NA24385','hgp','HLA1','HLA2','HLA3','HLA4','HLA5','HLA7','HLA10']

basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
def translate(seq):
    seq_r = []
    for character in seq:
        seq_r.append(basecomplement[character])
    seq_r_2 = "".join(seq_r)
    rseqn = seq_r_2[::-1]
    return rseqn


def Extract_SV_from_all_samples(sample_list,in_dir,out_dir,chr_start,chr_end):
    for chr_num_use in range(chr_start,chr_end + 1):
        print(chr_num_use)
        output_file = out_dir + "chr" + str(chr_num_use) + "_ins.txt"
        output_file_sorted = out_dir + "chr" + str(chr_num_use) + "_ins_sorted.txt"
        if chr_num_use == 23:
            chr_num_use = "X"
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
                            fw.writelines(sample_name + "\t" + chr_num + "\t" + str(start_) + "\t" +  str(end_) + "\t" + str(SV_len) + "\t" +   event_name + "\t" + GT + "\t" + contig_info + "\n")
        fw.close()

        sort_cmd = "cat " + output_file + " | sort -k3n > " + output_file_sorted 
        Popen(sort_cmd,shell=True).wait()
    print("done")


def Write_to_single_SV(save_data,fw,min_start,max_end):
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
    fw.writelines(chr_num + "\t" + str(min_start) + "\t" + str(max_end) + "\t" + ":".join(sample_num_list) + "\t" + ":".join(event_name_list) +"\t" + ":".join(SV_size_list) + "\t" +  ":".join(contig_info_list) + "\t" + ":".join(GT_list) + "\n")
    return fw


def Merge_SV(out_dir,chr_start,chr_end):
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
                    fw = Write_to_single_SV(save_data,fw,min_start,max_end)
                    merge_SV_dict[(min_start,max_end)] = save_data
                    save_start_list = []
                    save_data = []
                    save_start_list.append(_start)
                    save_data.append(data)
                    min_start = _start
                count += 1
        max_end = save_start_list[-1]
        fw = Write_to_single_SV(save_data,fw,min_start,max_end)
        merge_SV_dict[(min_start,max_end)] = save_data
        pickle.dump(merge_SV_dict, open(out_dir + "merge_INS_dict_chr" + str(chr_num_use)  + ".p", "wb"))
    return merge_SV_dict


def get_sam(fasta_file,sam_file,bam_file,chr_num,xin):
    align_cmd = "minimap2 -a /oak/stanford/groups/arend/Xin/AssemblyProj/reference_align_2/reference_files/ref_chr" + str(chr_num) + ".mmi " + fasta_file + " > " + sam_file
    sort_cmd = "samtools view -Sb " + sam_file + " | samtools sort  > " + bam_file
    idx_cmd = "samtools index " + bam_file
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

        pool = Pool(processes=2)
        pool.apply_async(get_sam,(hp1_file,sam1_file,bam1_file,chr_num,"xin"))
        pool.apply_async(get_sam,(hp2_file,sam2_file,bam2_file,chr_num,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()


def Delete_temp_files():
    _cmd = "rm target.fa target.mmi query.fa output.sam output.bam*"
    Popen(_cmd,shell=True).wait()


## target_seq: contig ref_seq
def Align_refseq_to_contig(ref_seq,target_seq):
    target_fa = open("target.fa","w")
    target_fa.writelines(">contig\n")
    target_fa.writelines(target_seq + "\n")
    target_fa.close()
    query_fa = open("query.fa","w")
    query_fa.writelines(">refseq\n")
    query_fa.writelines(ref_seq + "\n")
    query_fa.close()

    _cmd = "minimap2 -d " + " target.mmi " + " target.fa "
    Popen(_cmd,shell=True).wait()
    _cmd = "minimap2 -a target.mmi " + " query.fa " + " > " + " output.sam"
    Popen(_cmd,shell=True).wait()
    sam_to_bam_cmd = "samtools view -Sb " + "output.sam " + " | samtools sort  > " + "output.bam"
    Popen(sam_to_bam_cmd,shell=True).wait()
    index_cmd = "samtools index "  + "output.bam"
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


def Extract_fasta_for_SV(merge_SV_dict,out_dir,out_dir_chr,chr_num):
    ref_dict = pickle.load(open(out_dir_chr + "ref_seq_chr" + str(chr_num) + ".p","rb"))
    count_contig_0 = 0
    _count_miss = 0
    for key,value in merge_SV_dict.items():
        if len(value) > 1:
            _start = key[0]
            if _start == 5105385:  #5077846:  #5098515
                print(_start)
            _end = key[1]
            fw_ref = open(out_dir_chr + str(_start) + "_" + str(_end) + "_ref.fasta","w")
            fw_sv = open(out_dir_chr + str(_start) + "_" + str(_end) + "_sv.fasta","w")
            ref_seq = ref_dict[_start-500:_end+500]
            fw_ref.writelines(">human_ref" + "\n")
            fw_ref.writelines(ref_seq + "\n")
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
                event_name = one_value[5]
                chr_num = one_value[1]
                SV_size = one_value[4]
                ref_start = int(one_value[2])
                ref_end = int(one_value[3])
                _contig_1_dict = pickle.load(open(out_dir + sample_name + "/" + "contig_dict_" + sample_name + "_1_" + str(chr_num) + ".p" ,"rb"))
                _contig_2_dict = pickle.load(open(out_dir + sample_name + "/" + "contig_dict_" + sample_name + "_2_" + str(chr_num) + ".p","rb"))
                if contig_1_num != 0:
                    if contig_1_start >= 500 and (contig_1_end + 500) <= len(_contig_1_dict[contig_1_info]):
                        contig_1_seq = _contig_1_dict[contig_1_info][contig_1_start-500:contig_1_end+500]
                        fw_sv.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c1_" + GT + "_1_" + SV_size + "\n")
                        if contig_1_strand == "+":
                            fw_sv.writelines(contig_1_seq + "\n")
                        elif contig_1_strand == "-":
                            fw_sv.writelines(translate(contig_1_seq) + "\n")
                elif contig_1_num == 0:
                    count_contig_0 += 1
                    print(contig_1_num)
                    ref_seq_1 = ref_dict[ref_start-500:ref_end+500]
                    bam_file = out_dir + sample_name + "/" +  sample_name + "_" + chr_num + "_contig.1.bam"
                    sam_file_temp = "temp.sam"
                    bam_file_temp = "temp.bam"
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
                                    Align_refseq_to_contig(ref_seq_1,target_seq)
                                    sam_out = pysam.AlignmentFile("output.bam", "rb")
                                    for read_temp in sam_out.fetch():
                                        mapq_temp = int(read_temp.mapping_quality)
                                        if mapq_temp >= 20:
                                            align_start_temp = int(read_temp.pos)
                                            cigar_temp = read_temp.cigar
                                            total_use = Cal_extract_steps(cigar_temp)
                                            #### contig sequence extract
                                            if total_use >= 0.9*1000:
                                                target_seq_extract = target_seq[align_start_temp:align_start_temp+total_use]
                                                fw_ref.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c1_" + GT + "_0_"  + SV_size + "\n")
                                                fw_ref.writelines(translate(target_seq_extract) + "\n")
                                                Delete_temp_files()
                                                break
                                        else:
                                            _count_miss += 1
                                    sam_out.close()
                            else:
                                target_seq = read_seq
                                if len(target_seq) > len(ref_seq_1)*1.5:
                                    Align_refseq_to_contig(ref_seq_1,target_seq)
                                    sam_out = pysam.AlignmentFile("output.bam", "rb")
                                    for read_temp in sam_out.fetch():
                                        mapq_temp = int(read_temp.mapping_quality)
                                        if mapq_temp >= 20:
                                            align_start_temp = int(read_temp.pos)
                                            cigar_temp = read_temp.cigar
                                            total_use = Cal_extract_steps(cigar_temp)
                                            if total_use >= 0.9*1000:
                                                target_seq_extract = target_seq[align_start_temp:align_start_temp+total_use]
                                                fw_ref.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c1_" + GT + "_0_"  + SV_size + "\n")
                                                fw_ref.writelines(target_seq_extract + "\n")
                                                Delete_temp_files()
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
                        fw_sv.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c2_" + GT + "_1_" + SV_size + "\n")
                        if contig_2_strand == "+":
                            fw_sv.writelines(contig_2_seq + "\n")
                        elif contig_2_strand == "-":
                            fw_sv.writelines(translate(contig_2_seq) + "\n")
                elif contig_2_num == 0:
                    count_contig_0 += 1
                    ref_seq_2 = ref_dict[ref_start-500:ref_end+500]
                    print(contig_2_num)
                    bam_file = out_dir + sample_name + "/" + sample_name + "_" + chr_num + "_contig.2.bam"
                    sam_file_temp = "temp.sam"
                    bam_file_temp = "temp.bam"
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
                                    Align_refseq_to_contig(ref_seq_2,target_seq)
                                    sam_out = pysam.AlignmentFile("output.bam", "rb")
                                    for read_temp in sam_out.fetch():
                                        mapq_temp = int(read_temp.mapping_quality)
                                        if mapq_temp >= 20:
                                            align_start_temp = int(read_temp.pos)
                                            cigar_temp = read_temp.cigar
                                            total_use = Cal_extract_steps(cigar_temp)
                                            if total_use >= 0.9*1000:
                                                target_seq_extract = target_seq[align_start_temp:align_start_temp+total_use]
                                                fw_ref.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c2_" + GT + "_0_"  + SV_size + "\n")
                                                fw_ref.writelines(translate(target_seq_extract) + "\n")
                                                Delete_temp_files()
                                                break
                                        else:
                                            _count_miss += 1
                                    sam_out.close()
                            else:
                                target_seq = read_seq
                                if len(target_seq) > len(ref_seq_2)*1.5:
                                    Align_refseq_to_contig(ref_seq_2,target_seq)
                                    sam_out = pysam.AlignmentFile("output.bam", "rb")
                                    for read_temp in sam_out.fetch():
                                        mapq_temp = int(read_temp.mapping_quality)
                                        if mapq_temp >= 20:
                                            align_start_temp = int(read_temp.pos)
                                            cigar_temp = read_temp.cigar
                                            total_use = Cal_extract_steps(cigar_temp)
                                            if total_use >= 0.9*1000:
                                                target_seq_extract = target_seq[align_start_temp:align_start_temp+total_use]
                                                fw_ref.writelines(">" + sample_name + "_" + event_name.split("event")[1] + "_c2_" + GT + "_0_"  + SV_size + "\n")
                                                fw_ref.writelines(target_seq_extract + "\n")
                                                Delete_temp_files()
                                                break
                                        else:
                                            _count_miss += 1
                                    sam_out.close()
                        else:
                            _count_miss += 1

                    sam_file_curr.close()                

            fw_ref.close()
            fw_sv.close()
        
    print("count_contig_0: " + str(count_contig_0))
    print("_count_miss: " + str(_count_miss))





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
    Extract_SV_from_all_samples(sample_list,in_dir,out_dir,chr_start,chr_end)
    merge_SV_dict = Merge_SV(out_dir,chr_start,chr_end)

    """
    ### have done this for deletions
    for chr_num in range(chr_start,chr_end + 1):
        Align_hp_fasta_to_bam(assembly_dir,out_dir,chr_num)
    """
    for chr_num in range(chr_start,chr_end + 1):
        out_dir_chr = out_dir + "/" + "MSA_ins_results_chr" + str(chr_num) + "/"
        if os.path.exists(out_dir_chr):
            print("using existing output folder: " + out_dir_chr)
        else:
            os.makedirs(out_dir_chr)
        if chr_num == 23:
            extract_ref_chr(ref_file,"X",out_dir_chr)
            read_ref(out_dir_chr + "genome_ref_chrX.fasta",chr_num,out_dir_chr)
        else:
            extract_ref_chr(ref_file,chr_num,out_dir_chr)
            read_ref(out_dir_chr + "genome_ref_chr" + str(chr_num) + ".fasta",chr_num,out_dir_chr)

        for sample_name in sample_list:
            hp1_file = assembly_dir + sample_name + "/" + "Assembly_Contigs_files/" + "Aquila_Contig_chr" + str(chr_num) + "_hp1.fasta"  
            hp2_file = assembly_dir + sample_name + "/" + "Assembly_Contigs_files/" + "Aquila_Contig_chr" + str(chr_num) + "_hp2.fasta"  
            contig_dict_1 = extract_contig_fasta(hp1_file,sample_name + "_1_chr" + str(chr_num),out_dir + sample_name + "/")
            contig_dict_2 = extract_contig_fasta(hp2_file,sample_name + "_2_chr" + str(chr_num),out_dir + sample_name + "/")
        merge_SV_dict = pickle.load(open(out_dir + "merge_INS_dict_chr" + str(chr_num) + ".p", "rb"))
        Extract_fasta_for_SV(merge_SV_dict,out_dir,out_dir_chr,chr_num)
