import pdb
#pdb.set_trace()
import gzip
from collections import defaultdict
import pickle
import os
import sys
from subprocess import Popen
from argparse import ArgumentParser
import glob
from collections import Counter
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 


def Most_Common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]


def convert_vcf_to_fasta(vcf_file,output_fasta):
    f = open(vcf_file,"r")
    fw_fasta = open(output_fasta,"w")
    count = 0
    for line in f:
        data = line.rsplit()
        if data[0][0] != "#":
            ref = data[3]
            alt = data[4]
            GT = data[9].split(":")[0]
            count += 1
            event_name = data[2]
            if data[7] == "SVTYPE=INS":
                seq_used = data[4]
            elif data[7] == "SVTYPE=DEL":
                seq_used = data[3]
            fw_fasta.writelines(">" + event_name + "\n")
            fw_fasta.writelines(seq_used + "\n")
 
    fw_fasta.close()
    f.close()


def save_pickle_file(dict1,fname):
    for value in dict1:
        dict1[value] = dict(dict1[value])
    my_dict = dict(dict1)
    with open(fname, "wb") as f:
        pickle.dump(my_dict,f) 


def check_repeats_len(save_loci):
    repeats_len_dict = defaultdict(int)
    for loci in save_loci:
        _start = loci[0]
        _end = loci[1]
        for _step in range(_start, _end+1):
            repeats_len_dict[_step] = 1

    repeats_len = len(repeats_len_dict)
    return repeats_len


def Get_reads_match_tandem_repeats(tandem_repeats_file,human_ref_len):
    qname_match_repeats = defaultdict(int)
    f = open(tandem_repeats_file,"r")
    count = 0
    count_repeats = 0
    repeats_percent = -1
    for line in f:
        data = line.rsplit()
        if count >= 8:
            if data != []:
                if data[0] == "Sequence:":
                    if count_repeats > 1:
                        repeats_len = check_repeats_len(save_loci)
                        repeats_percent = float(repeats_len)/human_ref_len
                        print("repeats percent: " + str(repeats_percent))
                        qname_match_repeats[prev_qname] = repeats_percent
                        count_repeats = 0
                        save_loci = []
                    else:
                        count_repeats = 0
                        save_loci = []
                    prev_qname = data[1]
                    
                elif data[0] == "Parameters:":
                    count_repeats += 1
                else:
                    count_repeats += 1
                    _start = int(data[0])
                    _end = int(data[1])
                    save_loci.append([_start,_end])
        count += 1

    if count_repeats > 1:
        repeats_len = check_repeats_len(save_loci)
        repeats_percent = float(repeats_len)/human_ref_len
        print("repeats percent: " + str(repeats_percent))
        qname_match_repeats[prev_qname] = repeats_percent

    print("done~")
    return repeats_percent


def extract_SV_del_seq(vcf_file):
    f = open(vcf_file,"r")
    event_dict = defaultdict()
    for line in f:
        data = line.rsplit()
        if data[0][0] != "#":
            event_name = data[2]
            event_dict[event_name] = [int(data[1]),data[3]]
    return event_dict


def extract_SV_ins_seq(vcf_file):
    f = open(vcf_file,"r")
    event_dict = defaultdict()
    for line in f:
        data = line.rsplit()
        if data[0][0] != "#":
            event_name = data[2]
            event_dict[event_name] = [int(data[1]),data[4]]
    return event_dict


def remove_dash(val):
    new_val = []
    for one_val in val:
        if one_val != "-":
            new_val.append(one_val)
    return new_val


def Cal_SV_repeat_percentage_by_trf(in_dir,out_dir,chr_num):
    _flank_len = 2*10
    txt_files_all = glob.glob(in_dir +  "*_2.txt")
    count = 0
    ####trf_repeat_percent = defaultdict(float)
    SV_len_dict = defaultdict(int)
    fw_total = open(out_dir + "Consensus_seq_for_SV_chr" + str(chr_num) + ".fasta","w")
    for one_file in txt_files_all:
        fasta_file_name = one_file.split(".txt")[0] + ".fasta"
        sv_name = fasta_file_name.split("/")[-1]
        consensus_dict = defaultdict(list)
        with open(one_file,"r") as f:
            for line in f:
                data = line.rsplit()
                try:
                    use_seq = data[1]
                    for ii in range(len(use_seq)):
                        consensus_dict[ii].append(use_seq[ii])
                    flag = 1
                except:
                    flag = 0
        if flag == 1:
            consensus_string = ""
            for key,val in consensus_dict.items():
                val_new = remove_dash(val)
                if val_new != []:
                    consensus_string += Most_Common(val_new)

            sv_total_len = len(consensus_string)  
            fw = open(fasta_file_name,"w")
            fw.writelines(">" + sv_name + "\n")
            fw.writelines(consensus_string + "\n")
            fw.close()
            SV_len_dict[sv_name] = len(consensus_string) - _flank_len

            fw_total.writelines(">" + sv_name + "\n")
            fw_total.writelines(consensus_string + "\n")


        count += 1
      
        ####### run trf #######
        """
        if flag == 1:
            try:
                trf_cmd = code_path + "trf_tools/trf " + fasta_file_name + " 2 5 7 80 10 10 2000 -d -h" 
                Popen(trf_cmd,shell=True).wait()
            except:
                trf_cmd = "trf " + fasta_file_name + " 2 5 7 80 10 10 2000 -d -h"  
                Popen(trf_cmd,shell=True).wait()
            trf_file_name = fasta_file_name.split("/")[-1]
            try:
                repeats_percent = Get_reads_match_tandem_repeats(trf_file_name + ".2.5.7.80.10.10.2000.dat",sv_total_len)
                trf_repeat_percent[trf_file_name] = repeats_percent
                rm_cmd = "rm " + trf_file_name + ".2.5.7.80.10.10.2000.dat"
                Popen(rm_cmd,shell=True).wait()
            except:
                pass
       """
    fw_total.close()
    ####pickle.dump(trf_repeat_percent,open(in_dir + "trf_repeat_percent_chr" + str(chr_num) + ".p","wb"))
    pickle.dump(SV_len_dict,open(in_dir + "SV_len_dict_chr" + str(chr_num) + ".p","wb"))
    return SV_len_dict


"""
in_dir = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS_hg38/MSA_SV_results_chr21/MSA_SV_files/"
chr_num = 21
out_dir = "./"
Cal_SV_repeat_percentage_by_trf(in_dir,out_dir,chr_num)
"""
