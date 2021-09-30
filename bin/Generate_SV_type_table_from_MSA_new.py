import pdb
#pdb.set_trace()
import glob
import os
from collections import defaultdict
import pandas as pd
import pickle


def Generate_table(in_dir,output_file,SV_len_dict,SV_bk_dict,use_flag):
    break_flag = 0
    count_break = 0
    count = 0
    count_unsure = 0
    fw = open(output_file,"w")
    if use_flag == 1:
        fw.writelines("Name" + "\t" + "chr_num" + "\t" + "start" + "\t" + "end" + "\t" + "SV_length" + "\t" + "Num_of_Seqs" + "\t" + "Dash_To_N" + "\t" + "N_To_dash"  + "\t" + "Left_dash_break" + "\t" + "Right_dash_break"  + "\t" + "Num_of_Samples" + "\t" + "Num_of_Samples_with_dups" + "\n")
    else:
        fw.writelines("Name" + "\t" + "chr_num" + "\t" + "start" + "\t" + "end" + "\t" + "SV_length" + "\t" + "Num_of_Seqs" + "\t" + "Dash_To_N" + "\t" + "N_To_dash"  + "\t" + "Left_dash_break" + "\t" + "Right_dash_break"  + "\n")

    flank_len = 20
    msa_files_all = sorted(glob.glob(in_dir +  "*SV_Ref_bk_2.txt"),key=os.path.getsize,reverse=False)
    for one_file in msa_files_all:
        if os.path.getsize(one_file) > 0:
            sv_dict = defaultdict(str)
            count += 1
            file_name = one_file.split("/")[-1].split("_SV_Ref_bk_2.txt")[0]
            sv_len_list = []
            with open(one_file,"r") as f:
                for line in f:
                    data = line.rsplit()
                    _name = data[0]
                    if len(data) == 1:
                        break_flag = 1
                        break
                    _sv_string = data[1]
                    sv_dict[_name] = _sv_string
                    sv_len_list.append(len(_sv_string))
            if break_flag == 1:
                count_break += 1
                break_flag = 0
            else:
                if len(set(sv_len_list)) > 1:
                    count_unsure += 1
                sv_len = min(sv_len_list)
                dash_to_N = defaultdict(lambda: defaultdict(int))
                N_to_dash = defaultdict(lambda: defaultdict(int))
                for ii in range(1,sv_len):
                    for key, val in sv_dict.items():
                        prev_str = val[ii - 1]
                        cur_str = val[ii]
                        if prev_str == "-" and cur_str != "-":
                            dash_to_N[ii][prev_str,cur_str] += 1
                        if prev_str != "-" and cur_str == "-":
                            N_to_dash[ii][prev_str,cur_str] += 1
                
                count_dash_to_N = len(dash_to_N)
                count_N_to_dash = len(N_to_dash)
                num_of_seqs = len(sv_dict)
                if file_name in SV_bk_dict:
                    SV_start = SV_bk_dict[file_name][0] 
                    SV_end = SV_bk_dict[file_name][1]
                else:
                    SV_start = -1 
                    SV_end = -1
                if file_name + "_SV_Ref_bk_2.fasta" in SV_len_dict:
                    curr_sv_len = SV_len_dict[file_name + "_SV_Ref_bk_2.fasta"]
                else:
                    curr_sv_len = -1
                curr_chr_num = file_name.split("_")[0]
                fw.writelines(file_name + "\t" + curr_chr_num + "\t" + str(SV_start) + "\t" + str(SV_end) +  "\t" + str(curr_sv_len) + "\t" + str(num_of_seqs) +  "\t" + str(count_dash_to_N) + "\t" + str(count_N_to_dash) + "\t" )

                # Count the dinuc
                dinuc_dict = defaultdict(int)
                for key,val in sv_dict.items():
                    val_use = val.replace("-","")
                    for jj in range(len(val_use) - 1):
                        curr_dinuc = val_use[jj:jj+2]
                        dinuc_dict[curr_dinuc] += 1

                left_dash_break = 0
                right_dash_break = 0
                for key, val in sv_dict.items():
                    if val[0] == "-":
                        left_dash_break = 1
                        break
                    if val[-1] == "-":
                        right_dash_break = 1
                        break

                """
                for _dinuc in ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']:
                    fw.writelines(str(dinuc_dict[_dinuc]) + "\t")
                """
                fw.writelines(str(left_dash_break) + "\t" + str(right_dash_break) + "\t")

                #### check allele duplication for individuals ####
                if use_flag == 1:
                    sample_dict_1 = defaultdict(int)
                    sample_dict_2 = defaultdict(int)
                    with open(one_file,"r") as f :
                        for line in f:
                            data = line.rsplit()
                            if "human_ref" not in data[0]:
                                sample_info =  data[0].replace("*","").replace(">","").split("_")
                                if "c1" in sample_info:
                                    _use_idx = sample_info.index("c1")
                                elif "c2" in sample_info:
                                    _use_idx = sample_info.index("c2")
                                sample_name_hp = sample_info[0] + "_" + sample_info[_use_idx]
                                sample_dict_1[sample_name_hp] += 1
                                if "_c1" in sample_name_hp:
                                    sample_name = sample_name_hp.split("_c1")[0]
                                    sample_dict_2[sample_name] = 1
                                elif "c2" in sample_name_hp:
                                    sample_name = sample_name_hp.split("_c2")[0]
                                    sample_dict_2[sample_name] = 1

                    count_sample_dup = defaultdict(int)
                    for sample_name,val in sample_dict_1.items():
                        if val > 1:
                            count_sample_dup[sample_name.split("_")[0]] = 1
                    total_sample_dup = len(count_sample_dup)
                    total_sample_num = len(sample_dict_2)
                    fw.writelines("\t" + str(total_sample_num) + "\t" + str(total_sample_dup) + "\n")
                else:
                   fw.writelines("\n")
                

    fw.close()
    pd.read_csv(output_file, encoding='utf8', sep='\t', dtype='unicode').to_excel(output_file + '.xlsx', sheet_name='sheet1', index=False)
    print("count_unsure: " + str(count_unsure))
    print("count_break: " + str(count_break))
    print("count: " + str(count))
    print("done!")


"""
in_dir = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS_test/MSA_SV_results_chr21/HARP_derived_del/"
in_dir_2 = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS_test/MSA_SV_results_chr21/MSA_SV_files/"
output_file = "test.txt"
trf_repeats_percent = pickle.load(open(in_dir_2 + "trf_repeat_percent_chr21.p","rb"))
alu_dict = pickle.load(open(in_dir_2 + "Alu_chr21.p","rb"))
L1_dict = pickle.load(open(in_dir_2 + "L1_chr21.p","rb"))
Generate_table(in_dir,output_file,trf_repeats_percent,alu_dict,L1_dict,0)
"""







