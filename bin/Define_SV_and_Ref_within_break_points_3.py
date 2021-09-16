import pdb
#pdb.set_trace()
from collections import defaultdict
import glob
import os
import pickle
import re


def save_pickle_file(dict1,fname):
    for value in dict1:
        dict1[value] = dict(dict1[value])
    my_dict = dict(dict1)
    with open(fname, "wb") as f:
        pickle.dump(my_dict,f) 


def find_correct_idx(_idx_list,val_len):
    use_len = val_len/2
    if abs(_idx_list[0] - _idx_list[1]) <= 20 and _idx_list[0] >= 480:
        use_idx = _idx_list[0]
    else:
        use_idx = 100000 
        for _idx in _idx_list:
            if abs(_idx - use_len) < abs(use_idx - use_len):
                use_idx = _idx
    return use_idx
        

def Cal_step_1(val,left_step,right_step,ref_start,ref_end,add_more_flank):
    count_ = 0
    count_step = 0
    for one_bit in val:
        if one_bit != "-":
            count_ += 1
        count_step += 1
        if count_step == left_step:
            ref_start_new = ref_start - add_more_flank + count_
            break
    count_ = 0
    count_step = 0
    for one_bit in val[::-1]:
        if one_bit != "-":
            count_ += 1
        count_step += 1
        if count_step == (right_step)*(-1):
            ref_end_new = ref_end + add_more_flank - count_
            break
    return (ref_start_new, ref_end_new)


def Cal_step_2(val,left_step,ref_start,ref_end,add_more_flank):
    count_ = 0
    count_step = 0
    for one_bit in val:
        if one_bit != "-":
            count_ += 1
        count_step += 1
        if count_step == left_step:
            ref_start_new = ref_start - add_more_flank + count_
            break
    ref_end_new = ref_end + add_more_flank 
    return (ref_start_new, ref_end_new)


def Cal_step_3(val,right_step,ref_start,ref_end,add_more_flank):
    ref_start_new = ref_start - add_more_flank
    count_ = 0
    count_step = 0
    for one_bit in val[::-1]:
        if one_bit != "-":
            count_ += 1
        count_step += 1
        if count_step == (right_step)*(-1):
            ref_end_new = ref_end + add_more_flank - count_
            break
    return (ref_start_new, ref_end_new)


def find_dash_used_2(val_2):
    dash_used_2 = ""
    for i in range(len(val_2)):
        if val_2[i] == "-":
            dash_used_2 += val_2[i]
        elif val_2 != "-" and len(dash_used_2) > 1:
            break

    return dash_used_2


def Define_SV_within_bk(afa_file,SNP_dict,save_sv_coord_dict,linked_snp_folder):
    sv_name = afa_file.split("/")[-1].split(".afa")[0]
    ###output_afa = afa_file.split(".afa")[0] + "_SV_Ref_bk.afa"
    output_txt = afa_file.split(".afa")[0] + "_SV_Ref_bk.txt"
    output_txt_2 = afa_file.split(".afa")[0] + "_SV_Ref_bk_2.txt"
    #fw = open(output_afa,"w")
    fw_txt = open(output_txt,"w")
    fw_txt_2 = open(output_txt_2,"w")
    flanking_len = 400
    flanking_len_total = 500
    add_more_flank = 10
    dash_used = "-"*20
    seq_dict = defaultdict(str)
    count = 1
    with open(afa_file,"r") as f:
        try:
            for line in f:
                data = line.rsplit()
                if data[0][0] == ">":
                    seq_name = data[0]
                    if count > 1:
                        seq_dict[pre_seq_name] = seq_str
                    pre_seq_name = seq_name
                    seq_str = ""
                    count += 1
                else:
                    seq_str += data[0]
            seq_dict[pre_seq_name] = seq_str
        except:
            pass
    ################## Search Linked SNPs with each SV #####################
    count_step = 0
    count_len = 0
    human_ref = seq_dict[">human_ref"]
    human_ref_len = len(human_ref)
    ref_start = int(afa_file.split("/")[-1].split("_")[1])
    ref_end = int(afa_file.split("/")[-1].split("_")[2])
    chr_num = sv_name.split("_")[0]

    for ii in range(human_ref_len):
        ref_allele = human_ref[ii]
        if ref_allele != "-":
            for key,val in seq_dict.items():
                _allele = val[ii]
                if _allele != ref_allele and _allele != "-":
                    alt_allele = _allele
                    _pos = ref_start - 500 + count_len
                    SNP_dict[sv_name].append([chr_num,_pos,ref_allele,alt_allele])
                    break
            count_len += 1
        count_step += 1
        if count_len == flanking_len_total:
            break
    
    count_step = 0
    count_len = 0
    for ii in reversed(range(human_ref_len)):
        ref_allele = human_ref[ii]
        if ref_allele != "-":
            for key,val in seq_dict.items():
                _allele = val[ii]
                if _allele != ref_allele and _allele != "-":
                    alt_allele = _allele
                    _pos = ref_end + 500 - count_len - 1 
                    SNP_dict[sv_name].append([chr_num,_pos,ref_allele,alt_allele])
                    break
            count_len += 1
        count_step += 1
        if count_len == flanking_len_total:
            break
    ################# save linked snp for each sample ###########################
    SNP_dict_2 = defaultdict(lambda: defaultdict(list))
    Linked_SNP_dict = defaultdict(lambda: defaultdict(list))
    count_step = 0
    count_len = 0
     
    for ii in range(human_ref_len):
        ref_allele = human_ref[ii]
        if ref_allele != "-":
            for key,val in seq_dict.items():
                _allele = val[ii]
                if _allele == ref_allele:
                    _pos = ref_start - 500 + count_len
                    SNP_dict_2[(chr_num,_pos)][key] = ["ref",ref_allele]
                elif _allele != ref_allele and _allele != "-":
                    alt_allele = _allele
                    _pos = ref_start - 500 + count_len
                    SNP_dict_2[(chr_num,_pos)][key] = ["alt",alt_allele]
            count_len += 1
        count_step += 1
        if count_len == flanking_len_total:
            break
    
    count_step = 0
    count_len = 0
    for ii in reversed(range(human_ref_len)):
        ref_allele = human_ref[ii]
        if ref_allele != "-":
            for key,val in seq_dict.items():
                _allele = val[ii]
                if _allele == ref_allele:
                    _pos = ref_end + 500 - count_len - 1 
                    SNP_dict_2[(chr_num,_pos)][key] = ["ref",ref_allele]
                elif _allele != ref_allele and _allele != "-":
                    alt_allele = _allele
                    _pos = ref_end + 500 - count_len - 1 
                    SNP_dict_2[(chr_num,_pos)][key] = ["alt",alt_allele]
            count_len += 1
        count_step += 1
        if count_len == flanking_len_total:
            break
    
    for key, val_list in SNP_dict_2.items():
        ref_alt_list = []
        for sample_name,val in val_list.items():
            ref_alt_list.append(val[0])
        if "ref" in ref_alt_list and "alt" in ref_alt_list:
            for sample_name,val in val_list.items():
                Linked_SNP_dict[key][sample_name] = val
        
    ############################################################################
    count_flanking = 0
    count_left_step = 0
    for one_bit in seq_dict[">human_ref"]:
        if one_bit != "-":
            count_flanking += 1
        count_left_step += 1
        if count_flanking == flanking_len:
            #print(count_flanking,count_left_step)
            break
    
    count_flanking = 0
    count_right_step = 0
    for one_bit in seq_dict[">human_ref"][::-1]:
        if one_bit != "-":
            count_flanking += 1
        count_right_step += 1
        if count_flanking == flanking_len:
            #print(count_flanking,count_right_step)
            break
    save_dict = defaultdict(str) 
    for key, val in seq_dict.items():
        #print(key,len(val))
        ## just plot the SV or Ref
        ##seq_1 = val[count_left_step:]
        ##seq_2 = seq_1[:-1*count_right_step]
        
        ## plot the SV or Ref plus left and right 10bp flanking region
        seq_1 = val[count_left_step - add_more_flank:]
        seq_2 = seq_1[:-1*count_right_step + add_more_flank]
        save_dict[key] = seq_2

    #for key, val in save_dict.items():
        #print(key,len(val))
        #fw.writelines(key + "\n")
        #fw.writelines(val + "\n")
    #fw.close()
    
    all_keys = []
    for key,val in save_dict.items():
        all_keys.append(key)
    maxlen = len(max(all_keys, key=len))
    all_keys_new = [key.rjust(maxlen, '*') for key in all_keys]
    for key in all_keys_new:
        val = save_dict[key.replace("*","")]
        fw_txt.writelines(key[1:] + "\t")
        fw_txt.writelines(val + "\n")
    fw_txt.close()
    ###### double check the first dash "-" and only keep 10bp from the left, do the same for the right #####
    count = 1
    min_left_idx_2 = 1000
    min_right_idx_2 = 1000
    for key, val in seq_dict.items():
        val_2 = save_dict[key][110:-110]
        dash_used_2 = find_dash_used_2(val_2)
        if dash_used_2 != "":
            _idx_list = [m.start() for m in re.finditer(dash_used_2, val)]
            if len(_idx_list) > 1:
                val_len = len(val)
                _idx = find_correct_idx(_idx_list,val_len)
            else:
                if _idx_list == []:
                    _idx = 1000
                else:
                    if _idx_list[0] >= 480:
                        _idx = _idx_list[0]
                    else:
                        _idx = 1000
            if count == 1:
                if _idx == 1000:
                    min_left_idx_2 = len(val)
                else:
                    min_left_idx_2 = _idx
            elif count > 1 and min_left_idx_2 > _idx and _idx != 1000:
                min_left_idx_2 = _idx
            count += 1

    count = 1
    for key, val in seq_dict.items():
        val_2 = save_dict[key][110:-110][::-1]
        dash_used_2 = find_dash_used_2(val_2)
        if dash_used_2 != "":
            _idx_list = [m.start() for m in re.finditer(dash_used_2, val[::-1])]
            if len(_idx_list) > 1:
                val_len = len(val)
                _idx = find_correct_idx(_idx_list,val_len)
            else:
                if _idx_list == []:
                    _idx = 1000
                else:
                    if _idx_list[0] >= 480:
                        _idx = _idx_list[0]
                    else:
                        _idx = 1000
            if count == 1:
                if _idx == 1000:
                    min_right_idx_2 = len(val)
                else:
                    min_right_idx_2 = _idx
            elif count > 1 and min_right_idx_2 > _idx and _idx != 1000:
                min_right_idx_2 = _idx
            count += 1
    ###########################################################################
    count = 1
    for key, val in save_dict.items():
        _idx = val.find(dash_used)
        if count == 1:
            if _idx == -1:
                min_left_idx = len(val)
            else:
                min_left_idx = _idx
        elif count > 1 and min_left_idx > _idx and _idx != -1:
            min_left_idx = _idx
        count += 1

    count = 1
    for key, val in save_dict.items():
        _idx = val[::-1].find(dash_used)
        if count == 1:
            if _idx == -1:
                min_right_idx = len(val)
            else:
                min_right_idx = _idx
        elif count > 1 and min_right_idx > _idx and _idx != -1:
            min_right_idx = _idx
        count += 1
    ############################################################################
    if min_left_idx > add_more_flank and min_right_idx > add_more_flank:
        for key in all_keys_new:
            val = seq_dict[key.replace("*","")]
            ########### To calculate SV human ref coordinates #############
            if key.replace("*","") == ">human_ref":
                try:
                    left_step = min_left_idx_2 
                except:
                    pass
                right_step = (-1)*min_right_idx_2
                left_ref = val[:left_step]
                count_left_dash = left_ref.count("-")
                right_ref = val[right_step:]
                count_right_dash = right_ref.count("-")
                ref_start = int(afa_file.split("/")[-1].split("_")[1])
                ref_end = int(afa_file.split("/")[-1].split("_")[2])
                ref_start_new = ref_start - 500 + left_step - count_left_dash 
                ref_end_new = ref_end + 500 + right_step + count_right_dash 
                save_sv_coord_dict[sv_name] = [ref_start_new,ref_end_new]
            ##############################################################
        for key in all_keys_new:
            val = save_dict[key.replace("*","")]
            val_2 = val[min_left_idx - add_more_flank : (-1)*min_right_idx + add_more_flank]
            val_SV = val[min_left_idx: (-1)*min_right_idx]
            Linked_SNP_dict[(chr_num,ref_start_new,ref_end_new)][key.replace("*","")] = [val_SV]
            fw_txt_2.writelines(key[1:] + "\t")
            fw_txt_2.writelines(val_2 + "\n")
    elif min_left_idx > add_more_flank and min_right_idx <= add_more_flank:
        for key in all_keys_new:
            val = save_dict[key.replace("*","")]
            ########### To calculate SV human ref coordinates #############
            if key.replace("*","") == ">human_ref":
                left_step = min_left_idx_2 
                left_ref = val[:left_step]
                count_left_dash = left_ref.count("-")
                ref_start = int(afa_file.split("/")[-1].split("_")[1])
                ref_end = int(afa_file.split("/")[-1].split("_")[2])
                ref_start_new = ref_start - 500 + left_step - count_left_dash 
                ref_end_new = -1
                save_sv_coord_dict[sv_name] = [ref_start_new,ref_end_new]
            ##############################################################
        for key in all_keys_new:
            val = save_dict[key.replace("*","")]
            val_2 = val[min_left_idx - add_more_flank :]
            val_SV = val[min_left_idx: ]
            Linked_SNP_dict[(chr_num,ref_start_new,ref_end_new)][key.replace("*","")] = [val_SV]
            fw_txt_2.writelines(key[1:] + "\t")
            fw_txt_2.writelines(val_2 + "\n")
    elif min_left_idx <= add_more_flank and min_right_idx > add_more_flank:
        for key in all_keys_new:
            val = save_dict[key.replace("*","")]
            ########### To calculate SV human ref coordinates #############
            if key.replace("*","") == ">human_ref":
                right_step = (-1)*min_right_idx_2
                right_ref = val[right_step:]
                count_right_dash = right_ref.count("-")
                ref_start = int(afa_file.split("/")[-1].split("_")[1])
                ref_end = int(afa_file.split("/")[-1].split("_")[2])
                ref_start_new = -1 
                ref_end_new = ref_end + 500 + right_step + count_right_dash 
                save_sv_coord_dict[sv_name] = [ref_start_new,ref_end_new]
            ##############################################################
        for key in all_keys_new:
            val = save_dict[key.replace("*","")]
            val_2 = val[ : (-1)*min_right_idx + add_more_flank]
            val_SV = val[: (-1)*min_right_idx]
            Linked_SNP_dict[(chr_num,ref_start_new,ref_end_new)][key.replace("*","")] = [val_SV]
            fw_txt_2.writelines(key[1:] + "\t")
            fw_txt_2.writelines(val_2 + "\n")
    else:
        for key in all_keys_new:
            val = save_dict[key.replace("*","")]
            ########### To calculate SV human ref coordinates #############
            if key.replace("*","") == ">human_ref":
                ref_start = int(afa_file.split("/")[-1].split("_")[1])
                ref_end = int(afa_file.split("/")[-1].split("_")[2])
                ref_start_new = -1
                ref_end_new = -1
                save_sv_coord_dict[sv_name] = [ref_start_new,ref_end_new]
            ##############################################################
        for key in all_keys_new:
            val = save_dict[key.replace("*","")]
            val_SV = val
            Linked_SNP_dict[(chr_num,ref_start_new,ref_end_new)][key.replace("*","")] = [val_SV]
            fw_txt_2.writelines(key[1:] + "\t")
            fw_txt_2.writelines(val + "\n")
    fw_txt_2.close()
    ####################### Generate Four Different Haploytpes for Linked SNP and SV  ####################################################
    Linked_Haplotype_dict = defaultdict(lambda: defaultdict(list))
    for key, val_list in Linked_SNP_dict.items():
        for sample_name, val in val_list.items():
            Linked_Haplotype_dict[sample_name][key] = val
    fname = linked_snp_folder + sv_name + "_Linked_Haplotype_dict.p" 
    save_pickle_file(Linked_Haplotype_dict,fname)
    
    ######################################################################################################################################
    return (SNP_dict,save_sv_coord_dict)





def Define_SV_break_point(in_dir):
    linked_snp_folder = in_dir + "/Linked_SNP_and_SV" + "/"
    if os.path.exists(linked_snp_folder):
        print("using existing output folder: " + linked_snp_folder)
    else:
        os.makedirs(linked_snp_folder)
    afa_files_all = sorted(glob.glob(in_dir +  "*.afa"))
    SNP_dict = defaultdict(list)
    save_sv_coord_dict = defaultdict(list)
    for one_afa in afa_files_all:
        if "_1.afa" not in one_afa:
            SNP_dict,save_sv_coord_dict = Define_SV_within_bk(one_afa,SNP_dict,save_sv_coord_dict,linked_snp_folder)

    pickle.dump(SNP_dict,open(in_dir + "SNP_dict.p","wb"))
    pickle.dump(save_sv_coord_dict,open(in_dir + "SVs_breakpoints_dict.p","wb"))

    print("#: " + str(len(SNP_dict)))






"""

in_dir = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS_hg38/MSA_SV_results_chr21/MSA_SV_files_test/"
Define_SV_break_point(in_dir)
"""
