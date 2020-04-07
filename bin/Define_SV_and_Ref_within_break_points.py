import pdb
#pdb.set_trace()
from collections import defaultdict
import glob
import os

def Define_SV_within_bk(afa_file):
    #output_afa = afa_file.split(".afa")[0] + "_SV_Ref_bk.afa"
    output_txt = afa_file.split(".afa")[0] + "_SV_Ref_bk.txt"
    #fw = open(output_afa,"w")
    fw_txt = open(output_txt,"w")
    flanking_len = 500
    seq_dict = defaultdict(str)
    with open(afa_file,"r") as f:
        for line in f:
            data = line.rsplit()
            if data[0][0] == ">":
                seq_name = data[0]
            else:
                seq_dict[seq_name] += data[0]
    count_flanking = 0
    count_left_step = 0
    for one_bit in seq_dict[">human_ref"]:
        if one_bit != "-":
            count_flanking += 1
        count_left_step += 1
        if count_flanking == flanking_len:
            print(count_flanking,count_left_step)
            break
    
    count_flanking = 0
    count_right_step = 0
    for one_bit in seq_dict[">human_ref"][::-1]:
        if one_bit != "-":
            count_flanking += 1
        count_right_step += 1
        if count_flanking == flanking_len:
            print(count_flanking,count_right_step)
            break
    save_dict = defaultdict(str) 
    for key, val in seq_dict.items():
        print(key,len(val))
        ## just plot the SV or Ref
        ##seq_1 = val[count_left_step:]
        ##seq_2 = seq_1[:-1*count_right_step]
        
        ## plot the SV or Ref plus left and right 10bp flanking region
        seq_1 = val[count_left_step - 10:]
        seq_2 = seq_1[:-1*count_right_step + 10]
        save_dict[key] = seq_2

    for key, val in save_dict.items():
        print(key,len(val))
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
     

    return seq_dict





if __name__ == "__main__":
    in_dir  = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_Aquila/MSA_del_results_chr21" + "/"
    afa_files_all = sorted(glob.glob(in_dir +  "*.afa"),key=os.path.getsize,reverse=True)
    count = 1
    for one_afa in afa_files_all:
        Define_SV_within_bk(one_afa)
        print(count)
        count += 1
    print("done!")


