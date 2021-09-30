import pdb
#pdb.set_trace()
from collections import defaultdict
import pandas as pd


def Merge_table(input_table,output_table,derived_del_dict,derived_ins_dict,derived_complex_dict):
    fw = open(output_table,"w")
    count = 0
    with open(input_table,"r") as f:
        for line in f:
            data = line.rsplit()
            if count == 0:
                Name_idx = data.index("Name")
                fw.writelines("\t".join(data) + "\t" + "Ancestral_State" + "\n")
                total_fields_len = len(data) 
            else:
                _name = data[Name_idx]
                chr_num = _name.split("_")[0]
                chr_start = _name.split("_")[1]
                chr_end = _name.split("_")[2]
                total_len = len(data)
                if total_len > total_fields_len:
                    _diff = total_len - total_fields_len - 1
                    last_field = ""
                    for ii in range(total_fields_len - 1,total_len):
                        last_field += data[ii]
                    data_new =  data[:total_fields_len-1] 
                    data_new.append(last_field)
                else:
                    data_new = data
                
                if _name in derived_complex_dict:
                    fw.writelines("\t".join(data_new) + "\t" + "Com" + "\n")
                elif _name in derived_del_dict:
                    fw.writelines("\t".join(data_new) + "\t" + "Del" + "\n")
                elif _name in derived_ins_dict:
                    fw.writelines("\t".join(data_new) + "\t" + "Ins" + "\n")
                else:
                    fw.writelines("\t".join(data_new) + "\t" + "None" + "\n")
            count += 1
    fw.close()
    pd.read_csv(output_table, encoding='utf8', sep='\t', dtype='unicode').to_excel(output_table + '.xlsx', sheet_name='sheet1', index=False)


def read_table(input_table):
    count = 0
    _dict = defaultdict(int)
    with open(input_table,"r") as f:
        for line in f:
            data = line.rsplit()
            if count == 0:
                Name_idx = data.index("Name")
            else:
                _dict[data[Name_idx]] = 1
            count += 1
    return _dict 

            


"""
in_dir = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS_hg38/Final_tables/"
input_table = in_dir + "derived_del_msa_table_chr15.txt"
derived_del_dict = read_table(input_table)
input_table = in_dir + "derived_ins_msa_table_chr15.txt"
derived_ins_dict = read_table(input_table)
input_table = in_dir + "derived_complex_msa_table_chr15.txt"
derived_complex_dict = read_table(input_table)

#input_table = in_dir + "SV_msa_table_chr15_add_RepeatMasker.txt"
input_table = in_dir + "SV_msa_table_chr15.txt"
#output_table = in_dir + "SV_msa_table_chr15_add_RepeatMasker_add_AncestralState.txt"
output_table = in_dir + "SV_msa_table_chr15_add_AncestralState.txt"
Merge_table(input_table,output_table,derived_del_dict,derived_ins_dict,derived_complex_dict)
"""
