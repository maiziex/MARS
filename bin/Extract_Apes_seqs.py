import pdb
#pdb.set_trace()
from collections import defaultdict
import gzip
import pickle
import os
import glob
from subprocess import Popen


basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
def translate(seq):
    seq_r = []
    for character in seq:
        seq_r.append(basecomplement[character])
    seq_r_2 = "".join(seq_r)
    rseqn = seq_r_2[::-1]
    return rseqn


def get_match_num_revised_2(cigar):
    cigar_len = len(cigar)
    cigar_list = []
    num_string = ""
    for i in range(cigar_len):
        letter = cigar[i]
        if letter.isalpha():
            cigar_list.append(num_string)
            num_string = ""
            cigar_list.append(letter)
        else:
            num_string += letter

    indices_match = [i for i, x in enumerate(cigar_list) if x == "M"]
    indices_del = [i for i, x in enumerate(cigar_list) if x == "D"]
    match_list = []
    del_list = []
    for idx in indices_match:
        match_list.append(int(cigar_list[idx-1]))
    for idx in indices_del:
        del_list.append(int(cigar_list[idx-1]))

    return (match_list,del_list)


def Write_Apes_seq_for_multialign(derived_SV_dict,MSA_results,chr_num,out_dir,Ape_ref_dir,Ape_name):
    all_files_dict = defaultdict(int)
    fasta_files_all = sorted(glob.glob(MSA_results +  "*.fasta"),key=os.path.getsize,reverse=True)
    for one_file in fasta_files_all:
        one_file_name = one_file.split("/")[-1].split(".fasta")[0]
        all_files_dict[one_file_name] = 1
    count = 0
    use_derived_SV = defaultdict(list)
    for key, val in derived_SV_dict.items():
        _idx = key[::-1].index("_")
        _name = key[:(-1)*(_idx + 1)]
        if _name in all_files_dict:
            count += 1
            use_derived_SV[_name] = val

    use_ape_chr_num = "" 
    for event_name, val_list in use_derived_SV.items():
        if os.path.exists(out_dir + event_name + ".fasta"):
            fw = open(out_dir + event_name + ".fasta","a")
        else:
            fw = open(out_dir + event_name + ".fasta","w")
        for val in val_list:
            ape_reverse_flag = val[-1]
            ape_start = val[1]
            ape_end = val[2]
            ape_chr_num = val[0]
            if use_ape_chr_num != ape_chr_num:
                ape_ref_chr = pickle.load(open(Ape_ref_dir + Ape_name + "_chr" + ape_chr_num + ".p","rb"))
            
            fw.writelines(">" + Ape_name  +  "\n")
            if ape_reverse_flag == 1:
                ape_string = ape_ref_chr[ape_start:ape_end]
                ape_string_reverse = translate(ape_string)
                fw.writelines(ape_string_reverse + "\n")
            elif ape_reverse_flag == 0:
                ape_string = ape_ref_chr[ape_start:ape_end]
                fw.writelines(ape_string + "\n")

            use_ape_chr_num = ape_chr_num
        
        fw.close()


def Write_origin_SV_seq_for_multialign(out_dir,MSA_results):
    fasta_files_all = sorted(glob.glob(out_dir +  "*.fasta"),key=os.path.getsize,reverse=True)
    for one_file in fasta_files_all:
        one_file_name = one_file.split("/")[-1].split(".fasta")[0]
        one_fasta_origin = MSA_results + one_file_name + ".fasta"
        fw = open(one_file,"a")
        with open(one_fasta_origin,"r") as f:
            for line in f:
                fw.writelines(line)
        fw.close()


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


def Find_derived_complex_SV(out_dir_del,out_dir_ins,out_dir_complex):
    derived_del = []
    derived_ins = []
    fasta_files_all = sorted(glob.glob(out_dir_del +  "*.fasta"),key=os.path.getsize,reverse=True)
    for one_file in fasta_files_all:
        one_file_name = one_file.split("/")[-1].split(".fasta")[0]
        derived_del.append(one_file_name)

    fasta_files_all = sorted(glob.glob(out_dir_ins +  "*.fasta"),key=os.path.getsize,reverse=True)
    for one_file in fasta_files_all:
        one_file_name = one_file.split("/")[-1].split(".fasta")[0]
        derived_ins.append(one_file_name)

    derived_complex = set(derived_del).intersection(set(derived_ins))
    for one_file_name in derived_complex:
        fw = open(out_dir_complex + one_file_name + ".fasta","w")
        del_file = out_dir_del + one_file_name + ".fasta"
        ins_file = out_dir_ins + one_file_name + ".fasta"
        used_name = defaultdict(int)
        count = 0
        with open(del_file,"r") as f:
            for line in f:
                data = line.rsplit()
                if count%2 == 0:
                    _name = data[0]
                    fw.writelines(line)
                    used_name[_name] = 1
                elif count%2 == 1:
                    fw.writelines(line)
                count += 1
        count = 0
        flag = 0
        with open(ins_file,"r") as f:
            for line in f:
                data = line.rsplit()
                if count%2 == 0:
                    _name = data[0]
                    if _name not in used_name:
                        fw.writelines(line)
                        flag = 1
                elif count%2 == 1:
                    if flag == 1:
                        fw.writelines(line)
                        flag = 0
                count += 1

        fw.close()
        _cmd = "rm " + del_file
        #Popen(_cmd,shell=True).wait()
        _cmd = "rm " + ins_file
        #Popen(_cmd,shell=True).wait()
 
    print("extract derived complex sv done")


def Extract_derived_SV_with_Apes_seq(Ape_ref_list,in_dir,in_dir_2,chr_num):
    out_folder_derived_del = "HARP_derived_del/"
    out_folder_derived_ins = "HARP_derived_ins/"
    out_folder_derived_complex = "HARP_derived_complex/"
    Ape_ref_dir = in_dir + "Apes_ref/"
    if os.path.exists(Ape_ref_dir):
        print("using existing Apes_ref folder: " + Ape_ref_dir)
    else:
        os.makedirs(Ape_ref_dir)
    """
    ##### Move this block into MARS_step2.py
    for Ape_ref_file in Ape_ref_list:
        Ape_name = Ape_ref_file.split("/")[-1].split(".")[0]
        load_and_save_Apes_ref(Ape_ref_dir,Ape_ref_file,Ape_name)
    """
    MSA_results = in_dir + "MSA_SV_results_chr" + str(chr_num) + "/MSA_SV_files/"
    out_dir_del  = in_dir + "MSA_SV_results_chr" + str(chr_num) + "/" + out_folder_derived_del
    out_dir_ins  = in_dir + "MSA_SV_results_chr" + str(chr_num) + "/" + out_folder_derived_ins
    out_dir_complex  = in_dir + "MSA_SV_results_chr" + str(chr_num) + "/" + out_folder_derived_complex
    derived_del_dict =  defaultdict(list)
    derived_ins_dict =  defaultdict(list)
    derived_complex_dict =  defaultdict(list)
    ### for derived deletions ###
    for Ape_ref_file in Ape_ref_list:
        Ape_name = Ape_ref_file.split("/")[-1].split(".")[0]
        derived_SV_file = in_dir_2 + "derived_del_by_flanking_dict_chr" +str(chr_num) + "_" + Ape_name + ".p"
        derived_SV_dict = pickle.load(open(derived_SV_file,"rb")) 
        Write_Apes_seq_for_multialign(derived_SV_dict,MSA_results,chr_num,out_dir_del,Ape_ref_dir,Ape_name)
    ####Write_origin_SV_seq_for_multialign(out_dir_del,MSA_results) ## don't use it for profile alignment
    
    ### for derived insertions ###
    for Ape_ref_file in Ape_ref_list:
        Ape_name = Ape_ref_file.split("/")[-1].split(".")[0]
        derived_SV_file = in_dir_2 + "derived_ins_by_flanking_dict_chr" +str(chr_num) + "_" + Ape_name + ".p"
        derived_SV_dict = pickle.load(open(derived_SV_file,"rb")) 
        Write_Apes_seq_for_multialign(derived_SV_dict,MSA_results,chr_num,out_dir_ins,Ape_ref_dir,Ape_name)
    
    ####Write_origin_SV_seq_for_multialign(out_dir_ins,MSA_results) ## don't use it for profile alignment
    
    ### for derived complex SV
    Find_derived_complex_SV(out_dir_del,out_dir_ins,out_dir_complex)
    

"""
in_dir = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS/"
Ape_ref_list = ["/oak/stanford/groups/arend/Xin/AssemblyProj/reference_align_2/References_other_species/Gorilla_gorilla/Gorilla_gorilla_ref.fasta","/oak/stanford/groups/arend/Xin/AssemblyProj/reference_align_2/References_other_species/Pan_troglodytes/pan_troglodytes_ref.fasta","/oak/stanford/groups/arend/Xin/AssemblyProj/reference_align_2/References_other_species/Pongo_abelii/pongo_abelii_ref.fasta","/oak/stanford/groups/arend/Xin/AssemblyProj/reference_align_2/References_other_species/macaca_mulatta/macaca_mulatta_ref.fasta"]
Extract_derived_SV_with_Apes_seq(Ape_ref_list,in_dir,21)
"""


