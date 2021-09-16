import pdb
#pdb.set_trace()
import pickle
from collections import defaultdict
import csv
import glob


def Print_linked_SNPs_for_SV(out_dir,snp_dict,gnomad_dict,use_chr_num):
    fw = open(out_dir  + "SV_linked_with_gnomad_SNP_chr" + str(use_chr_num) + ".txt","w")
    snp_dict_rsid = defaultdict(list)
    #count = 0
    #count_err = 0
    count_rsid = 0
    count_rsid_corr = 0
    count_rsid_err = 0
    for key, val_list in snp_dict.items():
        for one_val in val_list:
            chr_num = one_val[0]
            _pos = one_val[1]
            ref = one_val[2]
            alt = one_val[3]
            ##cor_ref = ref_dict[_pos]
            """
            if ref == cor_ref:
                count += 1
            else:
                count_err += 1
            """
            if (chr_num,_pos + 1) in gnomad_dict:
                [gnomad_rsid,gnomad_ref,gnomad_alt]  = gnomad_dict[(chr_num,_pos + 1)]
                #print(gnomad_dict[(chr_num,_pos + 1)])
                #print(ref,alt)
                #print("--------------")
                count_rsid += 1
                if gnomad_ref == ref and gnomad_alt == alt:
                    count_rsid_corr += 1
                    snp_dict_rsid[key].append([chr_num,_pos,ref,alt,gnomad_rsid])
                else:
                    count_rsid_err += 1
                
    ##print(count,count_err,count_rsid)
    print(count_rsid_corr,count_rsid_err)
    #pickle.dump(snp_dict_rsid,open(out_dir + "snp_dict_rsid_chr" + str(use_chr_num) + ".p","wb"))  
    print("#: " + str(len(snp_dict)))
    print("#: " + str(len(snp_dict_rsid)))
    
    for key, val_list in snp_dict_rsid.items():
        fw.writelines(key + "\t" + str(len(val_list)) + "\t")
        for val in val_list:
            fw.writelines(val[4] + "\t" + str(val[1]) + "\t" + val[2] + "\t" + val[3] + "\t")
        fw.writelines("\n")
    
    fw.close()
    """
    with open(out_dir + "SV_linked_with_gnomad_SNP_chr" + str(use_chr_num) + ".txt", "r") as in_file:
        stripped = (line.strip() for line in in_file)
        lines = (line.split("\t") for line in stripped if line)
        with open(out_dir + "SV_linked_with_gnomad_SNP_chr" + str(use_chr_num) + ".csv", 'w') as out_file:
            writer = csv.writer(out_file)
            writer.writerows(lines)
    """

def Print_haplotypes_for_linked_SNP_and_SV(in_dir,gnomad_dict,output_file,output_file_2,chr_num,num_of_samples):
    fw = open(output_file,"w")
    fw_2 = open(output_file_2,"w")
    fw.writelines("Name" + "\t" + "sample_hap" + "\t" + "chr_num" + "\t" + "SV_start" + "\t" + "SV_end" + "\t" + "SNP_pos" + "\t" + "gnomad_rsid" + "\t" + "ref" + "\t" + "alt" + "\t" +  "Haplotype" + "\n")
    fw_2.writelines("Name" + "\t"  + "chr_num" + "\t" + "SV_start" + "\t" + "SV_end" + "\t" + "SNP_pos" + "\t" + "gnomad_rsid" + "\t" + "ref" + "\t" + "alt" + "\t" +  "count_00" + "\t" + "count_01" + "\t"  + "count_10" + "\t" + "count_11"  +  "\n")
    count_corr = 0
    count_wrong = 0
    count_error = 0
    pickle_files = sorted(glob.glob(in_dir + "Linked_SNP_and_SV/" + "*Linked_Haplotype_dict.p")) 
    for one_pickle in pickle_files:
        count_first = 1
        sv_name  = one_pickle.split("/")[-1].split("_Linked_Haplotype_dict.p")[0]
        ref_allele_dict = defaultdict(list)
        count_haplotype = defaultdict(lambda: defaultdict(int))
        linked_dict = pickle.load(open(one_pickle,"rb"))
        for sample_name, val_list in linked_dict.items():
            if sample_name == ">human_ref":
                for chr_pos_list,val in val_list.items():
                    if len(chr_pos_list) == 2:
                        chr_num = chr_pos_list[0]
                        _pos = chr_pos_list[1]
                        ref = val[1]
                        if (chr_num,_pos + 1) in gnomad_dict:
                            [gnomad_rsid,gnomad_ref,gnomad_alt]  = gnomad_dict[(chr_num,_pos + 1)]
                            if gnomad_ref == ref and len(ref) == 1 and len(gnomad_alt) == 1:
                                ref_allele_dict[_pos] = [ref,gnomad_rsid,gnomad_ref,gnomad_alt]
                            else:
                                count_wrong += 1
                    if len(chr_pos_list) == 3:
                        chr_num = chr_pos_list[0]
                        _start = chr_pos_list[1]
                        _end = chr_pos_list[2]
                        ref = val[0]
                        ref_allele_dict[(_start,_end)] = [ref]
        if len(ref_allele_dict) > 1:
            for sample_name, val_list in linked_dict.items():
                if sample_name != ">human_ref":
                    Haplotype_dict_snp = defaultdict(list)
                    Haplotype_dict_sv = defaultdict(str)
                    for chr_pos_list, val in val_list.items():
                        if len(chr_pos_list) == 2:
                            chr_num = chr_pos_list[0]
                            _pos = chr_pos_list[1]
                            _allele = val[1]
                            if _pos in ref_allele_dict:   ## in gnomad_dict
                                ref_allele = ref_allele_dict[_pos][0]
                                alt_allele = ref_allele_dict[_pos][3]
                                gnomad_rsid = ref_allele_dict[_pos][1]
                                gnomad_ref = ref_allele_dict[_pos][2]
                                if _allele == ref_allele:
                                    Haplotype_dict_snp[_pos] = ["0",gnomad_rsid,gnomad_ref,alt_allele]
                                elif _allele == alt_allele:
                                    Haplotype_dict_snp[_pos] = ["1",gnomad_rsid,gnomad_ref,alt_allele]
                        if len(chr_pos_list) == 3:
                            chr_num = chr_pos_list[0]
                            _start = chr_pos_list[1]
                            _end = chr_pos_list[2]
                            _allele = val[0]
                            if _start != -1 and _end != -1:
                                ref_allele = ref_allele_dict[(_start,_end)][0]
                                len_1 = len(_allele.replace("-",""))
                                len_2 = len(ref_allele.replace("-",""))
                                origin_hp = sample_name.split("_")[-2]
                                #if abs(len_1 - len_2)/(_end - _start) <= 0.1 and origin_hp == "0":
                                if ref_allele == _allele:
                                    Haplotype_dict_sv[(_start,_end)] = "0"
                                elif origin_hp == "1" and ref_allele != _allele:
                                    Haplotype_dict_sv[(_start,_end)] = "1"
                    for key_sv,val_sv in Haplotype_dict_sv.items():
                        for key_snp,val_snp in Haplotype_dict_snp.items():
                            _haplotype = val_sv + val_snp[0]
                            fw.writelines(sv_name + "\t" + sample_name + "\t" +  str(chr_num) + "\t" + str(key_sv[0]) + "\t" + str(key_sv[1]) + "\t" + str(key_snp) + "\t" + str(val_snp[1]) + "\t" + str(val_snp[2]) + "\t" + str(val_snp[3]) + "\t" +  _haplotype + "\n")
                            count_haplotype[(chr_num,key_sv[0],key_sv[1],key_snp,val_snp[1],val_snp[2],val_snp[3])][_haplotype] += 1
        if sv_name == "chr21_22745568_22745648_ins_del":
            print(sv_name)
        for key, val_list in count_haplotype.items():
            chr_num = key[0]
            sv_start = key[1]
            sv_end = key[2]
            snp_pos = key[3]
            snp_gnomad_rsid = key[4]
            snp_ref = key[5]
            snp_alt = key[6]
            count_00 = 0; count_01 = 0; count_10 = 0; count_11 = 0
            for _haplotype, _count in val_list.items():
                if _haplotype == "00":
                    count_00 = _count
                elif _haplotype == "01":
                    count_01 = _count
                elif _haplotype == "10":
                    count_10 = _count
                elif _haplotype == "11":
                    count_11 = _count
                sv_raw_start = int(sv_name.split("_")[1])
                sv_raw_end = int(sv_name.split("_")[2])
            total_haps = count_00 + count_01 + count_10 + count_11
            if count_00 <= num_of_samples and count_01 <= num_of_samples and count_10 <= num_of_samples and count_11 <= num_of_samples and total_haps <= 2*num_of_samples:
                if count_first == 1:
                    fw_2.writelines(sv_name + "\t" + str(chr_num) + "\t" + str(sv_start) + "\t" + str(sv_end) + "\t" + str(snp_pos) + "\t" + str(snp_gnomad_rsid) + "\t" + str(snp_ref) + "\t" + str(snp_alt) + "\t" +  str(count_00) + "\t" + str(count_01) +  "\t" + str(count_10) + "\t" + str(count_11) + "\n")
                    count_first += 1
                else:
                    #fw_2.writelines("" + "\t" + "" + "\t" + "" + "\t" + "" + "\t" + str(snp_pos) + "\t" + str(snp_gnomad_rsid) + "\t" + str(snp_ref) + "\t" + str(snp_alt) + "\t" +  str(count_00) + "\t" + str(count_01) +  "\t" + str(count_10) + "\t" + str(count_11) + "\n")
                    fw_2.writelines(sv_name + "\t" + str(chr_num) + "\t" + str(sv_start) + "\t" + str(sv_end) + "\t" + str(snp_pos) + "\t" + str(snp_gnomad_rsid) + "\t" + str(snp_ref) + "\t" + str(snp_alt) + "\t" +  str(count_00) + "\t" + str(count_01) +  "\t" + str(count_10) + "\t" + str(count_11) + "\n")
            else:
                count_error += 1
      
    fw.close()            
    fw_2.close()  
    print(count_error)



def Print_linked_SNPs_info(out_dir_chr,out_dir_table,gnomad_dict,use_chr_num,num_of_samples):
    snp_dict_file = out_dir_chr + "/SNP_dict.p"
    snp_dict = pickle.load(open(snp_dict_file,"rb"))

    #gnomad_dict_file = gnomad_dir + "/" + "gnomad_snp_dict_chr" + str(use_chr_num) + ".p"
    #gnomad_dict = pickle.load(open(gnomad_dict_file,"rb"))

    Print_linked_SNPs_for_SV(out_dir_table,snp_dict,gnomad_dict,use_chr_num)


    output_file = out_dir_table + "Haplotypes_by_Linked_SNP_and_SV_chr" + str(use_chr_num) + ".txt"
    output_file_2 = out_dir_table + "Haplotypes_by_Linked_SNP_and_SV_chr" + str(use_chr_num) + "_counts.txt"
    Print_haplotypes_for_linked_SNP_and_SV(out_dir_chr,gnomad_dict,output_file,output_file_2,use_chr_num,num_of_samples)

