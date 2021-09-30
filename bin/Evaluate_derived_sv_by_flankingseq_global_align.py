import pdb
#pdb.set_trace()
import pysam
from collections import defaultdict
import pickle


def get_match_num_revised_3(cigar):
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


def check_clipping_not_in_the_end(cigar):
    flag = 0
    if "H" not in cigar and "S" not in cigar:
        flag = 1
    else:
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
        if cigar_list[-1] == "S" or cigar_list[-1] == "H":
            flag = 0
        else:
            flag = 1
    return flag

    
def check_clipping_not_at_the_beginning(cigar):
    flag = 0
    if "H" not in cigar and "S" not in cigar:
        flag = 1
    else:
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
        if cigar_list[1] == "S" or cigar_list[1] == "H":
            flag = 0
        else:
            flag = 1
    return flag


def Evaluate_derived_sv(bam_file,in_dir,use_chr_num,Ape_name):
    sam_file = pysam.AlignmentFile(bam_file, "rb")
    flanking_seq_dict = defaultdict(lambda: defaultdict(list))
    align_score_list = []
    event_align_chr_dict = defaultdict()
    for read in sam_file.fetch():
        tags = read.get_tags()
        AS_field = [s for s in tags if "AS" in s]
        mapq = read.mapping_quality
        if AS_field != [] and mapq >= 20:
            align_score = int(AS_field[0][1])
            align_score_list.append(align_score)
            qname = read.qname
            start_pos = read.pos
            flanking_idx = int(qname.split("_")[0])
            event_chr = qname.split("_")[1].split("chr")[-1]
            event_name = qname[2:]
            event_size = int(qname.split("_")[-1])
            read_cigar = read.cigarstring
            read_chr = read.reference_name
            if read.is_reverse:
                read_reverse_flag = 1
            else:
                read_reverse_flag = 0
            if flanking_idx in flanking_seq_dict[event_name]:
                flanking_seq_dict[event_name][flanking_idx].append([read_chr, start_pos, event_size,read_reverse_flag, read_cigar,align_score])
            else:
                flanking_seq_dict[event_name][flanking_idx] = [[read_chr, start_pos, event_size,read_reverse_flag, read_cigar,align_score]]


    count = 0
    count_1 = 0
    count_2 = 0
    count_3 = 0
    derived_del = defaultdict(int)
    derived_ins = defaultdict(int)
    SV_ref_coord = defaultdict(list)
    event_cigar = defaultdict(list)
    for event_name, val in flanking_seq_dict.items():
        if len(val) == 2:
            flanking_seq_1_info_list = val[1]
            flanking_seq_2_info_list = val[2]
            for flanking_seq_1_info in flanking_seq_1_info_list:
                for flanking_seq_2_info in flanking_seq_2_info_list:
                    chr_num = flanking_seq_1_info[0]
                    chr_num_2 = flanking_seq_2_info[0]
                    if chr_num != chr_num_2:
                        break
                    reverse_flag_1 = flanking_seq_1_info[3]
                    reverse_flag_2 = flanking_seq_2_info[3]
                    align_score_1 = flanking_seq_1_info[-1]
                    align_score_2 = flanking_seq_2_info[-1]
                    if reverse_flag_1 == 0 and reverse_flag_2 == 0 and align_score_1 > 0 and align_score_2 > 0:
                        count_1 += 1
                        sv_size = flanking_seq_1_info[2]
                        cigar_1 = flanking_seq_1_info[4]
                        cigar_2 = flanking_seq_2_info[4]
                        flag_1 = check_clipping_not_in_the_end(cigar_1)
                        flag_2 = check_clipping_not_at_the_beginning(cigar_2)
                        if flag_1 == 1 and flag_2 == 1:
                            start_1 = flanking_seq_1_info[1]
                            match_list, del_list = get_match_num_revised_3(cigar_1)
                            total_num = sum(match_list) + sum(del_list)
                            end_1 = start_1 + total_num
                            
                            start_2 = flanking_seq_2_info[1]
                            match_list, del_list = get_match_num_revised_3(cigar_2)
                            end_2 = start_2 + sum(match_list) + sum(del_list)
                            event_cigar[event_name] = [cigar_1,cigar_2]
                            if abs(start_2 - end_1) <= 2:
                                if event_name in derived_ins:
                                    derived_ins[event_name].append([chr_num,start_1, end_2, reverse_flag_1]) 
                                else:
                                    derived_ins[event_name] = [[chr_num,start_1, end_2, reverse_flag_1]] 
                            elif float(start_2 - end_1)/sv_size >= 0.9 and float(start_2 - end_1)/sv_size <= 1.1:
                                if event_name in derived_del:
                                    derived_del[event_name].append([chr_num,start_1, end_2, reverse_flag_1]) 
                                else:
                                    derived_del[event_name] = [[chr_num,start_1, end_2, reverse_flag_1]]
                                    
                    elif reverse_flag_1 == 1 and reverse_flag_2 == 1 and align_score_1 > 0 and align_score_2 > 0 :
                        count_2 += 1
                        sv_size = flanking_seq_1_info[2]
                        cigar_1 = flanking_seq_1_info[4]
                        cigar_2 = flanking_seq_2_info[4]
                        flag_1 = check_clipping_not_in_the_end(cigar_2)
                        flag_2 = check_clipping_not_at_the_beginning(cigar_1)
                        if flag_1 == 1 and flag_2 == 1:
                            start_2 =  flanking_seq_2_info[1]  
                            match_list, del_list = get_match_num_revised_3(cigar_2)
                            total_num = sum(match_list) + sum(del_list)
                            end_2 = start_2 + total_num  # == end_1

                            start_1 = flanking_seq_1_info[1]   # == start_2
                            match_list, del_list = get_match_num_revised_3(cigar_1)
                            end_1 = start_1 + sum(match_list) + sum(del_list)
                            event_cigar[event_name] = [cigar_1,cigar_2]

                            if abs(start_1 - end_2) <= 2:
                                if event_name in derived_ins:
                                    derived_ins[event_name].append([chr_num,start_2, end_1,reverse_flag_1]) 
                                else:
                                    derived_ins[event_name] = [[chr_num,start_2, end_1,reverse_flag_1]] 
                            elif float(start_1 - end_2)/sv_size >= 0.9 and float(start_1 - end_2)/sv_size <= 1.1:
                                if event_name in derived_del:
                                    derived_del[event_name].append([chr_num,start_2, end_1, reverse_flag_1])
                                else:
                                    derived_del[event_name] = [[chr_num,start_2, end_1, reverse_flag_1]] 

                    else:
                        count_3 += 1
        else:
            count += 1
    pickle.dump(derived_del, open(in_dir + "derived_del_by_flanking_dict_chr" + str(use_chr_num) + "_" + Ape_name + ".p","wb"))
    pickle.dump(derived_ins, open(in_dir + "derived_ins_by_flanking_dict_chr" + str(use_chr_num) + "_" + Ape_name + ".p","wb"))

    print(len(derived_ins),len(derived_del))
    print(count_1,count_2,count_3)


def Evalute_derived_SV_by_flankingregion_global_align(in_dir,Ape_ref_list,chr_num):
    for Ape_ref_file in Ape_ref_list:
        Ape_name = Ape_ref_file.split("/")[-1].split(".")[0]
        input_bam = in_dir + "validated_SV_for_flankingseq_chr" + str(chr_num) + "_" + Ape_name + "_sorted.bam"
        Evaluate_derived_sv(input_bam,in_dir,chr_num,Ape_name)

