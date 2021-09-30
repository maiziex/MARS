import pdb
#pdb.set_trace()
import pickle
from collections import defaultdict
from subprocess import Popen
import os
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 

def get_max(sv_size_list):
    all_list = []
    for one in sv_size_list:
        sv_size = int(one)
        all_list.append(sv_size)
    return max(all_list)


def Extract_validated_SV(input_file,output_file,chr_num):
    f = open(input_file,"r")
    val_sv = defaultdict(list)
    for line in f:
        data = line.rsplit()
        cur_chr_num = data[0]
        if cur_chr_num == "chr" + str(chr_num):
            sv_type = data[1]
            _start = int(data[2])
            _end = int(data[3])

            if "del" in sv_type:
                sv_size = _end - _start
            elif "ins" in sv_type:
                sv_size_list = data[6].split(":")
                sv_size = get_max(sv_size_list)
                
            cur_name = data[0] + "_" + data[2] + "_"  + data[3] + "_" + data[1]
            _start_2 = _start - 500
            _end_2 = _end + 500
            val_sv[cur_name] = [_start, _end, _start_2, _end_2, sv_size]
    pickle.dump(val_sv, open(output_file,"wb"))


def write_validated_SV_fasta(input_file,ref_chr_dict,fw,chr_num):
    val_sv = pickle.load(open(input_file,"rb"))
    for key, val in val_sv.items():
        #print(val)
        _start = int(val[0])
        _end = int(val[1])
        _start_2 = int(val[2])
        _end_2 = int(val[3])
        sv_size = int(val[4])
        seq_1_flankingseq = ref_chr_dict[_start_2:_start]
        seq_2_flankingseq = ref_chr_dict[_end:_end_2]
        fw.writelines(">1" +  "_" + key + "_" + str(sv_size) + "\n")  
        fw.writelines(seq_1_flankingseq + "\n")
        fw.writelines(">2" + "_" + key + "_" + str(sv_size) + "\n")  
        fw.writelines(seq_2_flankingseq + "\n")

    return fw


def Extract_flankingregion_for_SV_to_align_to_Apes(in_dir,ref_dir,Ape_ref_list,chr_num):
    output_fasta = in_dir + "validated_SV_for_flankingseq_chr" + str(chr_num) + ".fasta" 
    output_pickle = in_dir + "validated_SV_for_flankingseq_chr" + str(chr_num) + ".p" 
    fw = open(output_fasta,"w")
    merged_sv_file = in_dir + "chr" + str(chr_num) + "_merged_sv.txt"
    ref_chr_file = ref_dir + "ref_seq_chr" + str(chr_num) + ".p"
    ref_chr_dict = pickle.load(open(ref_chr_file,"rb"))
    Extract_validated_SV(merged_sv_file,output_pickle,chr_num)
    fw = write_validated_SV_fasta(output_pickle,ref_chr_dict,fw,chr_num)
    fw.close()
   
    for Ape_ref_file in Ape_ref_list:
        Ape_name = Ape_ref_file.split("/")[-1].split(".")[0]
        output_sam = in_dir + "validated_SV_for_flankingseq_chr" + str(chr_num) + "_" + Ape_name + ".sam" 
        output_bam = in_dir + "validated_SV_for_flankingseq_chr" + str(chr_num) + "_" + Ape_name + ".bam" 
        output_bam_sorted = in_dir + "validated_SV_for_flankingseq_chr" + str(chr_num) + "_" + Ape_name + "_sorted.bam" 
        try:
            _cmd = "minimap2 -a " + Ape_ref_file + " " + output_fasta + "  >  " + output_sam
            Popen(_cmd,shell=True).wait()
        except:
            _cmd = code_path + "minimap2/" + "minimap2 -a " + Ape_ref_file + " " + output_fasta + "  >  " + output_sam
            Popen(_cmd,shell=True).wait()

        try: 
            _cmd = "samtools view -Sb  " + output_sam  + " > " +  output_bam
            Popen(_cmd,shell=True).wait()
        except:
            _cmd = code_path + "samtools/" + "samtools view -Sb  " + output_sam  + " > " +  output_bam
            Popen(_cmd,shell=True).wait()

        try:
            _cmd = "samtools sort " + output_bam + " -o " + output_bam_sorted
            Popen(_cmd,shell=True).wait()
        except:
            _cmd = code_path + "samtools/" + "samtools sort " + output_bam + " -o " + output_bam_sorted
            Popen(_cmd,shell=True).wait()

        try:
            _cmd = "samtools index " + output_bam_sorted
            Popen(_cmd,shell=True).wait()
        except:
            _cmd = code_path + "samtools/" + "samtools index " + output_bam_sorted
            Popen(_cmd,shell=True).wait()


        _cmd = "rm " + output_sam
        Popen(_cmd,shell=True).wait()
        
        _cmd = "rm " + output_bam
        Popen(_cmd,shell=True).wait()







