import pdb
#pdb.set_trace()
import glob
import pickle
from collections import defaultdict


def Check_linked_SNP_for_Ancestral_SV(in_dir,SNP_dict):
    SNP_dict_ancestral = defaultdict(list)
    fasta_files_all = sorted(glob.glob(in_dir +  "*.fasta"))
    for one_file in fasta_files_all:
        sv_name = one_file.split("/")[-1].split(".fasta")[0]
        if sv_name in SNP_dict:
            SNP_dict_ancestral[sv_name] = SNP_dict[sv_name]
    print("#: " + str(len(SNP_dict_ancestral)))

    pickle.dump(SNP_dict_ancestral,open(in_dir + "SNP_dict_derived.p","wb"))
    return SNP_dict_ancestral





in_dir_msa = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS/MSA_SV_results_chr21/MSA_SV_files/"
SNP_dict = pickle.load(open(in_dir_msa + "SNP_dict.p","rb"))

in_dir = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS/MSA_SV_results_chr21/Ancestral_with_Apes_Seqs_derived_del/"
Check_linked_SNP_for_Ancestral_SV(in_dir,SNP_dict)

in_dir = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS/MSA_SV_results_chr21/Ancestral_with_Apes_Seqs_derived_ins/"
Check_linked_SNP_for_Ancestral_SV(in_dir,SNP_dict)

in_dir = "/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS/MSA_SV_results_chr21/Ancestral_with_Apes_Seqs_derived_complex/"
Check_linked_SNP_for_Ancestral_SV(in_dir,SNP_dict)
