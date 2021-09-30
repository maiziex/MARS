from collections import defaultdict
from Bio import SeqIO
import glob
def longestConsecutive(num):
        # write your code here
        if num is None or len(num) == 0:
            return 0
        m = {}
        res = 0
        for i in num:
            if i not in m:
                l = 0
                r = 0
                if i - 1 in m:
                    l = m[i - 1]
                if i + 1 in m:
                    r = m[i + 1]
                m[i] = 1 + r + l
                m[i + r] = 1 + r + l
                m[i - l] = 1 + r + l
                res = max(res, m[i])
        return res
def generate_MSA_geno(MSAfile):
    hap_seq=[]
    for record in SeqIO.parse(MSAfile, "fasta"):
       # print(MSAfile)
        deletion=0
        insertion=0
        reference=0
        nameid=record.id
        seq=str(record.seq)
        if 'human_ref' in nameid:
            ref_seq=seq
        else:
            hap_seq.append(seq)
    for index in range(len(hap_seq)):
        ref_dash_list=[]
        hap_dash_list=[]
        for i in range(len(hap_seq[index])):
            if ref_seq[i]=='-' and hap_seq[index][i]!='-':
                ref_dash_list.append(i)
            elif ref_seq[i]!='-' and hap_seq[index][i]=='-':
                hap_dash_list.append(i)
        max_ref_dash=longestConsecutive(ref_dash_list)
        max_hap_dash=longestConsecutive(hap_dash_list)
        if max_ref_dash>=20:
            insertion+=1
        if max_hap_dash>=20:
            deletion+=1
        if max_ref_dash<20 and max_hap_dash<20:
            reference+=1
    return deletion,insertion,reference
def MSA_geno_all(MSA_path):
    MSA_geno=defaultdict(list)
    fileall=glob.glob(MSA_path+'/*.afa')
    for name in fileall:
        #print(name)
        keyid=name.split('/')[-1].split('.')[0]
        try:
            deletion,insertion,reference=generate_MSA_geno(name)
            #print(deletion)
            MSA_geno[keyid].append(deletion)
            MSA_geno[keyid].append(insertion)
            MSA_geno[keyid].append(reference)
        except:
            #print(name)
            MSA_geno[keyid].append('NA')
            MSA_geno[keyid].append('NA')
            MSA_geno[keyid].append('NA')
    return MSA_geno

def updatefile(original_file,new_file,MSA_path):
    MSA_geno=MSA_geno_all(MSA_path)
    #print(MSA_geno)
    outfile=open(new_file,'w')
    infile=open(original_file)
    count=0
    for line in infile:
        if count==0:
            outfile.write(line.strip('\n')+"\t#ins(MSA)\t#del(MSA)\t#ref_allele(MSA)\n")
        elif count!=0:
            outfile.write(line.strip('\n')+'\t')
            A=line.strip('\n').split('\t')
            if A[0] in MSA_geno.keys():
                del_ins_ref=MSA_geno[A[0]]
                outfile.write(str(del_ins_ref[1])+'\t')
                outfile.write(str(del_ins_ref[0])+'\t')
                outfile.write(str(del_ins_ref[2])+'\n')
            else:
                outfile.write('-\t')
                outfile.write('-\t')
                outfile.write('-\n')
        count+=1
    infile.close()
    outfile.close()

#updatefile("/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS_hg38/Final_tables/SV_msa_table_chr1.txt",'test','/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS_hg38/MSA_SV_results_chr1/MSA_SV_files')

#updatefile("/oak/stanford/groups/arend/Eric/aquallia_results/CCS_reads/svviz/chr21/SV_MSA_table_chr21_template","test1","/oak/stanford/groups/arend/Xin/AncestralProj/Results_MSA_MARS_hg38/MSA_SV_results_chr21/MSA_SV_files")
