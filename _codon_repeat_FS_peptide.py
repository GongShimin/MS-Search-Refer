#Construction of three codon repeats for
#each position of the transframe peptide reference
fr = open("CDS_species.fa","r")       #Use of CDS FASTA files
fw = open("reference_codon.fa","a")   #Output FASTA file of FS peptide reference
dic_cds={}
codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }
from collections import defaultdict
import re
def cut_text(text,lenth): 
    textArr = re.findall('.{'+str(lenth)+'}', text) 
    return textArr
for line in fr:
    if ">" in line:
        name = line.strip()
    else:
        seq = line.strip()[100:]
        if "codonrepeat" in seq :    #eg.TACTACTAC
            dic_cds[name]=seq  

dic_num  = defaultdict(list)
for k1,v1 in dic_cds.items():
        matches = re.finditer("codonrepeat", v1)  
        for match in matches:  
            l=match.start()
            if l%3==0:               #Ensure that in 0frame
                dic_num[k1].append(l)

for k,v in dic_num.items():
    for u in v:
        line1 = dic_cds[k]
        x = int(u)
        for count in range(3):
            start = count*3+x
            normal=line1[0:start]
            seq = line1[start-1:]
            seqence1=cut_text(normal,3)
            new_list = [codontable[a] for a in seqence1]
            aaseq="".join(new_list)
            seqence2=cut_text(seq,3)
            new_list2 = [codontable[a] for a in seqence2]
            aaseq2="".join(new_list2).split("*")[0]
            k1 = "K"+aaseq.split("K")[-1]
            r1 = "R"+aaseq.split("R")[-1]
            if len(k1) < len(r1):
                seq3 = k1+aaseq2
            else:
                seq3 = r1+aaseq2
            if len(seq3)<7:
                continue
            else:
                fw.write(k+":"+str(int(start/3))+":"+"CODON_MINUS"+"\n"+seq3+"\n")  

for k,v in dic_num.items():
    for u in v:
        line1 = dic_cds[k]
        x = int(u)
        for count in range(3):
            start = count*3+x
            normal=line1[0:start]
            seq = line1[start+1:]
            seqence1=cut_text(normal,3)
            new_list = [codontable[a] for a in seqence1]
            aaseq="".join(new_list)
            seqence2=cut_text(seq,3)
            new_list2 = [codontable[a] for a in seqence2]
            aaseq2="".join(new_list2).split("*")[0]
            k1 = "K"+aaseq.split("K")[-1]
            r1 = "R"+aaseq.split("R")[-1]
            if len(k1) < len(r1):
                seq3 = k1+aaseq2
            else:
                seq3 = r1+aaseq2
            if len(seq3)<7:
                continue
            else:
                fw.write(k+":"+str(int(start/3))+":"+"CODON_PLUS"+"\n"+seq3+"\n")
                