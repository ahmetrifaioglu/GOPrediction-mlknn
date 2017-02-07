import sys
import os
import subprocess
#/Users/trman/Dropbox/Academic/METU/PhD/PhDProject/VirEnvDrugProject/ccDeepFiles/GOTERMFiles

category = sys.argv[1]
path ="./.."
blastp = "./../bin/ncbi-blast-2.6.0+/bin/blastp"
mlknn_files ="./../mlknn"
blastdb_path = "%s/%s_mlknn_training.blastdb" %(mlknn_files,category)
training_fasta_path = "%s/%s_mlknn_training.fasta" %(mlknn_files,category)
test_fasta_path = "%s/%s_mlknn_test.fasta" %(mlknn_files,category)
fl_g_terms = open("%s/%sDeepFiles/GOTERMFiles/GOterms.txt" %(path,category),"r")
lst_fl_go_terms = fl_g_terms.read().split("\n")
fl_g_terms.close()

if "" in lst_fl_go_terms:
    lst_fl_go_terms.remove("")

lst_go_terms = []
for line in lst_fl_go_terms:
    go,count = line.split("\t")
    lst_go_terms.append(go)

#keys are trained GOs values are proteins associated with them
trained_gos_dict = dict()
#keys are proteins and values are GO terms annotations of them
train_prot_dict = dict()
#count dict is a dictionary that holds number of proteins that are associated with GO term has
count_dict = dict()
all_prots_set= set()
for go in lst_go_terms:
    fl_annots = open("%s/%sDeepFiles/Annots/train_%s_pos.ids" %(path,category,go),"r")
    lst_annots = fl_annots.read().split("\n")
    fl_annots.close()
    if "" in lst_annots:
        lst_annots.remove("")
    all_prots_set =set(lst_annots) | all_prots_set
    count_dict[go] = len(lst_annots)

    for prot in lst_annots:
        try:
            trained_gos_dict[go].add(prot)
        except:
            trained_gos_dict[go] = set()
            trained_gos_dict[go].add(prot)
        try:
            train_prot_dict[prot].add(go)
        except:
            train_prot_dict[prot] = set()
            train_prot_dict[prot].add(prot)
            
number_of_training_prots = len(all_prots_set)
#print training prot ids
#for prot in all_prots_set:
#    print(prot)

s=1.0
ph1l_dict = dict()
ph0l_dict = dict()
#computing prior probabilities
for g_term in lst_go_terms:
    #in the following for loop we just count number of proteins associated withe all trained GO terms.
    ph1l_dict[g_term] = (s+float(count_dict[g_term]))/(s*2+float(number_of_training_prots))
    ph0l_dict[g_term] = 1- ph1l_dict[g_term]

#subprocess.call([blastp,"-query",training_fasta_path,"-db",blastdb_path,"-outfmt","7","-out", "./../mlknn/%s_blast_train.out" %(category),"-evalue","50","-num_threads","4"])
k=20
training_blast_results=open("%s/%s_blast_train.out" %(mlknn_files,category),"r")
lst_training_blast_results= training_blast_results.read().split("\n")
training_blast_results.close()
lst_training_blast_results.remove("")
isDash=True
prot_train_knn_dict = dict()
#identify kNNs for training samples
for line in lst_training_blast_results:
    #print(line)
    if line.startswith("#"):
        isDash=True
    else:
        #Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        fields=line.split("\t")
        query_prot_id=fields[0].split("|")[1]
        target_prot_id = fields[1].split("|")[1]
        score = fields[11]
        #print(query_prot_id,target_prot_id,score)
        try:
            if len(prot_train_knn_dict[query_prot_id])<k:
                prot_train_knn_dict[query_prot_id].append([target_prot_id,float(score)])
        except:
            prot_train_knn_dict[query_prot_id] = []
            prot_train_knn_dict[query_prot_id].append([target_prot_id,float(score)])
prots_no_neighbours =set()
#c_dict is a dictinoary where keys are numbers between 0 and k and values are number of proteins whose neighbours are annotated by exactly key times
c_dict = dict()
cnot_dict = dict()
for go in lst_go_terms:
    for j in range(k+1):
        c_dict[j] = 0
        cnot_dict[j] = 0
    for prot in train_prot_dict.keys():
        #delta is number of neighbors that are annotated by the corresponding GO term
        delta = 0
        try:
            for neigh in prot_train_knn_dict[prot]:
                #print(neigh) 
                if len(train_prot_dict[neigh[0]]) and go in train_prot_dict[neigh[0]]:
                    delta +=1
        except:
            prots_no_neighbours.add(prot)
            pass
        if go in train_prot_dict[prot]:
            c_dict[delta] += 1
        else:
            cnot_dict[delta] +=1


#print(prots_no_neighbours)
#for key in prot_train_knn_dict.keys():
#    if len(prot_train_knn_dict[key])<k:
#        print(key,len(prot_train_knn_dict[key]))
