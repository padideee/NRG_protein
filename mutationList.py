import argparse
import subprocess
import sys
import re
import numpy
import pandas as pd
from pandas import DataFrame
from prettytable import PrettyTable
import itertools
import time
import progressbar


AMINO_DICT = {"ALA" : "A", "ARG" : "R", "ASN" : "N","ASP":"D","CYS":"C", "GLN":"Q",
"GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
"PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
AA_DICT=dict((v,k) for k,v in AMINO_DICT.items())

#the groups contain "p:polar" "h:hydrophobe" "-:negative charge" "+:positive charge" "l:large size" "s:small size" "a:aromatic"
groups={"p":["N","Q","R","K","D","E","H"] ,"h":["A","I","L","M","F","V","W","C"],
"-":["D","E"],"+":["K","R","H"],"l":["W","F","H","Y","R","K"],
"s":["G","A","D","C","S"],"a":["F","Y","W","H"]}
#dictionary of the closest substitution in terms of isoelectric point
dict_i={"A":"L","L":"I","R":"K","K":"R","N":"F","M":"S","D":"E","F":"N","C":"N","P":"I","Q":"Y",
"S":"Y","E":"D","T":"Y","Y":"T","G":"V","W":"V","H":"P","I":"L","V":"G"}
#dictionary of the closest substitution in terms of Chemical distance
dict_c={"A":"P","R":"K","N":"D","D":"N","C":"S","Q":"H","E":"Q","G":"P","H":"Q","I":"L","L":"I",
"K":"A","M":"V","F":"I","P":"G","S":"R","T":"P","W":"Y","Y":"I","V":"M"}
#dictionary of the closest substitution in terms of 1 PAM evolutionary distance
dict_m={"A":"S","R":"K","N":"S","D":"E","C":"S","Q":"E","E":"D","G":"A","H":"Q","I":"V","L":"V",
"K":"R","M":"L","F":"Y","P":"A","S":"A","T":"S","W":"R","Y":"F","V":"I"}

def clean(pdb_nc):
    #call the perl script to clean the pdb
    cmd=["perl", "clean_pdb.pl", "-i", pdb_nc, "-t","pdb_c.pdb"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    p.wait()

def repair():

    cmd=["./foldx", "--command=RepairPDB", "--pdb=pdb_c.pdb"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    p.wait()

def mutate_position(pos, toRead): #the raw position input, could contain requested mutations in form of AA or specificity names
	pos_input_list=pos.split(",")
	mutate_pos_list=[]
	specificity_dict={}
	for i in range(len(pos_input_list)):
		if pos_input_list[i]== "AUTO":
			#based on b_factors predict the most flexible sites
			cmd=["./build_encom","-i", toRead,"-cov","wt.cov", "-o", "wt.eigen"]
			p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			p.wait()

			lines=[]
			c_lines=[]
			with open("wt.cov") as f:
				for l in f:
					lines.append(l.rstrip())

			for i in range(len(lines)):
				if lines[i][0]=="C":
					c_lines.append(lines[i].split(" ")[5])

			val_ind=sorted( [(x,i) for (i,x) in enumerate(c_lines)], reverse=True )[:2]
			for i in range(len(val_ind)):
				mutate_pos_list.append(val_ind[i][1] +1)
				specificity_dict[val_ind[i][1] +1]="cm"

		elif "(" not in pos_input_list[i]:
			if "-" in pos_input_list[i]:
				dash_pos=pos_input_list[i].index("-")
				start_pos=int(pos_input_list[i][0:dash_pos])
				end_pos=int(pos_input_list[i][dash_pos+1:])
				while start_pos != end_pos:
					mutate_pos_list.append(start_pos)
					start_pos+=1
			elif "-" not in pos_input_list[i]:
				mutate_pos_list.append(int(pos_input_list[i]))
		elif "(" in pos_input_list[i]:
			paranthesis_pos=pos_input_list[i].index("(")
			specificity_str=pos_input_list[i][paranthesis_pos+1:-1]
			if "-" in pos_input_list[i]:
				dash_pos=pos_input_list[i].index("-")
				start_pos=int(pos_input_list[i][0:dash_pos])
				end_pos=int(pos_input_list[i][dash_pos+1:paranthesis_pos])
				while start_pos != end_pos:
					mutate_pos_list.append(start_pos)
					specificity_dict[start_pos]=specificity_str
					start_pos+=1
			elif "-" not in pos_input_list[i]:
				x=int(pos_input_list[i][0:paranthesis_pos])
				mutate_pos_list.append(x)
				specificity_dict[x]=specificity_str


	return(mutate_pos_list,specificity_dict)

def all_single_mut(toRead,mutate_pos,mutate_to): #build the individual_list.txt to write all the possible single mutations for the requested position
    f= open(toRead, "r")

    AminoAcid_list=[]
    chain_list=[]
    pos_list=[]

    for x in f:
        AminoAcid_list.append(x[17:21].rstrip())  #read the column for all AminoAcid name
        chain_list.append(x[21].rstrip()) #read the column for chainID
        pos_list.append(int(x[23:26].rstrip())) #read the columns for all position number

    i_pos = 0
    for x in range(len(pos_list)):
        if pos_list[x] == int(mutate_pos):
            i_pos=x
    i=int(i_pos)
    current_AminoAcid= AminoAcid_list[i] #the AminoAcid for the requested position
    current_chain= chain_list[i]


    Mut_info_list=[]
    wanted_AA=[]
    for i in range(len(mutate_to)):
        if mutate_to[i] in groups:
            wanted_AA.extend(groups[mutate_to[i]])
        elif mutate_to[i] == "i":
            wanted_AA.append(dict_i[AMINO_DICT[current_AminoAcid]])
        elif mutate_to[i] == "c":
            wanted_AA.append(dict_c[AMINO_DICT[current_AminoAcid]])
        elif mutate_to[i] == "m":
            wanted_AA.append(dict_m[AMINO_DICT[current_AminoAcid]])
        elif mutate_to[i] == "e":
            wanted_AA.append(dict_m[AMINO_DICT[current_AminoAcid]])
        else:
            wanted_AA.append(mutate_to[i])

    for i in range(len(wanted_AA)):
        if wanted_AA[i]== "x":
            if AMINO_DICT[current_AminoAcid] != wanted_AA[i]:
                if wanted_AA[i] != "x":
                    Mut_info= AMINO_DICT[current_AminoAcid]+current_chain+str(mutate_pos)+ wanted_AA[i]
                    Mut_info_list.append(Mut_info)
        else:
            if AMINO_DICT[current_AminoAcid] != wanted_AA[i]:
                Mut_info= AMINO_DICT[current_AminoAcid]+current_chain+str(mutate_pos)+ wanted_AA[i]
                Mut_info_list.append(Mut_info)
    Mut_info_list.append(AMINO_DICT[current_AminoAcid]+current_chain+str(mutate_pos)+ AMINO_DICT[current_AminoAcid])

    return (Mut_info_list,mutate_pos,current_chain,wanted_AA)


def file_generate(toRead,Mut_info_list,mutate_pos,current_chain,mut_to_AA):
	with open("individual_list.txt","w") as f:

		all_possibility=list(itertools.product(*Mut_info_list))

		for i in range(len(all_possibility)):
			for j in range(len(all_possibility[i])):
				f.write(all_possibility[i][j])
				if j < (len(all_possibility[i])-1):
					f.write(",")
			f.write(";\n")

	return all_possibility

def delta_foldingNRG():
	cmd=["./foldx","--command=BuildModel","--pdb=pdb_c_Repair.pdb","--mutant-file=individual_list.txt"]
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	p.wait()



def ENCoM_files(toRead,Mut_info_list,mutate_pos,current_chain,mut_to_AA):
	# create the normal modes and mode aplituted for WT
	cmd=["./build_encom","-i", toRead,"-cov","wt.cov", "-o", "wt.eigen"]
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	p.wait()

	# create the normal modes and mode amplituted for each Mut_file
	for i in range(len(Mut_info_list)):
		cmd=["./build_encom","-i", "pdb_c_Repair_"+str(i+1)+".pdb","-cov", str(i+1)+".cov", "-o", str(i+1)+".eigen"]
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		p.wait()



def compare_ENCoM(Mut_info_list,pdb_nc,data_table_txt):

	# the WT energy
	delta_energy_dict={}
	delta_energy_list=[]

	f=open("wt.cov", "r")
	lines=f.readlines()
	energy_WT=float(lines[1][7:])

	#each mutated file energy
	for	i in range(len(Mut_info_list)):
		f=open(str(i+1)+".cov", "r")
		lines=f.readlines()
		energy=float(lines[1][7:])
		delta_energy_dict[str(i+1)]=energy-energy_WT
		delta_energy_list.append(energy-energy_WT)

	with open(data_table_txt,"a+")as f:
		x=PrettyTable()
		x.add_column("name", Mut_info_list)
		x.add_column("ENCoM",delta_energy_list)


	name_DICT=dict(zip( delta_energy_list, Mut_info_list))
	delta_energy_list.sort(reverse=True)
	Mut_info_list_sorted=[]
	for i in range(len(delta_energy_list)):
		Mut_info_list_sorted.append(name_DICT[delta_energy_list[i]])

	with open(data_table_txt,"a+")as f:
		y=PrettyTable()
		y.add_column("name", Mut_info_list_sorted)
		y.add_column("ENCoM descending order",delta_energy_list)
		f.write(str(y)+"\n\n")

	return( delta_energy_dict,x)

def compare_FoldX(x,Mut_info_list,pdb_nc,data_table_txt):
	energy_dict={}
	delta_energy_list=[]
	f=open("Dif_pdb_c_Repair.fxout","r")
	lines=f.readlines()
	for i in range(len(Mut_info_list)):
		j=20
		while lines[9+i][j] != "\t":
			j=j+1
		energy=lines[9+i][19:j]
		energy=re.sub("\s+","",energy)
		energy_dict[str(i+1)]=float(energy)
		delta_energy_list.append(energy)

	with open (data_table_txt,"a+") as f:
		x.add_column("FoldX",delta_energy_list)


	name_DICT=dict(zip( delta_energy_list, Mut_info_list))
	delta_energy_list.sort()
	Mut_info_list_sorted=[]

	for i in range(len(delta_energy_list)):
		Mut_info_list_sorted.append(name_DICT[delta_energy_list[i]])

	with open(data_table_txt,"a+")as f:
		y=PrettyTable()
		y.add_column("name", Mut_info_list_sorted)
		y.add_column("FoldX ascending order",delta_energy_list)
		f.write(str(y)+"\n\n")

	return(energy_dict,x)


def protocol_coef(ENCoM,FoldX,x, Mut_info_list, data_table_txt):

	ENCoM_list=[]
	FoldX_list=[]
	with_coef=[]

	for i in range(len(Mut_info_list)):
		ENCoM_list.append(ENCoM[str(i+1)])
		FoldX_list.append(FoldX[str(i+1)])

	for i in range(len(Mut_info_list)):
		with_coef.append(-1.12*ENCoM_list[i]+0.38*FoldX_list[i])

	with open (data_table_txt,"a+") as f:
		x.add_column("-1.12*ENCoM + 0.38*FoldX",with_coef)
		f.write(str(x)+"\n\n")

		name_DICT=dict(zip( with_coef, Mut_info_list))
		with_coef.sort()
		Mut_info_list_sorted=[]

		for i in range(len(with_coef)):
			Mut_info_list_sorted.append(name_DICT[with_coef[i]])

		y=PrettyTable()
		y.add_column("name", Mut_info_list_sorted)
		y.add_column("predicted ddG ascending order",with_coef)
		f.write(str(y)+"\n\n")
		print(y)


def main():

	with open("help.txt") as f:
		lines=f.read()

	parser= argparse.ArgumentParser(description="the arguments.", add_help=False)
	parser.add_argument("-h","--help", action="help",help=print(lines))
	parser.add_argument("-i","--pdb_nc_input", action="store")
	parser.add_argument("-r","--pos", action="store")
	parser.add_argument("-o","--output_name", action="store")
	parser.add_argument("-v","--verbose", action="store")
	args=parser.parse_args()


	(mutate_pos_list, mutate_to_dict)=mutate_position(args.pos, "pdb_c.pdb")
	pdb_c= clean(args.pdb_nc_input)
	pdb_repaired=repair()


	individual_list=[]
	chain=[]
	mut_to_AA=[]
	for i in range(len(mutate_pos_list)):
		if mutate_pos_list[i] in mutate_to_dict:
			(Mut_info_list,pos,current_chain,wanted_AA)=all_single_mut("pdb_c.pdb",mutate_pos_list[i],mutate_to_dict[mutate_pos_list[i]])
		else:
			mutate_to="ARNDCQEGHILKMFPSTWYV"
			(Mut_info_list,pos,current_chain,wanted_AA)=all_single_mut("pdb_c.pdb",mutate_pos_list[i],mutate_to)
		individual_list.append(Mut_info_list)
		chain.append(current_chain)
		mut_to_AA.append(wanted_AA)

	id_list=file_generate("pdb_c.pdb",individual_list,mutate_pos_list,chain,mut_to_AA)
	with open("individual_list.txt") as f:
		lines=f.readlines()
		text_list=[]
		for i in range(len(lines)):
			text_list.append(lines[i].rsplit())

	delta_foldingNRG()
	ENCoM_files("pdb_c.pdb",text_list,mutate_pos_list,chain,mut_to_AA)
	(ENCoM,x)=compare_ENCoM(text_list,args.pdb_nc_input, args.output_name)
	(FoldX,y)=compare_FoldX(x,text_list,args.pdb_nc_input, args.output_name)
	protocol_coef(ENCoM,FoldX,y, text_list,args.output_name)


main()
