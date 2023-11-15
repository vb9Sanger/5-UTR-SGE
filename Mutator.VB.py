#Vanessa's Code to Create Variant Oligo Library for 5'UTR Pilot Study 
 
#region-specific changes possible:  pre-specified domains in input file 
#functions specify regions variants are applied

import sys
import re
import os
import random

Targetons_path_n_file = "/Users/vb9/Library/CloudStorage/GoogleDrive-vb9@sanger.ac.uk/My Drive/Python files/VBinput.txt"

OUTPUT_path = "/Users/vb9/Library/CloudStorage/GoogleDrive-vb9@sanger.ac.uk/My Drive/Python files"
#=========================================
def file_extractor():#extracts list of loci to mutate from text file
	found = False 
	ENTRIES = []#list of oligos from tg file
	with open(Targetons_path_n_file,"r") as f: #open as file object 
		for line_read in f.readlines():
			locus = re.split(r"\t",line_read)#split line on tabs
			ENTRIES.append({"tg_name":locus[0],"Primer1":locus[1],"pre-tg":locus[2],"Primer2":locus[3],"5'_intronic":locus[4],"SS_acceptor":locus[5],"UTR":locus[6],"CDS":locus[7],"SS_donor":locus[8],"3'_intronic":locus[9]}) #create entry of tg and features as a dictionary - corresponding to columns in VBinput.txt
		if ENTRIES != []:
			print("Successful extraction of loci")
			return(ENTRIES)
		else:
			print("Didn't find anything")
			return ([])
#===========================================			
def OUTPUT(OLIGO_LIST,name):#saves latest list to file. 
	global first_entry
	if first_entry: #boolean 
		mode = "w" #mode determines new file vs. appending existing file.
		first_entry = False
	else:
		mode = "a"
	with open(OUTPUT_path+"/"+name+".txt",mode) as f: #opens file in correct mode with .txt appended to name
		for oligo in OLIGO_LIST:
			f.write(oligo[0]+"\t"+oligo[1]+"\n") #writes tab delimited oligo info and newline "\n"
	return ([])
#=========================================
def MakeCoords(locusRegion):
	coords = list(locusRegion.split(","))
	coords[0] = int(coords[0])
	coords[1] = int(coords[1])
	
	return coords
#=========================================	
#=========================================				
def dels(locus,delsize,overlap,startPoint,endPoint,region): #create systematic set of deletions, only deleting full n regions - leaves regions < n
	pre_tg = locus["pre-tg"] #making local oligo variable as basis of modification
	library = [] #clear library as an empty list
	BaseCount = endPoint+1 #length of mutable pre targeton
	tg_name = locus["tg_name"]#used the dictionary
	step = delsize-overlap #step size to give required overlap
	for cursor in range (startPoint,BaseCount-(delsize-1),step): #run through oligo in step of size "step" 
		new_oligo = locus["Primer1"]+pre_tg[0:cursor]+"/"+pre_tg[cursor+delsize:]+locus["Primer2"]#assemble full tg using primers, slices of pre targeton sequence.  currently incudes backslash for testing output
		new_name = tg_name+"_"+region+"_"+str(cursor)+"_"+str(delsize)+"del"
		entry = [new_name,new_oligo] #latest oligo as a 2 item list
		library.append(entry)#add latest oligo to library list
	OUTPUT(library,tg_name)#save library to file
#==========================================
def INSERTION_SCAN(locus,insert,step,startPoint,endPoint,region): #create systematic set of insertions
	pre_tg = locus["pre-tg"]
	library = []
	BaseCount = endPoint+1
	tg_name = locus["tg_name"]

	names = {"gtcagacgagtcatcaa": "VANESSA", "atgattcatgagtgg": "MATTHEW", "tcagagatcatataagac": "SEBASTIAN", "gagctcatatcaagagtcat": "ELIZABETH"} #dictionary of names and corresponding sequence 
	insertions = {"tga": "STOP", "atgtga": "STARTSTOP", "gccatggc": "StrongKozak", "tttatgct": "WeakKozak", "ggggttgggggtgggtgggg": "RG4", "gggggaggttcgcctccccc": "StemLoop", "gcagtaagtaatacatgtaa": "SSD", "tcaccattatcgtttcagac": "SSA"} #dictionary of insertions and corresponding identifier 
	for cursor in range (startPoint,BaseCount+1,step):
		if insert in names: #separated simply for naming purposes 
			new_oligo = locus["Primer1"]+pre_tg[0:cursor]+insert+pre_tg[cursor:]+locus["Primer2"]
			new_name = tg_name+"_"+region+"_"+str(cursor)+"_"+insert+"_"+names[insert]+"_ins"
			entry = [new_name,new_oligo] #latest oligo as a 2 item list
			library.append(entry)#add latest oligo to library list
		elif insert in insertions:
			new_oligo = locus["Primer1"]+pre_tg[0:cursor]+insert+pre_tg[cursor:]+locus["Primer2"]
			new_name = tg_name+"_"+region+"_"+str(cursor)+"_"+insert+"_"+insertions[insert]+"_ins"
			entry = [new_name,new_oligo] #latest oligo as a 2 item list
			library.append(entry)#add latest oligo to library list
		else: #shouldn't be needed 
			new_oligo = locus["Primer1"]+pre_tg[0:cursor]+insert+pre_tg[cursor:]+locus["Primer2"]
			new_name = tg_name+"_"+region+"_"+str(cursor)+"_"+insert+"_ins"
			entry = [new_name,new_oligo] #latest oligo as a 2 item list
			library.append(entry)#add latest oligo to library list
	OUTPUT(library,tg_name)#save library to file
#=========================================
def SNV_SCAN(locus,startPoint,endPoint,region): #create systematic set of SNVs
	pre_tg = locus["pre-tg"]
	library = []
	BaseCount = endPoint+1
	bases = ["a","t","g","c"]
	tg_name = locus["tg_name"]
	for cursor in range (startPoint,BaseCount): #(start, stop (not included))
 		for base in bases:
 			if base != pre_tg[cursor].lower(): #check what base cursor is on 
 				if cursor != 0 and cursor != (BaseCount-1): #for bases not at beginning or end 
 					new_oligo = locus["Primer1"]+pre_tg[0:cursor]+base+pre_tg[cursor+1:]+locus["Primer2"] #assemble full tg using primers, slices of pre targeton sequence and alt base
 				elif cursor == 0:  #for first base 
 					new_oligo = locus["Primer1"]+base+pre_tg[cursor+1:]+locus["Primer2"] 
 				else: #for last base 
 					new_oligo = locus["Primer1"]+pre_tg[0:cursor]+base+pre_tg[cursor+1:]+locus["Primer2"]
 				new_name = tg_name+"_"+region+"_"+str(cursor)+"_"+base.upper()+"_SNV"#assemble an oligo name based on location and base change, with "_SNV" at the end		
 				entry = [new_name,new_oligo] #latest oligo as a 2 item list
 				library.append(entry)#add latest oligo to library list
	OUTPUT(library,tg_name)#save library to file
#=========================================
def uORF_SCAN(locus,START,STOP,step,startPoint,endPoint,region): #create set of uORF insertions- insert strong kozak sequence and a stop codon into endogenous sequence (20bp apart)
	pre_tg = locus["pre-tg"] #making local oligo variable as basis of modification
	library = [] #clear library as an empty list
	BaseCount = endPoint+1 #length of mutable pre targeton
	tg_name = locus["tg_name"]#used the dictionary
	for cursor in range (startPoint,BaseCount-19,step): #adjusted so that STOP insertion doesn't go beyond UTR region (BaseCount +1 -20)
		new_oligo = locus["Primer1"]+pre_tg[0:cursor]+START+pre_tg[cursor:cursor+20]+STOP+pre_tg[cursor+20:]+locus["Primer2"]#assemble full tg using primers, slices of pre targeton sequence.  currently incudes backslash for testing output
		new_name = tg_name+"_"+str(cursor)+"_"+"uORF"+"_ins"
		entry = [new_name,new_oligo] #latest oligo as a 2 item list
		library.append(entry)#add latest oligo to library list
	OUTPUT(library,tg_name)#save library to file
#===========================================
def SV40_large(locus,insert):#create oligos incorporating 270bp SV40 sequence 
	pre_tg = locus["pre-tg"]
	library = []
	tg_name = locus["tg_name"]
	SV40_L = {"tacgtagatccagacatgataagatacattgatgagtttggacaaaccacaactagaatgcagtgaaaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaaccattataagctgcaataaacaagttaacaacaacaattgcattcattttatgtttcaggttcagggggaggtgtgggaggttttttaattc": "LargeSV40_F", "gaattaaaaaacctcccacacctccccctgaacctgaaacataaaatgaatgcaattgttgttgttaacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatcttatcatgtctggatctacgta": "LargeSV40_R"}
	if locus["tg_name"] == "ENST00000282516.13_sg1" or "ENST00000282516.13_sg2": #NIPBL
		new_oligo = locus["Primer1"]+pre_tg[0:21]+insert+locus["Primer2"]
	if locus["tg_name"] == "ENST00000264010.10_sg3" or "ENST00000264010.10_sg4": #CTCF
		new_oligo = locus["Primer1"]+insert+pre_tg[207:]+locus["Primer2"]
	if locus["tg_name"] == "ENST00000426263.10_sg5" or "ENST00000426263.10_sg6": #SLC2A1
		new_oligo = locus["Primer1"]+insert+pre_tg[207:]+locus["Primer2"]
	if locus["tg_name"] == "ENST00000644876.2_sg7" or "ENST00000644876.2_sg8": #DDX3X
		new_oligo = locus["Primer1"]+insert+pre_tg[207:]+locus["Primer2"]
	if locus["tg_name"] == "ENST00000356577.10_sg9" or "ENST00000356577.10_sg10": #SON
		new_oligo = locus["Primer1"]+insert+pre_tg[207:]+locus["Primer2"]
	new_name = tg_name+"_"+SV40_L[insert]+"_replace"
	entry = [new_name,new_oligo] #latest oligo as a 2 item list
	library.append(entry)#add latest oligo to library list
	OUTPUT(library,tg_name)#save library to file
#=======
def SV40_small(locus,insert): #create oligos incorporating 120bp SV40 sequence 
	pre_tg = locus["pre-tg"]
	library = []
	tg_name = locus["tg_name"]
	SV40_S = {"taagatacattgatgagtttggacaaaccacaactagaatgcagtgaaaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaaccattataagctgcaataaacaagtt": "SmallSV40_F", "aacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatctta": "SmallSV40_R"}
	if locus["tg_name"] == "ENST00000282516.13_sg1" or "ENST00000282516.13_sg2": #NIPBL
		new_oligo = locus["Primer1"]+pre_tg[0:106]+insert+locus["Primer2"]
	if locus["tg_name"] == "ENST00000264010.10_sg3" or "ENST00000264010.10_sg4": #CTCF
		new_oligo = locus["Primer1"]+insert+pre_tg[122:]+locus["Primer2"]
	if locus["tg_name"] == "ENST00000426263.10_sg5" or "ENST00000426263.10_sg6": #SLC2A1
		new_oligo = locus["Primer1"]+pre_tg[0:40]+insert+pre_tg[163:]+locus["Primer2"]
	if locus["tg_name"] == "ENST00000644876.2_sg7" or "ENST00000644876.2_sg8": #DDX3X
		new_oligo = locus["Primer1"]+insert+pre_tg[122:]+locus["Primer2"]
	if locus["tg_name"] == "ENST00000356577.10_sg9" or "ENST00000356577.10_sg10": #SON
		new_oligo = locus["Primer1"]+insert+pre_tg[122:]+locus["Primer2"]
	new_name = tg_name+"_"+SV40_S[insert]+"_replace"
	entry = [new_name,new_oligo] #latest oligo as a 2 item list
	library.append(entry)#add latest oligo to library list
	OUTPUT(library,tg_name)#save library to file
#===========================================
def inversion(locus): #create single oligos that invert whole pretg sequence 
	pre_tg = locus["pre-tg"]
	library = []
	tg_name = locus["tg_name"]
	new_oligo = locus["Primer1"]+pre_tg[::-1]+locus["Primer2"]
	new_name = tg_name+"_"+"pretg_inversion"
	entry = [new_name,new_oligo] #latest oligo as a 2 item list
	library.append(entry)#add latest oligo to library list
	OUTPUT(library,tg_name)#save library to file
#===========================================
def permutation(locus): #create single oligos that are a permutation of corresponding pretg sequence 
	pre_tg = locus["pre-tg"]
	library = []
	tg_name = locus["tg_name"]
	permutation = "".join(random.sample(pre_tg, len(pre_tg))) 
	new_oligo = locus["Primer1"]+permutation+locus["Primer2"]
	new_name = tg_name+"_"+"pretg_permutation"
	entry = [new_name,new_oligo] #latest oligo as a 2 item list
	library.append(entry)#add latest oligo to library list
	OUTPUT(library,tg_name)#save library to file
#============================================
def mutify(Whole_list): #this launcher goes thru list of targetons and performs multiple mutagenesis steps
	global first_entry
	Kozak_step = 1 #frequency/spacing of START insertions with different strength kozak sequence
	TGA_step = 1 #frequency/spacing of TGA insertions 
	StartStop_step = 1 #frequency/spacing of Start/Stop insertions
	RG4_step = 15 #frequency/spacing of RG4 forming sequence insertions 
	Loop_step = 15 #frequency/spacing of stem loop forming sequence insertions 
	SSdonor_step = 10 #frequency/spacing of SS donor consensus sequence 
	SSacceptor_step = 10 #frequency/spacing of SS acceptor consensus sequence 
	names_step = 15 #frequency/spacing of Group member name insertions 
	uORF_step = 10 #frequency/spacing of uORF insertions #10bp overlap (20bp uORFs)
	for locus in Whole_list:
		first_entry = True #this flags new file creation vs append
		#=========================================
		############### DELETIONS ################
		coords = MakeCoords(locus["UTR"]) 
		if coords[1] < 200: #CTCF, DDX3X and SON
			dels(locus,3,2,coords[0],coords[1],"UTR") #3bp deletions -step of 1 (2bp overlap)
			if locus["SS_acceptor"] != "0": #identify is there is SS Acceptor sequence 
				coords = MakeCoords(locus["SS_acceptor"])
				dels(locus,3,2,coords[0],coords[1],"SSA") #also apply function here 
			if locus["SS_donor"] != "0":  #identify if there is SS donor sequence
				coords = MakeCoords(locus["SS_donor"]) 
				dels(locus,3,2,coords[0],coords[1],"SSD") #also apply function here 
		if coords[1] > 200: #NIPBL and SLC2A1
			dels(locus,3,1,coords[0],coords[1],"UTR") #3bp deletions (overlapping) - step of 2 (1bp overlap)
			if locus["SS_acceptor"] != "0": #identify is there is SS Acceptor sequence
				coords = MakeCoords(locus["SS_acceptor"])
				dels(locus,3,2,coords[0],coords[1],"SSA") #also apply function here
			if locus["SS_donor"] != "0": #identify if there is SS donor sequence
				coords = MakeCoords(locus["SS_donor"])
				dels(locus,3,2,coords[0],coords[1],"SSD") #also apply function here
 		#Always applied to UTR, applied to SS acceptor and/or SS donor if there is one 	
		coords = MakeCoords(locus["UTR"]) 
		dels(locus,10,5,coords[0],coords[1],"UTR") #10bp deletions (overlapping)
		#==========================================
		############### INSERTIONS ################
		coords = MakeCoords(locus["UTR"])
		if coords[1] < 200: #CTCF, DDX3X and SON
			INSERTION_SCAN(locus,"tga",TGA_step,coords[0],coords[1],"UTR")#add TGA STOP at every nth base
			if locus["CDS"] != "0": #identify if there is CDS
				coords = MakeCoords(locus["CDS"])
				INSERTION_SCAN(locus,"tga",TGA_step,coords[0]+1,coords[1],"CDS")
		if coords[1] > 200: #NIPBL and SLC2A1
			INSERTION_SCAN(locus,"tga",2,coords[0],coords[1],"UTR")#add TGA STOP at every nth base
			if locus["CDS"] != "0": #identify if there is CDS
				coords = MakeCoords(locus["CDS"])
				INSERTION_SCAN(locus,"tga",TGA_step,coords[0]+1,coords[1],"CDS")
		#Also applied to coding sequence if there is one
		coords = MakeCoords(locus["UTR"]) #default (always applied to UTR)
		if coords[1] < 200: #CTCF, DDX3X and SON
			INSERTION_SCAN(locus,"atgtga",StartStop_step,coords[0],coords[1],"UTR")#add STARTSTOP every nth base
		if coords[1] > 200: #NIPBL and SLC2A1
			INSERTION_SCAN(locus,"atgtga",2,coords[0],coords[1],"UTR")#add STARTSTOP every 2nd base
		coords = MakeCoords(locus["UTR"])
		if coords[1] < 200: #CTCF, DDX3X and SON
			INSERTION_SCAN(locus,"gccatggc",Kozak_step,coords[0],coords[1],"UTR")# add Strong Kozak every nth base
		if coords[1] > 200: #NIPBL and SLC2A1
			INSERTION_SCAN(locus,"gccatggc",2,coords[0],coords[1],"UTR")# add Strong Kozak every 2nd base
		coords = MakeCoords(locus["UTR"])
		if coords[1] < 200: #CTCF, DDX3X and SON
			INSERTION_SCAN(locus,"tttatgct",Kozak_step,coords[0],coords[1],"UTR")#add Weak Kozak every nth base
		if coords[1] > 200: #NIPBL and SLC2A1
			INSERTION_SCAN(locus,"tttatgct",2,coords[0],coords[1],"UTR")# add Weak Kozak every 2nd base 	
		coords = MakeCoords(locus["UTR"])
		INSERTION_SCAN(locus,"ggggttgggggtgggtgggg",RG4_step,coords[0],coords[1],"UTR")#RG4 forming sequence 
		INSERTION_SCAN(locus,"gggggaggttcgcctccccc",Loop_step,coords[0],coords[1],"UTR")#Stem loop forming sequence
		INSERTION_SCAN(locus,"tcaccattatcgtttcagac",SSacceptor_step,coords[0],coords[1],"UTR") #SS acceptor sequence
		INSERTION_SCAN(locus,"gcagtaagtaatacatgtaa",SSdonor_step,coords[0],coords[1],"UTR") #SS donor sequence
		INSERTION_SCAN(locus,"gtcagacgagtcatcaa",names_step,coords[0],coords[1],"UTR")#VANESSA
		INSERTION_SCAN(locus,"atgattcatgagtgg",names_step,coords[0],coords[1],"UTR")#MATTHEW
		INSERTION_SCAN(locus,"tcagagatcatataagac",names_step,coords[0],coords[1],"UTR")#SEBASTIAN
		INSERTION_SCAN(locus,"gagctcatatcaagagtcat",names_step,coords[0],coords[1],"UTR")#ELIZABETH
		# if statement correspond to functions that are applied to more than one region 
		#===========================================
		################ SNV Scan ##################
		coords = MakeCoords(locus["UTR"])
		SNV_SCAN(locus,coords[0],coords[1],"UTR")#run full SNV across entry
		if locus["SS_acceptor"] != "0": #identify if there is SS Acceptor sequence 
			coords = MakeCoords(locus["SS_acceptor"])
			SNV_SCAN(locus,coords[0],coords[1],"SSA")
		if locus["CDS"] != "0": #identify if there is CDS
			coords = MakeCoords(locus["CDS"])
			SNV_SCAN(locus,coords[0],coords[1],"CDS")
		if locus["SS_donor"] != "0": #identify if there is SS Donor sequence 
			coords = MakeCoords(locus["SS_donor"])
			SNV_SCAN(locus,coords[0],coords[1],"SSD")
		#Also applied to SS acceptor site / Coding sequence / SS donor site if there is one 
		#===========================================
		############### uORF Scan ################
		coords = MakeCoords(locus["UTR"])
		uORF_SCAN(locus,"gccatggc","tga",uORF_step,coords[0],coords[1],"UTR") #uses strong Kozak sequence 
		#===========================================
		################# SV40 #####################
		SV40_large(locus,"tacgtagatccagacatgataagatacattgatgagtttggacaaaccacaactagaatgcagtgaaaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaaccattataagctgcaataaacaagttaacaacaacaattgcattcattttatgtttcaggttcagggggaggtgtgggaggttttttaattc")
		SV40_large(locus,"gaattaaaaaacctcccacacctccccctgaacctgaaacataaaatgaatgcaattgttgttgttaacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatcttatcatgtctggatctacgta")
		SV40_small(locus,"taagatacattgatgagtttggacaaaccacaactagaatgcagtgaaaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaaccattataagctgcaataaacaagtt")
		SV40_small(locus,"aacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatctta")
		#===========================================
		################ Inversion ################
		inversion(locus)	
		#===========================================
		############## Permutation ################
		permutation(locus)	


#==============MAIN PROGRAM===========================
Oligo_list = file_extractor() #grab list of targetons to mutate

mutify(Oligo_list)#mutate the list!



