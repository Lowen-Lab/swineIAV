import os
import sys
import networkx as nx
import numpy as np
###################################   USER DEFINED VARIABLES   ###################################
project_dir = "./"

min_other_freq = 0.03
max_poly_prop = 0.2
max_N_prop = 0.1

ref_list = ['H1N1_HA.fa','H1N1_NA.fa','H1N1_NP.fa','H1N1_PB1.fa','H3N2_HA.fa','H3N2_NA.fa','H3N2_NP.fa','H3N2_PB1.fa','H1N1_M.fa','H1N1_NEP.fa','H1N1_PA.fa','H1N1_PB2.fa','H3N2_M.fa','H3N2_NEP.fa','H3N2_PA.fa','H3N2_PB2.fa']
H1N1_ref_list = ['H1N1_HA.fa','H1N1_NA.fa','H1N1_NP.fa','H1N1_PB1.fa','H1N1_M.fa','H1N1_NEP.fa','H1N1_PA.fa','H1N1_PB2.fa']
H3N2_ref_list = ['H3N2_HA.fa','H3N2_NA.fa','H3N2_NP.fa','H3N2_PB1.fa','H3N2_M.fa','H3N2_NEP.fa','H3N2_PA.fa','H3N2_PB2.fa']
write_ref_list = ['H1N1_HA.fa','H1N1_NA.fa','H3N2_HA.fa','H3N2_NA.fa','H3N2_NP.fa','H3N2_PB1.fa','H3N2_M.fa','H3N2_NEP.fa','H3N2_PA.fa','H3N2_PB2.fa']
write_ref_list_name_dict = {'H1N1_HA.fa':'H1N1_HA','H1N1_NA.fa':'H1N1_NA','H3N2_HA.fa':'H3N2_HA','H3N2_NA.fa':'H3N2_NA','H3N2_NP.fa':'H3N2_NP','H3N2_PB1.fa':'H3N2_PB1','H3N2_M.fa':'H3N2_M','H3N2_NEP.fa':'H3N2_NEP','H3N2_PA.fa':'H3N2_PA','H3N2_PB2.fa':'H3N2_PB2'}

############################################ FUNCTIONS ############################################

def make_nexus(seq_dict):
	outlines = "#NEXUS\n[ TITLE ]\nBEGIN TAXA;\n       DIMENSIONS NTAX="+str(len(seq_dict))+";\n       TAXLABELS\n"
	len_msa = 0
	for iso in seq_dict:
		outlines += "             "+iso+"\n"
		len_msa = len(seq_dict[iso])
	outlines += ";\nEND;\nBEGIN CHARACTERS;\n       DIMENSIONS NCHAR="+str(len_msa)+";\n       FORMAT MISSING=? GAP=- MATCHCHAR=. datatype=nucleotide;\nMATRIX\n\n"
	for iso in seq_dict:
		outlines += iso+"     "+seq_dict[iso].strip().replace("-","N").replace("N","?")+"\n"
	outlines += ";\nEND;\n"
	return outlines

def make_fasta_alignment(seq_dict):
	outlines = ''
	for iso in seq_dict:
		outlines += ">"+iso+"\n"+seq_dict[iso].strip()+"\n"
	return outlines
	

############################################   MAIN   #############################################

msa_infile_name = "swine_consensus_seqs.fasta"
msa_infile = open(project_dir+msa_infile_name,"r")
msa_dict = {}
segment_list = []
segment_accessions = {}
accession_list = []
for line in msa_infile:
	line = line.strip()
	if line[0] == ">":
		accession = line[1:len(line)].split("-")[0]
		segment = line[1:len(line)].split("-")[1]
		acc_seg = accession+"\t"+segment
		accession_list.append(accession)
		segment_list.append(segment)
	else:
		try:
			msa_dict[segment][accession] = line
			segment_accessions[segment].append(accession)
		except:
			msa_dict[segment] = {}
			msa_dict[segment][accession] = line
			segment_accessions[segment] = []
			segment_accessions[segment].append(accession)
msa_infile.close()
accession_list = list(set(accession_list))
segment_list = list(set(segment_list))


seq_poly_count_dict = {}
poly_mask_dict = {}
max_len_dict = {}
filename = "swine_consensus_seqs.poly_mask.txt"
infile = open(project_dir+filename,"r")
for line in infile:
	line = line.strip().split("\t")
	accession = line[0]
	segment = line[1]
	poly_list = line[2].split(",")
	try:
		curr_max_len = max_len_dict[segment]
		if len(poly_list) > curr_max_len:
			max_len_dict[segment] = len(poly_list)
	except:
		max_len_dict[segment] = len(poly_list)
	for loc in range(0,len(poly_list)):
		poly = float(poly_list[loc])
		if poly >= min_other_freq:
			try:
				poly_mask_dict[segment][accession][loc] = poly
			except:
				try:
					poly_mask_dict[segment][accession] = {}
					poly_mask_dict[segment][accession][loc] = poly
				except:
					poly_mask_dict[segment] = {}
					poly_mask_dict[segment][accession] = {}
					poly_mask_dict[segment][accession][loc] = poly
			try:
				seq_poly_count_dict[segment][accession] += 1
			except:
				try:
					seq_poly_count_dict[segment][accession] = 1
				except:
					seq_poly_count_dict[segment] = {}
					seq_poly_count_dict[segment][accession] = 1
infile.close()

mask_dict = {}
mask_info = open("segment_MSA_mask.txt","w")
for segment in segment_list:
	msa_info = open("segment_MSA_info."+segment+".txt","w")
	msa_info.write("loc\tN_count\tgap_count\tpoly_count\n")
	mask_dict[segment] = {}
	mask_info.write(segment)
	temp_accession_list = list(set(segment_accessions[segment]))
	for loc in range(0,max_len_dict[segment]):
		N_count = 0
		gap_count = 0
		poly_count = 0
		for i in range(0,len(temp_accession_list)):
			accession = temp_accession_list[i]
			nt = msa_dict[segment][accession][loc]
			if nt == "N":
				N_count +=1
			elif nt == "-":
				gap_count +=1
			else:
				try:
					poly = poly_mask_dict[segment][accession][loc]
					poly_count += 1
				except:
					pass
		N_prop = float(N_count)/float(len(temp_accession_list))
		gap_prop = float(gap_count)/float(len(temp_accession_list))
		poly_prop = float(poly_count)/float(len(temp_accession_list))
		keep = 0
		if N_prop+gap_prop > max_N_prop:
			keep += 1
		if poly_prop > max_poly_prop:
			keep += 1
		if keep > 0:
			mask_dict[segment][loc] = 0
			mask_info.write("\t0")
		else:
			mask_dict[segment][loc] = 1
			mask_info.write("\t1")
		msa_info.write(str(loc)+"\t"+str(N_count)+"\t"+str(gap_count)+"\t"+str(poly_count)+"\t"+str(N_prop+gap_prop)+"\t"+str(poly_prop)+"\t"+str(mask_dict[segment][loc])+"\n")
	mask_info.write("\n")
mask_info.close()
msa_info.close()

# segment_list = ["H3N2_PB1"]

low_poly_accessions = {}
for segment in segment_list:
	temp_accession_list = list(set(segment_accessions[segment]))
	low_poly_accessions[segment] = []
	for i in range(0,len(temp_accession_list)):
		accession = temp_accession_list[i]
		try:
			poly_count = seq_poly_count_dict[segment][accession]
		except:
			poly_count = 0
		if float(poly_count)/float(max_len_dict[segment]) <= 0.01:
			if float(len(msa_dict[segment][accession])-len(msa_dict[segment][accession].replace("-","").replace("N","")))/float(len(msa_dict[segment][accession])) <= 0.03:
				low_poly_accessions[segment].append(accession)

seq_out_dict = {}
for segment in segment_list:
	seq_out_dict[segment] = {}
	temp_accession_list = list(set(low_poly_accessions[segment]))
	for i in range(0,len(temp_accession_list)):
		accession = temp_accession_list[i]
		for loc in range(0,max_len_dict[segment]):
			if mask_dict[segment][loc] == 1:
				nt = msa_dict[segment][accession][loc]
				try:
					seq_out_dict[segment][accession] += nt
				except:
					seq_out_dict[segment][accession] = nt

	seq_outfile = open(segment+".masked.fasta","w")
	seq_out_lines = make_fasta_alignment(seq_out_dict[segment])
	seq_outfile.write(seq_out_lines)
	seq_outfile.close()

dSNP_outfile = open("dSNP_info.txt","w")
dist_dict = {}

n = nx.Graph()
node_hits = []
node_size_dict = {}

for segment in segment_list:
	print(segment)
	dist_outfile = open(segment+".dist.txt","w")
	hits = []
	dist_dict[segment] = {}
	temp_accession_list = list(set(low_poly_accessions[segment]))
	for i in range(0,len(temp_accession_list)):
		accession1 = temp_accession_list[i]
		dist_outfile.write("\t"+accession1)
	
	for i in range(0,len(temp_accession_list)):
		accession1 = temp_accession_list[i]
		dist_outfile.write(accession1)
		for j in range(0,len(temp_accession_list)):
			accession2 = temp_accession_list[j]
			dist = 0
			for loc in range(0,max_len_dict[segment]):
				if mask_dict[segment][loc] == 1:
					nt1 = msa_dict[segment][accession1][loc]
					nt2 = msa_dict[segment][accession2][loc]
					if nt1 != "N" and nt2 != "N" and nt1 != "-" and nt2 != "-":
						if nt1 != nt2:
							dist += 1
			dist_dict[segment][accession1+"\t"+accession2] = dist
			dist_dict[segment][accession2+"\t"+accession1] = dist
			dist_outfile.write("\t"+str(dist))
			if dist == 0:
				tup = (accession1,accession2)
				hits.append(tup)
		dist_outfile.write("\n")


	outfile = open(segment+".nodes.txt","w")
	g = nx.Graph()
	g.add_edges_from(hits)

	all_components = list(nx.connected_components(g))
	print(len(all_components))

	node_dict = {}
	# node_size_dict = {}
	local_node_sizes = []
	comp = 0
	for component in all_components:
		comp += 1
		component_nodes = list(component)
		# print(component_nodes)

		nodeID = "node_"+str(comp)
		outfile.write(nodeID +"\t")
		node_dict[nodeID] = []
		node_size_dict[segment+"_"+nodeID] = len(component_nodes)
		tup = (len(component_nodes),nodeID)
		local_node_sizes.append(tup)
		# compoent_list = []
		for i in range(0,len(component_nodes)):
			node_f =  component_nodes[i]
			# strain = component_nodes[i].split("_")
			# strain = strain[0]+"_"+strain[1]+"_"+strain[2]
			outfile.write(component_nodes[i] +"\t")
			node_dict[nodeID].append(component_nodes[i])
		outfile.write("\n")
	outfile.close()



	local_node_sizes = sorted(local_node_sizes, reverse=True)
	first_node = local_node_sizes[0][1]
	second_node = local_node_sizes[1][1]
	first_node_accessions = node_dict[first_node]
	second_node_accessions = node_dict[second_node]
	# main_node_accessions = first_node_accessions + second_node_accessions

	keep_column_list = []
	keep_column_dict = {}
	ref_nt_list = ['A','C','G','T']
	
	for loc in range(0,max_len_dict[segment]):
		if mask_dict[segment][loc] == 1:
			nt_count = {}
			base_list = []
			for accession in first_node_accessions:
				group = "first"
				nt = msa_dict[segment][accession][loc]
				base_list.append(nt)
				try:
					nt_count[group][nt] += 1
				except:
					try:
						nt_count[group][nt] = 1
					except:
						nt_count[group] = {'A':0,'C':0,'G':0,'T':0,'N':0,'-':0}
						nt_count[group][nt] = 1
			for accession in second_node_accessions:
				group = "second"
				nt = msa_dict[segment][accession][loc]
				base_list.append(nt)
				try:
					nt_count[group][nt] += 1
				except:
					try:
						nt_count[group][nt] = 1
					except:
						nt_count[group] = {'A':0,'C':0,'G':0,'T':0,'N':0,'-':0}
						nt_count[group][nt] = 1
			base_list = list(set(base_list))

			group_list = ['first','second']
			keep_column = False
			conflicts = 0
			if len(base_list) > 1:
				for group in group_list:
					poly_count = 0
					focus_nts = []
					for nt in base_list:
						if nt != "N" and nt != "-":
							try:
								count = nt_count[group][nt]
							except:
								count = 0
							if count >= 1:
								poly_count += 1
								focus_nts.append(nt)
					focus_nts = list(set(focus_nts))
					# if len(focus_nts) == 1:
					for focus_nt in focus_nts:
						other_found = 0
						for othergroup in group_list:
							if othergroup != group:
								try:
									count = nt_count[othergroup][focus_nt]
								except:
									count = 0
								if count >= 1:
									other_found += 1
						if other_found == 0:
							keep_column = True
							# print(focus_nts)
						else:
							conflicts += 1
			if keep_column == True and conflicts == 0:
				keep_column_list.append(loc)
				keep_column_dict[loc] = {}
				for group in group_list:
					keep_column_dict[loc][group] = []
					for nt in base_list:
						if nt!= "N" and nt!= "-":
							try:
								count = nt_count[group][nt]
							except:
								count = 0
							if count >= 1:
								keep_column_dict[loc][group].append(nt)
	keep_column_list = sorted(list(set(keep_column_list)))

	for z in range(0,len(keep_column_list)):
		loc = keep_column_list[z]
		dSNP_outfile.write(segment+"\t"+str(loc))
		for i in range(0,len(group_list)):
			group = group_list[i]
			nt_to_write = keep_column_dict[loc][group]
			for k in range(0,len(nt_to_write)):
				if k == 0:
					dSNP_outfile.write("\t"+nt_to_write[k])
				else:
					dSNP_outfile.write(","+nt_to_write[k])
		dSNP_outfile.write("\n")




	for nodeID1 in node_dict:
		ref_accession1 = node_dict[nodeID1][0]
		dist_list = []
		for nodeID2 in node_dict:
			ref_accession2 = node_dict[nodeID2][0]
			if nodeID1 != nodeID2:
				dist = dist_dict[segment][ref_accession1+"\t"+ref_accession2]
				tup = (dist,nodeID1,nodeID2)
				dist_list.append(tup)
		dist_list = sorted(dist_list)
		min_dist = dist_list[0][0]
		for i in range(0,len(dist_list)):
			ref1 = node_dict[dist_list[i][1]][0]
			ref2 = node_dict[dist_list[i][2]][0]
			dist = dist_dict[segment][ref1+"\t"+ref2]
			node_hits.append((segment+"_"+dist_list[i][1],segment+"_"+dist_list[i][2],{'weight':dist}))


n.add_edges_from(node_hits)
T=nx.minimum_spanning_tree(n)
edge_labels = nx.get_edge_attributes(T,'weight')
nx.set_node_attributes(T, node_size_dict, name='node_size')
nx.set_edge_attributes(T, edge_labels, "edge_weight")
nx.write_graphml(T,"all_segments.haplotype_network.graphml")
