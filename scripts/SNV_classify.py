import os
import sys
import numpy as np
import scipy.stats
###################################   USER DEFINED VARIABLES   ###################################
project_dir = "./"
trim_dir = project_dir+"IAV_trim/"
map_folder = 'rep_map'
map_dir = project_dir+map_folder+"/"
iSNV_global_naming = True

flagged_samples = []
output_suffix = "swineIAV"
accession_filename = "accessions.txt"
infile = open(project_dir+accession_filename,"r")
for line in infile:
	line = line.strip()
	flagged_samples.append(line)
infile.close

skip_translate = False
skip_writing_all_site_summaries = False

max_length_dict = {'H1N1_HA':1781,'H1N1_NA':1460,'H3N2_HA':1762,'H3N2_NA':1467,'H3N2_PB2':2341,'H3N2_PB1':2341,'H3N2_PA':2242,'H3N2_NP':1565,'H3N2_M':1027,'H3N2_NEP':890,'H1N1_PB2':2341,'H1N1_PB1':2341,'H1N1_PA':2242,'H1N1_NP':1565,'H1N1_M':1027,'H1N1_NEP':890}
ref_set_list = ['H3N2_PB2.fa-H1N1_PB2.fa','H3N2_PB1.fa-H1N1_PB1.fa','H3N2_PA.fa-H1N1_PA.fa','H3N2_NP.fa-H1N1_NP.fa','H3N2_M.fa-H1N1_M.fa','H3N2_NEP.fa-H1N1_NEP.fa','H3N2_HA.fa','H3N2_NA.fa','H1N1_HA.fa','H1N1_NA.fa']
ref_list = ['H1N1_HA.fa','H1N1_NA.fa','H3N2_HA.fa','H3N2_NA.fa','H3N2_NP.fa','H3N2_PB1.fa','H3N2_M.fa','H3N2_NEP.fa','H3N2_PA.fa','H3N2_PB2.fa']
internal_segment_list = ['H3N2_PB2.fa','H3N2_PB1.fa','H3N2_PA.fa','H3N2_NP.fa','H3N2_M.fa','H3N2_NEP.fa']

min_proportion = 0.03
min_qual = 37.0
min_avg_map_qual = 40.0
min_avg_read_loc = 30.0
max_avg_read_mismatch = 2.0
max_avg_read_indel = 1.0

min_major_map_qual = 43
max_major_read_mismatch = 1.5
max_major_read_indel = 0.5

min_cov_dict = {'run1':50,'run2':500}
end_skip_len = 20

amino_dict = {'A':['GCT','GCC','GCA','GCG'],
'R':['CGT','CGC','CGA','CGG','AGA','AGG'],
'N':['AAT','AAC'],
'D':['GAT','GAC'],
'C':['TGT','TGC'],
'Q':['CAA','CAG'],
'E':['GAA','GAG'],
'G':['GGT','GGC','GGA','GGG'],
'H':['CAT','CAC'],
'I':['ATT','ATC','ATA'],
'L':['TTA','TTG','CTT','CTC','CTA','CTG'],
'K':['AAA','AAG'],
'M':['ATG'],
'F':['TTT','TTC'],
'P':['CCT','CCC','CCA','CCG'],
'S':['TCT','TCC','TCA','TCG','AGT','AGC'],
'T':['ACT','ACC','ACA','ACG'],
'W':['TGG'],
'Y':['TAT','TAC'],
'V':['GTT','GTC','GTA','GTG'],
'*':['TAG','TGA','TAA']}

codon_to_aa_dict = {}
for aa in amino_dict:
	codons = amino_dict[aa]
	for codon in codons:
		codon_to_aa_dict[codon] = aa

samples_to_exclude = []
sites_to_exclude = {}

############################################ FUNCTIONS ############################################
def nt_diversity(nt_cov_list):
	N = 0.0
	D = 0.0
	if len(nt_cov_list) > 1:
		for nt in range(0,len(nt_cov_list)):
			nt_cov = int(float(nt_cov_list[nt]))
			if nt_cov > 0.0:
				D += nt_cov*(nt_cov-1.0)
				N += nt_cov
		if N > 0.0:
			D = ((N*(N-1.0)) - D)/((N*(N-1.0)))
	return D

def codon_finder(site, coding_start):
	codon_number = int(np.floor((site-coding_start)/3.0))
	codon_site = int((((float((site-coding_start))/3.0)-float((np.floor((site-coding_start)/3.0))))/0.333333+0.1))
	return codon_number,codon_site

def is_number(s):
	try:
		float(s)
	except ValueError:
		return False
	return True

def round_to_n_sig_figs(val,num_sig_figs):
	if is_number(val):
		val = float(val)
		if val == 0.0:
			return '0.0'
		num_sig_figs = int(num_sig_figs)
		if num_sig_figs == 0:
			num_sig_figs = 1
		sci_val = "{:.10e}".format(val)
		split_sci_val = sci_val.split("e")
		if len(split_sci_val) == 2:
			rounded_base_number = round(float(split_sci_val[0]),num_sig_figs-1)
			exponent = int(split_sci_val[1])
			if exponent == 0:
				val_out = str(rounded_base_number) + ((num_sig_figs)-1)*'0'
			elif exponent < 0:
				exponent*=-1
				val_out = '0.' + (exponent-1)*'0' + str(rounded_base_number).replace(".","")
				val_out = str(float(val_out))
			elif exponent > 0:
				val_out = str(rounded_base_number) +'e'+ (str(exponent))
			return val_out
		else:
			return val
	else:
		return val

def sub_site_count(seqin, coding_start, coding_stop,codon_to_aa_dict):
	nt_list = ['A','T','C','G']
	viable_codon_count = 0.0
	ns_count = 0.0
	s_count = 0.0
	for loc in range(coding_start,coding_stop+1):
		codon_number,codon_site = codon_finder(loc, coding_start)
		codon_bases = seqin[(codon_number*3+coding_start):(codon_number*3+coding_start)+3]
		if "-" not in codon_bases and "N" not in codon_bases and len(codon_bases)==3:
			viable_codon_count += 1.0
			for nt in nt_list:
				if nt != seqin[loc]:
					sub_codon = list(codon_bases)
					sub_codon[codon_site] = nt
					sub_codon = "".join(sub_codon)
					ref_aa = codon_to_aa_dict[codon_bases]
					sub_aa = codon_to_aa_dict[sub_codon]
					if ref_aa != sub_aa:
						ns_count += 0.3333
					elif ref_aa == sub_aa:
						s_count += 0.3333
	ns_prop = round(float(ns_count)/float(viable_codon_count),3)
	s_prop = round(float(s_count)/float(viable_codon_count),3)
	return ns_prop,s_prop,round(ns_count,2),round(s_count,2),int(viable_codon_count)

def translate_cds(seqin, coding_start, coding_stop,codon_to_aa_dict):
	aa_seq_out = ''
	for loc in range(coding_start,coding_stop+1,3):
		codon_number,codon_site = codon_finder(loc, coding_start)
		# print(str(loc)+"\t"+str(codon_number)+"\t"+str(codon_site))
		codon_bases = seqin[(codon_number*3+coding_start):(codon_number*3+coding_start)+3]
		if "-" not in codon_bases and "N" not in codon_bases and len(codon_bases)==3:
			aa = codon_to_aa_dict[codon_bases]
		else:
			aa = 'X'
		aa_seq_out += aa
	return aa_seq_out

###################################################################################################
### Load the group-wide consensus sequences for each segment 
consensus_seq_dict = {}
infile = open(project_dir+"swineIAV_consensus_seqs.single.fasta","r")
for line in infile:
	line = line.strip()
	if line[0] == ">":
		head = line[1:len(line)]
		group = head.split("-")[0]
		segment = head.split("-")[1]#+"_"+head.split("_")[1]
	else:
		try:
			consensus_seq_dict[group][segment] += line
		except:
			try:
				consensus_seq_dict[group][segment] = line
			except:
				consensus_seq_dict[group] = {}
				consensus_seq_dict[group][segment] = line

### Load consensus segment sequences for each sample 
genome_infile = open(project_dir+"swine_consensus_seqs.fasta","r")
seq_dict = {}
for line in genome_infile:
	line = line.strip()
	if line[0] == ">":
		accession = line[1:len(line)].split("-")[0]
		segment = line[1:len(line)].split("-")[1]
		acc_seg = accession+"-"+segment
	else:
		try:
			seq_dict[acc_seg] += line
		except:
			seq_dict[acc_seg] = line
genome_infile.close()

#load file that links the first sample sequenced from each pig
first_sample_dict = {}
filename = "pig_first_sample.txt"
infile = open(project_dir+filename,"r")
for line in infile:
	line = line.strip().split("\t")
	pigID = line[0]
	sampleID = line[2]
	first_sample_dict[pigID] = sampleID
infile.close()

#load file that links each accession name to the location of the reads
sample_run_dict = {}
filename = "accession_run.txt"
infile = open(project_dir+filename,"r")
for line in infile:
	line = line.strip().split("\t")
	sampleID = line[0]
	runID = line[1]
	# if runID == 'run2':
	if sampleID in flagged_samples:
		sample_run_dict[sampleID] = runID
	elif len(flagged_samples) == 0:
		sample_run_dict[sampleID] = runID
		output_suffix = "all"
		# flagged_samples.append(sampleID)
infile.close()

#load file that links samples to their metadata characteristics, such as pigID and time of sample
sample_dict = {}
pig_list = []
pig_to_samples = {}
pig_sample_list_dict = {}
sample_order_list = []
filename = "pig_to_sample.txt"
infile = open(project_dir+filename,"r")
for line in infile:
	line = line.strip().split("\t")
	pigID = line[0]
	sampleID = line[1]
	time = int(line[2])
	runID = ''
	try:
		runID = sample_run_dict[sampleID]
	except:
		runID = ''
	if runID != '':
		tup = (pigID,time)
		sample_dict[sampleID] = tup #(pigID,time)
		sample_dict[sampleID.split("_S")[0]] = tup
		sample_order_list.append(sampleID)
		pig_list.append(pigID)
		try:
			pig_to_samples[pigID][time] = sampleID
			pig_sample_list_dict[pigID].append(sampleID)
		except:
			pig_to_samples[pigID] = {}
			pig_to_samples[pigID][time] = sampleID
			pig_sample_list_dict[pigID] = []
			pig_sample_list_dict[pigID].append(sampleID)
infile.close()
pig_list = list(set(pig_list))

## load file that relates what genotype each segment is classified
group_dict = {}
group_list = []
pig_group_dict = {}
group_strain_lists = {}
infile = open(project_dir+"IAV_groups.txt","r")
for line in infile:
	line = line.strip().split("\t")
	segment = line[0]
	accession = line[1]
	group = line[2]
	if group == "H1N1" or group == "H3N2":
		group_list.append(group)
		try:
			group_dict[segment][accession] = group
		except:
			group_dict[segment] = {}
			group_dict[segment][accession] = group
		try:
			pigID = sample_dict[accession][0]
			try:
				prev_pig_group = pig_group_dict[pigID]
				if prev_pig_group != group:
					print(pigID+" attributed multiple groupIDs")
			except:				
				pig_group_dict[pigID] = group
		except:
			pass
infile.close()

#load file that indicates the location and nucleotide for SNPs that differentiate the groups in internal segments
dSNP_dict = {}
dSNP_group_nt_dict = {}
dSNP_sites = {}
segment_list = []
infile = open(project_dir+"swineIAV_dSNPs.txt","r")
first_line = True
for line in infile:
	line = line.strip().split("\t")
	if first_line == True:
		group_list = line
		first_line = False
	else:
		segment = line[0]
		segment_list.append(segment)
		loc = int(line[1])
		for i in range(2,len(line)):
			group = group_list[i-2]
			nt = line[i]
			try:
				dSNP_sites[segment].append(loc)
			except:
				dSNP_sites[segment] = []
				dSNP_sites[segment].append(loc)
			try:
				dSNP_dict[segment][loc][nt] = group
			except:
				try:
					dSNP_dict[segment][loc] = {}
					dSNP_dict[segment][loc][nt] = group
				except:
					dSNP_dict[segment] = {}
					dSNP_dict[segment][loc] = {}
					dSNP_dict[segment][loc][nt] = group
			try:
				dSNP_group_nt_dict[segment][loc][group] = nt
			except:
				try:
					dSNP_group_nt_dict[segment][loc] = {}
					dSNP_group_nt_dict[segment][loc][group] = nt
				except:
					dSNP_group_nt_dict[segment] = {}
					dSNP_group_nt_dict[segment][loc] = {}
					dSNP_group_nt_dict[segment][loc][group] = nt
infile.close()
segment_list = list(set(segment_list))
group_list = sorted(list(set(group_list))) #['H1N1','H3N2']


# load file that stores the location of CDS start and stop within the reference sequences that reads were aligned to
cds_loc_dict = {}
infile = open(project_dir+"ref_cds_loc.txt","r")
for line in infile:
	line = line.strip().split("\t")
	segment = line[1]
	locs = line[2]
	cds_loc_dict[segment] = locs
infile.close()

#Process consensus sequences to translate them, also to count the expectation for NS/S ratio
if skip_translate == False:
	outfile = open(project_dir+"S_NS_subcount."+output_suffix+".txt","w")
	seq_outfile = open(project_dir+"translated_seqs."+output_suffix+".faa","w")
	for acc_seg in seq_dict:
		accession = acc_seg.split("-")[0]
		segment = acc_seg.split("-")[1]
		seq = seq_dict[acc_seg]
		cds_start = int(cds_loc_dict[segment].split("-")[0])-1
		cds_stop = int(cds_loc_dict[segment].split("-")[1])-1
		ns_prop,s_prop,ns_count,s_count,viable_codon_count = sub_site_count(seq, cds_start, cds_stop,codon_to_aa_dict)
		outfile.write(segment+"\t"+accession+"\t"+str(ns_prop)+"\t"+str(s_prop)+"\t"+str(ns_count)+"\t"+str(s_count)+"\t"+str(viable_codon_count)+"\n")
		translated_seq = translate_cds(seq, cds_start, cds_stop,codon_to_aa_dict)
		seq_outfile.write(">"+accession+"_"+segment+"\n"+translated_seq+"\n")

	for seg_refs in ref_set_list:
		refs = seg_refs.split("-")
		for ref in refs:
			seq = ''
			segment = ref.split(".")[0]
			ref_segment = seg_refs.split("-")[0].split(".")[0]
			temp_infile = open(map_dir+ref,"r")
			for line in temp_infile:
				line = line.strip()
				if line[0] != ">":
					seq += line
			cds_start = int(cds_loc_dict[ref_segment].split("-")[0])-1
			cds_stop = int(cds_loc_dict[ref_segment].split("-")[1])-1
			ns_prop,s_prop,ns_count,s_count,viable_codon_count = sub_site_count(seq, cds_start, cds_stop,codon_to_aa_dict)
			outfile.write(segment+"\tref"+"\t"+str(ns_prop)+"\t"+str(s_prop)+"\t"+str(ns_count)+"\t"+str(s_count)+"\t"+str(viable_codon_count)+"\n")
			translated_seq = translate_cds(seq, cds_start, cds_stop,codon_to_aa_dict)
			seq_outfile.write(">ref_"+segment+"\n"+translated_seq+"\n")
	outfile.close()
	seq_outfile.close()
del skip_translate
	

site_count_all = {}
site_count_all_passcov = {}
group_count_dict = {}
dSNP_freq_dict = {}
poly_dict = {}
all_allele_info_dict = {}
poly_sum_dict = {}
S_NS_dict = {}
site_cov_dict = {}
sub_list = []
pig_iSNV_dict = {}
accession_list = []

info_tup_dict = {'cov':0,'prop':1,'qual':2,'map_qual':3,'read_len':4,'read_loc':5,'read_R_loc':6,'read_prop_loc':7,'mismatch':8,'indel':9}
base_tup_dict = {'A':0,'T':1,'C':2,'G':3}

iSNV_info = open(project_dir+"iSNV_info.self."+output_suffix+".txt","w")
all_sites_info = open(project_dir+"all_sites_info.self."+output_suffix+".txt","w")
dSNP_info = open(project_dir+"dSNP_info."+output_suffix+".txt","w")
iSNV_info.write("accession\tsegment\tloc\tcov\tmajor_nt\tminor_nt\tmajor_prop\tminor_prop\tmajor_qual\tminor_qual\tmajor_map_qual\tminor_map_qual\tmajor_read_len\tminor_read_len\tmajor_base_loc\tminor_base_loc\tmajor_base_R_loc\tminor_base_R_loc\tmajor_prop_loc\tminor_prop_loc\tmajor_mismatches\tminor_mismatches\tmajor_indels\tminor_indels\n")
site_count_outfile = open(project_dir+"site_count_info."+output_suffix+".txt","w")
for accession in sample_run_dict:
	skip = False
	try:
		pigID = sample_dict[accession][0]
		runID = sample_run_dict[accession]
	except:
		skip = True
	if skip == False:
		cov_thresh = min_cov_dict[runID]
		site_count_all[accession] = {}
		site_count_all_passcov[accession] = {}
		group_count_dict[accession] = {}
		dSNP_freq_dict[accession] = {}
		poly_dict[accession] = {}
		all_allele_info_dict[accession] = {}
		site_cov_dict[accession] = {}
		accession_list.append(accession)
		for seg_ref in ref_set_list:
			ref_filename = seg_ref.split("-")[0]
			segment = ref_filename.split(".f")[0]
			site_cov_dict[accession][segment] = {}
			all_allele_info_dict[accession][segment] = {}
			acc_seg = accession+"-"+segment
			substitutions_filename = map_dir+"substitutions/"+ref_filename.split(".f")[0]+"/"+accession+'.'+ref_filename.split(".f")[0]+".substitions.txt"
			if os.path.isfile(substitutions_filename) == True:
				sub_infile = open(substitutions_filename,"r")
				for sec_line in sub_infile:
					if len(sec_line) > 0:
						if sec_line[0] != "#" and sec_line[0] != "loc":
							sec_line = sec_line.strip().split("\t")
							loc = int(sec_line[0])
							cov = float(sec_line[1])
							A_prop,T_prop,C_prop,G_prop = float(sec_line[2]),float(sec_line[3]),float(sec_line[4]),float(sec_line[5])
							A_qual,T_qual,C_qual,G_qual = float(sec_line[6]),float(sec_line[7]),float(sec_line[8]),float(sec_line[9])
							A_map_qual,T_map_qual,C_map_qual,G_map_qual = float(sec_line[10]),float(sec_line[11]),float(sec_line[12]),float(sec_line[13])
							A_read_len,T_read_len,C_read_len,G_read_len = float(sec_line[14]),float(sec_line[15]),float(sec_line[16]),float(sec_line[17])
							A_read_loc,T_read_loc,C_read_loc,G_read_loc = float(sec_line[18]),float(sec_line[19]),float(sec_line[20]),float(sec_line[21])
							A_read_R_loc,T_read_R_loc,C_read_R_loc,G_read_R_loc = float(sec_line[22]),float(sec_line[23]),float(sec_line[24]),float(sec_line[25])
							A_read_prop_loc,T_read_prop_loc,C_read_prop_loc,G_read_prop_loc = float(sec_line[26]),float(sec_line[27]),float(sec_line[28]),float(sec_line[29])
							A_read_mismatch,T_read_mismatch,C_read_mismatch,G_read_mismatch = float(sec_line[30]),float(sec_line[31]),float(sec_line[32]),float(sec_line[33])
							A_read_indel,T_read_indel,C_read_indel,G_read_indel = float(sec_line[34]),float(sec_line[35]),float(sec_line[36]),float(sec_line[37])
							info_tup = (cov,(A_prop,T_prop,C_prop,G_prop),
								(A_qual,T_qual,C_qual,G_qual),
								(A_map_qual,T_map_qual,C_map_qual,G_map_qual),
								(A_read_len,T_read_len,C_read_len,G_read_len),
								(A_read_loc,T_read_loc,C_read_loc,G_read_loc),
								(A_read_R_loc,T_read_R_loc,C_read_R_loc,G_read_R_loc),
								(A_read_prop_loc,T_read_prop_loc,C_read_prop_loc,G_read_prop_loc),
								(A_read_mismatch,T_read_mismatch,C_read_mismatch,G_read_mismatch),
								(A_read_indel,T_read_indel,C_read_indel,G_read_indel))
							all_allele_info_dict[accession][segment][loc] = info_tup
							site_cov_dict[accession][segment][loc] = cov


###############################################################################################################################
##### Perform site specific relative mapping quality adjustment
avg_strain_map_qual_dict = {}
minor_prop_stat_dict = {}
major_read_loc_dict = {}
nt_list = ['A','T','C','G']

for seg_ref in ref_set_list:
	seg = seg_ref.split("-")[0]
	segment = seg.split(".")[0]
	avg_strain_map_qual_dict[segment] = {}
	minor_prop_stat_dict[segment] = {}
	for loc in range(20,max_length_dict[segment]-20):
		for a in range(0,len(accession_list)):
			accession = accession_list[a]
			runID = sample_run_dict[accession]
			cov_thresh = min_cov_dict[runID]
			nt_prop_list = []
			try:
				site_cov = all_allele_info_dict[accession][segment][loc][info_tup_dict['cov']]
			except:
				site_cov = 0
			if site_cov >= cov_thresh:
				for nt_num in range(0,4):
					try:
						nt_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][nt_num]
						# if nt_prop >= min_proportion:
						nt_prop_list.append((nt_prop,nt_list[nt_num]))
					except:
						pass
				if nt_prop_list != []:
					nt_prop_list = sorted(nt_prop_list,reverse=True)
					major_nt,major_prop = nt_prop_list[0][1],nt_prop_list[0][0]
					minor_nt,minor_prop = nt_prop_list[1][1],nt_prop_list[1][0]
					log_minor_prop = np.log10(minor_prop+0.0001)*-1
					major_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[major_nt]]
					minor_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[minor_nt]]
					major_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[major_nt]]
					minor_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[minor_nt]]
					try:
						major_read_loc_dict[segment][loc].append(major_read_loc)
					except:
						try:
							major_read_loc_dict[segment][loc] = []
							major_read_loc_dict[segment][loc].append(major_read_loc)
						except:
							major_read_loc_dict[segment] = {}
							major_read_loc_dict[segment][loc] = []
							major_read_loc_dict[segment][loc].append(major_read_loc)
					try:
						avg_strain_map_qual_dict[segment][accession].append(major_map_qual)
					except:
						avg_strain_map_qual_dict[segment][accession] = []
						avg_strain_map_qual_dict[segment][accession].append(major_map_qual)

					try:
						minor_prop_stat_dict[segment][accession].append(log_minor_prop)
					except:
						minor_prop_stat_dict[segment][accession] = []
						minor_prop_stat_dict[segment][accession].append(log_minor_prop)
for segment in avg_strain_map_qual_dict:
	for accession in avg_strain_map_qual_dict[segment]:
		if len(avg_strain_map_qual_dict[segment][accession]) >= int(max_length_dict[segment]*0.1):
			avg_map_qual = round(np.average(avg_strain_map_qual_dict[segment][accession]),1)
			avg_strain_map_qual_dict[segment][accession] = avg_map_qual
		else:
			avg_strain_map_qual_dict[segment][accession] = avg_map_qual

		if len(minor_prop_stat_dict[segment][accession]) >= int(max_length_dict[segment]*0.1):
			minor_prop_mean = round(np.mean(minor_prop_stat_dict[segment][accession]),4)
			minor_prop_std = round(np.std(minor_prop_stat_dict[segment][accession]),4)
			minor_prop_median = round(np.median(minor_prop_stat_dict[segment][accession]),4)
			minor_prop_stat_dict[segment][accession] = (minor_prop_mean,minor_prop_std,minor_prop_median)
		else:
			minor_prop_stat_dict[segment][accession] = ('nan','nan','nan')

##############################################################################################################################
#### Find iSNVs

nt_diversity_outfile = open(project_dir+"sample_nt_diversity.self."+output_suffix+".txt","w")
nt_diversity_outfile.write("accession\tsegment\tgroupID\tavg_cov\tseg_site_count\tseg_S_count\tseg_NS_count\tD_sub\tD_sub_S\tD_sub_NS\tpi_sub\tpi_S\tpi_NS\n")
sample_nt_diversity_dict = {}
avg_cov_dict = {}
for accession in all_allele_info_dict:
	avg_cov_dict[accession] = {}
	pigID = sample_dict[accession][0]
	runID = sample_run_dict[accession]
	cov_thresh = min_cov_dict[runID]
	if pigID not in pig_iSNV_dict:
		pig_iSNV_dict[pigID] = {}
	for segment in all_allele_info_dict[accession]:
		avg_cov_dict[accession][segment] = []
		acc_seg = accession+"-"+segment
		cds_start = int(cds_loc_dict[segment].split("-")[0])-1
		cds_stop = int(cds_loc_dict[segment].split("-")[1])-1
		segregating_site_count = 0
		sub_seg_site_count = 0
		seg_S_count = 0
		seg_NS_count = 0
		poly_S_count = 0
		poly_NS_count = 0
		D_all = 0.0
		D_sub = 0.0
		D_sub_S = 0.0
		D_sub_NS = 0.0
		if segment == "H1N1_HA" or segment == "H1N1_NA":
			groupID = "H1N1"
		elif segment == "H3N2_HA" or segment == "H3N2_NA":
			groupID = "H3N2"
		else:
			dSNP_freq_dict[accession][segment] = {}
			try:
				groupID = group_dict[segment.replace("NEP","NS")][accession]
			except:
				groupID = 'unknown'
		try:
			acc_seg_seq = seq_dict[acc_seg]
		except:
			acc_seg_seq = "-"*100000
		try:
			ref_seq = consensus_seq_dict[groupID][segment]
		except:
			ref_seq = "-"*100000
		site_count_all[accession][segment] = 0
		site_count_all_passcov[accession][segment] = 0

		try:
			minor_prop_stats = minor_prop_stat_dict[segment][accession]
			minor_prop_mean,minor_prop_std,minor_prop_median = minor_prop_stats[0],minor_prop_stats[1],minor_prop_stats[2]
		except:
			minor_prop_mean,minor_prop_std,minor_prop_median = 'nan','nan','nan'

		for loc in range(0,max_length_dict[segment]):
			skip_again = False
			if loc <= end_skip_len or loc >= (max_length_dict[segment]-end_skip_len):
				skip_again = True
				if cov > 0:
					site_count_all[accession][segment] += 1
			else:
				try:
					sites_to_exclude[segment+"-"+str(loc)]
					skip_again = True
				except:
					skip_again = False
			if skip_again == False:
				try:
					info_tup = all_allele_info_dict[accession][segment][loc]
				except:
					info_tup = ()
					avg_cov_dict[accession][segment].append(0)
				if info_tup != ():
					cov = info_tup[0]
					if cov >= cov_thresh:
						A_prop,T_prop,C_prop,G_prop = info_tup[1][0],info_tup[1][1],info_tup[1][2],info_tup[1][3]
						A_qual,T_qual,C_qual,G_qual = info_tup[2][0],info_tup[2][1],info_tup[2][2],info_tup[2][3]
						A_map_qual,T_map_qual,C_map_qual,G_map_qual = info_tup[3][0],info_tup[3][1],info_tup[3][2],info_tup[3][3]
						A_read_len,T_read_len,C_read_len,G_read_len = info_tup[4][0],info_tup[4][1],info_tup[4][2],info_tup[4][3]
						A_read_loc,T_read_loc,C_read_loc,G_read_loc = info_tup[5][0],info_tup[5][1],info_tup[5][2],info_tup[5][3]
						A_read_R_loc,T_read_R_loc,C_read_R_loc,G_read_R_loc = info_tup[6][0],info_tup[6][1],info_tup[6][2],info_tup[6][3]
						A_read_prop_loc,T_read_prop_loc,C_read_prop_loc,G_read_prop_loc = info_tup[7][0],info_tup[7][1],info_tup[7][2],info_tup[7][3]
						A_read_mismatch,T_read_mismatch,C_read_mismatch,G_read_mismatch = info_tup[8][0],info_tup[8][1],info_tup[8][2],info_tup[8][3]
						A_read_indel,T_read_indel,C_read_indel,G_read_indel = info_tup[9][0],info_tup[9][1],info_tup[9][2],info_tup[9][3]

						avg_cov_dict[accession][segment].append(cov)

						prop_nt_list = sorted([(A_prop,"A"),(T_prop,"T"),(C_prop,"C"),(G_prop,"G")],reverse=True)
						nt_count_list_full = [A_prop*cov,T_prop*cov,C_prop*cov,G_prop*cov]
						
						minor_nt = prop_nt_list[1][1]
						major_nt = prop_nt_list[0][1]

						minor_subID = segment+"-"+str(loc)+minor_nt
						major_subID = segment+"-"+str(loc)+major_nt

						minor_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[minor_nt]]
						minor_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[minor_nt]]
						minor_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[minor_nt]]
						minor_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[minor_nt]]
						minor_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[minor_nt]]
						minor_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[minor_nt]]
						minor_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[minor_nt]]
						minor_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][base_tup_dict[minor_nt]]
						
						major_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[major_nt]]
						major_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[major_nt]]
						major_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[major_nt]]
						major_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[major_nt]]
						major_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[major_nt]]
						major_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[major_nt]]
						major_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[major_nt]]
						major_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][base_tup_dict[major_nt]]
						
						major_qual_pass = False
						if major_qual >= min_qual and major_map_qual >= min_major_map_qual and major_read_loc >= min_avg_read_loc and major_read_mismatch <= max_major_read_mismatch and major_read_indel <= max_major_read_indel:
							major_qual_pass = True
						minor_qual_pass = False
						if minor_qual >= min_qual and minor_map_qual >= min_avg_map_qual and minor_read_loc >= min_avg_read_loc and minor_read_mismatch <= max_avg_read_mismatch and minor_read_indel <= max_avg_read_indel:
							minor_qual_pass = True

						if major_qual_pass == True:
							site_count_all_passcov[accession][segment] += 1
							if dSNP_loc_found == False:
								segregating_site_count += 1
								if minor_qual_pass == True:
									sub_seg_site_count += 1

						#Write out all dSNP allele info
						dSNP_loc_found = False
						try:
							dSNP_dict[segment][loc]
							dSNP_loc_found = True
						except:
							dSNP_loc_found = False
						if dSNP_loc_found == True:
							group1_nt = dSNP_group_nt_dict[segment][loc]['group1']
							group2_nt = dSNP_group_nt_dict[segment][loc]['group2']
							dSNP_info.write(accession+"\t"+segment+"\t"+str(loc)+"\t"+str(int(cov))+"\t"+group1_nt+"\t"+group2_nt)
							for val in range(1,10):
								group1_val = info_tup[val][base_tup_dict[group1_nt]]
								group2_val = info_tup[val][base_tup_dict[group2_nt]]
								dSNP_info.write("\t"+str(group1_val)+"\t"+str(group2_val))
							dSNP_info.write("\n")
						if minor_prop >= 0.001:
							all_sites_info.write(accession+"\t"+segment+"\t"+str(loc)+"\t"+str(int(cov))+"\t"+major_nt+"\t"+minor_nt+"\t"+str(dSNP_loc_found))
							for val in range(1,10):
								major_val = info_tup[val][base_tup_dict[major_nt]]
								all_sites_info.write("\t"+str(major_val))
							for val in range(1,10):
								minor_val = info_tup[val][base_tup_dict[minor_nt]]
								all_sites_info.write("\t"+str(minor_val))
							try:
								med_read_loc = np.median(major_read_loc_dict[segment][loc])
							except:
								med_read_loc = 'na'
							all_sites_info.write("\t"+str(med_read_loc)+"\t"+str(minor_qual_pass)+"\t"+str(major_qual_pass))
							all_sites_info.write("\n")

						
						if minor_prop >= min_proportion:
							iSNV_info.write(accession+"\t"+segment+"\t"+str(loc)+"\t"+str(int(cov))+"\t"+major_nt+"\t"+minor_nt+"\t"+str(dSNP_loc_found))
							for val in range(1,10):
								major_val = info_tup[val][base_tup_dict[major_nt]]
								iSNV_info.write("\t"+str(major_val))
							for val in range(1,10):
								minor_val = info_tup[val][base_tup_dict[minor_nt]]
								iSNV_info.write("\t"+str(minor_val))
							try:
								med_read_loc = np.median(major_read_loc_dict[segment][loc])
							except:
								med_read_loc = 'na'
							iSNV_info.write("\t"+str(med_read_loc)+"\t"+str(minor_qual_pass)+"\t"+str(major_qual_pass))
							iSNV_info.write("\n")

						# info_tup_dict = {'cov':0,'prop':1,'qual':2,'map_qual':3,'read_len':4,'read_loc':5,'read_R_loc':6,'read_prop_loc':7,'mismatch':8,'indel':9}
						nt_list = ['A','T','C','G']
						NS_counter = 0
						S_counter = 0
						pos_nt_list = []
						nt_count_list_sub = []
						S_nt_count_list_sub = []
						NS_nt_count_list_sub = []
						for nt_num in range(0,4):
							local_nt = nt_list[nt_num]
							local_sub_ID = segment+"-"+str(loc)+local_nt
							local_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][nt_num]
							local_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][nt_num]
							local_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][nt_num]
							# local_read_len = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_len']][nt_num]
							local_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][nt_num]
							# local_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][nt_num]
							# local_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][nt_num]
							local_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][nt_num]
							local_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][nt_num]

							local_qual_pass = False
							if local_qual >= min_qual and local_map_qual >= min_avg_map_qual and local_read_loc >= min_avg_read_loc and local_read_mismatch <= max_avg_read_mismatch and local_read_indel <= max_avg_read_indel:
								local_qual_pass = True

							#for counting total S and NS sites
							if loc >= cds_start and loc <= cds_stop and major_qual_pass == True and dSNP_loc_found == False:
								codon_number,codon_site = codon_finder(loc, cds_start)
								ref_acc_seg_codon = ref_seq[(codon_number*3+cds_start):(codon_number*3+cds_start)+3]
								ref_sub_codon = list(ref_acc_seg_codon)
								ref_sub_codon[codon_site] = local_nt
								ref_sub_codon = "".join(ref_sub_codon)
								if "-" not in ref_sub_codon and "N" not in ref_sub_codon and "-" not in ref_acc_seg_codon and "N" not in ref_acc_seg_codon:
									if codon_to_aa_dict[ref_acc_seg_codon] != codon_to_aa_dict[ref_sub_codon]:
										if local_nt != ref_acc_seg_codon[codon_site]:
											NS_counter += 0.3333334
										if local_qual_pass == True:#local_nt != major_nt and  and ref_sub_codon != ref_acc_seg_codon:
											NS_nt_count_list_sub.append(local_prop*cov)
									else:
										if local_nt != ref_acc_seg_codon[codon_site]:
											S_counter += 0.3333334
										if local_qual_pass == True:#local_nt != major_nt and  and ref_sub_codon != ref_acc_seg_codon:
											S_nt_count_list_sub.append(local_prop*cov)
							#Count dSNPs
							if "HA" not in segment and "NA" not in segment:
								if local_prop >= min_proportion and major_qual_pass == True and local_qual_pass == True and dSNP_loc_found == True:# and local_prop >= min_proportion and local_qual >= min_qual and local_map_qual >= min_avg_map_qual:
									try:
										group = dSNP_dict[segment][loc][local_nt]
									except:
										group = "other"
									try:
										group_count_dict[accession][segment][group] += 1
									except:
										try:
											group_count_dict[accession][segment][group] = 1
										except:
											group_count_dict[accession][segment] = {}
											group_count_dict[accession][segment][group] = 1
									try:
										dSNP_freq_dict[accession][segment][loc][group] = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[local_nt]]
									except:
										dSNP_freq_dict[accession][segment][loc] = {}
										dSNP_freq_dict[accession][segment][loc][group] = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[local_nt]]
							if local_qual_pass == True and major_qual_pass == True and dSNP_loc_found == False and local_prop >= min_proportion:
								sub_list.append(local_sub_ID)
								nt_count_list_sub.append(local_prop*cov)
								pos_nt_list.append(local_nt)
								poly_dict[accession][local_sub_ID] = local_prop
								if iSNV_global_naming == False:
									try:
										poly_sum_dict[pigID][local_sub_ID] += local_prop #for naming sub_labels from subID
									except:
										try:
											poly_sum_dict[pigID][local_sub_ID] = local_prop
										except:
											poly_sum_dict[pigID] = {}
											poly_sum_dict[pigID][local_sub_ID] = local_prop
								elif iSNV_global_naming == True:
									try:
										poly_sum_dict[local_sub_ID] += local_prop #for naming sub_labels from subID
									except:
										poly_sum_dict[local_sub_ID] = local_prop

								try:
									pig_iSNV_dict[pigID][local_sub_ID] = ''
								except:
									pig_iSNV_dict[pigID] = {}
									pig_iSNV_dict[pigID][local_sub_ID] = ''

								#For sites that are within the CDS, collect info on whether a substitution is NS or S
							if loc >= cds_start and loc <= cds_stop:# and minor_prop >= min_proportion:# and local_nt != major_nt:# and minor_prop >= min_proportion and minor_qual_pass == True:
								codon_number,codon_site = codon_finder(loc, cds_start)
								
								ref_acc_seg_codon = ref_seq[(codon_number*3+cds_start):(codon_number*3+cds_start)+3]
								ref_sub_codon = list(ref_acc_seg_codon)
								ref_sub_codon[codon_site] = local_nt
								ref_sub_codon = "".join(ref_sub_codon)
								
								# print(accession+"\t"+segment+"\t"+str(loc)+"\t"+sub_codon+"\t"+ref_sub_codon)
								if "-" not in ref_sub_codon and "N" not in ref_sub_codon and "-" not in ref_acc_seg_codon and "N" not in ref_acc_seg_codon:# and ref_sub_codon != ref_acc_seg_codon:
									local_ns_type = "nan"
									if codon_to_aa_dict[ref_acc_seg_codon] != codon_to_aa_dict[ref_sub_codon]:
										S_NS_dict[local_sub_ID] = "NS"
									else:
										S_NS_dict[local_sub_ID] = "S"
						pos_nt_list = list(set(pos_nt_list))
						seg_S_count += S_counter
						seg_NS_count += NS_counter
						if len(pos_nt_list)>1 and minor_qual_pass == True and major_qual_pass == True and minor_prop >= min_proportion:
							D_sub += nt_diversity(nt_count_list_sub)
							D_sub_S += nt_diversity(S_nt_count_list_sub)
							D_sub_NS += nt_diversity(NS_nt_count_list_sub)

		site_count_outfile.write(accession+"\t"+segment+"\t"+str(site_count_all[accession][segment])+"\t"+str(site_count_all_passcov[accession][segment])+"\n")
		if D_sub > 0 and sub_seg_site_count > 0:
			pi_subs = D_sub/float(sub_seg_site_count)
		else:
			pi_subs = 0
		if D_sub_S > 0 and seg_S_count > 0:
			pi_S_subs = D_sub_S/float(round(seg_S_count,2))
		else:
			pi_S_subs = 0
		if D_sub_NS > 0 and seg_NS_count > 0:
			pi_NS_subs = D_sub_NS/float(round(seg_NS_count,2))
		else:
			pi_NS_subs = 0
		try:
			avg_cov = round(np.average(avg_cov_dict[accession][segment]),2)
		except:
			avg_cov = 0.0
		nt_diversity_outfile.write(accession+"\t"+segment+"\t"+groupID+"\t"+str(avg_cov)+"\t"+
			str(sub_seg_site_count)+"\t"+str(round(seg_S_count,2))+"\t"+str(round(seg_NS_count,2))+"\t"+
			str(round_to_n_sig_figs(D_sub,3))+"\t"+str(round_to_n_sig_figs(D_sub_S,3))+"\t"+str(round_to_n_sig_figs(D_sub_NS,3))+"\t"+
			str(round_to_n_sig_figs(pi_subs,3))+"\t"+str(round_to_n_sig_figs(pi_S_subs,3))+"\t"+str(round_to_n_sig_figs(pi_NS_subs,3))+"\n")
		try:
			div_tup_old = sample_nt_diversity_dict[accession]
			div_tup_new = (D_sub+div_tup_old[0],D_sub_S+div_tup_old[1],D_sub_NS+div_tup_old[2],sub_seg_site_count+div_tup_old[3],seg_S_count+div_tup_old[4],seg_NS_count+div_tup_old[5])
			sample_nt_diversity_dict[accession] = div_tup_new
		except:
			div_tup_new = (D_sub,D_sub_S,D_sub_NS,sub_seg_site_count,seg_S_count,seg_NS_count)
			sample_nt_diversity_dict[accession] = div_tup_new
site_count_outfile.close()
sub_infile.close()
nt_diversity_outfile.close()
sub_list = sorted(list(set(sub_list)))
iSNV_info.close()
all_sites_info.close()
dSNP_info.close()
print("len(S_NS_dict): "+str(len(S_NS_dict)))

nt_diversity_outfile = open(project_dir+"sample_nt_diversity.all_segments."+output_suffix+".txt","w")
nt_diversity_outfile.write("pigID\ttime\taccession\tgroupID\tseg_site_count\tseg_S_count\tseg_NS_count\tD_sub\tD_sub_S\tD_sub_NS\tpi_sub\tpi_S\tpi_NS\n")
for accession in sample_nt_diversity_dict:
	div_tup = sample_nt_diversity_dict[accession]
	D_sub = div_tup[0]
	D_sub_S = div_tup[1]
	D_sub_NS = div_tup[2]
	sub_seg_site_count = div_tup[3]
	seg_S_count = div_tup[4]
	seg_NS_count = div_tup[5]
	if D_sub > 0 and sub_seg_site_count > 0:
		pi_subs = D_sub/float(sub_seg_site_count)
	else:
		pi_subs = 0
	if D_sub_S > 0 and seg_S_count > 0:
		pi_S_subs = D_sub_S/float(round(seg_S_count,2))
	else:
		pi_S_subs = 0
	if D_sub_NS > 0 and seg_NS_count > 0:
		pi_NS_subs = D_sub_NS/float(round(seg_NS_count,2))
	else:
		pi_NS_subs = 0

	info_tup = sample_dict[accession]
	pigID = info_tup[0]
	time = info_tup[1]
	try:
		groupID = pig_group_dict[pigID]
	except:
		groupID = "nan"
	nt_diversity_outfile.write(pigID+"\t"+str(time)+"\t"+accession+"\t"+groupID+"\t"+
		str(sub_seg_site_count)+"\t"+str(round(seg_S_count,2))+"\t"+str(round(seg_NS_count,2))+"\t"+
		str(round_to_n_sig_figs(D_sub,3))+"\t"+str(round_to_n_sig_figs(D_sub_S,3))+"\t"+str(round_to_n_sig_figs(D_sub_NS,3))+"\t"+
		str(round_to_n_sig_figs(pi_subs,3))+"\t"+str(round_to_n_sig_figs(pi_S_subs,3))+"\t"+str(round_to_n_sig_figs(pi_NS_subs,3))+"\n")
nt_diversity_outfile.close()


# info_tup_dict = {'cov':0,'prop':1,'qual':2,'map_qual':3,'read_len':4,'read_loc':5,'read_R_loc':6,'read_prop_loc':7,'mismatch':8,'indel':9}
# base_tup_dict = {'A':0,'T':1,'C':2,'G':3}
dSNP_freq_outfile = open(project_dir+"all_dSNP_freq."+output_suffix+".txt","w")
dSNP_freq_outfile.write("accession\tsegment\tloc\tcov\t"+group_list[0]+"\t"+group_list[1]+"\n")
for accession in accession_list:
	for seg in internal_segment_list:
		segment = seg.split(".")[0]
		dSNP_locs = list(set(dSNP_sites[segment])) 
		for l in range(0,len(dSNP_locs)):
			loc = dSNP_locs[l]
			g1_nt = dSNP_group_nt_dict[segment][loc][group_list[0]]
			g2_nt = dSNP_group_nt_dict[segment][loc][group_list[1]]
			try:
				site_cov = all_allele_info_dict[accession][segment][loc][info_tup_dict['cov']]
				g1_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[g1_nt]]
				g2_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[g2_nt]]
				g1_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[g1_nt]]
				g2_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[g2_nt]]
				g1_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[g1_nt]]
				g2_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[g2_nt]]
				g1_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[g1_nt]]
				g2_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[g2_nt]]
				g1_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[g1_nt]]
				g2_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[g2_nt]]
			except:
				site_cov = 0

			runID = sample_run_dict[accession]
			if site_cov >= min_cov_dict[runID]:
				dSNP_freq_outfile.write(accession+"\t"+segment+"\t"+str(loc)+"\t"+str(int(site_cov))+"\t"+str(g1_prop)+"\t"+str(g2_prop)+"\t"+str(g1_qual)+"\t"+str(g2_qual)+"\t"+str(g1_map_qual)+"\t"+str(g2_map_qual)+"\t"+str(g1_read_loc)+"\t"+str(g2_read_loc)+"\t"+str(g1_read_R_loc)+"\t"+str(g2_read_R_loc)+"\n")
dSNP_freq_outfile.close()

if skip_writing_all_site_summaries == False:
	nt_list = ['A','T','C','G']
	num_allele = 1
	for num_allele in range(0,2):
		if num_allele == 1:
			cov_outfile = open(project_dir+"all_cov."+output_suffix+".txt","w")
		prop_outfile = open(project_dir+"all_sites.prop."+output_suffix+"."+str(num_allele)+".txt","w")
		readloc_outfile = open(project_dir+"all_sites.read_loc."+output_suffix+"."+str(num_allele)+".txt","w")
		qual_outfile = open(project_dir+"all_sites.qual."+output_suffix+"."+str(num_allele)+".txt","w")
		mapqual_outfile = open(project_dir+"all_sites.map_qual."+output_suffix+"."+str(num_allele)+".txt","w")
		mismatch_outfile = open(project_dir+"all_sites.mismatch."+output_suffix+"."+str(num_allele)+".txt","w")
		indel_outfile = open(project_dir+"all_sites.indel."+output_suffix+"."+str(num_allele)+".txt","w")
		if num_allele == 1:
			cov_outfile.write("segment\tloc")
		prop_outfile.write("segment\tloc")
		readloc_outfile.write("segment\tloc")
		qual_outfile.write("segment\tloc")
		mapqual_outfile.write("segment\tloc")
		mismatch_outfile.write("segment\tloc")
		indel_outfile.write("segment\tloc")

		for a in range(0,len(accession_list)):
			accession = accession_list[a]
			if num_allele == 1:
				cov_outfile.write("\t"+accession)
			prop_outfile.write("\t"+accession)
			readloc_outfile.write("\t"+accession)
			qual_outfile.write("\t"+accession)
			mapqual_outfile.write("\t"+accession)
			mismatch_outfile.write("\t"+accession)
			indel_outfile.write("\t"+accession)
		if num_allele == 1:
			cov_outfile.write("\n")
		prop_outfile.write("\n")
		readloc_outfile.write("\n")
		qual_outfile.write("\n")
		mapqual_outfile.write("\n")
		mismatch_outfile.write("\n")
		indel_outfile.write("\n")
		for seg_ref in ref_set_list:
			seg = seg_ref.split("-")[0]
			segment = seg.split(".")[0]
			for loc in range(0,max_length_dict[segment]):
				if num_allele == 1:
					cov_outfile.write(segment+"\t"+str(loc))
				prop_outfile.write(segment+"\t"+str(loc))
				readloc_outfile.write(segment+"\t"+str(loc))
				qual_outfile.write(segment+"\t"+str(loc))
				mapqual_outfile.write(segment+"\t"+str(loc))
				mismatch_outfile.write(segment+"\t"+str(loc))
				indel_outfile.write(segment+"\t"+str(loc))
				for a in range(0,len(accession_list)):
					accession = accession_list[a]
					runID = sample_run_dict[accession]
					cov_thresh = min_cov_dict[runID]
					try:
						minor_prop_stats = minor_prop_stat_dict[segment][accession]
						minor_prop_mean,minor_prop_std,minor_prop_median = minor_prop_stats[0],minor_prop_stats[1],minor_prop_stats[2]
					except:
						minor_prop_mean,minor_prop_std,minor_prop_median = 'nan','nan','nan'

					nt_prop_list = []
					for nt_num in range(0,4):
						try:
							nt_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][nt_num]
							nt_prop_list.append((nt_prop,nt_list[nt_num]))
						except:
							pass

					cov = 0
					allele_freq = ''
					read_loc = ''
					map_qual = ''
					rel_map_qual = ''
					mismat = ''
					indel = ''
					if nt_prop_list != []:
						nt_prop_list = sorted(nt_prop_list,reverse=True)
						major_nt = nt_prop_list[0][1]
						minor_nt = nt_prop_list[num_allele][1]
						try:
							cov = all_allele_info_dict[accession][segment][loc][info_tup_dict['cov']]
							if cov >= cov_thresh:
								allele_freq = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[minor_nt]]
								read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[minor_nt]]
								qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[minor_nt]]
								map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[minor_nt]]
								mismat = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[minor_nt]]
								indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[minor_nt]]
							else:
								allele_freq = ''
								read_loc = ''
								qual = ''
								map_qual = ''
								rel_map_qual = ''
								mismat = ''
								indel = ''
						except:
							pass
					if num_allele == 1:
						cov_outfile.write("\t"+str(cov))
					prop_outfile.write("\t"+str(allele_freq))
					readloc_outfile.write("\t"+str(read_loc))
					qual_outfile.write("\t"+str(qual))
					mapqual_outfile.write("\t"+str(map_qual))
					mismatch_outfile.write("\t"+str(mismat))
					indel_outfile.write("\t"+str(indel))
				if num_allele == 1:
					cov_outfile.write("\n")
				prop_outfile.write("\n")
				readloc_outfile.write("\n")
				qual_outfile.write("\n")
				mapqual_outfile.write("\n")
				mismatch_outfile.write("\n")
				indel_outfile.write("\n")
		if num_allele == 1:
			cov_outfile.close()
		prop_outfile.close()
		readloc_outfile.close()
		qual_outfile.close()
		mapqual_outfile.close()
		mismatch_outfile.close()
		indel_outfile.close()





debug = open("sub_labels."+output_suffix+".txt","w")
sub_loc_prop_dict = {}
sub_loc_min_prop_dict = {}
sub_ID_dict = {}
sub_label_dict = {}
sub_label_list = []
non_dichotomous = 0
poly_in_count = 0
sub_count = 0
sub_fail_count = 0

if iSNV_global_naming == False:
	for pigID in poly_sum_dict:
		# if pigID in flagged_pigs:
		sub_loc_prop_dict[pigID] = {}
		sub_loc_min_prop_dict[pigID] = {}
		sub_ID_dict[pigID] = {}
		sub_label_dict[pigID] = {}
		for sub_ID in poly_sum_dict[pigID]:
			poly_in_count += 1
			nt_prop_sum = poly_sum_dict[pigID][sub_ID]
			segment = sub_ID.split("-")[0]
			sub_loc = sub_ID.split("-")[1][0:-1]
			# print(pigID+"\t"+sub_ID+"\t"+segment+"\t"+str(sub_loc))
			seg_loc = segment+"-"+sub_loc
			nt = sub_ID.split("-")[1][-1]
			tup = (nt_prop_sum,nt)
			try:
				sub_loc_prop_dict[pigID][seg_loc].append(tup)
			except:
				sub_loc_prop_dict[pigID][seg_loc] = []
				sub_loc_prop_dict[pigID][seg_loc].append(tup)

			try:
				pig_iSNV_dict[pigID][sub_ID]
				try:
					sub_loc_min_prop_dict[pigID][seg_loc].append(tup)
				except:
					sub_loc_min_prop_dict[pigID][seg_loc] = []
					sub_loc_min_prop_dict[pigID][seg_loc].append(tup)
			except:
				pass

		for seg_loc in sub_loc_prop_dict[pigID]:
			temp_sub_list = sorted(sub_loc_prop_dict[pigID][seg_loc], reverse=True)
			try:
				num_minprop_alleles = len(sub_loc_min_prop_dict[pigID][seg_loc])
			except:
				num_minprop_alleles = 0
			if num_minprop_alleles == 2:
				segment = seg_loc.split("-")[0]
				sub_loc = seg_loc.split("-")[1]
				major_subID = segment+"-"+sub_loc+temp_sub_list[0][1]
				minor_subID = segment+"-"+sub_loc+temp_sub_list[1][1]
				sub_label = segment+"-"+temp_sub_list[0][1]+sub_loc+temp_sub_list[1][1]
				sub_ID_dict[pigID][major_subID] = sub_label
				sub_ID_dict[pigID][minor_subID] = sub_label
				sub_label_dict[pigID][sub_label] = minor_subID
				debug.write(pigID+"\t"+major_subID+"\t"+minor_subID+"\t"+sub_label+"\n")
				sub_label_list.append(sub_label)
				sub_count += 1
			elif num_minprop_alleles > 2:
				non_dichotomous += 1
				segment = seg_loc.split("-")[0]
				sub_loc = seg_loc.split("-")[1]
			elif num_minprop_alleles == 1:
				sub_fail_count += 1
elif iSNV_global_naming == True:
	for sub_ID in poly_sum_dict:
		poly_in_count += 1
		nt_prop_sum = poly_sum_dict[sub_ID]
		segment = sub_ID.split("-")[0]
		sub_loc = sub_ID.split("-")[1][0:-1]
		seg_loc = segment+"-"+sub_loc
		nt = sub_ID.split("-")[1][-1]
		tup = (nt_prop_sum,nt)
		try:
			sub_loc_prop_dict[seg_loc].append(tup)
		except:
			sub_loc_prop_dict[seg_loc] = []
			sub_loc_prop_dict[seg_loc].append(tup)
		for pigID in pig_iSNV_dict:
			try:
				pig_iSNV_dict[pigID][sub_ID]
				try:
					sub_loc_min_prop_dict[pigID][seg_loc].append(nt)
				except:
					try:
						sub_loc_min_prop_dict[pigID][seg_loc] = []
						sub_loc_min_prop_dict[pigID][seg_loc].append(nt)
					except:
						sub_loc_min_prop_dict[pigID] = {}
						sub_loc_min_prop_dict[pigID][seg_loc] = []
						sub_loc_min_prop_dict[pigID][seg_loc].append(nt)
			except:
				pass

	for seg_loc in sub_loc_prop_dict:
		temp_sub_list = sorted(sub_loc_prop_dict[seg_loc], reverse=True)
		segment = seg_loc.split("-")[0]
		sub_loc = seg_loc.split("-")[1]
		major_nt = temp_sub_list[0][1]
		major_subID = segment+"-"+sub_loc+major_nt
		for pigID in sub_loc_min_prop_dict:
			try:
				num_minprop_alleles = len(sub_loc_min_prop_dict[pigID][seg_loc])
			except:
				num_minprop_alleles = 0
			if num_minprop_alleles == 2:
				if sub_loc_min_prop_dict[pigID][seg_loc][0] == major_nt:
					minor_nt = sub_loc_min_prop_dict[pigID][seg_loc][1]
				elif sub_loc_min_prop_dict[pigID][seg_loc][1] == major_nt:
					minor_nt = sub_loc_min_prop_dict[pigID][seg_loc][0]
				minor_subID = segment+"-"+sub_loc+minor_nt
				sub_label = segment+"-"+temp_sub_list[0][1]+sub_loc+temp_sub_list[1][1]
				try:
					sub_ID_dict[pigID][major_subID] = sub_label
				except:
					sub_ID_dict[pigID] = {}
					sub_ID_dict[pigID][major_subID] = sub_label
				try:
					sub_ID_dict[pigID][minor_subID] = sub_label
				except:
					sub_ID_dict[pigID] = {}
					sub_ID_dict[pigID][minor_subID] = sub_label
				try:
					sub_label_dict[pigID][sub_label] = minor_subID
				except:
					sub_label_dict[pigID] = {}
					sub_label_dict[pigID][sub_label] = minor_subID
				try:
					S_NS_minor = S_NS_dict[minor_subID]
				except:
					S_NS_minor = 'nan'
				try:
					S_NS_major = S_NS_dict[major_subID]
				except:
					S_NS_major = 'nan'
				debug.write(pigID+"\t"+segment+"\t"+sub_loc+"\t"+major_nt+"\t"+minor_nt+"\t"+major_subID+"\t"+minor_subID+"\t"+sub_label+"\t"+S_NS_major+"\t"+S_NS_minor+"\t"+"\n")
				sub_label_list.append(sub_label)
				sub_count += 1
			elif num_minprop_alleles > 2:
				non_dichotomous += 1
				# segment = seg_loc.split("-")[0]
				# sub_loc = seg_loc.split("-")[1]
				# debug.write(pigID+"\t"+seg_loc+"\t"+segment+"\t"+str(sub_loc)+"\n"+str(temp_sub_list)+"\n")
			elif num_minprop_alleles == 1:
				sub_fail_count += 1
				# if sub_fail_count == 10:
				# 	print(pigID+"\t"+seg_loc+"\t"+segment+"\t"+str(sub_loc)+"\t"+str(temp_sub_list)+"\t"+str(sub_loc_min_prop_dict[pigID][seg_loc]))
debug.close()

sub_label_list = list(set(sub_label_list))
print("non_dichotomous: "+str(non_dichotomous))
print("poly_in_count: "+str(poly_in_count))
print("sub_count: "+str(sub_count))
print("no_sub_sites: "+str(sub_fail_count))
print("len(sub_label_list): "+str(len(sub_label_list)))


group_list.append("other")
group_list = sorted(list(set(group_list)))
outfile = open(project_dir+"swineIAV_group_count."+output_suffix+".self.txt","w")
outfile.write("\t")
for i in range(0,len(group_list)):
	group = group_list[i]
	outfile.write("\t"+group)
outfile.write("\tpoly_count\tsub_count\n")

for a in range(0,len(accession_list)):
	accession = accession_list[a]
	for s in range(0,len(segment_list)):
		segment = segment_list[s]
		outfile.write(accession+"\t"+segment)
		for i in range(0,len(group_list)):
			group = group_list[i]
			try:
				count = group_count_dict[accession][segment][group]
			except:
				count = 0
			outfile.write("\t"+str(count))
		try:
			segregating_site_count = site_count_all_passcov[accession][segment]
		except:
			segregating_site_count = 0
		outfile.write("\t"+str(segregating_site_count))
		try:
			sub_count = len(poly_dict[accession])
		except:
			sub_count = 0
		try:
			avg_cov = round(np.average(avg_cov_dict[accession][segment]),2)
		except:
			avg_cov = 0.0
		outfile.write("\t"+str(sub_count)+"\t"+str(avg_cov)+"\n")
outfile.close()


outfile = open(project_dir+"swineIAV_group_count_table."+output_suffix+".self.txt","w")
for i in range(0,len(group_list)):
	group = group_list[i]
	for s in range(0,len(internal_segment_list)):
		segment = internal_segment_list[s].split(".")[0]
		outfile.write("\t"+group+"-"+segment)
outfile.write("\n")
for a in range(0,len(accession_list)):
	accession = accession_list[a]
	outfile.write(accession)
	for i in range(0,len(group_list)):
		group = group_list[i]
		for s in range(0,len(internal_segment_list)):
			segment = internal_segment_list[s].split(".")[0]
			try:
				count = group_count_dict[accession][segment][group]
			except:
				count = 0
			outfile.write("\t"+str(count))
	outfile.write("\n")
outfile.close()

outfile = open(project_dir+"swineIAV_avg_cov_table."+output_suffix+".self.txt","w")
for s in range(0,len(ref_list)):
	segment = ref_list[s].split(".")[0]
	outfile.write("\t"+segment)
outfile.write("\n")
for a in range(0,len(accession_list)):
	accession = accession_list[a]
	outfile.write(accession)
	for s in range(0,len(ref_list)):
		segment = ref_list[s].split(".")[0]
		try:
			avg_cov = round(np.average(avg_cov_dict[accession][segment]),2)
		except:
			avg_cov = 0.0
		outfile.write("\t"+str(avg_cov))
	outfile.write("\n")
outfile.close()

temp_sub_list = []
for i in range(0,len(accession_list)):
	accession = accession_list[i]
	runID = sample_run_dict[accession]
	cov_thresh = min_cov_dict[runID]
	pigID = sample_dict[accession][0]
	# for j in range(0,len(sub_list)):
	for sub_ID in pig_iSNV_dict[pigID]:
		# sub_ID = sub_list[j]
		segment = sub_ID.split("-")[0]
		loc = sub_ID.split("-")[1]
		loc = int(loc[0:len(loc)-1])
		sub_nt = sub_ID[-1]
		try:
			site_cov = all_allele_info_dict[accession][segment][loc][info_tup_dict['cov']]
			sub_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[sub_nt]]
			sub_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[sub_nt]]
			sub_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[sub_nt]]
			sub_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[sub_nt]]
			# sub_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[sub_nt]]
			sub_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][base_tup_dict[sub_nt]]
			sub_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[sub_nt]]
			sub_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[sub_nt]]
			# if sub_prop >= min_proportion and sub_qual >= min_qual and sub_map_qual >= min_avg_map_qual and sub_read_loc >= min_avg_read_loc and sub_read_prop_loc >= min_avg_read_prop_loc and sub_read_prop_loc <= (1.0-min_avg_read_prop_loc):
			if site_cov >= cov_thresh:
				if sub_prop >= min_proportion and sub_qual >= min_qual and sub_map_qual >= min_avg_map_qual and sub_read_loc >= min_avg_read_loc and sub_read_mismatch <= max_avg_read_mismatch and sub_read_indel <= max_avg_read_indel:
					temp_sub_list.append(sub_ID)
		except:
			pass

temp_sub_list = list(set(temp_sub_list))
outfile = open(project_dir+"swineIAV_iSNV_table."+output_suffix+".subID_self.txt","w")
for j in range(0,len(temp_sub_list)):
	outfile.write("\t"+temp_sub_list[j])
outfile.write("\n")
for i in range(0,len(accession_list)):
	accession = accession_list[i]
	runID = sample_run_dict[accession]
	cov_thresh = min_cov_dict[runID]
	pigID = sample_dict[accession][0]
	outfile.write(accession)
	for j in range(0,len(temp_sub_list)):
		sub_ID = temp_sub_list[j]
		segment = sub_ID.split("-")[0]
		loc = sub_ID.split("-")[1]
		loc = int(loc[0:len(loc)-1])
		sub_nt = sub_ID.split("-")[1][-1]
		try:
			site_cov = all_allele_info_dict[accession][segment][loc][info_tup_dict['cov']]
		except:
			site_cov = 0
		if site_cov >= cov_thresh:
			try:
				sub_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[sub_nt]]
				sub_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[sub_nt]]
				sub_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[sub_nt]]
				sub_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[sub_nt]]
				sub_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[sub_nt]]
				sub_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][base_tup_dict[sub_nt]]
				sub_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[sub_nt]]
				sub_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[sub_nt]]
				qual_pass = False
				if sub_qual >= min_qual and sub_map_qual >= min_avg_map_qual and sub_read_loc >= min_avg_read_loc and sub_read_mismatch <= max_avg_read_mismatch and sub_read_indel <= max_avg_read_indel:
					qual_pass = True
				
				if sub_prop < min_proportion:# and qual_pass == True:
					sub_prop = 0.0
				elif qual_pass == True:
					sub_prop = round_to_n_sig_figs(sub_prop,3)
				if qual_pass == False:
					sub_prop = 'na'		
			except:
				sub_prop = 'na'
		else:
			sub_prop = 'na'
		outfile.write("\t"+str(sub_prop))
	outfile.write("\n")
outfile.close()



debug_outfile = open(project_dir+"debug."+output_suffix+".missed.txt","w")
count1 = 0
count2 = 0
outfile = open(project_dir+"swineIAV_iSNV_table."+output_suffix+".subLabel_self.txt","w")
for j in range(0,len(sub_label_list)):
	outfile.write("\t"+sub_label_list[j])
outfile.write("\n")
for i in range(0,len(accession_list)):
	accession = accession_list[i]
	runID = sample_run_dict[accession]
	cov_thresh = min_cov_dict[runID]
	pigID = sample_dict[accession][0] #(pigID,time)
	time = sample_dict[accession][1] #(pigID,time)
	outfile.write(accession)
	for j in range(0,len(sub_label_list)):
		sub_label = sub_label_list[j]
		try:
			minor_subID = sub_label_dict[pigID][sub_label]
		except:
			minor_subID = "na"
		segment = sub_label.split("-")[0]
		loc = sub_label.split("-")[1]
		loc = int(loc[1:-1])

		try:
			site_cov = all_allele_info_dict[accession][segment][loc][info_tup_dict['cov']]
		except:
			site_cov = 0
		if site_cov < cov_thresh:
			minor_sub_prop = "nan"
		else:
			if minor_subID != "na":
				minor_subID_nt = minor_subID[-1]
				major_subID_nt = sub_label.split("-")[1][0]
				try:
					minor_sub_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[minor_subID_nt]]
					major_sub_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[major_subID_nt]]

					minor_sub_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[minor_subID_nt]]
					minor_sub_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[minor_subID_nt]]
					minor_sub_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[minor_subID_nt]]
					minor_sub_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[minor_subID_nt]]
					minor_sub_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[minor_subID_nt]]
					# minor_sub_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[minor_subID_nt]]
					# minor_sub_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][base_tup_dict[minor_subID_nt]]

					major_sub_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[major_subID_nt]]
					major_sub_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[major_subID_nt]]
					major_sub_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[major_subID_nt]]
					major_sub_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[major_subID_nt]]
					major_sub_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[major_subID_nt]]
					major_sub_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[major_subID_nt]]
					# major_sub_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[major_subID_nt]]
					# major_sub_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][base_tup_dict[major_subID_nt]]
					qual_pass = False
					if minor_sub_qual >= min_qual and minor_sub_map_qual >= min_avg_map_qual and minor_sub_read_loc >= min_avg_read_loc and minor_sub_read_mismatch <= max_avg_read_mismatch and minor_sub_read_indel <= max_avg_read_indel:
						if major_sub_qual >= min_qual and major_sub_map_qual >= min_major_map_qual and major_sub_read_loc >= min_avg_read_loc and major_sub_read_mismatch <= max_major_read_mismatch and major_sub_read_indel <= max_major_read_indel:
							qual_pass = True
					if minor_sub_prop < min_proportion:# or minor_sub_qual < min_qual or minor_sub_map_qual < min_avg_map_qual or minor_sub_read_loc < min_avg_read_loc or minor_sub_read_prop_loc < min_avg_read_prop_loc or minor_sub_read_prop_loc > (1.0-min_avg_read_prop_loc) or minor_sub_read_mismatch > max_mismatch_count:#rel_map_qual < min_rel_map_qual or 
						minor_sub_prop = 0.0
					elif qual_pass == True:
						minor_sub_prop = round_to_n_sig_figs(minor_sub_prop,3)
					if qual_pass == False:
						minor_sub_prop = "nan"
				except:
					minor_sub_prop = "nan"
					count1 += 1
			else:
				minor_sub_prop = "nan"
				count2 += 1
				ref_nt = sub_label.split("-")[1][0]
				sub_nt = sub_label.split("-")[1][-1]
				segment = sub_label.split("-")[0]
				subloc_nt = sub_label.split("-")[1][-1]
				subloc = int(sub_label.split("-")[1][1:-1])
				dSNP_loc_found = "False"
				try:
					dSNP_dict[segment][subloc]
					dSNP_loc_found = 'True'
				except:
					dSNP_loc_found = 'False'
				try:
					minor_sub_prop = all_allele_info_dict[accession][segment][subloc][info_tup_dict['prop']][base_tup_dict[sub_nt]]
					minor_sub_qual = all_allele_info_dict[accession][segment][subloc][info_tup_dict['qual']][base_tup_dict[sub_nt]]
					minor_sub_map_qual = all_allele_info_dict[accession][segment][subloc][info_tup_dict['map_qual']][base_tup_dict[sub_nt]]
					minor_sub_read_loc = all_allele_info_dict[accession][segment][subloc][info_tup_dict['read_loc']][base_tup_dict[sub_nt]]
					major_sub_prop = all_allele_info_dict[accession][segment][subloc][info_tup_dict['prop']][base_tup_dict[ref_nt]]
					major_sub_qual = all_allele_info_dict[accession][segment][subloc][info_tup_dict['qual']][base_tup_dict[ref_nt]]
					major_sub_map_qual = all_allele_info_dict[accession][segment][subloc][info_tup_dict['map_qual']][base_tup_dict[ref_nt]]
					major_sub_read_loc = all_allele_info_dict[accession][segment][subloc][info_tup_dict['read_loc']][base_tup_dict[ref_nt]]
					major_sub_mismatch = all_allele_info_dict[accession][segment][subloc][info_tup_dict['mismatch']][base_tup_dict[ref_nt]]
					major_sub_indel = all_allele_info_dict[accession][segment][subloc][info_tup_dict['indel']][base_tup_dict[ref_nt]]
					if minor_sub_prop < min_proportion:
						minor_sub_prop = 0.0
					if minor_sub_qual < min_qual or minor_sub_map_qual < min_avg_map_qual or minor_sub_read_loc < min_avg_read_loc or major_sub_map_qual < min_major_map_qual or major_sub_mismatch > max_major_read_mismatch or major_sub_indel > max_major_read_indel:
						minor_sub_prop = "nan"
				except:
					minor_sub_prop = "nan"
					debug_outfile.write(accession+"\t"+pigID+"\t"+str(time)+"\t"+segment+"\t"+str(subloc)+"\t"+sub_label+"\t"+dSNP_loc_found+"\t"+str(site_cov)+"\t"+
						str(major_sub_prop)+"\t"+str(major_sub_qual)+"\t"+str(major_sub_map_qual)+"\t"+str(major_sub_read_loc)+"\t"+
						str(minor_sub_prop)+"\t"+str(minor_sub_qual)+"\t"+str(minor_sub_map_qual)+"\t"+str(minor_sub_read_loc)+"\n")
		outfile.write("\t"+str(minor_sub_prop))
	outfile.write("\n")
outfile.close()
debug_outfile.close()

print("no entry in all_allele_info_dict count: "+str(count1))
print("minor_subID=='na' count: "+str(count2))


sublabel_S_NS_dict = {}
for pigID in sub_label_dict:
	for sub_label in sub_label_dict[pigID]:
		minor_subID = sub_label_dict[pigID][sub_label]
		try:
			sub_type = S_NS_dict[minor_subID]
			sublabel_S_NS_dict[sub_label] = sub_type
		except:
			pass
print("len(sublabel_S_NS_dict): "+str(len(sublabel_S_NS_dict)))

print(pig_to_samples)
print(sub_label_dict)
outfile = open(project_dir+"swineIAV_iSNV_table.reformat."+output_suffix+".txt","w")
outfile.write("pigID\tsegment\tref_nt\tiSNV_loc\tiSNV_nt")
for time in range(1,8):
	outfile.write("\td"+str(time-1)+"_tot\td"+str(time-1)+"_sub")
outfile.write("\tNS-S\tsubtype\n")
for pigID in sub_label_dict:
	for sub_label in sub_label_dict[pigID]:
		minor_subID = sub_label_dict[pigID][sub_label]
		ref_nt = sub_label.split("-")[1][0]
		sub_nt = sub_label.split("-")[1][-1]
		loc = int(sub_label.split("-")[1][1:-1])
		segment = sub_label.split("-")[0]
		dSNP_loc = False
		try:
			dSNP_dict[segment][loc]
			dSNP_loc = True
		except:
			dSNP_loc = False
		if dSNP_loc == False:
			out_line = pigID+"\t"+segment+"\t"+ref_nt+"\t"+str(loc)+"\t"+sub_nt
			alt_out_line = pigID+"\t"+segment+"\t"+sub_nt+"\t"+str(loc)+"\t"+ref_nt
			first_minor_nt = ''
			subcount = 0
			for time in range(1,8):
				try:
					accession = pig_to_samples[pigID][time-1]
				except:
					accession = ''
				try:
					runID = sample_run_dict[accession]
				except:
					runID = ''
				try:
					total_site_cov = int(site_cov_dict[accession][segment][loc])
				except:
					total_site_cov = "nan"
				try:
					sub_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[sub_nt]]
					sub_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[sub_nt]]
					sub_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[sub_nt]]
					sub_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[sub_nt]]
					sub_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[sub_nt]]
					sub_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[sub_nt]]
					sub_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[sub_nt]]
					sub_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][base_tup_dict[sub_nt]]
					
					major_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[ref_nt]]
					major_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[ref_nt]]
					major_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[ref_nt]]
					major_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[ref_nt]]
					major_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[ref_nt]]
					major_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[ref_nt]]
					major_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[ref_nt]]
					major_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][base_tup_dict[ref_nt]]
					
					qual_pass = False
					if sub_qual >= min_qual and sub_map_qual >= min_avg_map_qual and sub_read_loc >= min_avg_read_loc and sub_read_mismatch <= max_avg_read_mismatch and sub_read_indel <= max_avg_read_indel:
						if major_qual >= min_qual and major_map_qual >= min_major_map_qual and major_read_loc >= min_avg_read_loc and major_read_mismatch <= max_major_read_mismatch and major_read_indel <= max_major_read_indel:
							qual_pass = True
							if first_minor_nt == '':
								if sub_prop <= major_prop:
									first_minor_nt = sub_nt
								else:
									first_minor_nt = ref_nt
					if qual_pass == False:
						sub_prop = 'nan'
				except:
					sub_prop = 'nan'
					major_prop = 'nan'
					qual_pass = False
				if runID != "":
					if total_site_cov != "nan" and sub_prop != 'nan':
						if total_site_cov >= min_cov_dict[runID]:
							iSNV_site_cov = int(sub_prop*total_site_cov)
							alt_iSNV_site_cov = int(major_prop*total_site_cov)
							if sub_prop >= min_proportion and major_prop >= min_proportion:
								subcount += 1
					else:
						iSNV_site_cov = 'nan'
						alt_iSNV_site_cov = 'nan'
				else:
					total_site_cov = 'nan'
					iSNV_site_cov = 'nan'
					alt_iSNV_site_cov = 'nan'
					count7 += 1
				out_line += "\t"+str(total_site_cov)+"\t"+str(iSNV_site_cov)
				alt_out_line += "\t"+str(total_site_cov)+"\t"+str(alt_iSNV_site_cov)
			if subcount >=1:
				try:
					groupID = pig_group_dict[pigID]
				except:
					groupID = "nan"
					print("noGroupID\t"+pigID+"\t"+str(time-1)+"\t"+segment+"\t"+accession)
				try:
					sub_type = sublabel_S_NS_dict[sub_label]
				except:
					sub_type = "nan"
				if first_minor_nt == sub_nt:
					outfile.write(out_line+"\t"+sub_type+"\t"+groupID+"\n")
				elif first_minor_nt == ref_nt:
					outfile.write(alt_out_line+"\t"+sub_type+"\t"+groupID+"\n")
outfile.close()


outfile = open(project_dir+"S_NS_subs.self."+output_suffix+".txt","w")
outfile2 = open(project_dir+"all_subs.self."+output_suffix+".txt","w")
for pigID in sub_label_dict:
	pig_accession_list = []
	for time in range(1,8):
		try:
			accession = pig_to_samples[pigID][time-1]
			pig_accession_list.append(accession)
		except:
			accession = ''
	pig_accession_list = list(set(pig_accession_list))
	for sub_label in sub_label_dict[pigID]:
		minor_subID = sub_label_dict[pigID][sub_label]
		segment = sub_label.split("-")[0]
		ref_nt = sub_label.split("-")[1][0]
		sub_nt = sub_label.split("-")[1][-1]
		loc = int(sub_label.split("-")[1][1:-1])
		for accession in pig_accession_list:
			runID = sample_run_dict[accession]
			cov_thresh = min_cov_dict[runID]
			try:
				sub_type = sublabel_S_NS_dict[sub_label]
			except:
				sub_type = 'na'
			try:
				site_cov = int(site_cov_dict[accession][segment][loc])
			except:
				site_cov = 0
			if site_cov >= cov_thresh:
				try:
					sub_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[sub_nt]]
					sub_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[sub_nt]]
					sub_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[sub_nt]]
					sub_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[sub_nt]]
					sub_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[sub_nt]]
					sub_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[sub_nt]]
					sub_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[sub_nt]]
					sub_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][base_tup_dict[sub_nt]]

					major_prop = all_allele_info_dict[accession][segment][loc][info_tup_dict['prop']][base_tup_dict[ref_nt]]
					major_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['qual']][base_tup_dict[ref_nt]]
					major_map_qual = all_allele_info_dict[accession][segment][loc][info_tup_dict['map_qual']][base_tup_dict[ref_nt]]
					major_read_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_loc']][base_tup_dict[ref_nt]]
					major_read_R_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_R_loc']][base_tup_dict[ref_nt]]
					major_read_mismatch = all_allele_info_dict[accession][segment][loc][info_tup_dict['mismatch']][base_tup_dict[ref_nt]]
					major_read_indel = all_allele_info_dict[accession][segment][loc][info_tup_dict['indel']][base_tup_dict[ref_nt]]
					major_read_prop_loc = all_allele_info_dict[accession][segment][loc][info_tup_dict['read_prop_loc']][base_tup_dict[ref_nt]]
					
					qual_pass = False
					if sub_qual >= min_qual and sub_map_qual >= min_avg_map_qual and sub_read_loc >= min_avg_read_loc and sub_read_mismatch <= max_avg_read_mismatch and sub_read_indel <= max_avg_read_indel:
						if major_qual >= min_qual and major_map_qual >= min_major_map_qual and major_read_loc >= min_avg_read_loc and major_read_mismatch <= max_major_read_mismatch and major_read_indel <= max_major_read_indel:
							qual_pass = True

					if qual_pass == False:
						sub_prop = 'na'
					elif sub_prop < min_proportion:
						sub_prop = 0
				except:
					sub_prop = 'na'
					sub_qual = 'na'
					sub_map_qual = 'na'
					sub_read_loc = 'na'
					sub_read_R_loc = 'na'
					sub_read_prop_loc = 'na'
					rel_map_qual = 'na'
			else:
				sub_prop = 'na'
				sub_qual = 'na'
				sub_map_qual = 'na'
				sub_read_loc = 'na'
				sub_read_R_loc = 'na'
				sub_read_prop_loc = 'na'
				rel_map_qual = 'na'
			if sub_type != 'na':
				outfile.write(pigID+"\t"+accession+"\t"+sub_label+"\t"+sub_type+"\t"+str(sub_prop)+"\t"+str(site_cov)+"\t"+str(sub_qual)+"\t"+str(sub_map_qual)+"\t"+str(sub_read_loc)+"\t"+str(sub_read_R_loc)+"\n")#+"\t"+str(sub_read_prop_loc)+"\t"+str(rel_map_qual)+"\n")
			outfile2.write(pigID+"\t"+accession+"\t"+sub_label+"\t"+sub_type+"\t"+str(sub_prop)+"\t"+str(site_cov)+"\t"+str(sub_qual)+"\t"+str(sub_map_qual)+"\t"+str(sub_read_loc)+"\t"+str(sub_read_R_loc)+"\n")#+"\t"+str(sub_read_prop_loc)+"\t"+str(rel_map_qual)+"\n")
outfile.close()
outfile2.close()


