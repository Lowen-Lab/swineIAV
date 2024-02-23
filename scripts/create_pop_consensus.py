###################################   USER DEFINED VARIABLES   ###################################
project_dir = "./"
group_filename = "IAV_groups.all.txt"
outfile_name = "swineIAV_consensus_seqs.single.fasta"
outfile = open(project_dir+outfile_name,"w")

max_other_count = 1

ref_set_list = ['H3N2_PB2.fa-H1N1_PB2.fa','H3N2_PB1.fa-H1N1_PB1.fa','H3N2_PA.fa-H1N1_PA.fa','H3N2_NP.fa-H1N1_NP.fa','H3N2_M.fa-H1N1_M.fa','H3N2_NEP.fa-H1N1_NEP.fa','H3N2_HA.fa','H3N2_NA.fa','H1N1_HA.fa','H1N1_NA.fa']


group_dict = {}
group_list = []
group_strain_lists = {}
infile = open(project_dir+group_filename,"r")
for line in infile:
	line = line.strip().replace("group1","H1N1").replace("group2","H3N2").replace("H3N2_NS","H3N2_NEP").split("\t")
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
infile.close()
group_list = list(set(group_list))
print(group_dict)
first = True

for seg_refs in ref_set_list:
	MSA_infile_name = seg_refs.split("-")[0].replace('.fa','.consensus.fa')
	segment = MSA_infile_name.split(".")[0]
	print(segment+"\t"+MSA_infile_name)
	accession_list = []
	msa_infile = open(project_dir+""+MSA_infile_name,"r")
	msa_dict = {}
	skip = True
	for line in msa_infile:
		line = line.strip()
		if line[0] == ">":
			accession = line[1:len(line)]
			try:
				group_dict[segment][accession]
				msa_dict[accession] = ''
				skip = False
				accession_list.append(accession)
			except:
				skip = True
		elif skip == False:
			try:
				msa_dict[accession] += line
			except:
				msa_dict[accession] = line
	msa_infile.close()

	len_align = len(msa_dict[accession_list[2]])


	keep_column_list = []
	outseq_dict = {}
	ref_nt_list = ['A','C','G','T']
	for column in range(0,len_align):
		nt_count = {}
		base_list = []
		for accession in msa_dict:
			group = group_dict[segment][accession]
			nt = msa_dict[accession][column]
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


		for group in group_list:
			loc_nt_list = []
			for nt in ref_nt_list:
				try:
					tup = (nt_count[group][nt],nt)
					loc_nt_list.append(tup)
				except:
					pass
			if len(loc_nt_list) > 0:
				loc_nt_list = sorted(loc_nt_list, reverse=True)
				all_base_count = (nt_count[group]['A']+nt_count[group]['C']+nt_count[group]['G']+nt_count[group]['T']+nt_count[group]['N']+nt_count[group]['-'])
				if nt_count[group]['-']/all_base_count <= 0.25:
					if loc_nt_list[0][0]/all_base_count >= 0.6:
						try:
							outseq_dict[group] += loc_nt_list[0][1]
						except:
							outseq_dict[group] = loc_nt_list[0][1]
					else:
						try:
							outseq_dict[group] += "N"
						except:
							outseq_dict[group] = "N"
				elif nt_count[group]['-'] > 0:
					try:
						outseq_dict[group] += "-"
					except:
						outseq_dict[group] = "-"

	for group in group_list:
		try:
			outseq = outseq_dict[group]
			outfile.write(">"+group+"-"+segment+"\n"+outseq+"\n")
		except:
			pass
outfile.close()

