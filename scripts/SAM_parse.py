import os
import sys
import math
import random
import numpy as np
import time

##################################  USER DEFINED VARIABLES  #######################################
dedupe_reads = True
parallel_process = True
parallel_max_cpu = -2
sampleID_file_name = "accessions.txt"

project_dir = "./"
input_read_dir =project_dir+"raw_data/"
ref_dir = project_dir+"refs/"
trim_dir = project_dir+"trim/"
map_folder = 'rep_map'
map_dir = project_dir+map_folder+"/"
qual_dir = map_dir+"qual/"
sam_dir = map_dir+"sam/"
sub_dir = map_dir+"substitutions/"
temp_dir = map_dir+"temp/"


ref_set_list = ['H3N2_PB2.fa-H1N1_PB2.fa','H3N2_PB1.fa-H1N1_PB1.fa','H3N2_PA.fa-H1N1_PA.fa','H3N2_NP.fa-H1N1_NP.fa','H3N2_M.fa-H1N1_M.fa','H3N2_NEP.fa-H1N1_NEP.fa','H3N2_HA.fa','H3N2_NA.fa','H1N1_HA.fa','H1N1_NA.fa']

min_proportion = 0.03 # for read allele summary files only, if enabled

###################################################################################################
####
#Set up parallel processing if enabled
if parallel_process == True:
	import multiprocessing
	from joblib import Parallel, delayed
	max_cores = multiprocessing.cpu_count()
	if parallel_max_cpu == 0:
		num_cores = max_cores
	elif parallel_max_cpu > 0:
		num_cores = min(parallel_max_cpu,max_cores)
	elif parallel_max_cpu < 0:
		num_cores = max_cores+parallel_max_cpu
else:
	num_cores = 1
####
############################################ FUNCTIONS ############################################
def run_command(command,mode="quiet"):
	'''
	The quiet function on this command redirects the bash shell output (and error messages) to dev/null, 
	which is quiet, but will prevent some programs from working properly. such as samtools. So, if you
	are debugging or using an affected program, send commands with the "verbose" flag
	'''
	run_status = False
	if mode == "verbose":
		return_code = os.system(command)
	elif mode == "quiet":
		return_code = os.system(command+" >/dev/null 2>&1")
	if return_code == 0:
		run_status = True
	return run_status

def run_command(command):
	run_status = False
	# return_code = os.system(command)
	return_code = os.system(command+" >/dev/null 2>&1")
	# subprocess.check_output(command,stderr=subprocess.STDOUT)
	if return_code == 0:
		run_status = True
	return run_status


def check_if_reads_mapped(sampleID,ref_set,sam_dir,read_count_dict):
	reads_mapped = True
	# for ref_set in ref_set_list:
	refs = ref_set.split("-")
	primary_ref = ref_set.split("-")[0].split(".f")[0]
	for ref_filename in refs:
		refID = ref_filename.split('.f')[0]
		if os.path.isfile(sam_dir+refID+"/"+sampleID+"."+refID+'.sam') == False:
			try:
				read_count = read_count_dict[refID]
			except:
				read_count = 0
			if read_count >= 1000:
				reads_mapped = False
	return reads_mapped


def check_if_processing_finished(sampleID,ref_set,sub_dir):
	sam_processed = True
	# for ref_set in ref_set_list:
	refs = ref_set.split("-")
	for ref_filename in refs:
		refID = ref_filename.split(".f")[0]
		out_substitutions_filename = sub_dir+refID+"/"+sampleID+"."+refID+".substitions.txt"
		if os.path.isfile(out_substitutions_filename) == False:
			sam_processed = False
	return sam_processed


def load_read_counts(sampleID,trim_dir):
	read_count_dict = {}
	for ref_set in ref_set_list:
		refs = ref_set.split("-")
		primary_ref = ref_set.split("-")[0].split(".f")[0]
		for ref_filename in refs:
			segment = ref_filename.split(".f")[0]
			read_count_filename = trim_dir+segment+"/"+sampleID+'.'+segment+'.read_count.txt'
			if os.path.isfile(read_count_filename) == True:
				read_count_file = open(read_count_filename,'r')
				for line in read_count_file:
					line = line.strip()
					read_count_dict[segment] = int(line)
				read_count_file.close()
	return read_count_dict


def dedupe_samfile(sam_filename_in,temp_sam_filename,sam_filename_out):
	sam_infile = open(sam_filename_in,"r")
	read_map_strings = {}
	for sam_line in sam_infile:
		if len(sam_line) > 0:
			if sam_line[0] != "@":
				sam_line = sam_line.strip().split("\t")
				read_seq = sam_line[9]
				if read_seq !="*":
					readID = sam_line[0].split(" ")[0]
					map_site = sam_line[3]
					map_qual = sam_line[4]
					cigar = sam_line[5]
					map_string = map_site+"_"+map_qual+"_"+cigar
					try:
						read_map_strings[readID] += " "+map_string
					except:
						read_map_strings[readID] = " "+map_string
	sam_infile.close()
	map_string_dedupe = {}
	for readID in read_map_strings:
		map_string = read_map_strings[readID]
		if len(map_string.split(" ")[0]) == 2:
			map_string_R = map_string.split(" ")[0]+" "+map_string.split(" ")[1]
			try:
				one_readID = map_string_dedupe[map_string]
			except:
				try:
					one_readID = map_string_dedupe[map_string_R]
				except:
					map_string_dedupe[map_string] = readID
		else:
			try:
				one_readID = map_string_dedupe[map_string]
			except:
				map_string_dedupe[map_string] = readID
	del read_map_strings
	pass_readID_dict = {}
	for map_string in map_string_dedupe:
		readID = map_string_dedupe[map_string]
		pass_readID_dict[readID] = ''
	del map_string_dedupe
	dedupe_sam_outfile = open(temp_sam_filename,"w")
	sam_infile = open(sam_filename_in,"r")
	for sam_line in sam_infile:
		if len(sam_line) > 0:
			sam_line = sam_line.strip()#.split("\t")
			if sam_line[0] != "@":
				split_line = sam_line.split("\t")
				readID = split_line[0].split(" ")[0]
				read_seq = split_line[9]
				if read_seq !="*":
					try:
						pass_readID_dict[readID]
						dedupe_sam_outfile.write(sam_line+"\n")
					except:
						pass
			else:
				dedupe_sam_outfile.write(sam_line+"\n")
	sam_infile.close()
	dedupe_sam_outfile.close()
	os.rename(temp_sam_filename,sam_filename_out)
	return len(pass_readID_dict)


def SAM_parse(map_position,cigar_string,quality_string,read_seq,mapping_qual): #mapping position is indexed to 1, not zero
	cigar_dict = {'M':'imperfect_match','D':'insertion_in_ref','I':'insertion_in_read','S':'soft_clip_reference','H':'soft_clip_read','=':'perfect_match','X':'mismatch'}
	read_key = ''
	ref_site_map = []
	match_count_dict = {}
	cur_position = int(map_position)-2
	cur_num = ''
	cur_operator = ''
	all_valid_operators = True
	for i in range(0,len(cigar_string)):
		character = cigar_string[i]
		try:
			cigar_dict[character]
			cur_operator = character

			cur_num = int(cur_num)
			try:
				match_count_dict[cur_operator] += cur_num
			except:
				match_count_dict[cur_operator] = cur_num
			for num in range(0,cur_num):
				if cur_operator == "D":
					cur_position += 1
				elif cur_operator == "I":
					cur_position += 0
					ref_site_map.append('')
					read_key += cur_operator
				elif cur_operator == "S":
					cur_position += 0
					ref_site_map.append('')
					read_key += cur_operator
				elif cur_operator == "H":
					all_valid_operators = False
				elif cur_operator == "=": #perfect match
					cur_position += 1
					ref_site_map.append(cur_position)
					read_key += cur_operator
				elif cur_operator == "M": #imperfect match
					cur_position += 1
					ref_site_map.append(cur_position)
					read_key += cur_operator
				elif cur_operator == "X": #mismatch
					cur_position += 1
					ref_site_map.append(cur_position)
					read_key += cur_operator
				else:
					all_valid_operators = False
			cur_num = ''
			cur_operator = ''
		except:
			cur_num += character
	base_dict_fun = {}
	read_end_site = 0
	if all_valid_operators == True:
		for num in range(0,len(read_key)):
			operator = read_key[num]
			ref_loc = ref_site_map[num]
			if ref_loc != '':
				ref_loc = int(ref_loc)
				try:
					qual_char = quality_string[num]
				except:
					sys.exit("Index ")
				qual_score = ord(qual_char)-33
				read_nt = read_seq[num]
				read_length = len(read_seq)
				nt_prop_location = round(num/read_length,2)
				R_end_dist = read_length-num
				tup = (read_nt,qual_score,mapping_qual,read_length,num,R_end_dist,nt_prop_location)
				base_dict_fun[ref_loc] = tup

	return base_dict_fun,match_count_dict

###########################################################

def main_pipeline(sample_ref_tup,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,temp_dir,sub_dir,command_prefix,dedupe_reads,min_proportion):
	max_indel_count = 10
	max_mismatch_count = 15
	sampleID = sample_ref_tup[0]
	ref_set = sample_ref_tup[1]
	refs = ref_set.split("-")
	primary_refID = refs[0].split(".f")[0]

	read_count_dict = load_read_counts(sampleID,trim_dir)
	reads_mapped = check_if_reads_mapped(sampleID,ref_set,sam_dir,read_count_dict)
	if reads_mapped == False:
		return None
	
	map_files_processed = check_if_processing_finished(sampleID,ref_set,sub_dir)
	if map_files_processed == True:
		return None
	
	print(sampleID+" "+primary_refID)
	for ref_filename in refs:
		read_map_info_dict = {}
		map_qual_dict = {}
		indel_dict = {}
		mismatch_dict = {}
		refID = ref_filename.split(".f")[0]
		sam_filename = sam_dir+refID+"/"+sampleID+"."+refID+".sam"
		cont = True
		if dedupe_reads  == False:
			sam_filename_to_open = sam_filename
			if os.path.isfile(sam_filename) == False:
				cont = False
		elif dedupe_reads == True:
			if os.path.isfile(sam_filename) == True:
				temp_sam_filename = temp_dir+sampleID+"."+ref_filename.split(".f")[0]+".dedupe.sam"
				dedupe_samfilename = sam_dir+refID+"/"+sampleID+"."+ref_filename.split(".f")[0]+".dedupe.sam"
				num_reads_after_dedupe = dedupe_samfile(sam_filename,temp_sam_filename,dedupe_samfilename)
				sam_filename_to_open = dedupe_samfilename
			else:
				cont = False
		if cont == True:
			sam_infile = open(sam_filename_to_open,"r")
			for sam_line in sam_infile:
				if len(sam_line) > 0:
					if sam_line[0] != "@":
						sam_line = sam_line.strip().split("\t")
						map_qual = int(sam_line[4])
						read_seq = sam_line[9]
						if read_seq !="*":
							readID = sam_line[0].split(" ")[0]
							map_site = int(sam_line[3])
							cigar = sam_line[5] #{'M':'imperfect_match','D':'insertion_in_ref','I':'insertion_in_read','S':'soft_clip_reference','H':'soft_clip_read','=':'perfect_match','X':'mismatch'}
							read_seq_qual = sam_line[10]
							read_base_dict,CIGAR_operator_count_dict = SAM_parse(map_site,cigar,read_seq_qual,read_seq,map_qual)
							if len(read_base_dict) > 0: #per base in read: (read_nt,qual_score,mapping_qual,read_length,num,R_end_dist,nt_prop_location)
								try:
									in_count = CIGAR_operator_count_dict['I']
								except:
									in_count = 0
								try:
									del_count = CIGAR_operator_count_dict['D']
								except:
									del_count = 0
								indel_count = in_count + del_count
								try:
									mismatch_count = CIGAR_operator_count_dict['X']
								except:
									mismatch_count = 0
								if indel_count <= max_indel_count and mismatch_count <= max_mismatch_count:
									keep_reading = False
									try:
										prev_map_qual = read_map_info_dict[readID+"-"+read_seq[0:5]]#map_qual_dict[readID+"-"+read_seq[0:6]]
										if map_qual > prev_map_qual: #if higher quality mapping, keep reading to replace the entry for this read
											keep_reading = True
										else:
											keep_reading = False
									except: #if you haven't encountered the read ID-read_seq[0:6] combo, keep reading to store seq info
										keep_reading = True
									if keep_reading == True:
										read_map_info_dict[readID+"-"+read_seq[0:5]] = (map_qual,mismatch_count,indel_count,read_base_dict)

			sam_infile.close()

			if dedupe_reads == True:
				os.remove(dedupe_samfilename)

		if verbose == True:
			print(sampleID+"\t"+ref_set+"\tlen(read_map_info_dict):\t"+str(len(read_map_info_dict)))

		min_base_qual_for_counting = 20
		min_map_qual_for_counting = 25
		min_loc_for_counting = 10


		segment_isnv_link_count_dict = {}
		max_loc = 0
		sub_info_dict = {}
		allele_count_dict = {}
		count_dict = {}
		total_read_count = 0
		for readID_flag in read_map_info_dict:
			read_map_info_tup = read_map_info_dict[readID_flag]
			map_qual = read_map_info_dict[readID_flag][0]#read_map_info_dict[readID_flag][0]
			mismatch_count = read_map_info_dict[readID_flag][1]
			indel_count = read_map_info_dict[readID_flag][2]
			read_base_info_dict = read_map_info_dict[readID_flag][3]
			read_subIDs = []
			if map_qual >= min_map_qual_for_counting:
				total_read_count += 1
			for loc in read_base_info_dict:
				base = read_base_info_dict[loc][0]
				base_qual = read_base_info_dict[loc][1]
				map_qual = read_base_info_dict[loc][2]
				read_len = read_base_info_dict[loc][3]
				loc_in_read = read_base_info_dict[loc][4]
				loc_from_R_end = read_base_info_dict[loc][5]
				prop_loc = read_base_info_dict[loc][6]
				
				min_end_dist = min(loc_in_read,loc_from_R_end)
				max_end_dist = max(loc_in_read,loc_from_R_end)
				if min_end_dist >= min_loc_for_counting and base_qual >= min_base_qual_for_counting:
					loc_first = False
					try:
						count_dict[loc] += 1
					except:
						count_dict[loc] = 1
						sub_info_dict[loc] = {}
						allele_count_dict[loc] = {}
						loc_first = True
						if loc > max_loc:
							max_loc = loc
					base_first = False
					try:
						allele_count_dict[loc][base] += 1
					except:
						allele_count_dict[loc][base] = 1
						base_first = True
					if base_first == False:
						if skip_qual == True:
							sub_info_dict[loc][base][0] += base_qual
							sub_info_dict[loc][base][1] += map_qual
							sub_info_dict[loc][base][2] += read_len
							sub_info_dict[loc][base][3] += min_end_dist
							sub_info_dict[loc][base][4] += max_end_dist
							sub_info_dict[loc][base][5] += prop_loc
							sub_info_dict[loc][base][6] += mismatch_count
							sub_info_dict[loc][base][7] += indel_count
						elif skip_qual == False:
							temp_list1,temp_list2,temp_list3,temp_list4,temp_list5,temp_list6,temp_list7,temp_list8 = sub_info_dict[loc][base][0],sub_info_dict[loc][base][1],sub_info_dict[loc][base][2],sub_info_dict[loc][base][3],sub_info_dict[loc][base][4],sub_info_dict[loc][base][5],sub_info_dict[loc][base][6],sub_info_dict[loc][base][7]
							temp_list1.append(base_qual) # base_qual
							temp_list2.append(map_qual) # map_qual
							temp_list3.append(read_len) # read_len
							temp_list4.append(min_end_dist) # min_end_dist
							temp_list5.append(max_end_dist) # max_end_dist
							temp_list6.append(prop_loc) # prop_loc
							temp_list7.append(mismatch_count) # mismatch_count
							temp_list8.append(indel_count) # indel_count
							tup_out = (temp_list1,temp_list2,temp_list3,temp_list4,temp_list5,temp_list6,temp_list7,temp_list8)
							sub_info_dict[loc][base] = tup_out
					elif base_first == True:
						if skip_qual == True:
							tup_out = [base_qual,map_qual,read_len,min_end_dist,max_end_dist,prop_loc,mismatch_count,indel_count]
						elif skip_qual == False:
							tup_out = ([base_qual],[map_qual],[read_len],[loc_in_read],[loc_from_R_end],[prop_loc],[mismatch_count],[indel_count])
						sub_info_dict[loc][base] = tup_out
						if skip_linked_SNV_counting == False:
							subID = str(loc)+base
							read_subIDs.append(subID)
			read_map_info_dict[readID_flag] = ''
			if skip_linked_SNV_counting == False:
				read_subIDs = list(set(read_subIDs))
				if len(read_subIDs)>1:
					for num1 in range(0,len(read_subIDs)):
						for num2 in range(num1,len(read_subIDs)):
							if num1 != num2:
								sub1 = read_subIDs[num1]
								sub2 = read_subIDs[num2]
								loc1 = int(sub1[0:-1])
								loc2 = int(sub2[0:-1])
								if loc1 < loc2:
									pair = sub1+"_"+sub2
								elif loc1 > loc2:
									pair = sub2+"_"+sub1
								try:
									segment_isnv_link_count_dict[pair] += 1
								except:
									segment_isnv_link_count_dict[pair] = 1
		if verbose == True:
			print(sampleID+"\t"+ref_set+"\tlen(sub_info_dict):\t"+str(len(sub_info_dict)))

		del read_map_info_dict

		if skip_qual == False:
			site_loc_outfilename = temp_dir+sampleID+"."+refID+".read_loc_hist.txt"
			read_len_outfilename = temp_dir+sampleID+"."+refID+".read_len_hist.txt"
			prop_loc_outfilename = temp_dir+sampleID+"."+refID+".read_prop_loc_hist.txt"
			indel_outfilename = temp_dir+sampleID+"."+refID+".indel_count_hist.txt"
			mismatch_outfilename = temp_dir+sampleID+"."+refID+".mismatch_count_hist.txt"
			map_qual_outfilename = temp_dir+sampleID+"."+refID+".base_map_qual.txt"

			out_site_loc_outfilename = qual_dir+refID+"/"+sampleID+".read_loc_hist.txt"
			out_read_len_outfilename = qual_dir+refID+"/"+sampleID+".read_len_hist.txt"
			out_prop_loc_outfilename = qual_dir+refID+"/"+sampleID+".read_prop_loc_hist.txt"
			out_indel_outfilename = qual_dir+refID+"/"+sampleID+".indel_count_hist.txt"
			out_mismatch_outfilename = qual_dir+refID+"/"+sampleID+".mismatch_count_hist.txt"
			out_map_qual_outfilename = qual_dir+refID+"/"+sampleID+".base_map_qual.txt"
			
			site_loc_outfile = open(site_loc_outfilename,"w")
			site_loc_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			for count_val in range(0,300):
				site_loc_outfile.write("\t"+str(count_val))
			site_loc_outfile.write("\n")

			read_len_outfile = open(read_len_outfilename,"w")
			read_len_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			for count_val in range(0,300):
				read_len_outfile.write("\t"+str(count_val))
			read_len_outfile.write("\n")

			prop_loc_outfile = open(prop_loc_outfilename,"w")
			prop_loc_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			step_size = 0.01
			z = -1*step_size
			while z < 1.0:
				z += step_size
				prop_loc_outfile.write("\t"+str(round(z,3)))
			prop_loc_outfile.write("\n")

			indel_outfile = open(indel_outfilename,"w")
			indel_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			for count_val in range(0,20):
				indel_outfile.write("\t"+str(count_val))
			indel_outfile.write("\t>=20\n")

			mismatch_outfile = open(mismatch_outfilename,"w")
			mismatch_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			for count_val in range(0,40):
				mismatch_outfile.write("\t"+str(count_val))
			mismatch_outfile.write("\t>=40\n")

			map_qual_outfile = open(map_qual_outfilename,"w")
			map_qual_outfile.write("loc\tbase\tfreq\tavg_qual\tavg_map_qual\tavg_base_loc\tavg_prop_loc\tavg_mismatch\tavg_indel")
			for count_val in range(0,45):
				map_qual_outfile.write("\t"+str(count_val))
			map_qual_outfile.write("\n")

		temp_substitutions_filename = temp_dir+sampleID+"."+refID+".substitions.txt"
		out_substitutions_filename = sub_dir+refID+"/"+sampleID+"."+refID+".substitions.txt"

		substitutions_outfile = open(temp_substitutions_filename,"w")
		substitutions_outfile.write("#loc\t"+str(total_read_count))
		info_list = ['base_freq','avg_qual','map_qual','read_len','min_end_dist','max_end_dist','prop_loc','mismatches','indels']
		info_string = ''
		base_string = ''
		for category in info_list:
			for i in range(0,len(bases)):
				base = bases[i]
				info_string += "\t"+category
				base_string += "\t"+base
		substitutions_outfile.write(info_string+"\n#loc\tdepth"+base_string+"\n")

		for loc in range(0,max_loc):
			try:
				loc_count = count_dict[loc]
			except:
				loc_count = 0

			substitutions_outfile.write(str(loc)+"\t"+str(loc_count))
			base_freq_string = ''
			avg_qual_string = ''
			map_qual_string = ''
			read_len_string = ''
			base_loc_string = ''
			base_R_loc_string = ''
			prop_loc_string = ''
			mismatch_count_string = ''
			indel_count_string = ''

			poly_count = 0

			for i in range(0,len(bases)):
				base = bases[i]
				try:
					base_count = allele_count_dict[loc][base]
					# base_count = len(sub_info_dict[loc][base][0])
				except:
					base_count = 0
				if base_count > 0:
					#sub_info_dict[loc][base]([base_qual],[map_qual],[read_len],[loc_in_read],[loc_from_R_end],[prop_loc],[mismatch_count],[indel_count])
					#sub_info_dict[loc][base]([    0    ],[    1   ],[    2   ],[     3     ],[      4       ],[    5   ],[      6       ],[     7     ])
					if skip_qual == True:
						prop_base = round((base_count/loc_count),4)
						avg_qual = round(sub_info_dict[loc][base][0]/base_count,1)
						avg_map_qual = round(sub_info_dict[loc][base][1]/base_count,1)
						avg_read_len = round(sub_info_dict[loc][base][2]/base_count,1)
						avg_base_loc = round(sub_info_dict[loc][base][3]/base_count,1)
						avg_base_loc_from_R = round(sub_info_dict[loc][base][4]/base_count,1)
						avg_prop_loc = round(sub_info_dict[loc][base][5]/base_count,2)
						avg_mismatch_count = round(sub_info_dict[loc][base][6]/base_count,1)
						avg_indel_count = round(sub_info_dict[loc][base][7]/base_count,1)
					elif skip_qual == False:
						prop_base = round((base_count/loc_count),4)
						avg_qual = round(np.average(sub_info_dict[loc][base][0]),1)
						avg_map_qual = round(np.average(sub_info_dict[loc][base][1]),1)
						avg_read_len = round(np.average(sub_info_dict[loc][base][2]),1)
						avg_base_loc = round(np.average(sub_info_dict[loc][base][3]),1)
						avg_base_loc_from_R = round(np.average(sub_info_dict[loc][base][4]),1)
						avg_prop_loc = round(np.average(sub_info_dict[loc][base][5]),2)
						avg_mismatch_count = round(np.average(sub_info_dict[loc][base][6]),1)
						avg_indel_count = round(np.average(sub_info_dict[loc][base][7]),1)
				else:
					prop_base = 0
					avg_qual = 0
					avg_map_qual = 0
					avg_read_len = 0
					avg_base_loc = 0
					avg_base_loc_from_R = 0
					avg_prop_loc = 0
					avg_mismatch_count = 0
					avg_indel_count = 0

				base_freq_string +="\t"+str(prop_base)
				avg_qual_string +="\t"+str(avg_qual)
				map_qual_string +="\t"+str(avg_map_qual)
				read_len_string +="\t"+str(avg_read_len)
				base_loc_string +="\t"+str(avg_base_loc)
				base_R_loc_string +="\t"+str(avg_base_loc_from_R)
				prop_loc_string +="\t"+str(avg_prop_loc)
				mismatch_count_string +="\t"+str(avg_mismatch_count)
				indel_count_string +="\t"+str(avg_indel_count)

			substitutions_outfile.write(base_freq_string+avg_qual_string+map_qual_string+read_len_string+base_loc_string+base_R_loc_string+prop_loc_string+mismatch_count_string+indel_count_string+"\n")
		substitutions_outfile.close()
		if os.path.isfile(out_substitutions_filename) == True:
			os.remove(out_substitutions_filename)
		os.rename(temp_substitutions_filename,out_substitutions_filename)

		if verbose == True:
			print(sampleID+"\t"+ref_set+"\tFinished writing substitions file")

		if skip_qual == False:
			for loc in range(0,max_loc):
				sam_infile.close()
				try:
					loc_count = count_dict[loc]
				except:
					loc_count = 0

				for i in range(0,len(bases)):
					base = bases[i]
					try:
						prop_base = round(float(len(sub_info_dict[loc][base][0]))/float(loc_count),4)
					except:
						prop_base = 0
					if prop_base >= min_proportion:
						try:
							avg_qual = round(np.average(sub_info_dict[loc][base][0]),1)
							avg_map_qual = round(np.average(sub_info_dict[loc][base][1]),1)
							avg_read_len = round(np.average(sub_info_dict[loc][base][2]),1)
							avg_base_loc = round(np.average(sub_info_dict[loc][base][3]),1)
							avg_base_loc_from_R = round(np.average(sub_info_dict[loc][base][4]),1)
							avg_prop_loc = round(np.average(sub_info_dict[loc][base][5]),1)
							avg_mismatch_count = round(np.average(sub_info_dict[loc][base][6]),1)
							avg_indel_count = round(np.average(sub_info_dict[loc][base][7]),1)
						except:
							avg_qual = 0
							avg_map_qual = 0
							avg_read_len = 0
							avg_base_loc = 0
							avg_base_loc_from_R = 0
							avg_prop_loc = 0
							avg_mismatch_count = 0
							avg_indel_count = 0

						try:
							loc_bincount = np.bincount(sub_info_dict[loc][base][3],minlength=300)
						except:
							loc_bincount =[]
							for z in range(0,300):
								loc_bincount.append(0) 
						site_loc_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						for count_val in range(0,300):
							site_loc_outfile.write("\t"+str(loc_bincount[count_val]))
						site_loc_outfile.write("\n")
						del loc_bincount

						try:
							len_bincount = np.bincount(sub_info_dict[loc][base][2],minlength=300)
						except:
							len_bincount =[]
							for z in range(0,300):
								len_bincount.append(0) 
						read_len_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						for count_val in range(0,300):
							read_len_outfile.write("\t"+str(len_bincount[count_val]))
						read_len_outfile.write("\n")
						del len_bincount
						
						step_size = 0.01
						prop_loc_bincount = {}
						prop_loc_list = sub_info_dict[loc][base][5]
						z = -1*step_size
						while z < 1.0:
							z += step_size
							counter = 0
							for prop_val in prop_loc_list:
								if prop_val >= z and prop_val < z+step_size:
									counter += 1
							prop_loc_bincount[z] = counter
						prop_loc_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						z = -1*step_size
						while z < 1.0:
							z += step_size
							prop_loc_outfile.write("\t"+str(prop_loc_bincount[z]))
						prop_loc_outfile.write("\n")
						del prop_loc_bincount


						try:
							qual_bincount = np.bincount(sub_info_dict[loc][base][1],minlength=45)
						except:
							qual_bincount = []
							for z in range(0,45):
								qual_bincount.append(0) 
						map_qual_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						for count_val in range(0,45):
							map_qual_outfile.write("\t"+str(qual_bincount[count_val]))
						map_qual_outfile.write("\n")
						del qual_bincount


						try:
							mismatch_bincount = np.bincount(sub_info_dict[loc][base][6],minlength=200)
						except:
							mismatch_bincount =[]
							for z in range(0,200):
								mismatch_bincount.append(0) 
						mismatch_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						for count_val in range(0,40):
							mismatch_outfile.write("\t"+str(mismatch_bincount[count_val]))
						greater_than_count = 0
						for count_val in range(40,200):
							try:
								greater_than_count += mismatch_bincount[count_val]
							except:
								pass
						mismatch_outfile.write("\t"+str(greater_than_count))
						mismatch_outfile.write("\n")
						del mismatch_bincount

						try:
							indel_bincount = np.bincount(sub_info_dict[loc][base][7],minlength=200)
						except:
							indel_bincount =[]
							for z in range(0,200):
								indel_bincount.append(0) 
						indel_outfile.write(str(loc)+"\t"+base+"\t"+str(prop_base)+"\t"+str(avg_qual)+"\t"+str(avg_map_qual)+"\t"+str(avg_base_loc)+"\t"+str(avg_prop_loc)+"\t"+str(avg_mismatch_count)+"\t"+str(avg_indel_count))
						for count_val in range(0,20):
							indel_outfile.write("\t"+str(indel_bincount[count_val]))
						greater_than_count = 0
						for count_val in range(20,200):
							try:
								greater_than_count += indel_bincount[count_val]
							except:
								pass
						indel_outfile.write("\t"+str(greater_than_count))
						indel_outfile.write("\n")
						del indel_bincount
			site_loc_outfile.close()
			map_qual_outfile.close()
			indel_outfile.close()
			mismatch_outfile.close()
			prop_loc_outfile.close()
			read_len_outfile.close()
			os.rename(site_loc_outfilename,out_site_loc_outfilename)
			os.rename(read_len_outfilename,out_read_len_outfilename)
			os.rename(prop_loc_outfilename,out_prop_loc_outfilename)
			os.rename(indel_outfilename,out_indel_outfilename)
			os.rename(mismatch_outfilename,out_mismatch_outfilename)
			os.rename(map_qual_outfilename,out_map_qual_outfilename)
		del sub_info_dict

#################################################################################################################################################################################

skip_qual = True
if skip_qual == False:
	if not os.path.exists(qual_dir):
		os.makedirs(qual_dir)

process_remain_list = []
infile = open(project_dir+sampleID_file_name,"r")
for line in infile:
	line = line.strip().split("\t")
	sampleID = line[0]
	for ref_set in ref_set_list:
		subfile_found = True
		refs = ref_set.split("-")
		primary_ref = ref_set.split("-")[0].split(".f")[0]
		for ref_filename in refs:
			ref_ID = ref_filename.split(".f")[0]
			if os.path.isfile(sub_dir+ref_ID+"/"+sampleID+"."+ref_ID+'.substitutions.txt') == False:
				subfile_found = False
		if subfile_found == False:
			sample_ref_tup = (sampleID,ref_set)
			process_remain_list.append(sample_ref_tup)
infile.close()
process_remain_list = sorted(list(set(process_remain_list)))

processed_list = Parallel(n_jobs=num_cores)(delayed(main_pipeline)(sample_ref_tup,ref_set_list,input_read_dir,ref_dir,trim_dir,map_dir,qual_dir,sam_dir,temp_dir,sub_dir,command_prefix,dedupe_reads,min_proportion) for sample_ref_tup in process_remain_list)

