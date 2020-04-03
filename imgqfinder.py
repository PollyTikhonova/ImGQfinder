import sys, os
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import argparse
from multiprocessing import Pool

import pickle

import numpy as np

class QuadruplexFinder(object):

	def __init__(self, fasta_file, output_path = '', 
		GC='G', L=7, q=4, nquadruplets=4, mdef=1, tetdef=1, len_bulge=1, max_bulge = 1, 
		bulge_priority=False, repeats=False, verbose=False, nthreads=1):
	#  parse arg 
		self.fasta_file = fasta_file
		self.output_path = output_path 
		self.GC = GC
		self.L = L
		self.q = q
		self.nquadruplets = nquadruplets
		self.mdef = mdef
		self.tetdef = tetdef
		self.repeats = repeats
		self.verbose = verbose
		self.len_bulge = len_bulge
		self.max_bulge = max_bulge
		self.bulge_priority = bulge_priority
		self.nthreads = nthreads

	def load_fasta(self):
		sequences = []
		for record in SeqIO.parse(self.fasta_file, "fasta"):
			sequences.append((record.seq, record.id))
		return sequences

	def find_quadruplets_without_bulges(self, fasta):
		quadruplets = []
		stack = [self.QuadrupletDetector(nuc) for nuc in fasta[:self.q]]
		current_state = sum(stack)
		if current_state >= self.q - self.tetdef:
			quadruplets.append((0, self.q - current_state, self.q))
		for i in tqdm(range(self.q, len(fasta)), desc='Qadrupleting', disable = self.verbose):
			stack.append(self.QuadrupletDetector(fasta[i]))
			current_state = current_state + stack[-1] - stack.pop(0)
			if current_state >= self.q - self.tetdef:
				quadruplets.append((i-self.q+1, self.q - current_state, self.q))
		return quadruplets

	def QuadrupletDetector(self, quadr):
		if self.repeats:
			quadr = quadr.upper()
		return 1 if quadr == self.GC.upper() else 0

	def find_quadruplets_wrapper(self, data):
		return self.find_quadruplets(**data)

	def find_quadruplets(self, fasta, shift=0, tqdm_keep_silence=None):
		'''
		bulge_stack - a set of numbers - amounts how many non-G nucleotides was before + 1
		'''
		tqdm_keep_silence = self.verbose if tqdm_keep_silence is None else tqdm_keep_silence
		quadruplets = []
		quadruplets_sequences = []
		open_bulge = 0
		bulge_stack = []
		sequence_stack = ''
		bulge_current_state = 0
		bulge_current_num_state = 0
		bulge_num_state = 0
		n_bulges = 0

		def add_bulge(nuc):
			nonlocal open_bulge, bulge_current_num_state, bulge_current_state, bulge_num_state, bulge_stack, n_bulges 
			if self.QuadrupletDetector(nuc):
					bulge_stack.append(open_bulge+1)
					if len(bulge_stack) == 1:
						bulge_stack[0] = 1
					open_bulge = 0
					if bulge_current_num_state < self.q:
						bulge_current_num_state += 1
						bulge_current_state += bulge_stack[-1]
						if bulge_stack[-1] != 1:
							n_bulges += 1
					else:
						bulge_num_state += 1
			else:
				open_bulge += 1


		def remove_bulge(nuc):
			nonlocal bulge_num_state, bulge_current_state, bulge_current_num_state, bulge_stack, n_bulges
			if self.QuadrupletDetector(nuc):
				if bulge_num_state > 0:
					bulge_current_state += bulge_stack[bulge_current_num_state]
					bulge_num_state -= 1
					if bulge_stack[bulge_current_num_state] != 1:
						n_bulges += 1
				else:
					bulge_current_num_state -= 1
				bulge_current_state -= bulge_stack.pop(0) 
				if len(bulge_stack) > 0:
					pop = bulge_stack.pop(0)
					if pop != 1:
						n_bulges -= 1
					bulge_current_state -= pop - 1
					bulge_stack.insert(0, 1)

		for i, nuc in enumerate(fasta[:(self.q+self.len_bulge)]):
			add_bulge(nuc)
			sequence_stack = sequence_stack+nuc
			
		if ((bulge_current_num_state == self.q) & (n_bulges <= self.max_bulge) &
			(self.QuadrupletDetector(fasta[0])) & (self.QuadrupletDetector(fasta[bulge_current_state-1]))):
				quadruplets.append((0+shift, n_bulges, bulge_current_state))
				quadruplets_sequences.append(sequence_stack[:bulge_current_state])

		stack = [self.QuadrupletDetector(nuc) for nuc in fasta[:self.q]]
		current_state = sum(stack)
		if ((current_state >= self.q - self.tetdef) & (current_state < self.q) &
			(self.QuadrupletDetector(fasta[0])) & (self.QuadrupletDetector(fasta[self.q-1]))):
			quadruplets.append((0+shift, self.q - current_state, self.q))
			quadruplets_sequences.append(sequence_stack[:self.q])
		for i in tqdm(range(self.q, len(fasta)), desc='Quadrupleting', disable = tqdm_keep_silence):
			remove_bulge(fasta[i-self.q])

			i_bulge = i + self.len_bulge
			if i_bulge < len(fasta):
				add_bulge(fasta[i_bulge])
				sequence_stack = sequence_stack+fasta[i_bulge]

			stack.append(self.QuadrupletDetector(fasta[i]))
			current_state = current_state + stack[-1] - stack.pop(0)
			sequence_stack = sequence_stack[1:] 

			if self.QuadrupletDetector(fasta[i-self.q+1]):
				if ((bulge_current_num_state == self.q) & (n_bulges <= self.max_bulge) & 
					(self.QuadrupletDetector(fasta[i-self.q+bulge_current_state]))):
					quadruplets.append((i-self.q+1+shift, n_bulges, bulge_current_state))
					quadruplets_sequences.append(sequence_stack[:bulge_current_state])
				if ((current_state >= self.q - self.tetdef) & (current_state < self.q) & 
					(self.QuadrupletDetector(fasta[i]))):
					quadruplets.append((i-self.q+1+shift, self.q - current_state, self.q))
					quadruplets_sequences.append(sequence_stack[:self.q])
		return quadruplets, quadruplets_sequences

	def find_quadruplets_in_parallel(self, fasta):
		pool = Pool(processes=self.nthreads)
		minimal_chunk_length = self.q + self.len_bulge
		base_chunk_length = len(fasta) // self.nthreads
		if base_chunk_length < minimal_chunk_length:
			base_chunk_length = minimal_chunk_length
		fasta_chunks_starts = list(range(0, len(fasta), base_chunk_length))
		if len(fasta) % base_chunk_length != 0:
			fasta_chunks_starts = fasta_chunks_starts[:-1]		
		fasta_chunks_ends = fasta_chunks_starts[1:] + [len(fasta)-minimal_chunk_length]
		quadruplets_list = pool.map(self.find_quadruplets_wrapper, ({'fasta':fasta[start:(end+minimal_chunk_length)],
																	 'shift':start, 
																	 'tqdm_keep_silence':None if silence_ind==len(fasta_chunks_starts)-1 else True}
							   for silence_ind, (start, end) in enumerate(zip(fasta_chunks_starts, fasta_chunks_ends))))
		pool.close()
		pool.join()
		quadruplets = []
		quadruplets_sequences = []
		quadruplets_list_ = []
		quadruplets_seq_list_ = []
		for quad, quad_seq in quadruplets_list:
			if len(quad) != 0:
				quadruplets_list_.append(quad)
				quadruplets_seq_list_.append(quad_seq)
		del quadruplets_list

		for quadruplet_now, quadruplet_next, quadruplet_seq_now, quadruplet_seq_next in zip(
												   quadruplets_list_[:-1], quadruplets_list_[1:], 
												   quadruplets_seq_list_[:-1], quadruplets_seq_list_[1:]):
			first_next_quad = quadruplet_next[0]
			num_quad_now = -1
			while (first_next_quad == quadruplet_now[num_quad_now]) or (first_next_quad[0] <= quadruplet_now[num_quad_now][0]):
				num_quad_now -= 1
			num_quad_now += 1
			if num_quad_now != 0:
				quadruplet_now = quadruplet_now[:num_quad_now]
				quadruplet_seq_now = quadruplet_seq_now[:num_quad_now]
			quadruplets.extend(quadruplet_now)
			quadruplets_sequences.extend(quadruplet_seq_now)
		quadruplets_sequences.extend(quadruplet_seq_next)
		quadruplets.extend(quadruplet_next)
		del quadruplets_list_
		del quadruplets_seq_list_
		return quadruplets, quadruplets_sequences

	def find_quadruplexes_wrapper(self, data):
		return self.find_quadruplexes(**data)

	def find_quadruplexes(self, quadruplets, tqdm_keep_silence=None):
		'''
			quadruplex: [[Q1-Start,	Q1-Defects,	Q1-Length]]*self.nquadruplets
		'''
		tqdm_keep_silence = self.verbose if tqdm_keep_silence is None else tqdm_keep_silence
		total_wrongs = 0  #number of the quadruplets with defect
		wrongNum = 0

		def check_conditions():
			nonlocal total_wrongs, wrongNum
			if i == 0:
				total_wrongs = 0
				wrongNum = 0
			elif (quadruplets[k][0] - quadruplets[quadruplex_set[i-1]][0] <= quadruplets[quadruplex_set[i-1]][2]):
				return 'too close'
			elif (quadruplets[k][0] - (quadruplets[quadruplex_set[i-1]][0] + quadruplets[quadruplex_set[i-1]][2]) > self.L):
				return 'too far'
			if quadruplets[k][1] != 0:
				wrongNum = i+1
				total_wrongs += 1
			if total_wrongs > self.mdef:
				total_wrongs -= 1
				return False
			else:
				return True

		def revert_wrongs():
			nonlocal total_wrongs, wrongNum
			if (i >= 0) & (quadruplets[quadruplex_set[i]][1] != 0):
				total_wrongs -= 1
				if wrongNum == i+1:
					for j in range(i):
						if quadruplets[quadruplex_set[j]][1] != 0:
							wrongNum == j+1
							break
				if wrongNum == i+1:
					wrongNum = 0

		quadruplexes = []
		quadruplex_set = list(range(-1, self.nquadruplets))

		i = 0
		k = quadruplex_set[i]

		with tqdm(desc='Qadruplexing', total=len(quadruplets), disable = tqdm_keep_silence) as pbar:
			while i >= 0:
				k = quadruplex_set[i]+1	
				if i == 0:
					pbar.update(1)
				if i == self.nquadruplets:
					quadruplex = tuple([quadruplets[qu] for qu in quadruplex_set[:-1]] + [total_wrongs])
					quadruplexes.append(list(quadruplex))
					i -= 1
					revert_wrongs()
				elif k >= len(quadruplets) - self.nquadruplets + 1 + i:
					i -= 1
					revert_wrongs()
				else:
					status = check_conditions()
					if status == True:
						quadruplex_set[i] = k
						i += 1
						quadruplex_set[i] = quadruplex_set[i-1]	
					elif status == 'too far':
						i -= 1
						revert_wrongs()
					else:
						quadruplex_set[i] = k
			pbar.update(len(quadruplets) - pbar.n)
		return quadruplexes

	def group_quadruplexes(self, quadruplexes):
		groups = []
		q1 = 0
		q2 = 1
		with tqdm(desc='Grouping', 	total=len(quadruplexes)-1, disable = self.verbose) as pbar:
			while q1 < len(quadruplexes)-1:
				while q2 < len(quadruplexes):
					pbar.update(1)
					tetrads_lehgth_q1 = sum([quadruplexes[q1][i][2]+quadruplexes[q1][i][1] for i in range(self.nquadruplets)])
					tetrads_lehgth_q2 = sum([quadruplexes[q2][i][2]+quadruplexes[q2][i][1] for i in range(self.nquadruplets)])
					general_length_q1 = quadruplexes[q1][self.nquadruplets - 1][0] + quadruplexes[q1][self.nquadruplets - 1][2] - 1 - quadruplexes[q1][0][0]
					general_length_q2 = quadruplexes[q2][self.nquadruplets - 1][0] + quadruplexes[q2][self.nquadruplets - 1][2] - 1 - quadruplexes[q2][0][0]
					if (quadruplexes[q2][0][0] > quadruplexes[q1][self.nquadruplets - 1][0] + quadruplexes[q1][self.nquadruplets - 1][2] - 1):
						groups.append(quadruplexes[q1])
						q1 = q2
						if (q2 == len(quadruplexes)-1):
							groups.append(quadruplexes[q2])
							q1 = len(quadruplexes)
					elif ((tetrads_lehgth_q2 < tetrads_lehgth_q1) & (not self.bulge_priority) or
						  (tetrads_lehgth_q2 >= tetrads_lehgth_q1) & (self.bulge_priority) or
						  (general_length_q2 < general_length_q1) & (not self.bulge_priority) or
						  (general_length_q2 < general_length_q1) & (self.bulge_priority)):
							q1 = q2
							if (q2 == len(quadruplexes)-1):
								groups.append(quadruplexes[q2])
								q1 = len(quadruplexes)
					elif (q2 == len(quadruplexes)-1):
								groups.append(quadruplexes[q1])
								q1 = len(quadruplexes)
					q2 += 1
		return groups

	def find_quadruplexes_in_parallel(self, quadruplets):
		pool = Pool(processes=self.nthreads)		
		minimal_chunk_length = (self.q + self.len_bulge + self.L)*(self.nquadruplets)-self.L
		if len(quadruplets) > self.nthreads:
			base_chunk_length = len(quadruplets) // self.nthreads
		else:
			base_chunk_length = 1
		quadruplets_chunks_starts = list(range(0, len(quadruplets), base_chunk_length))
		if len(quadruplets) % base_chunk_length != 0:
			quadruplets_chunks_starts = quadruplets_chunks_starts[:-1]
		quadruplets_chunks_ends = []
		for start_tmp in quadruplets_chunks_starts[1:]:
			end_ind = start_tmp
			end_val = quadruplets[start_tmp][0]
			tmp_end_val = quadruplets[end_ind][0]
			while (end_ind < len(quadruplets)) and (tmp_end_val - end_val <= minimal_chunk_length):
				end_ind += 1
				tmp_end_val = quadruplets[end_ind][0]
			quadruplets_chunks_ends.append(end_ind-1)
		quadruplets_chunks_ends.append(len(quadruplets))

		quadruplexes_list = pool.map(self.find_quadruplexes_wrapper, ({'quadruplets':quadruplets[start:end],
					 'tqdm_keep_silence':None if silence_ind==len(quadruplets_chunks_starts)-1 else True}
							   for silence_ind, (start, end) in enumerate(zip(quadruplets_chunks_starts, quadruplets_chunks_ends))))
		pool.close()
		pool.join()
		quadruplexes_list_ = []
		for quad in quadruplexes_list:
			if len(quad) != 0:
				quadruplexes_list_.append(quad)
		del quadruplexes_list
		quadruplexes = []
		for quadruplex_now, quadruplex_next in zip(quadruplexes_list_[:-1], quadruplexes_list_[1:]):
			first_next_quad = quadruplex_next[0]
			num_quad_now = -1
			while first_next_quad[0][0] <= quadruplex_now[num_quad_now][0][0]:
				if (first_next_quad == quadruplex_now[num_quad_now]) or (first_next_quad[0][0] <= quadruplex_now[num_quad_now][0][0]):
					num_quad_now -= 1
			num_quad_now += 1
			if num_quad_now != 0:
				quadruplex_now = quadruplex_now[:num_quad_now]
			quadruplexes.extend(quadruplex_now)
		quadruplexes.extend(quadruplex_next)
		del quadruplexes_list_
		return quadruplexes

	def group_to_ranges(self, groups, fasta_id):
		ranges = []
		for group in tqdm(groups, desc='Converting to ranges', disable = self.verbose):
			start = group[0][0]
			end = group[self.nquadruplets-1][0]+group[self.nquadruplets-1][2]-1
			ranges.append((fasta_id, start, end))
		return ranges

	def prepare_quadruplets_toprint(self, quadruplets, quadruplets_sequences, tqdm_keep_silence=None):
		tqdm_keep_silence =  self.verbose if tqdm_keep_silence is None else tqdm_keep_silence
		quadruplets_toprint = []
		for quad, seq in tqdm(list(zip(quadruplets, quadruplets_sequences)), 
							  desc='Postprocessing quadruplets', disable=tqdm_keep_silence):
			quad = list(quad)
			quad.append(seq)
			quadruplets_toprint.append(quad)
		return quadruplets_toprint

	def prepare_quadruplexes_toprint(self, quadruplexes, fasta_di, tqdm_keep_silence=None):
		tqdm_keep_silence = self.verbose if tqdm_keep_silence is None else tqdm_keep_silence
		quadruplexes_toprint = []
		[(shift, fasta)] = fasta_di.items()
		for quadruplex in tqdm(quadruplexes, desc='Postprocessing quadruplexes', disable=tqdm_keep_silence):
			seq = ''
			quadruplex_toprint = []
			for qu1, qu2 in zip(quadruplex[:-2], quadruplex[1:-1]):
				seq = seq + fasta[(qu1[0]-shift):(qu1[0]+qu1[2]-shift)].upper()+\
							fasta[(qu1[0]+qu1[2]-shift):(qu2[0]-shift)].lower()
				quadruplex_toprint.extend(list(qu1))
			quadruplex_toprint.extend(list(qu2))
			quadruplex_toprint.append(quadruplex[-1])
			seq = seq+fasta[(qu2[0]-shift):(qu2[0]+qu2[2]-shift)].upper()
			quadruplex_toprint.append(seq)
			quadruplexes_toprint.append(tuple(quadruplex_toprint))
		return quadruplexes_toprint

	def prepare_groups_toprint(self, groups, fasta, tqdm_keep_silence=None):
		tqdm_keep_silence = self.verbose if tqdm_keep_silence is None else tqdm_keep_silence
		groups_toprint = []
		for group in tqdm(groups, desc='Postprocessing groups', disable=tqdm_keep_silence):
			seq = ''
			group_toprint = []
			for qu1, qu2 in zip(group[:-2], group[1:-1]):
				seq = seq + fasta[qu1[0]:(qu1[0]+qu1[2])].upper()+fasta[(qu1[0]+qu1[2]):qu2[0]].lower()
				group_toprint.extend(qu1)
			group_toprint.extend(qu2)
			group_toprint.append(group[-1])
			seq = seq+fasta[qu2[0]:(qu2[0]+qu2[2])].upper()
			group_toprint.append(seq)
			groups_toprint.append(tuple(group_toprint))
		return groups_toprint

	def split_args_for_prepare_quadruplexes_toprint(self, quadruplexes, fasta, n):
		quad_len = len(quadruplexes) // n   
		minimal_chunk_length = (self.q + self.len_bulge + self.L)*(self.nquadruplets)
		parts = list(range(0, len(quadruplexes), quad_len))[1:]
		if len(quadruplexes) % n != 0:
			parts = parts[:-1]
		quadruplexes_parts = [quadruplexes[start:end] for start, end in zip(
														  [0]+parts, parts+[len(quadruplexes)])]
		fasta_parts_coordinates = [(quadruplex_set[0][0][0], quadruplex_set[-1][-2][0]+minimal_chunk_length) 
								   for quadruplex_set in quadruplexes_parts]
		fasta_parts = [{start:fasta[start:end]} for start, end in fasta_parts_coordinates]
		show_status = [True]*(len(quadruplexes_parts)-1)+[None]
		return list(zip(quadruplexes_parts, fasta_parts, show_status))

	def postprocess_wrapper(self, kwargs):
		'''
		args: {'args': args, 'func':function}
		'''
		return kwargs['func'](*kwargs['args'])

	def save_tables(self, df, columns, fasta_id, which_table):
		n = len(columns)
		with open('{}/{}_{}'.format(self.output_path, fasta_id, which_table), 'w') as f:
			if n == 4:
				f.write('{}\t{}\t{}\t{}\n'.format(columns[0], columns[1], columns[2], columns[3]))
				for row in df:
					f.write('{}\t{}\t{}\t{}\n'.format(row[0], row[1], row[2], row[3]))
			elif n == 14:
				f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
												  columns[0], columns[1], columns[2], columns[3], 
												  columns[4], columns[5], columns[6], columns[7], 
												  columns[8], columns[9], columns[10], columns[11], 
												  columns[12], columns[13]))
				for row in df:
					f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
												  row[0], row[1], row[2], row[3], 
												  row[4], row[5], row[6], row[7], 
												  row[8], row[9], row[10], row[11], 
												  row[12], row[13]))
			elif n==3:
				for row in df:
					f.write('{}\t{}\t{}\n'.format(row[0], row[1], row[2]))

	def save_tables_wrapper(self, args):
		return self.save_tables(*args)
	
	def description_file(self):
		all_members = self.__dict__.keys()
		columns_description =  '''\n\nColumns Description

	Qudruplets File:
		Start: an index of a quadruplet begining
		Number of Defects: a total number of the mismatches or a number of bulges
		Length: a length of a quadruplet
		Sequence: quadruplet sequences if not suppressed

	Quadruplex & Group Files:
		Qi-Start: an index of a quadruplet begining
		Qi-Defects: a total number of the mismatches or a number of bulges
		Qi-Length: a length of a quadruplet
		Defective: a number of quadruplets with defects (mismatches or bulges)
		Sequence: a sequence of a quadruplex with loops if not suppressed, quadruplets are uppercase
		
	Ranges File: bed-format
		Fasta ID: fasta id
		Start: an index of a quadruplex begining
		End: an index of the end of the quadruplex
		'''
		description_file = 'Parametres\n'+'\n'.join(['\t%s = %s'%(item, self.__dict__[item]) for item in all_members if (not item.startswith("_")) & ('nquadruplets' not in item)]) + columns_description
		with open('%s/description.txt'%(self.output_path), 'w') as f:
			f.write(description_file)
	
	def run(self, print_sequences=True, output_type = [4]):
		print('Loading %s'%self.fasta_file)
		sequences = self.load_fasta()
		print('This fasta file contains %d sequences.'%len(sequences))		
		for fasta, fasta_id in sequences:
			quadruplets, quadruplets_sequences, quadruplexes, groups, ranges, quadruplets_toprint, quadruplexes_toprint, groups_toprint = [[]]*8
			print('Processing %s:'%fasta_id)
			if (self.nthreads == 1) or (self.len_bulge == 0) or (self.max_bulge == 0):
				quadruplets, quadruplets_sequences = self.find_quadruplets(fasta) if ((self.len_bulge > 0) or (self.max_bulge != 0)) else self.find_quadruplets_without_bulges(fasta)
			else:
				quadruplets, quadruplets_sequences = self.find_quadruplets_in_parallel(fasta)
			if output_type[-1] > 0:
				if (self.nthreads == 1):
					quadruplexes = self.find_quadruplexes(quadruplets)
				else:
					quadruplexes = self.find_quadruplexes_in_parallel(quadruplets)
			if output_type[-1] > 1:
				groups = self.group_quadruplexes(quadruplexes)
			if output_type[-1] > 2:
				ranges = self.group_to_ranges(groups, fasta_id)

			columns_set1 = ['Start', 'Number of Defects', 'Length']
			columns_set2 = []
			[columns_set2.extend(['Q%d-Start'%i, 'Q%d-Defects'%i, 'Q%d-Length'%i]) for i in range(1, self.nquadruplets+1)]
			columns_set2.extend(['Defective'])

			if output_type[0] < 3:
				k = sum([0 if ind_num in output_type else 1 for ind_num in [0, 2]])
				if print_sequences:
					if self.nthreads > 1: 
						pool = Pool(processes=self.nthreads)
						args_quadruplexes = []
						n = 1
						if 1 in output_type:
							if self.nthreads - k > 2:
								n = self.nthreads
								args_quadruplexes = self.split_args_for_prepare_quadruplexes_toprint(quadruplexes, fasta, n)
							elif self.nthreads - k > 1:
								n = self.nthreads - k 
								args_quadruplexes = self.split_args_for_prepare_quadruplexes_toprint(quadruplexes, fasta, n)
							else:
								n = 1
								args_quadruplexes = {0:fasta}
						args_dict = {0: [(quadruplets, quadruplets_sequences)], 
									 1: args_quadruplexes,
									 2: [(groups, fasta)]}
						func_dict = {0: [self.prepare_quadruplets_toprint], 
									 1: [self.prepare_quadruplexes_toprint]*n, 
									 2: [self.prepare_groups_toprint]}
						results_inds_dict = {0: [0],
											 1: [1]*n,
											 2: [2]
						}
						args_all = []
						functions = []
						results_inds = []
						for output_ind in output_type:
							if output_ind < 3:
								functions.extend(func_dict[output_ind])
								args_all.extend(args_dict[output_ind])
								results_inds.extend(results_inds_dict[output_ind])
						uni, inds, counts = np.unique(results_inds, return_index=True, return_counts=True)
						slice_dict = {}
						for un, ind, count in zip(uni, inds, counts):
						    slice_dict[un] = (ind,(ind+count))

						results = pool.map(self.postprocess_wrapper, ({'func':func, 'args': args}
													   for func, args in zip(functions, args_all)))
						if 0 in slice_dict.keys():
							[quadruplets_toprint] = results[slice_dict[0][0]:slice_dict[0][1]]
						if 1 in slice_dict.keys():
							quadruplexes_toprint_all = results[slice_dict[1][0]:slice_dict[1][1]]
							quadruplexes_toprint = []
							[quadruplexes_toprint.extend(quad) for quad in quadruplexes_toprint_all];
						if 2 in slice_dict.keys():
							[groups_toprint] = results[slice_dict[2][0]:slice_dict[2][1]]

						pool.close()
						pool.join()
					else:
						if 0 in output_type:
							quadruplets_toprint = self.prepare_quadruplets_toprint(quadruplets, quadruplets_sequences)
						if 1 in output_type:
							quadruplexes_toprint = self.prepare_quadruplexes_toprint(quadruplexes, {0:fasta})
						if 2 in output_type:
							groups_toprint = self.prepare_groups_toprint(groups, fasta)

					columns_set1.extend(['Sequence'])
					columns_set2.extend(['Sequence'])
				else:
					if 0 in output_type:
						quadruplets_toprint = quadruplets
					if 1 in output_type:
						quadruplexes_toprint = []
						for quadruplex in quadruplexes:
							seq = ''
							quadruplex_toprint = []
							for qu1 in quadruplex[:-1]:
								quadruplex_toprint.extend(qu1)
							quadruplex_toprint.append(quadruplex[-1])
							quadruplexes_toprint.append(tuple(quadruplex_toprint))
					if 2 in output_type:
						groups_toprint = []
						for group in groups:
							seq = ''
							group_toprint = []
							for qu1 in group[:-1]:
								group_toprint.extend(qu1)
							group_toprint.append(group[-1])
							groups_toprint.append(tuple(group_toprint))
					
			for i in tqdm(range(1), desc='Saving tables', disable=self.verbose):
				pool = Pool(processes=self.nthreads)
				data = np.array([(quadruplets_toprint, columns_set1, fasta_id, 'quadruplets.csv'), 
						(quadruplexes_toprint, columns_set2, fasta_id, 'quadruplexes.csv'), 
						(groups_toprint, columns_set2, fasta_id, 'groups.csv'), 
						(ranges, ['Fasta ID', 'Start', 'End'], fasta_id, 'ranges.bed')])[output_type]
				pool.map(self.save_tables_wrapper, data)
				pool.close()
				pool.join()

		self.description_file()
		print('Finished')
# Disable prints
def blockPrint():
	sys.stdout = open(os.devnull, 'w')

def main():
	parser = argparse.ArgumentParser(prog='ImGQFinder', description='The tool for finding G-, C- quadruplexes. The output positions are represented in a zero based counting.')
	parser.add_argument('-i', '--input', help='Assembly scaffolds/contigs or full genomes, required.', required=True)
	parser.add_argument('-o', '--output', default='', help='Name/path of a folder for output files. Saves to the current folder if not provided.')
	parser.add_argument('-GC', default='G', help='Quad type, G- or C-. By default, G.')
	parser.add_argument('-L', default=7, help='Maximum loop length. By default, 7.')
	parser.add_argument('-q', default=4, help="The length of a quadruplet.") # the length of a tetrad
	parser.add_argument('-nq', '--nquadruplets', default=4, help=argparse.SUPPRESS) # 'Number of quadruplets. By default, 4.'
	parser.add_argument('-mdef', default=1, help='Allowed number of defective tetrads. By default, 1.')
	parser.add_argument('-bulgelen', default=1, help='Total length of bulges in one quadruplet. By default, 1.')
	parser.add_argument('-maxbulge', default=1, help='Maximum number of bulges per quadruplet. By default, 1.')
	parser.add_argument('-bp', '--bulge_priority', action='store_true', help='By default, quadrouplexes with shorter bulge or without them are preferable while grouping. This behaviour can be changed with this parameter.')
	parser.add_argument('-tetdef', default=1, help='Allowed number of defective nucleotides in tetrads. By default, 1.')
	parser.add_argument('-ns', '--no-sequences', action='store_true', help='Not to include sequences to the output.')
	parser.add_argument('-r', '--repeats', action='store_true', help='To include soft-masked genome areas. By default, not included.')
	parser.add_argument('-v', '--verbose', action='store_true', help='Show the status of procesing or not. By default print stages info.')
	parser.add_argument('--nthreads', default=1, help='Number of kernels to use.')
	parser.add_argument('--output_type', default=['all'], nargs='+', help='List the numbers of file types you need the tool to generate or write all if you want all files. All - is the default. 0 - quadruplets, 1 - quadruplexes, 2 - groups, 3 - ranges. For example, --output_type 1 2 will generate only 2 files: quadruplexes and groups.')

	args = parser.parse_args()
	if not os.path.isdir(args.output):
		os.mkdir(args.output)
	args.output = os.path.dirname(args.output)
	if args.verbose:
		blockPrint()

	args.output_type = [atype.lower() for atype in args.output_type]
	output_type_dict = {'all':4, '0':0, '1':1, '2':2, '3':3, '4':4}
	output_type_dict_report = {'all':'all', '0':'quadruplets', '1':'quadruplexes', '2':'groups', '3':'ranges', '4':'all'}
	output_type = sorted([output_type_dict[user_type] for user_type in args.output_type if user_type in list(output_type_dict.keys())])
	output_type_report = [output_type_dict_report[user_type] for user_type in args.output_type if user_type in list(output_type_dict.keys())]
	if output_type[-1] == 4:
		output_type = [0, 1, 2, 3]
		output_type_report = [output_type_dict_report[user_type] for user_type in ['0', '1', '2', '3']]
	if 'all' in output_type:
		output_type = ['quadruplets', 'quadruplexes', 'groups', 'ranges']
	if len(output_type) == 1:
		print('The ImGQfinder will generate %s file.'%(output_type_report[0]))
	else:
		print('The ImGQfinder will generate %s and %s files.'%(', '.join(output_type_report[:-1]), output_type_report[-1]))
	
	if int(args.mdef) < int(args.tetdef):
		print('Warning:  The allowed number of defective nucleotides (-tetdef) is more than the number of nucleotides (-mdef).', end='\n\n')
	finder = QuadruplexFinder(args.input, output_path = args.output, verbose = args.verbose, repeats=args.repeats, 
		GC=args.GC, L=int(args.L) , q=int(args.q), nquadruplets=int(args.nquadruplets), mdef=int(args.mdef), tetdef=int(args.tetdef),
		len_bulge = int(args.bulgelen), max_bulge = int(args.maxbulge), bulge_priority = args.bulge_priority, nthreads = int(args.nthreads))
	finder.run(print_sequences= not args.no_sequences, output_type = output_type)

if __name__ == '__main__':
	main()
