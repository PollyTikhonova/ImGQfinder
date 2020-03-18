import sys, os
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import argparse

import numpy as np

class QuadruplexFinder(object):

	def __init__(self, fasta_file, output_path = '', 
		GC='G', L=7, q=4, nquadruplets=4, mdef=1, tetdef=1, len_bulge=1, max_bulge = 1, 
		bulge_priority=False, repeats=False, verbose=False):
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

	def find_quadruplets(self, fasta):
		quadruplets = []
		open_bulge = 0
		bulge_stack = []
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
					# print('+1', n_bulges, bulge_stack)
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
				# print('-1', n_bulges, bulge_stack)

		for i, nuc in enumerate(fasta[:(self.q+self.len_bulge)]):
			add_bulge(nuc)

		if (bulge_current_num_state == self.q) & (n_bulges <= self.max_bulge):
				quadruplets.append((0, n_bulges, bulge_current_state))
				# quadruplets.append((0, 1 if bulge_current_state-self.q > 0 else 0, bulge_current_state))
				 
		stack = [self.QuadrupletDetector(nuc) for nuc in fasta[:self.q]]
		current_state = sum(stack)
		if current_state >= self.q - self.tetdef:
			quadruplets.append((0, self.q - current_state))
		for i in tqdm(range(self.q, len(fasta)), desc='Qadrupleting', disable = self.verbose):
			remove_bulge(fasta[i-self.q])

			i_bulge = i + self.len_bulge
			if i_bulge < len(fasta):
				add_bulge(fasta[i_bulge])

			# if sum(np.array(bulge_stack) > 1) != n_bulges:
			# 	raise Exception('o-o')

			stack.append(self.QuadrupletDetector(fasta[i]))
			current_state = current_state + stack[-1] - stack.pop(0)

			if self.QuadrupletDetector(fasta[i-self.q+1]):
				# print(bulge_current_state)
				if (bulge_current_num_state == self.q) & (n_bulges <= self.max_bulge):
					# quadruplets.append((i-self.q+1, 1 if bulge_current_state-self.q > 0 else 0, bulge_current_state))
					quadruplets.append((i-self.q+1, n_bulges, bulge_current_state))
				if (current_state >= self.q - self.tetdef) & (current_state < self.q) & (self.QuadrupletDetector(fasta[i])):
					quadruplets.append((i-self.q+1, self.q - current_state, self.q))

		return quadruplets


	def find_quadruplexes(self, quadruplets):
		'''
			quadruplex: [[Q1-Start,	Q1-Defects,	Q1-Length]]*self.nquadruplets
		'''
		total_wrongs = 0  #number of the quadruplet with defect
		wrongNum = 0

		def check_conditions():
			nonlocal total_wrongs, wrongNum
			if i == 0:
				total_wrongs = 0
				wrongNum = 0
			elif (quadruplets[k][0] - quadruplets[quadruplex_set[i-1]][0] <= quadruplets[quadruplex_set[i-1]][2]):
				return 'too close'
			# elif (quadruplets[k][0] - quadruplets[quadruplex_set[i-1]][0] <= self.q):
			# 	return 'too close'
			# elif (quadruplets[k][0] - (quadruplets[quadruplex_set[i-1]][0]+ self.q) > self.L):
			# 	return 'too far'
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
		
		with tqdm(desc='Qadruplexing', total=len(quadruplets), disable = self.verbose) as pbar:
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

	def group_to_ranges(self, groups):
		ranges = []
		for group in groups:
			start = group[0][0]
			end = group[self.nquadruplets-1][0]+group[self.nquadruplets-1][2]-1
			ranges.append((start, end))
			return ranges

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
		'''
		description_file = 'Parametres\n'+'\n'.join(['\t%s = %s'%(item, self.__dict__[item]) for item in all_members if (not item.startswith("_")) & ('nquadruplets' not in item)]) + columns_description
		with open('%s/description.txt'%(self.output_path), 'w') as f:
			f.write(description_file)


	def run(self, print_sequences=True):
		print('Loading %s'%self.fasta_file)
		sequences = self.load_fasta()
		print('This fasta file contains %d sequences.'%len(sequences))
		for fasta, fasta_id in sequences:
			print('Processing %s:'%fasta_id)
			quadruplets = self.find_quadruplets(fasta) if ((self.len_bulge > 0) or (self.max_bulge == 0)) else self.find_quadruplets_without_bulges(fasta)
			quadruplexes = self.find_quadruplexes(quadruplets)
			groups = self.group_quadruplexes(quadruplexes)
			ranges = self.group_to_ranges(groups)

			columns_set1 = ['Start', 'Number of Defects', 'Length']
			columns_set2 = []
			[columns_set2.extend(['Q%d-Start'%i, 'Q%d-Defects'%i, 'Q%d-Length'%i]) for i in range(1, self.q+1)]
			# [columns_set2.extend(['Q%d Start'%i]) for i in range(1, self.q+1)]
			columns_set2.extend(['Defective'])

			if print_sequences:
				quadruplets_toprint = []
				for quadruplet in quadruplets:
					quadruplet = list(quadruplet)
					quadruplet.append(fasta[quadruplet[0]:quadruplet[0]+quadruplet[2]])
					quadruplets_toprint.append(tuple(quadruplet))
				# quadruplets = quadruplets_new
				quadruplexes_toprint = []
				for quadruplex in quadruplexes:
					seq = ''
					quadruplex_toprint = []
					for qu1, qu2 in zip(quadruplex[:-2], quadruplex[1:-1]):
						seq = seq + fasta[qu1[0]:(qu1[0]+qu1[2])].upper()+fasta[(qu1[0]+qu1[2]):qu2[0]].lower()
						quadruplex_toprint.extend(list(qu1))
					quadruplex_toprint.extend(list(qu2))
					quadruplex_toprint.append(quadruplex[-1])
					seq = seq+fasta[qu2[0]:(qu2[0]+qu2[2])].upper()
					quadruplex_toprint.append(seq)
					quadruplexes_toprint.append(tuple(quadruplex_toprint))
				groups_toprint = []
				for group in groups:
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

				columns_set1.extend(['Sequence'])
				columns_set2.extend(['Sequence'])
			else:
				quadruplets_toprint = quadruplets
				quadruplexes_toprint = []
				for quadruplex in quadruplexes:
					seq = ''
					quadruplex_toprint = []
					for qu1 in quadruplex[:-1]:
						quadruplex_toprint.extend(qu1)
					quadruplex_toprint.append(quadruplex[-1])
					quadruplexes_toprint.append(tuple(quadruplex_toprint))
				groups_toprint = []
				for group in groups:
					seq = ''
					group_toprint = []
					for qu1 in group[:-1]:
						group_toprint.extend(qu1)
					group_toprint.append(group[-1])
					groups_toprint.append(tuple(group_toprint))
			# print(columns_set1)
			pd.DataFrame(quadruplets_toprint, columns=columns_set1).to_csv('%s/%s_quadruplets.csv'%(self.output_path, fasta_id), index=False, sep='\t')
			pd.DataFrame(quadruplexes_toprint, columns=columns_set2).to_csv('%s/%s_quadruplexes.csv'%(self.output_path, fasta_id), index=False, sep='\t')
			pd.DataFrame(groups_toprint, columns=columns_set2).to_csv('%s/%s_groups.csv'%(self.output_path, fasta_id), index=False, sep='\t')
			pd.DataFrame(ranges, columns=['start', 'end']).to_csv('%s/%s_ranges.csv'%(self.output_path, fasta_id), index=False, sep='\t')

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
    parser.add_argument('-ns', '--no-sequences', action='store_true', help='Not to include sequences to the output')
    parser.add_argument('-r', '--repeats', action='store_true', help='To include soft-masked genome areas. By default, not included.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show the status of procesing or not. By default print stages info')

    args = parser.parse_args()
    if not os.path.isdir(args.output):
    	os.mkdir(args.output)
    if args.output.endswith('\\') or args.output.endswith('/'):
    	args.output = args.output.replace('\\', '').replace('/', '')
    if args.verbose:
    	blockPrint()
    if int(args.mdef) < int(args.tetdef):
    	print('Warning:  The allowed number of defective nucleotides (-tetdef) is more than the number of nucleotides (-mdef).', end='\n\n')
    finder = QuadruplexFinder(args.input, output_path = args.output, verbose = args.verbose, repeats=args.repeats, 
    	GC=args.GC, L=int(args.L) , q=int(args.q), nquadruplets=int(args.nquadruplets), mdef=int(args.mdef), tetdef=int(args.tetdef),
    	len_bulge = int(args.bulgelen), max_bulge = int(args.maxbulge), bulge_priority = args.bulge_priority)
    finder.run(print_sequences= not args.no_sequences)

if __name__ == '__main__':
    main()
