'''
Created on Feb 9, 2018

@author: mmp
'''
import unittest, os
from CodonsTable import CodonsTable

class Test(unittest.TestCase):

	
	def test_count_codons(self):
		codons_table = CodonsTable()
		codons_table.count_codons(os.path.join('.', 'tests', 'test_fasta_file.fasta'))
		self.assertEqual(3, codons_table.lst_data[0][21])
		self.assertEqual(1, codons_table.lst_data[21][21])
		self.assertEqual(2, codons_table.lst_data[0][63])
		
	def test_count_codons_by_sequence(self):
		codons_table = CodonsTable()
		self.assertEqual(0, codons_table.lst_data[0][21])
		codons_table.count_codons_by_sequence('AAACCC')
		self.assertEqual(1, codons_table.lst_data[0][21])
		codons_table.count_codons_by_sequence('AAACCC')
		self.assertEqual(2, codons_table.lst_data[0][21])
		codons_table.count_codons_by_sequence('AAACCCCCC')
		self.assertEqual(3, codons_table.lst_data[0][21])
		self.assertEqual(1, codons_table.lst_data[21][21])
		
		self.assertEqual(0, codons_table.lst_data[0][63])
		codons_table.count_codons_by_sequence('AAATTT')
		self.assertEqual(1, codons_table.lst_data[0][63])
		codons_table.count_codons_by_sequence('AAATTTC')
		self.assertEqual(2, codons_table.lst_data[0][63])
		
	def test_count_pair_codons(self):
		codons_table = CodonsTable()
		self.assertEqual(0, codons_table.lst_data[0][0])
		codons_table.count_pair_codons('A', 'A', 'A', 'A', 'A', 'A')
		self.assertEqual(1, codons_table.lst_data[0][0])
		
		self.assertEqual(0, codons_table.lst_data[63][0])
		codons_table.count_pair_codons('T', 'T', 'T', 'A', 'A', 'A')
		self.assertEqual(1, codons_table.lst_data[63][0])
		
		self.assertEqual(0, codons_table.lst_data[63][63])
		codons_table.count_pair_codons('T', 'T', 'T', 'T', 'T', 'T')
		self.assertEqual(1, codons_table.lst_data[63][63])
		codons_table.count_pair_codons('T', 'T', 'T', 'T', 'T', 'T')
		self.assertEqual(2, codons_table.lst_data[63][63])
		
		codons_table.count_pair_codons('T', 'T', 'T', 'A', 'c', 'd')
		
		
	def test_get_pos_codon(self):
		codons_table = CodonsTable()
		self.assertEqual(0, codons_table.get_pos_codon('A', 'A', 'A'))
		self.assertEqual(-1, codons_table.get_pos_codon('A', 'A', 'a'))
		self.assertEqual(-1, codons_table.get_pos_codon('A', 'a', 'C'))
		self.assertEqual(-1, codons_table.get_pos_codon('z', 'A', 'G'))
		self.assertEqual(63, codons_table.get_pos_codon('T', 'T', 'T'))
		self.assertEqual(3, codons_table.get_pos_codon('A', 'A', 'T'))
		self.assertEqual(15, codons_table.get_pos_codon('A', 'T', 'T'))
		self.assertEqual(21, codons_table.get_pos_codon('C', 'C', 'C'))

	def test_get_codon_by_number(self):
		codons_table = CodonsTable()
		self.assertEqual('AAA', codons_table.get_codon_by_number(0))
		self.assertEqual('AAA', codons_table.get_codon_by_number(64))
		self.assertEqual('AAC', codons_table.get_codon_by_number(1))
		self.assertEqual('TTT', codons_table.get_codon_by_number(63))
		self.assertEqual('CCC', codons_table.get_codon_by_number(21))
		self.assertEqual('ATT', codons_table.get_codon_by_number(15))
		
		
		
if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.test_count_codons_by_sequence']
	unittest.main()