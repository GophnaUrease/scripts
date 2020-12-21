'''
Utilities for reading a MEGAN gene-centric assembler output file
'''

def get_coverages(header):
  '''
  returns a list of the average coverage of the contigs in a scaffold in a MEGAN gene-centric assembler output fasta file
  '''
  coverage_strings = [substring for substring in header.replace('[', ' ').replace(']', ' ').split() if substring.find('avCoverage=') != -1]
  return [float(coverage_string.split('=')[1]) for coverage_string in coverage_strings]
  
def get_normalized_coverage(header):
  '''
  returns the average coverage of a nucleotide * contig length of a contig in a MEGAN gene-centric assembler output fasta file
  '''
  coverage = float(header.split('avCoverage=')[1])
  length = float(header.split('length=')[1].split(' ')[0])
  return coverage * length

def get_name(header):
  '''
  get the name of a gontig by its header in a MEGAN gene-centric assembler output fasta file
  '''
  return header.split()[0][1:]
