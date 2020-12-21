'''
utililties for reading and writing fasta files to / from a {header : sequence} dictionary
'''

def read_file(filename):
  '''
  read a fasta file, returns a {header : sequence} dictionary
  '''
  sequence = {}
  with open(filename, 'r') as f:
    head = ''
    for line in f:
      if line[0] == '>':
        key = line
        sequence[key] = ''
      else:
        sequence[key] += line[:-1]
  return sequence
  
def write_file(filename, sequence):
  '''
  creates a fasta file with the name filename containing the sequences in the sequence dictionary ({header : sequence} format)
  '''
  with open(filename, 'w') as f:
    for head in sequence:
      f.write(head)
      f.write(sequence[head] + '\n')
      