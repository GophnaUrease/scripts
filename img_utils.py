'''
Utilitiese for reading files from IMG
'''

def get_strain(genome, info):
  '''
  get the strain of a genome from a given IMG genome info file
  '''
  return head[head.find('[') + 1 : head.rfind(']')]

def get_taxonomy(genome, info, level):
  '''
  get the taxonomy at a given level of a genome from a given IMG genome info file
  '''
  if genome not in info:
    return "unknown"
  return info[genome][level]
  
def get_complete_taxonomy(genome, info):
  '''
  get the taxonomy at all levels of a given genome from a given IMG genome info file
  '''
  taxonomy_levels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
  return {taxonomy_level : info[genome][taxonomy_level] for taxonomy_level in taxonomy_levels}

def get_taxonomy_table(genome_info, level):
  '''
  get a complete taxonomy table from IMG genome info file
  '''
  return {level : get_complete_taxonomy(genome, info) for genome in genome_info}
  
def get_gene(head):
  '''
  get the gene ID of a gene via its header in a fasta file from IMG
  '''
  return head.split()[0][1:]
  
def get_genome(head, gene_info):
  '''
  get the genome ID of a gene via its header in a fasta file from IMG, and the gene_info dictionary containing information about this gene
  '''
  gene = get_gene(head)
  return gene_info[get_gene(head)]['Genome ID']

def read_info(filename):
  '''
  read a IMG file as a dictionary of dictionaries (each entry will have a dictionary of field name : value)
  '''
  info = {}
  with open(filename, "r") as info_file:
    key = info_file.readline().split('\t')
    for line in info_file:
      line = line.split('\t')
      id = line[0]
      info[id] = {}
      for i in xrange(len(line)):
        info[id][key[i]] = line[i]
  return info