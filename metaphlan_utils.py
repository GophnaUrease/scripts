'''
Utilities for interpretation of the metaphlan relative abundance files
'''

def read_metaphlan(filename):
  '''
  read the metaphlan file relative abundance file and return a dictionary containing the relative abundance data
  '''
  species_data = {}
  with open(filename, 'r') as metaphlan:
    samples = metaphlan.readline().rstrip().split('\t')
    for line in metaphlan:
      data = line.rstrip().split('\t')
      species_data[data[0]] = data[1:]
  relative_abundances = {}
  for i in range(1, len(samples)):
    relative_abundances[samples[i]] = {}
    for species in species_data:
      if species not in relative_abundances[samples[i]]:
        relative_abundances[samples[i]][species] = 0.0
      relative_abundances[samples[i]][species] += float(species_data[species][i-1])
  return relative_abundances
  
def get_genus_frmt1(metaphlan_name):
  '''
  get the taxon name in the "genus_species" format (e.g. Streptococcus_thermophilus) and return the genus name
  '''
  return metaphlan_name.split('_')[0]
  
def get_species_frmt1(metaphlan_name):
  '''
  get the taxon name in the "genus_species" format (e.g. Streptococcus_thermophilus) and return the species name
  '''
  return metaphlan_name.split('_')[1]
    
def get_genus_frmt2(metaphlan_name):
  '''
  get the taxon name in the "g__genus|s__genus_species" format (g__Streptococcus|s__Streptococcus_thermophilus) and return the genus name
  '''
  return metaphlan_name.split('|')[0][3:]
  
def get_species_frmt2(metaphlan_name):
  '''
  get the taxon name in the "g__genus|s__genus_species" format (g__Streptococcus|s__Streptococcus_thermophilus) and return the species name
  '''
  return metaphlan_name.split('|')[1][3:].split('_')[1]