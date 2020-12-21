'''
analysis pipeline of the urease gene data, normalization by the gyra gene as a marker
'''

################################################# imports #################################################

import sys
import os

import fasta_utils
import img_utils
import gene_analysis_utils
import metaphlan_utils
import ureasa_tools

################################################ constants ################################################

# general data #
UREASA_DATABASE = # a fasta file containing all the UreA genes (DNA nucleotides)
GYRA_DATABASE = # a fasta file containing all the GyrA genes (DNA nucleotides)
IMG_GENOME_INFO_FILE = # a xls file downloaded from IMG with the data about the genomes
                       # of the all the species the UreA and GyrA genes were taken from
UREASA_GENE_INFO_FILE = # a xls file downloaded from IMG with the data about all the UreA genes
GYRA_GENE_INFO_FILE = # a xls file downloaded from IMG with the data about all the GyrA genes
UREASA_GENE_LENGTH = 1719 # average lengths
GYRA_GENE_LENGTH = 2571

# PRISM data #
PRISM_UREASA_FOLDER = # the folder containing all the data regarding UreA in each sample (contigs, blastn output)
PRISM_GYRA_FOLDER = # the folder with all the data regarding GyrA in each sample (contigs, blastn output)
PRISM_OUTPUT_FOLDER = # the folder that the output files of the samples will be written to
PRISM_RUN_TABLE = # the run-table of PRISM, NLIBD and LLDeep samples from SRA
PRISM_META_DATA = # the metadata of the samples from PRISM, NLIBD and LLDeep
PRISM_METAPHLAN2 = # the metaphlan2 output of the samples from PRISM
PRISM_SAMPLES = ['SRR' + str(i) for i in xrange(6468499, 6468719)]

# NLIBD and LLDEEP data #
DUTCH_OUTPUT_FOLDER = # the folder that the output files of the samples from NLIBD and LLDeep will be written to
DUTCH_METAPHLAN2 = # the metaphlam2 output of the samples from NLIBD and LLDeep

# 1000IBD data #
C_1000IBD_UREASA_FOLDER = # the folder containing all the data regarding UreA in each sample (contigs, blastn output)
C_1000IBD_GYRA_FOLDER = # the folder with all the data regarding GyrA in each sample (contigs, blastn output)
C_1000IBD_OUTPUT_FOLDER = # the folder that the output files of the samples will be written to
C_1000IBD_metadata = # the metadata of the samples from 1000IBD
C_1000IBD_METAPHLAN2 = # the metaphlan2 output of the samples from 1000IBD
C_1000IBD_SAMPLES = set(filename.split('_')[0] for filename in reduce(lambda l1, l2: l1 + l2, [filenames for (dirpath, dirnames, filenames) in os.walk(C_1000IBD_UREASA_FOLDER)]))
C_1000IBD_SAMPLES.remove('UC')
C_1000IBD_SAMPLES.remove('CD')
C_1000IBD_SAMPLES.remove('IBDU')

# pouch data #
POUCH_UREASA_FOLDER = # the folder containing all the data regarding UreA in each sample (contigs, blastn output)
POUCH_GYRA_FOLDER = # the folder with all the data regarding GyrA in each sample (contigs, blastn output)
POUCH_OUTPUT_FOLDER = # the folder that the output files of the samples will be written to
POUCH_META_DATA = # the metadata of the samples from patients with a pouch
POUCH_METAPHLAN2 = # the metaphlan2 output of the samples from patients with a pouch
POUCH_SAMPLES = set(filename.split('_')[0] for filename in reduce(lambda l1, l2: l1 + l2, [filenames for (dirpath, dirnames, filenames) in os.walk(POUCH_UREASA_FOLDER)]))

# output configuration #
COMBINED_DATA_FOLDER =  # the folder that the output files of all the samples will be written to
FREQUENCY_DATA_FILENAME = 'frequency_data'
SAMPLE_DATA_FILENAME = 'sample_data'

# cuttoffs #
FREQ_CUTOFF = 0 # minimum relative abundance of a taxon in at least SAMPLES_CUTOFF samples to pass filtering (0-100)
SAMPLES_CUTOFF = 0 # minimum number of samples the taxon appears in with relative abundance of at least FREQ_CUTOFF to pass filtering
QUALITY_THRESHOLD = 40.0 # minimum total coverage of marker gene in the sample to pass filtering
UREASE_MIN_COVEREGE_REQUEST = 0.7 # minimum % of nucleotides in a urease gene from a taxon covered by cotings to include it in the sample

# streptococci data #
UREASA_POSITIVE_STREPTOCOCCI = [species for species in ureasa_tools.get_ureasa_positive_species(UREASA_DATABASE, UREASA_GENE_INFO_FILE, IMG_GENOME_INFO_FILE) if species.startswith('Streptococcus')]
DEFAULT_STREPTOCOCCI_RATIOS = {'Streptococcus' : {'Streptococcus thermophilus': 0.5, 'Streptococcus non-thermophilus': 0.5}} # shouldn't be used because it was validated all samples that pass our quality threshold have the metaphlan data

################################################ functions ################################################

# metadata organization #

def get_split_streptococci_abudnace(filename, get_genus, get_species):
  '''
  Parameters: filename - the name of the file with the methaphlan2 relative abundance table, 
         get_genus and get_species - functions that get the genus and sepcies names from the taxon name in the metaphlan2 table
  Return value: a dictionary than contains for each sample the relative abundance of streptococcus thermophilus and of other urease producing streptococci,
                among all the urease producing streptococci
  '''
  relative_abundances = metaphlan_utils.read_metaphlan(filename)
  ratio_dict = {}
  for sample in relative_abundances:
    thermophilus, non_thermophilus = 0.0, 0.0
    for metaphlan_name in relative_abundances[sample]:
      genus = get_genus(metaphlan_name)
      if genus == 'Streptococcus':
        is_thermophilus = get_species(metaphlan_name).startswith("thermoph")
        if is_thermophilus:
          thermophilus += relative_abundances[sample][metaphlan_name]
        elif genus + ' ' + get_species(metaphlan_name) in UREASA_POSITIVE_STREPTOCOCCI: # if it's a urease producing streptococcus other than streptococcus thermophilus
          non_thermophilus += relative_abundances[sample][metaphlan_name]
    total_urease_positive_streptococci = thermophilus + non_thermophilus
    if total_urease_positive_streptococci != 0:
      thermophilus, non_thermophilus = thermophilus / total_urease_positive_streptococci, non_thermophilus / total_urease_positive_streptococci
    ratio_dict[sample] = {'Streptococcus' : {'Streptococcus thermophilus': thermophilus, 'Streptococcus non-thermophilus': non_thermophilus}}
  return ratio_dict

def get_pouch_metadata():
  '''
  Returns a dictionary, the keys are the sample IDs of the pouch samples, the values are: (grouped_phenotype, phenotype, calprotectin level, calprotectiv level, streptococcus ratio, patient ID).
  The calprotectin appears twice for compatibility with the gerenal metadata format we use here for the main IBDs.
  '''  
  
  def _categorize_calprotectin(calprotectin):
    '''
    gets the calprotectin level in the sample (mcg / g) and returns the clinical interpretation - low (<=300) or high (>300)
    '''
    if calprotectin in {"NA", " NA"}:
      return "NA"
    if float(calprotectin) <= 300:
      return "low (<=300)"
    return "high (>300)"
    
  def _get_grouped_phenotype(phenotype):
    if phenotype == "Normal" or phenotype == "FAP" or phenotype == "Ultra Normal":
      return "Normal"
    if phenotype == "Chronic" or phenotype == "CLDP":
      return "Chronic+CLDP"
    return phenotype
    
  metadata = {}
  strepococci_abundances = get_split_streptococci_abudnace(POUCH_METAPHLAN2, metaphlan_utils.get_genus_frmt1, metaphlan_utils.get_species_frmt1)
  with open(POUCH_META_DATA, "r") as metadata_file:
    for line in metadata_file:
      sample_data = line.split(',')
      if sample_data[1] != 'sample_name':
        sample_id, patient_id, phenotype, calprotectin  = sample_data[1], sample_data[2], sample_data[3], sample_data[9]
        streptococci_abundance = DEFAULT_STREPTOCOCCI_RATIOS
        if sample_id in strepococci_abundances:
          streptococci_abundance = strepococci_abundances[sample_id]
        else:
          print sample_id
        metadata[sample_id] = (_get_grouped_phenotype(phenotype), phenotype, calprotectin, _categorize_calprotectin(calprotectin), streptococci_abundance, patient_id)
  return metadata

def get_C_1000IBD_metadata():
  '''
  Returns a dictionary, the keys are the sample IDs of the C_1000IBD samples, the values are: (disease, disease+location, calprotectin level, calprotectin category (low / intermediate / high), streptococcus ratio).
  '''   
  def _categorize_calprotectin(calprotectin):
    '''
    gets the calprotectin level in the sample (mcg / g) and returns the common clinical interpretation - low (<50), intermediate (50-200) or high (>200)
    '''
    if calprotectin == "NA":
      return "NA"
    if calprotectin == 'under40' or int(calprotectin) < 50:
      return "low (<50)"
    if int(calprotectin) <= 200:
      return "intermediate (50-200)"
    return "high (>200)"

  def _get_location(phenotype, location):
    '''
    gets the phenotype and its location, and returns the location string in the common format used
    '''
    if phenotype != "CD":
      return phenotype
    if location == "L1terminalileum":
      return "CD (ileal)"
    if location == "L2colon":
      return "CD (Colon)"
    if location == "L3ileocolon":
      return "CD (Ileocolon)"
    return "CD (Unknown locatoin)"

  metadata = {}
  strepococci_abundances = get_split_streptococci_abudnace(C_1000IBD_METAPHLAN2, metaphlan_utils.get_genus_frmt2, metaphlan_utils.get_species_frmt2)
  with open(C_1000IBD_metadata, "r") as metadata_file:
    for line in metadata_file:
      sample_data = line.split(',')
      if sample_data[0] != '1000IBDID':
        sample_id, phenotype, location, calprotectin = sample_data[0], sample_data[5], sample_data[6], sample_data[4]
        metaphlan_sample_id = 'profiled_' + sample_id
        streptococci_abundance = DEFAULT_STREPTOCOCCI_RATIOS
        if metaphlan_sample_id in strepococci_abundances:
          streptococci_abundance = strepococci_abundances[metaphlan_sample_id]
        else:
          print metaphlan_sample_id
        metadata[sample_id] = (phenotype, _get_location(phenotype, location), calprotectin, _categorize_calprotectin(calprotectin), streptococci_abundance)
  return metadata
            
def get_prism_metadata():
  '''
  Returns 2 dictionaries, the first is for PRISM (American), the second is for dutch (NLIBD and LLDEEP) samples.
  The keys are the sample SRRs, the values are: (disease, disease+location, calprotectin level, calprotectin category (low / intermediate / high), streptococcus ratio).
  '''
  def _categorize_calprotectin(calprotectin):
    '''
    gets the calprotectin level in the sample (mcg / g) and returns the common clinical interpretation - low (<50), intermediate (50-200) or high (>200)
    '''
    if calprotectin == "#N/A":
      return "NA"
    if int(calprotectin) < 50:
      return "low (<50)"
    if int(calprotectin) <= 200:
      return "intermediate (50-200)"
    return "high (>200)"
    
  def _match_sra():
    '''
    retuns a sra : srr dictionary of all the samples
    '''
    srr_to_sra = {}
    with open(PRISM_RUN_TABLE, "r") as f:
      for line in f:
        srr, sra = line.split('\t')[11], line.split('\t')[6]
        srr_to_sra[srr] = sra
    return srr_to_sra
      
  def _get_metadata_by_sra():
    '''
    Returns a dictionary, the keys are the sample SRAs, the values are: (disease, disease+location, calprotectin level, calprotectiv category (low / medium / high), streptococcus ratio, geographic origin (American - true, Dutch - false)).  
    '''
    def _get_location(phenotype, location):
      '''
      gets the phenotype and its location, and returns the location string in the common format used
      '''
      if phenotype != 'CD':
        return phenotype
      if location == 'NA':
        return 'CD (unknown location)'
      return 'CD ' + location.split(' ')[1]
      
    with open(PRISM_META_DATA, "r") as metadata_file:
      for i, line in enumerate(metadata_file):
        if i == 1:
          cohort = line[:-2]
        if i == 2:
          SRAs = line[:-2]
        if i == 4:
          phenotype = line[:-2]
        if i == 5:
          calprotectin = line[:-2]
        if i == 10:
          location = line[:-2]
    SRAs = SRAs.split(',')[1:]
    phenotype = phenotype.split(',')[1:]
    calprotectin = calprotectin.split(',')[1:]
    location = location.split(',')[1:]
    cohort = cohort.split(',')[1:]
    metadata = {}
    strepococci_abundances = get_split_streptococci_abudnace(PRISM_METAPHLAN2, metaphlan_utils.get_genus_frmt2, metaphlan_utils.get_species_frmt2)
    dutch_strepococci_abundances = get_split_streptococci_abudnace(DUTCH_METAPHLAN2, metaphlan_utils.get_genus_frmt2, metaphlan_utils.get_species_frmt2)
    strepococci_abundances.update(dutch_strepococci_abundances)
    for i in xrange(len(SRAs)):
      if SRAs[i] not in strepococci_abundances:
        print SRAs[i]
        strepococci_abundances[SRAs[i]] = DEFAULT_STREPTOCOCCI_RATIOS
    for i in xrange(len(SRAs)):
      metadata[SRAs[i]] = (phenotype[i], _get_location(phenotype[i], location[i]), calprotectin[i], _categorize_calprotectin(calprotectin[i]), strepococci_abundances[SRAs[i]], cohort[i].startswith("PRISM"))
    return metadata
    
  SRR_TO_SRA = _match_sra()
  SRA_TO_STATUS = _get_metadata_by_sra()
  return {sample : SRA_TO_STATUS[SRR_TO_SRA[sample]] for sample in PRISM_SAMPLES if SRA_TO_STATUS[SRR_TO_SRA[sample]][5]}, {sample : SRA_TO_STATUS[SRR_TO_SRA[sample]] for sample in PRISM_SAMPLES if not SRA_TO_STATUS[SRR_TO_SRA[sample]][5]}
         
# taxonomic assignment #
         
def get_taxonomy_data(genome_info_file, gene_info_file, database_file):
  '''
  Input: database (fasta of sequences), information of genes and genomes (to map gene IDs to the bacterial strains they come from)
  Output: 2 dictionaries, key is gene ID, value is taxonomy for calculations in the first, and the group for dominating urease producer analysis in the second
  '''  
   
  def _get_customized_taxonomy(head, gene_info, genome_info):
    '''
    Get the taxon name from the fasta header:
    - Usually genus name
    - Blautia and Ruminococcus are combined, and so Citrobacter and Enterobacter
    - Species name in chosen genera where only one species appears in the results
    '''
    specify_species = ("Bifidobacterium", "Fusicatenibacter", "Lactobacillus", "Eisenbergiella") # genera where species name is returned
    genus = genome_info[img_utils.get_genome(head, gene_info)]['Genus']
    if genus in specify_species:
      return img_utils.get_taxonomy(img_utils.get_genome(head, gene_info), genome_info, 'Species')
    if genus == "Blautia" or genus == "Ruminococcus":
      return "Blautia or Ruminococcus"
    if genus == "Citrobacter" or genus == "Enterobacter":
      return "Citrobacter or Enterobacter"
    return genus
  
  def _get_group(head, gene_info, genome_info):
    '''
    Get the urease producing group from the fasta header by our classification
    '''
    species = img_utils.get_taxonomy(img_utils.get_genome(head, gene_info), genome_info, 'Species')
    if species == 'Streptococcus thermophilus' or species == 'Bifidobacterium longum':
      return 'Probiotic'
    order = img_utils.get_taxonomy(img_utils.get_genome(head, gene_info), genome_info, 'Order')
    if order == 'Bacteroidales' or order == 'Clostridiales':
      return 'Bacteroidales and Clostridiales'
    clas = img_utils.get_taxonomy(img_utils.get_genome(head, gene_info), genome_info, 'Class')
    phylum = img_utils.get_taxonomy(img_utils.get_genome(head, gene_info), genome_info, 'Phylum')
    if clas == 'Bacilli' or phylum == 'Proteobacteria':
      return 'Bacilli and Proteobacteria'
    return clas
      
  genome_info = img_utils.read_info(genome_info_file)
  gene_info = img_utils.read_info(gene_info_file)
  database = fasta_utils.read_file(database_file)
  taxonomy_dict = {}
  group_dict = {}
  for head in database.keys():
    taxon, group = _get_customized_taxonomy(head, gene_info, genome_info), _get_group(head, gene_info, genome_info)
    taxonomy_dict[img_utils.get_gene(head)] = taxon
    if group not in group_dict:
      group_dict[group] = set([])
    group_dict[group].add(taxon)
  return taxonomy_dict, group_dict

# organization of data after the analysis #

def organize_analyzed_data(pouch_metadata):
  '''
  Combine the dutch data, and all the ibd data
  Add patient ID to the pouch data
  '''
  # combine the dutch frequency data
  with open(C_1000IBD_OUTPUT_FOLDER + FREQUENCY_DATA_FILENAME, 'a') as fa:
    with open(DUTCH_OUTPUT_FOLDER + FREQUENCY_DATA_FILENAME, 'r') as fr:
      fa.write(fr.read())
  # combine the dutch sample data
  with open(C_1000IBD_OUTPUT_FOLDER + SAMPLE_DATA_FILENAME, 'a') as fa:
    with open(DUTCH_OUTPUT_FOLDER + SAMPLE_DATA_FILENAME, 'r') as fr:
      fa.write(fr.read())
      
  # combine all the ibd frequency data
  with open(COMBINED_DATA_FOLDER + FREQUENCY_DATA_FILENAME, 'w') as fw:
    with open(C_1000IBD_OUTPUT_FOLDER + FREQUENCY_DATA_FILENAME, 'r') as fr:
      fw.write(fr.read())
    with open(PRISM_OUTPUT_FOLDER + FREQUENCY_DATA_FILENAME, 'r') as fr:
      fw.write(fr.read())
  # combine all the ibd frequency data
  with open(COMBINED_DATA_FOLDER + SAMPLE_DATA_FILENAME, 'w') as fw:
    with open(C_1000IBD_OUTPUT_FOLDER + SAMPLE_DATA_FILENAME, 'r') as fr:
      fw.write(fr.read())
    with open(PRISM_OUTPUT_FOLDER + SAMPLE_DATA_FILENAME, 'r') as fr:
      fw.write(fr.read())
      
  # add patient ID to pouch frequency data 
  with open(POUCH_OUTPUT_FOLDER + FREQUENCY_DATA_FILENAME + '_tmp', 'w') as fw:
    with open(POUCH_OUTPUT_FOLDER + FREQUENCY_DATA_FILENAME, 'r') as fr:
      for line in fr:
        fw.write(line[:-1] + '\t' + pouch_metadata[line.split('\t')[4]][5] + '\n')
  os.remove(POUCH_OUTPUT_FOLDER + FREQUENCY_DATA_FILENAME)
  os.rename(POUCH_OUTPUT_FOLDER + FREQUENCY_DATA_FILENAME + '_tmp', POUCH_OUTPUT_FOLDER + FREQUENCY_DATA_FILENAME)   
  # add patient ID to pouch sample data   
  with open(POUCH_OUTPUT_FOLDER + SAMPLE_DATA_FILENAME + '_tmp', 'w') as fw:
    with open(POUCH_OUTPUT_FOLDER + SAMPLE_DATA_FILENAME, 'r') as fr:
      for line in fr:
        fw.write(line[:-1] + '\t' + pouch_metadata[line.split('\t')[2]][5] + '\n')
  os.remove(POUCH_OUTPUT_FOLDER + SAMPLE_DATA_FILENAME)
  os.rename(POUCH_OUTPUT_FOLDER + SAMPLE_DATA_FILENAME + '_tmp', POUCH_OUTPUT_FOLDER + SAMPLE_DATA_FILENAME) 
  
# main #
  
def main():
  taxonomy, groups = get_taxonomy_data(IMG_GENOME_INFO_FILE, UREASA_GENE_INFO_FILE, UREASA_DATABASE)
  taxonomy2, _ = get_taxonomy_data(IMG_GENOME_INFO_FILE, GYRA_GENE_INFO_FILE, GYRA_DATABASE)
  taxonomy.update(taxonomy2)
  groups['Probiotic'].add('Streptococcus thermophilus')
  groups['Bacilli and Proteobacteria'].add('Streptococcus non-thermophilus')
  LEFT = {'Bacteroidales and Clostridiales'} # to show on the left of the graph (i.e. low index for ordering)
  RIGHT = {'Bacilli and Proteobacteria'} # to show on the right of the graph (i.e. high index for ordering)
  group_data = (groups, (LEFT, RIGHT))
  
  # PRISM american cohort
  prism_metadata, dutch_metadata = get_prism_metadata()
  gene_analysis_utils.analyze((UREASA_GENE_LENGTH, UREASE_MIN_COVEREGE_REQUEST, GYRA_GENE_LENGTH, None), taxonomy, prism_metadata, (PRISM_UREASA_FOLDER, '_ureasa'), (PRISM_GYRA_FOLDER, '_gyra'), group_data, (PRISM_OUTPUT_FOLDER, FREQUENCY_DATA_FILENAME, SAMPLE_DATA_FILENAME), (FREQ_CUTOFF, SAMPLES_CUTOFF, QUALITY_THRESHOLD))
  # Dutch cohorts
  gene_analysis_utils.analyze((UREASA_GENE_LENGTH, UREASE_MIN_COVEREGE_REQUEST, GYRA_GENE_LENGTH, None), taxonomy, dutch_metadata, (PRISM_UREASA_FOLDER, '_ureasa'), (PRISM_GYRA_FOLDER, '_gyra'), group_data, (DUTCH_OUTPUT_FOLDER, FREQUENCY_DATA_FILENAME, SAMPLE_DATA_FILENAME), (FREQ_CUTOFF, SAMPLES_CUTOFF, QUALITY_THRESHOLD))
  # C_1000IBD dutch cohort
  C_1000IBD_metadata = get_C_1000IBD_metadata()
  gene_analysis_utils.analyze((UREASA_GENE_LENGTH, UREASE_MIN_COVEREGE_REQUEST, GYRA_GENE_LENGTH, None), taxonomy, C_1000IBD_metadata, (C_1000IBD_UREASA_FOLDER, '_ureasa'), (C_1000IBD_GYRA_FOLDER, '_gyra'), group_data, (C_1000IBD_OUTPUT_FOLDER, FREQUENCY_DATA_FILENAME, SAMPLE_DATA_FILENAME), (FREQ_CUTOFF, SAMPLES_CUTOFF, QUALITY_THRESHOLD))
  # Our Pouch Data
  pouch_metadata = get_pouch_metadata()
  gene_analysis_utils.analyze((UREASA_GENE_LENGTH, UREASE_MIN_COVEREGE_REQUEST, GYRA_GENE_LENGTH, None), taxonomy, pouch_metadata, (POUCH_UREASA_FOLDER, '_ureasa'), (POUCH_GYRA_FOLDER, '_gyra'), group_data, (POUCH_OUTPUT_FOLDER, FREQUENCY_DATA_FILENAME, SAMPLE_DATA_FILENAME), (FREQ_CUTOFF, SAMPLES_CUTOFF, QUALITY_THRESHOLD))

  organize_analyzed_data(pouch_metadata)
  
if __name__ == "__main__":
	sys.exit(main())
