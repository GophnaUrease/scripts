'''
general analysis functions for the abundance of single copy genes using a universal marker
'''
import fasta_utils
import gc_assembler_utils
import blast_utils
import math

#START_SIGN should be smaller than END_SIGN, this is what the capitalization is for
START_SIGN = 'S'
END_SIGN = 'e'

def build_abundance_table(filename, taxonomy, gene_length, perc_coverage = None):
  '''
  Input: file name (for blastn file of gene matches, and assembly file of contigs for coverage), taxonomy dictionary (gene IDs to taxonomy),
         gene_length, and percenteage of the nucleotides covered in the gene for it to be considered present in the sample - or None if such an analysis is
         not required
  Output: 2 dictionaries, for every taxonomic unit, its relative abundance (first dictionary) and absolute coverage (second dictionary)
  '''
  def _get_best_matches(blastn_output, taxonomy_dict):
    '''
    match each assembled contig to its best match (by taxonomy), returns for each query the best taxonomic match, it's beginning (nucleotide position) and end
    '''
    best_matches = {}
    for query in blastn_output:
      classified_hits = [hit for hit in blastn_output[query] if taxonomy_dict[hit['subject acc.ver']] != "unclassified"]
      if len(classified_hits) != 0:
        max_score = max(hit['bit score'] for hit in classified_hits)
        candidates = set([taxonomy_dict[hit['subject acc.ver']] for hit in classified_hits if hit['bit score'] == max_score]) # all matches
        start, end =  [(hit['s. start'], hit['s. end']) for hit in classified_hits if hit['bit score'] == max_score][0]
        best_matches[query] = (' or '.join(candidates), int(start), int(end))
    return best_matches
  
  def _calculate_abundances(matches, coverage, gene_length, perc_coverage):
    '''
    input: matches - matches the contig ID to its best match (taxonomy, start, end), 
    coverage - matches the conting ID to its normalized coverage, 
    perc_coverage - percenteage of the nucleotides covered in the gene for it to be considered present in the sample - or None if such an analysis is
         not required
    returns 2 dictionaries, relative abundance of each taxon and absolute coverage of each taxon
    '''
    def _normalize(dictionary):
      '''
      input: coverage of each taxon
      output: relative abundance of each taxon
      '''
      total_sum = sum(dictionary.values())
      if total_sum == 0:
        return {}
      return {key : dictionary[key] * 100.0 / total_sum for key in dictionary}
    
    def _get_taxon_coverages(matches, coverage, gene_length):
      '''
      swipe line algorithm, from the matches and coverages of each contig calculates for each taxon the % of nucleotides of its gene that are covered, and its normalized coverage (by all the contigs covering it)
      '''
      def _get_segment_coverage(segments):
        '''
        the swipe line algorithm itself, gets the start and end of all the contigs covering the gene from a specific taxon, and calculates the % of its nucleotides that are covered by them
        '''
        covered = 0
        last_open = 0
        count_open = 0
        for segment in sorted(segments):
          if segment[1] == START_SIGN:
            if count_open == 0:
              last_open = segment[0]
            count_open += 1
          else:
            count_open -= 1
            if count_open == 0:
              covered += segment[0] - last_open + 1
        return covered
      
      taxons = set([taxon for taxon, _, _ in matches.values()])
      segments = {taxon : [] for taxon in taxons}
      for match in matches:
        taxon, start, end = matches[match]
        segments[taxon].append((start, START_SIGN, coverage[match]))
        segments[taxon].append((end, END_SIGN, coverage[match]))
      taxon_coverages = {}

      for taxon in taxons:
        taxon_coverages[taxon] = (float(_get_segment_coverage(segments[taxon])) / gene_length, sum(cov for _, _, cov in segments[taxon]) / 2)
      return taxon_coverages
      
    taxon_coverages = _get_taxon_coverages(matches, coverage, gene_length) 
    if perc_coverage:
      abundances = {taxon : taxon_coverages[taxon][1] for taxon in taxon_coverages if taxon_coverages[taxon][0] > perc_coverage}
    else:
      abundances = {taxon : taxon_coverages[taxon][1] for taxon in taxon_coverages}
    return _normalize(abundances), abundances
  
  try:
    assembly_output = fasta_utils.read_file(filename + '.fasta')
    blastn_output = blast_utils.read_outfmt7_file(filename + '.blastn')
  except IOError:
    assembly_output = {}
    blastn_output = {}
  coverage_dict = {gc_assembler_utils.get_name(head) : gc_assembler_utils.get_normalized_coverage(head) / gene_length for head in assembly_output.keys()} # contig ID : coverage (normalized by length)
  taxonomy_matches = _get_best_matches(blastn_output, taxonomy) # contig ID : taxonomy, start, end
  abundance_dict, coverage_dict = _calculate_abundances(taxonomy_matches, coverage_dict, gene_length, perc_coverage)
  
  return abundance_dict, coverage_dict
  
def filter_and_complete_frequencies(gene_frequencies, marker_coverage, filter_cutoffs):
  '''
  input: relative abundances of the analyzed gene, the absolute total coverages of the marker gene, the cutoffs for filtering:
  (minimum relative abundance of a taxon in at least the given number of samples, minimum number of samples the taxon appears in with at least the given relative abundance, minimum total coverage of marker gene in the sample)
  '''
  def _filter_taxons(frequencies, filter_cutoffs):
    '''
    frequencies - relative abundance of each taxon
    filter_cuttoffs = (minimum relative abundance of a taxon in at least the given number of samples, minimum number of samples the taxon appears in with at least the given relative abundance, minimum total coverage of marker gene in the sample)
    '''
    occurences = {}
    taxons = set([])
    for sample in frequencies:
      for taxon in frequencies[sample]:
        taxons.add(taxon)
        if taxon not in occurences:
            occurences[taxon] = 0
        if frequencies[sample][taxon] >= filter_cutoffs[0]:
          occurences[taxon] += 1
    return set(taxon for taxon in taxons if occurences[taxon] >= filter_cutoffs[1])
    
  def _filter_and_complete(frequencies, filtered_taxons):
    '''
    input: relative abundance of each taxon, a set of the taxons that passed the filtering
    updates the frequencies: for each sample add a 0.0 for each taxon that passed the filtering but is not found in the sample, remove all taxa that didn't pass the filtering and add 'Other' to the frequency table with their combined relative abundance
    '''
    for sample in frequencies:
      other_freq = 0.0
      taxons = frequencies[sample].keys()
      for taxon in taxons:
        if taxon not in filtered_taxons:
          other_freq += frequencies[sample].pop(taxon)
      frequencies[sample]['Other'] = other_freq
      for taxon in filtered_taxons:
        if taxon not in frequencies[sample]:
          frequencies[sample][taxon] = 0.0
    
  filtered_taxons = _filter_taxons(gene_frequencies, filter_cutoffs)
  _filter_and_complete(gene_frequencies, filtered_taxons)
  bad_quality = [sample for sample in marker_coverage if marker_coverage[sample] < filter_cutoffs[2]]

  for sample in bad_quality:
    marker_coverage.pop(sample)
    gene_frequencies.pop(sample)
      
def calculate_load(gene_coverages, marker_coverage):
  '''
  calculate gene load - i.e. the percentage of the gene producers in the community
  '''
  loads = {}
  for sample in marker_coverage:
    if marker_coverage[sample] != 0:
      loads[sample] = sum(gene_coverages[sample].values()) * 100.0 / marker_coverage[sample]
    else:
      loads[sample] = 0 
  return loads

def output_frequencies_data(output_folder, frequencies, metadata):
  '''
  output the relative abundance of the gene for each sample to the given folder
  '''
  for sample in frequencies:
    file_name = metadata[sample][0] + '_' + sample
    data = reversed(sorted([(frequencies[sample][taxon], taxon) for taxon in frequencies[sample]]))
    with open(output_folder + file_name, 'w') as output_file:
      for line in data:
        output_file.write(line[1] + '\t' + str(line[0]) + '%\n')
   
def output_all_data(output_names, gene_frequencies, metadata, separator, indices, loads, categories):
  '''
  output the data for all samples
  output_names: the output folder, the name of the file to create with the relative abundance of all taxa in all sam[les, the name of the file to create with the sample data (for each sample: the main disease, the specific disease, the sample ID, the index, the load (% producers), the calprotectin value, the calprotectin category (low / intermediate / high), the main producers category  
  '''
  output_folder, frequency_filename, sample_filename = output_names
  with open(output_folder + frequency_filename, 'w') as frequencies_file:
    for sample_id in gene_frequencies:
      sample_frequencies, main_disease, disease, index = gene_frequencies[sample_id], metadata[sample_id][0], metadata[sample_id][1], str(indices[sample_id])
      for taxon in sample_frequencies:
          frequencies_file.write(separator.join([main_disease, disease, taxon, str(sample_frequencies[taxon]), sample_id, index]) + '\n')
  with open(output_folder + sample_filename, 'w') as f: 
    for sample_id in gene_frequencies:
      main_disease, disease, calprotectin, calprotectin_category, index, load, category = metadata[sample_id][0], metadata[sample_id][1], metadata[sample_id][2], metadata[sample_id][3], str(indices[sample_id]), str(loads[sample_id]), categories[sample_id]
      f.write(separator.join([main_disease, disease, sample_id, index, load, calprotectin, calprotectin_category, category]) + '\n')

def calculate_index(gene_frequencies, marker_coverage, group_data):
  '''
  calculate the index for each sample by the given groups.
  gene_frequencies: relative abundance of the analyzed gene
  marker_coverage: total coverage of each gene
  group_data: the taxons belonging to each group for the dominating group analysis, and the groups that are counted for the index ("good" for smaller indices, "bad" for higher ones)
  '''
  INF = 1000000
  def _calc_load(group_taxa, groups, frequency):
    '''
    calculates the relative abundance of gene producers from some given groups
    '''
    taxa = set.union(*[group_taxa[group] for group in groups])
    return sum([frequency[taxon] for taxon in taxa if taxon in frequency])
    
  index = {}
  for sample in gene_frequencies:
    bad_load = _calc_load(group_data[0], group_data[1][1], gene_frequencies[sample])
    good_load = _calc_load(group_data[0], group_data[1][0], gene_frequencies[sample])
    depth = marker_coverage[sample] + 1
    if sum(gene_frequencies[sample].values()) == 0:
      bad_load = INF # to put all gene-nagative samples last
    if good_load == 0:
      good_load = 2.0 / depth
    if bad_load == 0:
      bad_load = 3.0 / depth
    index[sample] = math.log(bad_load/good_load)
  return index
 
def categorize(urease_frequencies, groups):
  '''
  categorizes the samples by the dominating group (>50.0 %) of gene producers
  '''
  def _calc_frequencies(taxa, frequency):
    return sum([frequency[taxon] for taxon in taxa if taxon in frequency])

  categories = {}
  for sample in urease_frequencies:
    if sum(urease_frequencies[sample].values()) == 0:
      categories[sample] = "No producers" # no gene producers at all
    else:
      for group in groups:
        if _calc_frequencies(groups[group], urease_frequencies[sample]) > 50.0:
          categories[sample] = group + " dominated"
      if sample not in categories:
        categories[sample] = "No dominating group"
  return categories
   
def split_taxa(gene_frequencies, gene_coverages, metadata):
  '''
  split the relative abundance of the given taxa (in the metadata) by their ratios in the sample (given in the metadata)
  '''
  for sample in metadata:
    split_dict = metadata[sample][4]
    for taxon in split_dict:
      if taxon in gene_frequencies[sample]:
        freq = gene_frequencies[sample].pop(taxon)
        for son_taxon in split_dict[taxon]:
          gene_frequencies[sample][son_taxon] = freq * split_dict[taxon][son_taxon]
        freq = gene_coverages[sample].pop(taxon)
        for son_taxon in split_dict[taxon]:
          gene_coverages[sample][son_taxon] = freq * split_dict[taxon][son_taxon]
         
def analyze(gene_lengths, taxonomy_classification, metadata, gene_folder, marker_folder, groups, output_folder, filter_cutoffs=(0,0,40)):
  '''
  analyses and saves the relative abundances of a single copy gene, and additional data from assembled contigs and blastn matches to a reference database of a single copy gene and a universal marker
  takes the data of the gene for the analysis and the marker from the given folders.
  Parameters:
  gene_lengths: a tuple, (the length of the gene for the analysis, the min % of nucleotides requested for it to be covered / None if no such requirement exist, length of merker, min % nucleotides that sohuld be covered / None)
  taxonomy_classification: a dictionary that matches gene IDs in the blastn matches to the taxonomy of this gene
  metadata: the metadata for each sample: a dictionary of tuple for each sample: main phenotype, detailed phenotype, calprotectin level, calprotectin category, ration data for splitting the taxonomy post analysis (in cases where the matches cannot be differentiated but there is additional data that the taxon can be splitted by)
  gene_folder: the folder where all the data files (assembled contig and their matches) of the gene to be analyzed is located
  marker_folder: the folder where all the data files (assembled contig and their matches) of the marker gene is located
  groups: the taxons belonging to each group for the dominating group analysis, and the groups that are counted for the index ("good" for smaller indices, "bad" for higher ones).
  output_folder: a tuple of 3 values: output folder, the name of the file to create with the relative abundance frequency data, and the name of the file to create with the sample data
  filter_cutoffs: (minimum relative abundance of a gene from a taxon in at least the given number of samples, minimum number of samples the gene from a taxon appears in with at least the given relative abundance, minimum total coverage of marker gene in the sample)
  '''

  gene_length, gene_cov_cutoff, marker_length, marker_cov_cutoff = gene_lengths
  gene_frequencies = {}
  gene_coverages = {}
  marker_coverages = {}
  for sample in metadata:
    _, marker_coverages[sample] = build_abundance_table(marker_folder[0] + sample + marker_folder[1], taxonomy_classification, marker_length, marker_cov_cutoff)
    gene_frequencies[sample], gene_coverages[sample] = build_abundance_table(gene_folder[0] + sample + gene_folder[1], taxonomy_classification, gene_length, gene_cov_cutoff)
  marker_coverage = {sample : sum(marker_coverages[sample].values()) for sample in marker_coverages}
  split_taxa(gene_frequencies, gene_coverages, metadata)
  
  output_frequencies_data(gene_folder[0], gene_frequencies, metadata)
  filter_and_complete_frequencies(gene_frequencies, marker_coverage, filter_cutoffs)
  output_all_data(output_folder, gene_frequencies, metadata, '\t', calculate_index(gene_frequencies, marker_coverage, groups), calculate_load(gene_coverages, marker_coverage), categorize(gene_frequencies, groups[0]))
 