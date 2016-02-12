# run in bash
# perl /software/TargetScan/6.2/targetscan_60.pl miRNA.txt msa.aln targetscan_out.txt

# parse and reformat table with python

# define functions
def read_infile(infile):
    I = open(infile)
    miRNA_sites = {}
    # structure will be: miRNA_sites[site_id] = {miRNA_name, gene_name, species, # number of species, site_type}
    header = I.readline()
    for line in I:
	if '\tx\t' in line:
	    species = line.replace('\n', '').split('\t')[11].split(' ')
	    ID = line.split('\t')[0].split('_')
	    chromosome = ID[0]
	    start = int(ID[1])
	    end = int(ID[2])
	    gene_name = ID[3]
	    gene_id = '%s_%s' %(ID[4], ID[5])
	    block = ID[6].split('.')[-1]
	    utr_length = ID[7]
	    strand = ID[-1]
	    miRNA_sites[int(line.split('\t')[7])] = {'miRNA_name' : line.split('\t')[1], 'chr' : chromosome, 'start': start, 'end':end, 'gene_id': gene_id, 'block': block, 'utr_length': utr_length, 'strand': strand ,'gene_name' : gene_name, 'species' : species, 'number_of_species': len(species), 'site_type': line.split('\t')[8], 'MSA_start': int(line.split('\t')[3]), 'UTR_start': int(line.split('\t')[5])}
    I.close()
    return(miRNA_sites)

def write_stats(outfile, miRNA_sites):
    o = open(outfile, 'w')
    o.write('miRNA_name\tgene_id\tgene_name\tchr\texon_start\texon_end\tblock_in_exon\tstrand\testimated_utr_length\tMSA_site_start\tUTR_site_start\tsite_type\tsite_id\tnum_species\tspecies\n')
    for site in miRNA_sites:
	o.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %( miRNA_sites[site]['miRNA_name'],  miRNA_sites[site]['gene_id'],  miRNA_sites[site]['gene_name'],  miRNA_sites[site]['chr'],  miRNA_sites[site]['start'],  miRNA_sites[site]['end'],  miRNA_sites[site]['block'],  miRNA_sites[site]['strand'],  miRNA_sites[site]['utr_length'],  miRNA_sites[site]['MSA_start'],  miRNA_sites[site]['UTR_start'],  miRNA_sites[site]['site_type'], site,  miRNA_sites[site]['number_of_species'], ','.join(miRNA_sites[site]['species'])))
    o.close()
    return

# run script
miRNA_hits = read_infile('targetscan_out.txt')
write_stats('targetscan_out.summarized.txt', miRNA_hits)
