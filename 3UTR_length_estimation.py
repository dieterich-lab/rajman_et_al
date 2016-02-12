# run in bash

# java -Xmx2048m -jar ~/isoSCM/IsoSCM-2.0.7.jar assemble -bam /data/projects/external/Gerhar_Schratt/Marek_Rajman/rat_hippocampal_neurons/EMBLE/RNAseq/isoSCM_utrs/bam/Control.sorted.bam -base control -s reverse_forward -t 3

# sed -i 's/;type /-/g' control.coverage.gtf
# awk '{print $1 "\t" $4 "\t" $5 "\t" $11 "\t" $10 "\t" $7}' control.coverage.gtf > control.coverage.bed
# sed -i 's/"//g' control.coverage.bed
# sed -i 's/;locus_id//g' control.coverage.bed

# sort control.coverage.bed > control.coverage.sorted.bed
# intersectBed -s -wo -a control.coverage.sorted.bed -b rn5.RefSeq.genes.bed > control.coverage.genes.sorted.bed
# intersectBed -s -wao -a control.coverage.genes.sorted.bed -b start_codon.rn5.RefSeq.sorted.bed > start_codon.control.coverage.sorted.bed 
# intersectBed -s -wao -a control.coverage.genes.sorted.bed -b stop_codon.rn5.RefSeq.sorted.bed > stop_codon.control.coverage.sorted.bed 



# run in python

# define input files
stop = 'stop_codon.control.coverage.sorted.bed' 
start = 'start_codon.control.coverage.sorted.bed' 
ref_names = 'refGene_id2name.txt'


# define function
def read_intersected_bed(Infile):
    I = open(Infile)
    utrs = {}
    while True:
	Line = I.readline()
	if not Line:break
	L = Line.split('\t')
	if L[3].startswith('locus'):
  	    # safe entries as utrs[locus_id] = {(chr, start, end, gene_id) : {'gene_id': 'codon_id' , 'start/stop_codon': , 'strand': , 'coverage': }}
	    if not L[3].split('-')[0] in utrs:
		utrs[L[3].split('-')[0]] = {}
	    #print(L)
	    if not L[9] in utrs[L[3].split('-')[0]]:
		utrs[L[3].split('-')[0]][L[9]] = {}
	    if L[16] == '.' or L[9] == L[16]:
		utrs[L[3].split('-')[0]][L[9]][(L[0], int(L[1]), int(L[2]))] = {'gene_id': L[16]  , 'codon': (int(L[14]), int(L[15])), 'strand': L[5], 'coverage': float(L[4]), 'exon_type' : L[3].split('-')}
    I.close()
    return(utrs)

def read_names_file(Infile):
    I = open(Infile)
    Names = {}
    while True:
	Line = I.readline()
	if not Line:break
	Names[Line.replace('\n', '').split('\t')[1]] = Line.split('\t')[0] 
    I.close()
    return(Names)


def find_left(gene): 
    left_end = 0
    discarded = []
    considered = []
    gene_name = ''
    coord = []
    utr_start = 0
    codon = False
    for lola in sorted(gene):
	if not codon and gene[lola]['codon'] == (-1,-1) and gene[lola]['coverage'] >= 10 :
	    coord += [lola]
	elif not gene[lola]['codon'] == (-1,-1):
	    codon = True
	    gene_name = gene[lola]['gene_id']
	    utr_start = gene[lola]['codon'][0]
	    if gene[lola]['coverage'] >= 10: 
		coord += [lola]
    if len(coord) > 0:
	for i, lola in enumerate(coord):
	    for j, forrest in enumerate(coord):
		if i < j:
		    if lola[1] == forrest[1] or lola[2] == forrest[2]:
			if gene[lola]['coverage'] >= gene[forrest]['coverage']:
			    discarded += [forrest]
			else:
			    discarded += [lola]
	    if not lola in discarded:
		considered += [lola]
	if not coord[-1] in discarded:
	    considered += [lola]
	considered = set(considered)
	for lola in sorted(considered):
	    left_end += lola[2]- lola[1]
	if utr_start > 0:
	    left_end -= sorted(considered)[-1][2] - utr_start
    else:
	g = gene.keys()
	for lola in g:
	    if not gene[lola]['coverage'] >= 10:
		del gene[lola]
	if len(gene) > 0:
	    left_end = sorted(gene)[0][2]- sorted(gene)[0][1]
	    considered = (sorted(gene)[0])
	else:
	    left_end = 0
	    considered = ('NA', 'NA', 'NA')
	print('no stop codon')
    return(left_end, gene_name, considered, utr_start)



def find_right(gene): 
    right_end = 0
    discarded = []
    considered = []
    gene_name = ''
    coord = []
    utr_start = 0
    codon = True
    for lola in sorted(gene):
	if not codon and gene[lola]['codon'] == (-1,-1) and gene[lola]['coverage'] >= 10 :
	    coord += [lola]
	elif not gene[lola]['codon'] == (-1,-1):
	    codon = False
	    gene_name = gene[lola]['gene_id']
	    utr_start = gene[lola]['codon'][1]
	    if gene[lola]['coverage'] >= 10: 
		coord += [lola]
    if len(coord) > 0:
	for i, lola in enumerate(coord):
	    for j, forrest in enumerate(coord):
		if i < j:
		    if lola[1] == forrest[1] or lola[2] == forrest[2]:
			if gene[lola]['coverage'] >= gene[forrest]['coverage']:
			    discarded += [forrest]
			else:
			    discarded += [lola]
	    if not lola in discarded:
		considered += [lola]
	if not coord[-1] in discarded:
	    considered += [lola]
	considered = set(considered)
	for lola in sorted(considered):
	    right_end += lola[2]- lola[1]
	if utr_start > 0:
	    right_end -= utr_start - sorted(considered)[0][1]
    else:
	g = gene.keys()
	for lola in g:
	    if not gene[lola]['coverage'] >= 10:
		del gene[lola]
	if len(gene) > 0:
	    right_end = sorted(gene)[-1][2]- sorted(gene)[-1][1]
	    considered = (sorted(gene)[-1])
	else:
	    left_end = 0
	    considered = ('NA', 'NA', 'NA')
	print('no_stop_codon')
    return(right_end, gene_name, considered, utr_start)



def define_utr(exons, end):
    utrs = {}
    for lola in exons:
	utrs[lola] = {}
	for forrest in exons[lola]:
	    if len(exons[lola][forrest]) > 0:
		utr = 'NA'
		strand = exons[lola][forrest][exons[lola][forrest].keys()[0]]['strand']
		if strand == '+' and end == 'start':
		    length, gene_name, coord, utr_start = find_left(exons[lola][forrest])
		    utr = '5p'
		elif strand == '+' and end == 'stop':
		    length, gene_name, coord, utr_start = find_right(exons[lola][forrest])
		    utr = '3p'
		elif strand == '-' and end == 'start':
		    utr = '5p'
		    length, gene_name, coord, utr_start = find_right(exons[lola][forrest])
		elif strand == '-' and end == 'stop':
		    utr = '3p'
		    length, gene_name, coord, utr_start = find_left(exons[lola][forrest])
		else:
		    print('Fatal error, either no end or no strand')
		    length = 0
		    gene_name = 'NA'
		    coord = set([('NA', 0, 0)])
		gene_name = forrest
		utrs[lola][forrest] = {'utr': utr ,'length':length, 'gene_name' : gene_name,'strand': strand, 'coord': coord, 'utr_start' : utr_start}
    return(utrs)

def write_output_table(Outfile, utrs, ref):
    O = open(Outfile, 'w')
    O.write('gene_id\tgene_name\tstrand\tutr\tlength\tutr_start\tcoordinates\n')
    for lola in utrs:
	for forrest in utrs[lola]:
	    if forrest in ref:
		O.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(utrs[lola][forrest]['gene_name'], ref[forrest], utrs[lola][forrest]['strand'], utrs[lola][forrest]['utr'], utrs[lola][forrest]['length'],utrs[lola][forrest]['utr_start'] ,list(utrs[lola][forrest]['coord'])))
    O.close()
    return


# run script
STOP = read_intersected_bed(stop)
START = read_intersected_bed(start)
REF = read_names_file(ref_names)

utr3 = define_utr(STOP, 'stop')
utr5 = define_utr(START, 'start')

write_output_table('control.3utr.txt', utr3, REF)
write_output_table('control.5utr.txt', utr5, REF)

