import os

# define functions

def read_starbase_bed(infile):
    i = open(infile, 'r')
    rbp = {}
    while True:
        Line = i.readline()
        if not Line:break
        L = Line.split('\t')
	chromosome = L[0]
	start = int(L[1])
	end = int(L[2])
	score = L[4]
	strand = L[5].replace('\n','')
	cluster = L[3]
	if chromosome in rbp:
	    rbp[chromosome][(start, end)] = [strand, score, cluster] 
	else:
	    rbp[chromosome] = {(start,end):[strand, score, cluster]}
    i.close()
    print('reading starbase bed done')
    return(rbp)

def read_bed_regions(file_name):
    F = open(file_name, 'r')
    Dict = {}
    counter = 0
    while True:
	Line = F.readline()
	if not Line:break
	else:
	    counter += 1
	    if Line.startswith('chr'):
		L = Line.replace('\n', '').split('\t') 
		if L[0] in Dict:
		    if 'NM' in L[3]:
		        Dict[L[0]][(int(L[1]), int(L[2]))]=['NM_%s' %(L[3].split('_')[1]), L[3].split('_')[2]]
		    else:
		        Dict[L[0]][(int(L[1]), int(L[2]))]=['NR_%s' %(L[3].split('_')[1]),  L[3].split('_')[2]]
		else:
		    if 'NM' in L[3]:
		        Dict[L[0]]= {(int(L[1]), int(L[2])):['NM_%s' %(L[3].split('_')[1]), L[3].split('_')[2]]}
		    else:
		        Dict[L[0]]= {(int(L[1]), int(L[2])):['NR_%s' %(L[3].split('_')[1]),  L[3].split('_')[2]]}
		    
    print('generated dictionary')
    print(counter)
    F.close()
    return(Dict)


def read_bed_regions_with_length(file_name):
    F = open(file_name, 'r')
    Dict = {}
    counter = 0
    while True:
	Line = F.readline()
	if not Line:break
	else:
	    counter += 1
	    if Line.startswith('chr'):
		L = Line.replace('\n', '').split('\t') 
		if 'NM' in L[3]:
		    if 'NM_%s' %(L[3].split('_')[1]) in Dict:
			Dict['NM_%s' %(L[3].split('_')[1])][L[3].split('_')[2]] += [int(L[2])-int(L[1])]
		    else:
		        Dict['NM_%s' %(L[3].split('_')[1])] = {'utr3':[], 'utr5':[], 'cds':[]}
		        Dict['NM_%s' %(L[3].split('_')[1])][L[3].split('_')[2]] += [int(L[2])-int(L[1])]
		else:
		    if 'NR_%s' %(L[3].split('_')[1]) in Dict:
			Dict['NR_%s' %(L[3].split('_')[1])][L[3].split('_')[2]] += [int(L[2])-int(L[1])]
		    else:
		        Dict['NR_%s' %(L[3].split('_')[1])] = {'utr3':[], 'utr5':[], 'cds':[]}
		        Dict['NR_%s' %(L[3].split('_')[1])][L[3].split('_')[2]] += [int(L[2])-int(L[1])]
    print('generated dictionary')
    print(counter)
    F.close()
    return(Dict)



def read_bed_names(file_name):
    F = open(file_name, 'r')
    Dict = {}
    counter = 0
    while True:
	Line = F.readline()
	if not Line:break
	else:
	    counter += 1
	    if Line.startswith('N'):
		L = Line.replace('\n', '').split('\t')
		Dict[L[0]] = L[5]
    print('generated dictionary')
    print(counter)
    F.close()
    return(Dict)


def read_silac(file_name):
    i = open(file_name)
    silac_genes = {}
    i.readline()
    while True:
        Line = i.readline()
        if not Line:break
        L = Line.split('\t')
        silac_genes[L[2]] = [L[0], L[1]]
    print('Done reading SILAC targets')
    return(silac_genes)


def annotate(bed_regions, bed_names, RBP):
    Annotated = {}
    chroms = set(RBP.keys()) & set(bed_regions.keys())
    for lola in chroms:
        Annotated[lola]={}
    for lola in chroms:
        for region in RBP[lola]:
            for forrest in bed_regions[lola]:
	        S1, E1 = forrest
	        S2, E2 = region
	        # Case 1
	        if not (S2, E2) in Annotated[lola]:
	            if S1 <= S2 and E1 >= E2:
		        Annotated[lola][(S2, E2)] =  [bed_names[bed_regions[lola][forrest][0]], bed_regions[lola][forrest][1]]
	            # Case 2
	            elif S1 >= S2 and S1 <= E2:
		        Annotated[lola][(S2, E2)] = [bed_names[bed_regions[lola][forrest][0]], bed_regions[lola][forrest][1]]
	            # Case 3
	            elif S2 <= E1 and E1 <= E2:
		        Annotated[lola][(S2, E2)] = [bed_names[bed_regions[lola][forrest][0]], bed_regions[lola][forrest][1],]
	            # Case 4
	            elif S1 >= S2 and E1 <= E2: #
		        Annotated[lola][(S2, E2)] = [bed_names[bed_regions[lola][forrest][0]], bed_regions[lola][forrest][1],]
	print(lola)
    print('Done annotating')
    return(Annotated)


def write_output(outfile, RBP, RBP_annotated, silac):
    o = open(outfile, 'w')
    o.write('chr\tstart\tend\tcluster\tscore\tstrand\tgene\tregion\tsilac_id\tsilac_FC\n')
    KK = RBP.keys()
    KK.sort()
    for lola in KK:
        Keys = RBP[lola].keys()
        Keys.sort()
        for forrest in Keys:
	    o.write('%s\t%s\t%s\t%s\t%s\t%s\t' %(lola, forrest[0], forrest[1], RBP[lola][forrest][2], RBP[lola][forrest][1], RBP[lola][forrest][0]))
	    if lola in RBP_annotated:
		if forrest in RBP_annotated[lola]:
		    o.write('%s\t%s\t' %(RBP_annotated[lola][forrest][0], RBP_annotated[lola][forrest][1]))
		    if RBP_annotated[lola][forrest][0] in silac.keys():
			o.write('%s\t%s\n' %(silac[RBP_annotated[lola][forrest][0]][0], silac[RBP_annotated[lola][forrest][0]][1]))
		    else:
			o.write('NA\tNA\n')
		else:
		    o.write('NA\tNA\tNA\tNA\n')
	    else:
	        o.write('NA\tNA\tNA\tNA\n')
    o.close()
    print('writing output done')
    return

def write_region_length(bed_regions_length, bed_names, silac, outputfile):
    o = open(outputfile, 'w')
    o.write('gene\tcds\tutr3\tutr5\n')
    for transcript in bed_names:
	if bed_names[transcript] in silac:
	    o.write('%s\t%s\t%s\t%s\n' %(bed_names[transcript], sum(bed_regions_length[transcript]['cds']), sum(bed_regions_length[transcript]['utr3']),sum(bed_regions_length[transcript]['utr5'])))  
    o.close()
    return




# run script

# define paths
RBP_path = '<PATH>/CLIP_targetlists'
out_path = '<PATH>/analysis'
RBPs = os.listdir(RBP_path)

# read input files
bed_regions = read_bed_regions('mm9.bed')
bed_names = read_bed_names('mm9.genenames.bed')
silac = read_silac('silac.txt')

# important for down stream analsysi to adjust number of motifs by gene length
bed_regions_length = read_bed_regions_with_length('mm9.bed')
write_region_length(bed_regions_length, bed_names, silac, 'silac_region_length.txt')

# iterate over public datasets
for protein in RBPs:
    targets = read_starbase_bed('%s/%s'%(RBP_path,protein))
    annotated_targets = annotate(bed_regions,bed_names,targets)
    write_output('%s/%s_annotated_.txt' %(out_path, protein.replace('.bed','')), targets,annotated_targets, silac)



  