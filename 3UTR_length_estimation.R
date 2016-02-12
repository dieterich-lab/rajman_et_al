# utr length analysis
utr3_control = read.table('control.3utr.txt', as.is = T, header = T)
utr5_control = read.table('control.5utr.txt', as.is = T, header = T)
utr3_ptx = read.table('ptx.3utr.txt', as.is = T, header = T)
utr5_ptx = read.table('ptx.5utr.txt', as.is = T, header = T)

ptx = merge(utr5_ptx, utr3_ptx, by = 'gene_id', all = T)
utr3 =  merge(utr3_control, utr3_ptx, by = 'gene_id', all = T)
utr5 =  merge(utr5_control, utr5_ptx, by = 'gene_id', all = T)
control =  merge(utr5_control, utr3_control, by = 'gene_id', all = T)

All = merge(control, ptx, by = 'gene_id', all = T)

X = subset(All, gene_id %in% names(subset(table(All$gene_id), table(All$gene_id)>1)))

All_easy = subset(All, !gene_id %in% X$gene_id)


write.table('All_genes_with_one_start_stop_codon.txt', sep = '\t', quote = F, row.names = F)

gene_names = All$gene_name.x.x
gene_names[is.na(gene_names)] <- All$gene_name.y.x[is.na(gene_names)]
gene_names[is.na(gene_names)] <- All$gene_name.y.y[is.na(gene_names)]
gene_names[is.na(gene_names)] <- All$gene_name.x.y[is.na(gene_names)]


gene_strand = All$strand.x.x
gene_strand[is.na(gene_strand)] <- All$strand.y.x[is.na(gene_strand)]
gene_strand[is.na(gene_strand)] <- All$strand.y.y[is.na(gene_strand)]
gene_strand[is.na(gene_strand)] <- All$strand.x.y[is.na(gene_strand)]

start_codon = All$utr_start.x.x
start_codon[is.na(start_codon)] <- All$utr_start.x.y[is.na(start_codon)]

stop_codon = All$utr_start.y.x
stop_codon[is.na(stop_codon)]<- All$utr_start.y.y[is.na(stop_codon)]

XX <-data.frame(gene_id = All$gene_id, gene_name = gene_names, strand = gene_strand, start_codon = start_codon, stop_codon = stop_codon, length.control.5p =All$length.x.x , length.control.3p= All$length.y.x, length.ptx.5p=All$length.x.y , length.ptx.3p = All$length.y.y)

XX <- XX[!duplicated(XX),]

write.table(XX, 'gene_info.txt', sep = '\t', quote = F, row.names = F)

# manually deleted duplicated columns

utr <- read.table('gene_info.txt', as.is = T, header = T)
regulation <-read.table('regulation_high_confidence.csv', as.is = T, header = T)
silac <-read.table('silac_high_confidence.csv', as.is = T, header = T)

diff_utr5 <- utr$length.control.5p - utr$length.ptx.5p
diff_utr3 <- utr$length.control.3p - utr$length.ptx.3p

utr <- cbind(utr, diff_utr5, diff_utr3)

regulation_names = toupper(regulation$RN5.gene.name)
utr_names = toupper(utr$gene_name)
silac_names = toupper(silac$RN5.gene.name)

utr = cbind(utr, utr_names)
regulation = cbind(regulation, regulation_names)
silac = cbind(silac, silac_names)

silac_utr <- merge(silac, utr, by.x = 'silac_names', by.y = 'utr_names', all.x = T)
regulation_utr <- merge(regulation, utr, by.x = 'regulation_names', by.y = 'utr_names', all.x = T)

write.table(silac_utr, 'silac_length.txt', quote = F, sep = '\t', row.names = F)
write.table(regulation_utr, 'regulation_length.txt', quote = F, sep = '\t', row.names = F)


# proceed with manual checking, if ok
need_to_check_these_5p = subset(utr, abs(diff_utr5) > 20)
need_to_check_these_3p = subset(utr, abs(diff_utr3) > 20)

write.table(need_to_check_these_3p, paste(main_folder, structure[4], 'needs_validation_3p', sep = '/'), sep = '\t', quote = F, row.names = F)
write.table(need_to_check_these_5p, paste(main_folder, structure[4], 'needs_validation_5p', sep = '/'), sep = '\t', quote = F, row.names = F)



