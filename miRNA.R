# script to analyse Targetscan reformatted tables with R
table = c('targetscan_out.summarized.txt')

bio_tables = c('silac.high_confidence.csv', 'regulation.high_confidence.csv')

main_sample = 'control.3p'
main_mirna = 'FDR < 0.05'
main_bio_tables = c('DEP (high confidence)', 'DEG (high confidence)')



for(j in 1:length(bio_tables)){
    print(j)
    miRNA = read.table(table, as.is = T, header = T)
    trans = read.table(bio_tables[j], as.is = T, header = T)

    all = merge(trans, miRNA, by = 'gene_id')

    conservation = rep(NA, length(all[,1]))
    regulation = rep(NA, length(all[,1]))
    print('ok')
    regulation <- with(all,ifelse(all[,4] < 0 ,1 ,0)) # 1 is down regulation and 0 is up regulation
    conservation <- with(all,ifelse(num_species <= 3 , 1 ,conservation))
    conservation <- with(all,ifelse((num_species > 3  & num_species < 10), 2 ,conservation))
    conservation <- with(all,ifelse(num_species >= 10 , 3 ,conservation))
    print('ok')
    all = cbind(all, regulation, conservation)
    
    # generate summary table rows = genes, columns = miRNAs
    miRNA_names <- names(table(all$miRNA_name))
    genes <- unique(all$gene_name)
    con1 <- matrix(rep(NA, (length(miRNA_names) * length(genes))), ncol = length(miRNA_names))
    con2 <-  matrix(rep(NA, (length(miRNA_names) *  length(genes))), ncol = length(miRNA_names))
    con3 <-  matrix(rep(NA, (length(miRNA_names) *  length(genes))), ncol = length(miRNA_names))
    total <- matrix(rep(NA, (length(miRNA_names) * length(genes))), ncol = length(miRNA_names))
    
    for(x in 1:length(miRNA_names)){
      for(y in 1:length(genes)){
	con1[y,x] = length(subset(all, miRNA_name == miRNA_names[x] & gene_name == genes[y] & conservation == 1)[,1])
	con2[y,x] = length(subset(all, miRNA_name == miRNA_names[x] & gene_name == genes[y] & conservation == 2)[,1])
	con3[y,x] = length(subset(all, miRNA_name == miRNA_names[x] & gene_name == genes[y] & conservation == 3)[,1])
	total[y,x] = length(subset(all, miRNA_name == miRNA_names[x] & gene_name == genes[y])[,1])
      }
    }
    
    con1 <- as.data.frame(con1)
    con2 <- as.data.frame(con2)
    con3 <- as.data.frame(con3)
    total<- as.data.frame(total)
    
    row.names(con1) <- genes
    names(con1) <- miRNA_names
    row.names(con2) <- genes
    names(con2) <- miRNA_names
    row.names(con3) <- genes
    names(con3) <- miRNA_names
    row.names(total) <- genes
    names(total) <- miRNA_names

    write.table(con1, paste('miRNA_summary_table', main_sample, main_mirna, main_bio_tables[j], 'poorly_conserved.txt', sep = '.'), sep = '\t', quote = F, row.names = T)
    write.table(con2, paste('miRNA_summary_table', main_sample, main_mirna, main_bio_tables[j], 'conserved.txt', sep = '.'), sep = '\t', quote = F, row.names = T)
    write.table(con3, paste('miRNA_summary_table', main_sample, main_mirna, main_bio_tables[j], 'highly_conserved.txt', sep = '.'), sep = '\t', quote = F, row.names = T)
    write.table(total, paste('miRNA_summary_table', main_sample, main_mirna, main_bio_tables[j], 'txt', sep = '.'), sep = '\t', quote = F, row.names = T)

    
    print(head(all))
    conservation_matrix = matrix(rep(NA, 6), ncol = 3)
      for(x in 1:3){
	conservation_matrix[1,x] = length(subset(all, conservation == x & regulation == 0)[,1])
	conservation_matrix[2,x] = length(subset(all, conservation == x & regulation == 1)[,1])
      }
    up = subset(all, regulation == 0)
    down = subset(all, regulation == 1)

    miRNAs_up = names(table(up$miRNA_name))
    miRNAs_down = names(table(down$miRNA_name))
 
    # filling up the up matrix
    UP = matrix(rep(NA, length(miRNAs_up) * 3), nrow = 3)
    if(!is.null(miRNAs_up)){
      for(x in 1:length(miRNAs_up)){
	UP[1,x] = length(subset(up, miRNA_name == miRNAs_up[x] & conservation == 1)[,1])
	UP[2,x] = length(subset(up, miRNA_name == miRNAs_up[x] & conservation == 2)[,1])
	UP[3,x] = length(subset(up, miRNA_name == miRNAs_up[x] & conservation == 3)[,1])
      }
    } 
    # filling up the down matrix
    DOWN = matrix(rep(NA, length(miRNAs_down) * 3), nrow = 3)
    if(!is.null(miRNAs_down)){
      for(x in 1:length(miRNAs_down)){
	DOWN[1,x] = length(subset(down, miRNA_name == miRNAs_down[x] & conservation == 1)[,1])
	DOWN[2,x] = length(subset(down, miRNA_name == miRNAs_down[x] & conservation == 2)[,1])
	DOWN[3,x] = length(subset(down, miRNA_name == miRNAs_down[x] & conservation == 3)[,1])
      }
    } 

    pdf(paste('miRNA_seeds', main_sample, main_mirna, main_bio_tables[j], 'pdf', sep = '.'))
      par(mar = c(5,5,10,5))
      bp = barplot(conservation_matrix, beside = T, main = paste( main_bio_tables[j], '\nNumber of miRNA seeds\nfrom', main_mirna, 'in', main_sample), ylab = 'number of seeds', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2, col = c('firebrick2', 'dodgerblue2') )
      mtext( at = colMeans(bp), c('poorly\nconserved', '\nconserved', 'highly\nconserved'), cex = 1.7, side = 1, line = 3)
      legend('topright', c('up-regulated', 'down-regulated'), fill = c('firebrick2', 'dodgerblue2'), bty = 'n', cex = 1.7)
      if(!is.null(miRNAs_up)){
	par(mar = c(12,5,5,1))
	barplot(UP, names = miRNAs_up, las = 2, ylim = c(0, max(colSums(UP)) + 5), main = 'miRNA seeds in up regulated genes', ylab = 'number of seeds', cex.axis = 1.7, cex.lab = 1.7, cex.names = 1.3, col = c('dodgerblue3', 'firebrick2', 'gold1'))
	legend('topleft', c('highly conserved', 'conserved', 'poorly conserved'), fill = c('gold1', 'firebrick2', 'dodgerblue3'), bty = 'n', cex = 1.3)
      }
      if(!is.null(miRNAs_down)){
	par(mar = c(12,5,5,1))
	barplot(DOWN, names = miRNAs_down, las = 2, ylim = c(0, max(colSums(DOWN)) + 5), main = 'miRNA seeds in down regualted genes', ylab = 'number of seeds', cex.axis = 1.7, cex.lab = 1.7, cex.names = 1.3, col = c('dodgerblue3', 'firebrick2', 'gold1'))
	legend('topleft', c('highly conserved', 'conserved', 'poorly conserved'), fill = c('gold1', 'firebrick2', 'dodgerblue3'), bty = 'n', cex = 1.3)
      }
    dev.off()
    print('end')
  write.table(all, paste('miRNA_seeds_table', main_sample, main_mirna, main_bio_tables[j], 'txt', sep = '.'), sep = '\t', quote = F, row.names = F)
  }

