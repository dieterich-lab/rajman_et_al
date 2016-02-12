# Rscript to run EISA following online examples
# load edgeR library
library(edgeR)
library(locfit)

# functions
# collapse multiple transcript entries
collapse_transcripts <-function(D){
  genes = table(D$id)
  collapsed = aggregate(D[,5] ~ id, data = D, FUN = sum)
  for(i in 6:length(D)){
    temp <- aggregate(D[,i] ~ id, data = D, FUN = sum)
    collapsed = merge(collapsed, temp, by = 'id')  
  }
  names(collapsed) <- names(D[,4:length(D)])
  return(collapsed)
}

#call edgeR
call_edgeR<-function(counts,group){

  y <- DGEList(counts=counts,group=group)
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  et <- exactTest(y, pair = c('WT', 'PTX'))
  topTags(et)

  return(et)
}

# files
gene_lists = c('regulation_high_confidence.csv', 'silac_high_confidence.csv')
titles = c('DEG (high confidence)', 'SILAC (high_confidence)')
names = c('differentially_expressed_transcripts_high_confidence.txt', 'differentially_expressed_proteins_high_confidence.txt')
indeces = c(4, 7)
conditions = c('Vehicle', 'Vehicle', 'Vehicle', 'PTX', 'PTX', 'PTX')

for(g in 1:length(gene_lists)){
  # read in tables
  introns = read.table('introns.txt', as.is = T, header = T) # generated using bedtools coverage
  exons = read.table('exons.txt', as.is = T, header = T) # generated using bedtools coverage
  gene_list = read.table(gene_lists[g], header = T, as.is = T)

  # merge introns / exons from one gene
  E = collapse_transcripts(exons)
  I = collapse_transcripts(introns)

  # merge and separate Exons and Introns to have the same sorting
  A = merge(E, I, by = 'id', all = T)
  A[is.na(A)] <- 0
  # split data framge again
  Ex = A[,2:length(E)]
  In = A[,(length(E)+1):length(A)]

  # add row names again
  row.names(Ex) <- A$id
  row.names(In) <- A$id

  # normalize like in paper
  E.norm <- data.frame(t(mean(colSums(Ex))*t(Ex)/colSums(Ex))) 
  I.norm <- data.frame(t(mean(colSums(In))*t(In)/colSums(In)))

  # select genes which have enough coverage and are in a given gene list
  genes.sel <- rownames(Ex) %in% gene_list$ID & rowMeans(log2(E.norm+8))>=5 & rowMeans(log2(I.norm+8))>=5 # or 

  # for exons/introns individually
  Exons = subset(Ex, rownames(Ex) %in% names(subset(genes.sel, genes.sel)))
  Introns = In[genes.sel,]

  I.edgeR = call_edgeR(Introns, conditions)
  E.edgeR = call_edgeR(Exons, conditions)

  I.fc = data.frame(ID = row.names(I.edgeR$table), FC = I.edgeR$table$logFC)
  E.fc = data.frame(ID = row.names(E.edgeR$table), FC = E.edgeR$table$logFC)

  output_table = data.frame(ID = row.names(E.edgeR$table), E.edgeR$table, I.edgeR$table)
  output_table = merge(gene_list, output_table, by = 'ID')
  write.table(output_table, names[g], sep = '\t', quote = F)

}
