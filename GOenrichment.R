library(biomaRt)
library(CellPlot)
library(plyr)
library(multtest)
library(annotate)
library(topGO)

# define functions
cuff2topGO <- function (d, g, id = "external_gene_name")
{
   # input:
   # d: cuffdiff table, for single pair, as returned by read.cuffdiff
   # g: biomaRt:::getBM output
   #    attr:
   # "external_gene_name","go_id","name_1006","definition_1006","namespace_1003"
   #    mart: useDataset("celegans_gene_ensembl", mart = useMart("ensembl"))
   # id: the name of the attribute in g which has the same ids as in
   # 'gene' column in d
   # create a list of character vectors with goIDs for each gene
   # list names are gene ids, attribute 'stat' are the q values
   g <- dlply(g, id, function(x)x$go_id)
   j <- sapply(g, length) < 1 & sapply(g, function(x) {nchar(x[1])}) < 1
   g <- g[!j]
   i <- sort(intersect(d$gene, names(g)))
   i <- setdiff(i, "-")
   g <- g[i]
   q <- merge.data.frame(data.frame(gene=i,stringsAsFactors=F),
                         d[,c("gene","p.adj")],
                         sort=F, all.x=T, all.y=F, by.x="gene",by.y="gene")
   idup <- q$gene[duplicated(q$gene)]
   qdup <- subset(q, gene %in% idup)
   gdup <- g[names(g)%in%idup]
   g <- g[!names(g)%in%idup]
   q <- subset(q, !gene %in% idup)
   attr(g, "q") <- q
   attr(g, "gdup") <- gdup
   attr(g, "qdup") <- qdup
   attr(g, "i") <- i
   attr(g, "idup") <- idup
   return(g)
}

create_artificial_padj <- function(D, genelist){
  # for genes in gene_list, generate an artificial 'padj' value which is 'significant' and try as normal function
  gene.select = rep(1, length(D[,1]))
  for(g in 1:length(D[,1])){
  gene.select[g] <- D$gene[g] %in% genelist}
  D <- cbind(D, gene.select)
  return(D)
  }

# run script
all_gene_lists = c('regulation.txt', 'silac.txt')
titles = c('Differentially expressed genes', 'Differentially expressed proteins')

mart = useMart(biomart = "ensembl", dataset = 'rnorvegicus_gene_ensembl')
go_gene_map = getBM(attributes = c("external_gene_name","go_id"), mart = mart)
go_gene_map$external_gene_name <- toupper(go_gene_map$external_gene_name)


for(i in 1:length(all_gene_lists)){
  tiff(gsub('.txt', '.tiff', all_gene_lists[i]))
  gene_list = read.table(all_gene_lists[i], as.is = T, header = F)
  gene_list$V1 <- toupper(gene_list$V1)

  d = read.table('gene_exp.diff', as.is = T, header = T)
  d$gene <- toupper(d$gene)
  d = create_artificial_padj(d, gene_list$V1)

  goterms <- cuff2topGO(d, go_gene_map)

  godata <- new(
    "topGOdata", ontology = "BP", description = 'golub',
    allGenes = setNames(d$gene.select, d$gene), nodeSize = 25,
    geneSelectionFun = function (allScore) { allScore },
    annotationFun = annFUN.gene2GO, gene2GO = goterms)

  golubstat <- go.enrich(godata, d)

  x <- subset(golubstat, p<=.05 & significant>4 & !duplicated(genes))
  x <- head(x, 20)
  sym.plot( ticksize = 5, x = setNames(x$loge, x$term), cells = x$deg.log2fc, 
          x.annotated = x$annotated, y.mar = c(0.1,0), x.mar = c(.47, 0), key.n=7, 
          key.lab = all_gene_lists[i], cex = 1.6, axis.cex=0.8, group.cex=0.7, 
          grid.lwd=3, main = titles[i])

  dev.off()
}