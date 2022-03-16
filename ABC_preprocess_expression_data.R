library(biomaRt)
library(dplyr)

RNA = read.table("/home/darragh/ABC/ENCFF292FVY_HMEC.tsv", sep = "\t", header = TRUE)
RNA = RNA[650:nrow(RNA), c("gene_id", "TPM")]
RNA$gene_id = gsub("\\..*","",RNA$gene_id)
head(RNA)

#convert ensemble ids to gene ids using biomart
mart = useMart('ensembl')
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
gene <- getBM( attributes = c("ensembl_gene_id","external_gene_name"),values = RNA$gene_id, mart = ensembl)

#match RNA ids to biomart gene ids
id <- match(RNA$gene_id , gene$ensembl_gene_id)
#Create a new symbol column with the biomart gene ids
RNA$Symbol <- gene$external_gene_name[id]

RNA = RNA %>% dplyr::select(Symbol, everything() )

#I only want to keep the entries wit valid gene names.
RNA = RNA[!(is.na(RNA$Symbol) | RNA$Symbol==""), ]

RNA = RNA[, c("Symbol", "TPM")]

write.table(RNA, "/home/darragh/ABC/HMEC_RNA_abc.bed",row.names = F, col.names =F, sep = "\t", quote = FALSE, )


Genes = read.table("/home/darragh/hg19_protein_coding_minimal.gtf", sep = "\t", header = F)
Genes$V5 = gsub("\\;","",Genes$V5)
Genes$V6 = 0
Genes = Genes[, c(1,2,3,5,6,4)]
write.table(Genes, "/home/darragh/ABC/hg19_gene_bounds_abc.bed",row.names = F, col.names =F, sep = "\t", quote = FALSE, )
