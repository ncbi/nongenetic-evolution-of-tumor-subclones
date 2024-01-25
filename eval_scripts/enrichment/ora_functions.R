library(tidyr)
library(dplyr)
library(stringr)
library(clusterProfiler)

##########
### Helper functions to run over-representation enrichment analysis.
##########

kegg_enrich <- function(query_genes, ref_genes, save_file) {
    query_entrez <- bitr(query_genes$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    ref_entrez <- bitr(ref_genes$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    ego_kegg <- enrichKEGG(
                        gene = query_entrez$ENTREZID,
                        universe = ref_entrez$ENTREZID,
                        keyType = "ncbi-geneid",
                        organism = "mmu",
                        pvalueCutoff = .05,
                        pAdjustMethod = "fdr",
                        qvalueCutoff = .05
    )
    ego_kegg_df <- data.frame(ego_kegg)
    if (nrow(ego_kegg_df) != 0) {
        for (i in 1:nrow(ego_kegg_df)){
            entrez_ids <- str_split(ego_kegg_df[i, "geneID"], "/")
            gene_names <- bitr(entrez_ids[[1]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
            ego_kegg_df[i, paste("geneName", sep="")] <- paste(unique(gene_names$SYMBOL), collapse="/")
        }
    }
    write.csv(ego_kegg_df, save_file, row.names = FALSE)
}

go_enrich <- function(query_genes, ref_genes, save_file) {
    go_all <- enrichGO(
                    gene = query_genes$gene_name,
                    universe = ref_genes$gene_name,
                    OrgDb = "org.Mm.eg.db",
                    keyType = "SYMBOL",
                    ont = "ALL",
                    pvalueCutoff = .05,
                    pAdjustMethod = "fdr",
                    qvalueCutoff = .05,
                    readable = TRUE
                )
    write.csv(go_all, save_file, row.names = FALSE)
}