kegg.res <- fread("KEGG_plot_top15.xls",sep='\t',header=T)


geneID<-kegg.res$geneID
geneID_all <- unlist(apply(as.matrix(geneID),1,function(x) unlist(strsplit(x,'/'))))
x <- new("enrichResult",
         result         = kegg.res,
         pvalueCutoff   = 1.0,
         pAdjustMethod  = "BH",
         qvalueCutoff   = 1.0,
         gene           = geneID_all,
         universe       = 'Unknown',
         geneSets       = list(),
         organism       = "hsa",
         keytype        = "Unknown",
         readable       = FALSE
)

dotplot(x,showCategory = 15)
ggsave('C_vs_Q_DGEs.Upregulated_KEGG_selected_top15_dotplot.pdf',width = 8,height = 6)
