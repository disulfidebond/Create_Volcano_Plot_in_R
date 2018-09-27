# start with data frame
# if using data from DESEQ2, copy data from DESEQ2 object to data frame:
df_results <- as.data.frame(results)
# or
rlog_data_matrix <- assay(rlog_transformed_data)
df_rlog_data_matrix <- as.data.frame(rlog_data_matrix)

# if not using a Biological data set with genes, or row names will not be used as part of heatmap, skip this step
df_names <- rownames(df_results)
df_results$genes <- df_names

# the following steps assume the dataframe will have the format:
# colN,colN,padj,log2FoldChange,colN,colN,...,genes,log10,colorflags,expression
# where colN is any column that will not be used but may be included
# genes is a list of gene names, log10 is a list of log10 transformed p values, colorflags is used to create expression,
# and expression is used to create the gradient and labels for the gradient
df_results$log10 <- -log10(df_results$padj)
# IMPORTANT: this will set values outside a filter of log10(pval) > 1.309 to NA, and they will NOT be included in the chart.
# you may want to comment out this next line
df_results$GeneNames <- ifelse(df_results$log10 > 1.309, df_results$genes, NA)
df_results$colorflags <- ifelse(df_results$padj < 0.05, ifelse(df_results$log2FoldChange > 1, 1, ifelse(df_results$log2FoldChange < -1, -1, 0)), 0)
df_results$expression <- ifelse(df_results$colorflags == -1, "down", ifelse(df_results$colorflags == 0, "notsignificant", ifelse(df_results$colorflags == 1, "up", 0)))

# OPTIONAL: decrease the number of results
# df_results <- df_results[1:1216,1:11]

p <- ggplot(df_results, aes(x=log2FoldChange, y=log10, colour=expression)) + geom_point(shape=15, size=1)
# aes: set x-axis to log2foldChange, y-axis to log10; you may log-transform the values here with y=log10(padj)
# shape: change the shape of the point, size: change the size
p <- p + geom_vline(xintercept = -1) + geom_vline(xintercept = 1) + geom_hline(yintercept = 1.3) + xlab("Log2 Fold Change") + scale_x_continuous(limits=c(-10,10)) + ylab("-log10 p-value") # + geom_text(aes(y=log10, label=heatmap_df_trtrg_trt$GeneNames), size=2, hjust=0)
# vline: create a vertical line, hline: create a horizontal line, xlab/ylab: label axes
# scale_x_continuous: scale x axis, use with caution, usually has unexpected results
# geom_text() set the text for the labels with automatic labeling, use y or x axis as a basis, and adjust by y=log10 + N, size == size of text
# (hjust=1) == right justify ; (hjust=0) == left justify ; (vjust=0) == set v to baseline of text ; (vjust=1) == set v to top of text ; use caution if adjusting x/y values and justify together

# alternatively, labels can be labeled individually, basically by selecting a point in the Cartesian plane that matches to the point on the graph, then adding a label there.
# adjust the label by adding/subtracting x and y values
p + annotate("text", x=2.4, y=1.47, label="Label1", size=2, hjust=0) + annotate("text", x=2.50, y=1.56, label="Label2", size=2, hjust=0)
# note that using annotate() and geom_text() is not a good idea, as one will overwrite the other
# running this command or typing "p" will display the graph in RStudio
# to print the graph, use 
ggsave("myimage010101.png", width = 16, height = 16, units = "cm")
# NOTE THE FOLLOWING:
# # this will save the last graph created
# # it will overwrite any existing file without asking
# # it will autodetect the file type from the name provided, i.e. "myimage010101.png" will be a png image, but PNG, PDF, and SVG files are all possible
# # it will NOT save multi-page plots
