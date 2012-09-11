
#chrom =
#cell_type =
#MBstart = 
#MBend = 
#out_folder =
#region.ID = 
put_loc = '/home/jkb4y/ubs/work/Achilleas/hg19/yank_put.txt'
#subset.by.region<- function(df,region.info, expansion=0){
subset.by.region<- function(df,region, expansion=0){
	#MBstart = as.numeric(region.info[[3]]) - expansion
	MBstart = as.numeric(region$region_start) - expansion
	MBend = as.numeric(region$region_end) + expansion
	#MBend = as.numeric(region.info[[4]]) + expansion
	chrom =as.numeric(region$gene_chr)
	#chrom =as.numeric(region.info[[2]])
	BPstart <- MBstart  * 1e6
	BPend <- MBend  * 1e6
	sub_df <- as.numeric(df$CHR) == chrom &
	 as.numeric(df$BP) >= BPstart &
	  as.numeric(df$BP) <= BPend;
	region_results <- df[sub_df, ]
	#print(head(region_results))
	return(region_results)
}

subset.by.gene <- function(df, gene){
	sub_df <- df$GENE == gene
	gene_df <- df[sub_df, ]
	#gene_results <- subset(df, sub_df)
	return(gene_df)
}

unique.gene.list <- function(df){
	genes <- as.character(df$GENE)
	gene.list <- unique(genes)
	#print(gene.list)
	return(gene.list)
}

#write.genelist <- function(table_loc, region.info, out_folder,yank_loc){
write.genelist <- function(table_loc, region, out_folder,yank_loc){
	#print(region.info)
	#region.ID = as.character(region.info[[1]])
	region.ID = as.character(region$region_id)
	
	out_name = paste0(region.ID,"_genes.list")
	region_folder = file.path(out_folder, region.ID)
	out_loc <- file.path(region_folder,out_name, fsep = .Platform$file.sep)

	results <- read.table(table_loc, T,strip.white = TRUE)
	gene_region_results <- subset.by.region(results,region)
	#gene_region_results <- subset.by.region(results,region.info)
	sig_gene_results <- as.numeric(gene_region_results$EMP1) <= 5e-3
	sig_region_results <-gene_region_results[sig_gene_results,]
	##sig_region_results <-subset(gene_region_results, sig_gene_results)
	snp_region_results <- subset.by.region(results, region,expansion=.5)
	#snp_region_results <- subset.by.region(results, region.info,expansion=.5)
	gene.list <- unique.gene.list(sig_region_results)
	
	#genes <- as.character(gene_region_results$GENE)
	#gene.list <- unique(genes)
	#print(gene.list)
	if (!(length(gene.list) == 0)){
		if (! file.exists(region_folder)){
		dir.create(region_folder)
		}
	}
	#print(out_folder)
	if ( file.exists(out_loc)){
	unlink(out_loc)
	}
	lapply(gene.list, write, out_loc, append=TRUE, ncolumns=1)
	lapply(gene.list, function(x){write.gene.max(x, snp_region_results, region_folder)})
	yank_thing <- adply(gene.list, 1,function(x){retrieve.yank.lines(x, snp_region_results, yank_loc,region)})
		#yank_thing <- adply(gene.list, 1,function(x){retrieve.yank.lines(x, snp_region_results, yank_loc,region.info)})
	#yank_thing <- lapply(gene.list, function(x){retrieve.yank.lines(x, snp_region_results, yank_loc,region.info)})
	print(yank_thing)

	return(yank_thing)

}


write.title <- function(yank.loc, df){
names <- colnames(df)
names <- append(names, c("REGION","START","END","\n"))
cat(names, file = yank.loc, sep = "\t", fill = FALSE, labels = NULL, append = FALSE)
}


do.stuff <-function(region, table_loc, out_folder,yank_loc){
	#region.info <- read.region.info(region)
	#yank_thing <- write.genelist(table_loc, region.info, out_folder,yank_loc)
	yank_thing <- write.genelist(table_loc, region, out_folder,yank_loc)
	return(yank_thing)
}

read.region.info <- function(region){
	print(region)
	#chrom <- as.numeric(region[1])
	chrom <- as.numeric(region$gene_chr)
	#MBstart <- as.numeric(region[3])
	MBstart <- as.numeric(region$region_start)
	#MBend <- as.numeric(region[4])
	MBend <- as.numeric(region$region_end)
	region.ID <-as.character(region$region_id)
	region.info = c(region.ID,chrom,MBstart,MBend)
	#print(region.info)
	return(region.info)
}



write.gene.max <- function(gene, df, out_folder){
	ref_name = paste0(gene,"_ref.txt")
	#print(ref_name)
	ref_loc <-file.path(out_folder, ref_name, fsep = .Platform$file.sep)
	gene_results <- subset.by.gene(df, gene)
	#sub_df <- df$GENE == gene
	#gene_results <- subset(df, sub_df)
	ref_snp <- as.character( gene_results$EXTREME_LZ[ order(gene_results$EMP1)[1] ] );
	sink(ref_loc)
	cat(ref_snp)
	sink()

}

#include.region.cols <- function(df, region.info){
include.region.cols <- function(df, region){
	df$REGION <- as.character(region$region_id)
	#df$REGION <- region.info[1]
	#df$START <- region.info[3]
	df$START <- as.character(region$region_start)
	#df$END <- region.info[4]
	df$END <- as.character(region$region_end)
	return(df)
}

write.yank.lines <- function(df,yank.loc){
write.table(df, file = yank.loc, append = TRUE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

}
bind.yank.lines <- function(total_df, input_df){
total_df <- rbind(total_df, input_df)
}


#retrieve.yank.lines <- function(gene, df,yank.loc, region.info){
retrieve.yank.lines <- function(gene, df,yank.loc, region){
	gene_results <- subset.by.gene(df, gene)
	#sub_df <- df$GENE == gene
	#gene_results <- subset(df, sub_df)
	gene_results <- include.region.cols(gene_results, region)
	#gene_results <- include.region.cols(gene_results, region.info)
	#gene_results$REGION <-region.ID
	min.p = as.numeric(gene_results$EMP1[order(gene_results$EMP1)[1]])
	#min.p = as.numeric(ordered$EMP1[1])
	#print(min.p)
	yank_sub_results <- as.numeric(gene_results$EMP1) == min.p
	yank_region_results <- gene_results[yank_sub_results,]
	#yank_region_results <-subset(gene_results, yank_sub_results)
	#print(yank_region_results)
	write.yank.lines(yank_region_results,yank.loc)
	return(yank_region_results)

}

main <-function(cell.type=NULL,region_loc=NULL,out_folder=NULL){
library(plyr)
folder_loc = "/home/jkb4y/ubs/work/data/Achilleas/cis-eQTL"
file_name = paste0(cell.type,"_cis_eqtls_permwBP_JB.txt")
table_loc <- file.path(folder_loc, file_name, fsep = .Platform$file.sep)

regions <- read.table(region_loc,T,strip.white = TRUE)
yank_name = paste0(cell.type,"_yank.tbl")
yank_second = paste0(cell.type,"_yank2.tbl")
yank_loc <- file.path(out_folder, yank_name)
yank2_loc <- file.path(out_folder,yank_second)
out_folder <-file.path(out_folder, cell.type, "GeneLists")
df <- read.table(table_loc,T)
write.title(yank_loc,df)
df <- NULL
thing <- adply(regions, 1, function(x){do.stuff(x, table_loc=table_loc,out_folder=out_folder, yank_loc=yank_loc)})
	myvars <- c("REGION", "GENE","SNP_AP","SNP_IM","CHR","BP","EMP1")
	newdata <- thing[myvars]
	write.table(newdata, file = yank2_loc, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)

#thing <- apply(regions, 1, function(x){do.stuff(x, table_loc=table_loc,out_folder=out_folder, yank_loc=yank_loc)})
print(thing)
}
