

subset.by.region<- function(df,region, expansion=0){
# Subsets a dataframe using region information
# Args:
#   df: dataframe of genetic data
#	region: vector from region dataframe
#	expansion: number of MBs to expand the region on either side. Default is 0.
# Returns:
#	subsetted dataframe
	MBstart = as.numeric(region$region_start) - expansion
	MBend = as.numeric(region$region_end) + expansion
	chrom =as.numeric(region$chr)
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
#subset data frame by a specific gene in the GENE column
	sub_df <- df$GENE == gene
	gene_df <- df[sub_df, ]
	#gene_results <- subset(df, sub_df)
	return(gene_df)
}
unique.exon.list <- function(df){
# create a list genes in the dataframe
	exons <- as.character(df$EXON)
	exon.list <- unique(exons)
	#print(gene.list)
	#print(gene.list)
	return(exon.list)
}


unique.gene.list <- function(df){
# create a list genes in the dataframe
	genes <- as.character(df$GENE)
	gene.list <- unique(genes)
	#print(gene.list)
	#print(gene.list)
	return(gene.list)
}

write.genelist <- function(results, region, out_folder,pmax, expansion=0,exon){
#write.genelist <- function(table_loc, region, out_folder,pmax, perm){
#write.genelist <- function(table_loc, region, out_folder,yank_loc, pmax, perm){
	region.ID = as.character(region$region_id)
	#print(region.ID)
	
	out_name = paste0(region.ID,"_genes.list")
	region_folder = file.path(out_folder, region.ID)
	out_loc <- file.path(region_folder,out_name, fsep = .Platform$file.sep)

	#results <- read.table(table_loc, T,strip.white = TRUE)
	#if (perm){results$P = results$EMP1}
	gene_region_results <- subset.by.region(results,region)
	sig_gene_results <- as.numeric(gene_region_results$P) <= as.numeric(pmax)
	
	sig_region_results <-gene_region_results[sig_gene_results,]
	snp_region_results <- subset.by.region(results, region,expansion=expansion)
	gene.list <- unique.gene.list(sig_region_results)
	#genes <- as.character(gene_region_results$GENE)
	#gene.list <- unique(genes)
	#print(gene.list)
	if (!(length(gene.list) == 0)){
		if (! file.exists(region_folder)){
		dir.create(region_folder)
		}
	}
	if ( file.exists(out_loc)){
	unlink(out_loc)
	}
	lapply(gene.list, write, out_loc, append=TRUE, ncolumns=1)   #########put me back!!
	lapply(gene.list, function(x){write.gene.max(x, snp_region_results, region_folder,exon)})    #########put me back!!
	yank_thing <- adply(gene.list, 1,function(x){retrieve.yank.lines(x, snp_region_results, region)})
	#yank_thing <- adply(gene.list, 1,function(x){retrieve.yank.lines(x, snp_region_results, yank_loc,region)})
	#yank_thing <- lapply(gene.list, function(x){retrieve.yank.lines(x, snp_region_results, yank_loc,region.info)})
	#print(yank_thing)

	return(yank_thing)

}


write.title <- function(yank.loc, df){
names <- colnames(df)
names <- append(names, c("REGION","START","END","\n"))
cat(names, file = yank.loc, sep = "\t", fill = FALSE, labels = NULL, append = FALSE)
}

#do.stuff <-function(region, table_loc, out_folder,yank_loc, pmax, perm){
#do.stuff <-function(region, table_loc, out_folder,pmax, perm){
do.stuff <-function(region, table.df, out_folder,pmax, expansion,exon){
	yank_thing <- write.genelist(table.df, region, out_folder, pmax=pmax, expansion=expansion, exon=exon)
	#yank_thing <- write.genelist(table_loc, region, out_folder, pmax=pmax,perm=perm)
	#yank_thing <- write.genelist(table_loc, region, out_folder,yank_loc, pmax=pmax,perm=perm)
	return(yank_thing)
}

read.region.info <- function(region){
	print(region)
	chrom <- as.numeric(region$chr)
	MBstart <- as.numeric(region$region_start)
	MBend <- as.numeric(region$region_end)
	region.ID <-as.character(region$region_id)
	region.info = c(region.ID,chrom,MBstart,MBend)
	return(region.info)
}



write.gene.max <- function(gene, df, out_folder,exon){
	ref_file = paste0(gene,"_ref.txt")
	#print(ref_name)
	ref_loc <-file.path(out_folder, ref_file, fsep = .Platform$file.sep)
	gene_results <- subset.by.gene(df, gene)
	if (exon){
		exon_file = paste0(gene, "_exons.txt")
		exon_loc <- file.path(out_folder, exon_file, fsep = .Platform$file.sep)
		exon.list <- unique.exon.list(gene_results)
		lapply(exon.list, write, exon_loc, append=TRUE, ncolumns=1)
		}
	#sub_df <- df$GENE == gene
	#gene_results <- subset(df, sub_df)
	top = order(gene_results$P)[1]
	#ref_name = as.character(gene_results$SNP[top])
	ref_snp <- as.character(gene_results$EXTREME_LZ[ top ] );
	ref_name <- as.character( gene_results$SNP_LZ_JB[ top ] );
	if ( is.na(ref_name) |!(substr(ref_name,start=1,stop=2) == 'rs')){
	ref_name = as.character(gene_results$SNP_IM[ top ] )
	}
	print(ref_name)
	write_me = paste(ref_snp, ref_name, sep='\t')
	#print(write_me)
	sink(ref_loc)
	cat(paste(ref_snp, ref_name, sep='\t', append = FALSE))
	sink()

}


include.region.cols <- function(df, region){
	df$REGION <- as.character(region$region_id)
	df$START <- as.character(region$region_start)
	df$END <- as.character(region$region_end)
	return(df)
}

write.yank.lines <- function(df,yank.loc){
write.table(df, file = yank.loc, append = TRUE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

}

retrieve.yank.lines <- function(gene, df, region){
#retrieve.yank.lines <- function(gene, df,yank.loc, region){
	gene_results <- subset.by.gene(df, gene)
	#sub_df <- df$GENE == gene
	#gene_results <- subset(df, sub_df)
	gene_results <- include.region.cols(gene_results, region)
	min.p = as.numeric(gene_results$P[order(gene_results$P)[1]])
	yank_sub_results <- as.numeric(gene_results$P) == min.p
	yank_region_results <- gene_results[yank_sub_results,]
	#yank_region_results <-subset(gene_results, yank_sub_results)
	#print(yank_region_results)
	#write.yank.lines(yank_region_results,yank.loc)
	return(yank_region_results)

}

read_and_filter <- function(table_loc, filter, perm, cis=TRUE, exon=FALSE){
	df <- read.delim(table_loc, T, stringsAsFactors=FALSE, sep='\t')
	names(df)[names(df)=='SNP']<-'SNP_AP'
	print(head(df))
	#if (! cis & ! exon){
	#	df <- df[df$TRANS.CIS == 'trans',]
	#}
	if (perm){
		df$P = df$EMP1
	}
	sub <- as.numeric(df$P) <= as.numeric(filter)
	df <-df[sub,]
	return(df)

}

region_yank <-function(cell.type=NULL,region_loc=NULL,perm=FALSE, pmax = 1, exon=FALSE, filter=1e-4, cis=TRUE, expansion=0){
	library(plyr)
	folder_loc = "/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Feb2013_pcaCorrected"
	#folder_loc = "/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Jan2013/data"
	#out_basefolder = "/home/jkb4y/ubs/work/results/Achilleas/eQTLs_Jan2013"
	out_basefolder = "/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected"
	columns <- c("REGION", "GENE","SNP_AP","SNP_IM","CHR","BP","P")
	print(columns)
	perm_flag = ""
	if (perm == TRUE){
		perm_flag = '_perm'
		columns <- c(columns[-which(columns %in% c("P"))], c("EMP1"))
		print(columns)
		}
	if (cis){
	cisflag = 'cis'
	#transcriptflag = 'transcript'
	#file_name = paste0(cell.type,"_cis_transcript_eqtls",perm_flag,"_all_new_noNA_JB.txt")
	#file_name = paste0(cell.type,"_cis_transcript_eqtls",perm_flag,"_filtered_1e-4_JB.txt")
	if (exon){
		transcriptflag = 'exon'
		#file_name = paste0(cell.type,"_exon_cis_eqtl",perm_flag,"_results_noNA_JB.txt")
		file_name = paste0("exonCis", cell.type,".txt")
		columns <- c(columns, c("EXON"))
		print(columns)
		}
	else{
		transcriptflag = 'transcript'
		file_name = paste0('transcriptCis',cell.type,'.txt')
		}
	}
	else{
		cisflag = 'trans'
		if (exon){
		#file_name = paste0(cell.type,'_all_data_JB.txt')
		#if (exon){
		transcriptflag = 'exon'
		file_name = paste0("exonTrans",cell.type,"_filtered_","1e-04",".txt")
		}
		else{ 
		transcriptflag = 'transcript'
		file_name = paste0('transcriptTrans',cell.type,'_filtered_1e-04.txt')
		}
		}
	out_basefolder = file.path(out_basefolder, transcriptflag, cisflag)
	yank_name = paste0(cell.type,'_',cisflag,perm_flag,"_",transcriptflag,"_R_yank.tbl")
	print(yank_name)
	#perm_flag = ""
	#if (perm == TRUE){
	#perm_flag = '_perm'
	#}
	#if (! cis){
	#folder_loc = "/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Jan2013"
	#file_name = paste0(cell.type,'_all_data_JB.txt')
	#out_basefolder = "/home/jkb4y/ubs/work/results/Achilleas/transcript_eqtl_results_Oct12"
	#yank_name = paste0(cell.type,"_trans_R_yank.tbl")
	
	#}
	#file_name = paste0(cell.type,"_cis_transcript_eqtls",perm_flag,"_all_new_noNA_JB.txt")
	table_loc <- file.path(folder_loc, file_name, fsep = .Platform$file.sep)
	print(table_loc)
	table.df <- read_and_filter(table_loc, filter, perm, cis, exon)
	
	#print(names(table.df))
	regions <- read.table(region_loc,T,strip.white = TRUE)
	#yank_name = paste0(cell.type,perm_flag,"_",transflag,"_R_yank.tbl")
	#yank_second = paste0(cell.type,perm_flag,"_yank2.tbl")
	#yank_loc <- file.path(out_folder, yank_name)
	#yank2_loc <- file.path(out_folder,yank_second)
	cell_folder = file.path(out_basefolder, paste0(cell.type, perm_flag))
	print(cell_folder)
	expansion_folder = file.path(cell_folder, paste0("expansion_",as.character(expansion),"MB"))
	print(expansion_folder)
	#out_folder <-file.path(out_basefolder, paste0(cell.type,perm_flag), "GeneLists")
	out_folder <- file.path(expansion_folder, "GeneLists")
	yank_folder <-file.path(cell_folder, "RegionYank")
	yank_loc <- file.path(yank_folder, yank_name)
	if (! file.exists(cell_folder)){
		dir.create(cell_folder)
		}
	if (! file.exists(expansion_folder)){
		dir.create(expansion_folder)
		}
	if (! file.exists(out_folder)){
		dir.create(out_folder)
		}
	if (! file.exists(yank_folder)){
		dir.create(yank_folder)
		}
	#df <- read.table(table_loc,T)
	#write.title(yank_loc,df)
	#df <- NULL
	thing <- adply(regions, 1, function(x){do.stuff(x, table.df = table.df,out_folder=out_folder, pmax=pmax, expansion=expansion, exon=exon)})
	
	#thing <- adply(regions, 1, function(x){do.stuff(x, table_loc=table_loc,out_folder=out_folder, pmax=pmax,perm=perm)})
	#thing <- adply(regions, 1, function(x){do.stuff(x, table_loc=table_loc,out_folder=out_folder, yank_loc=yank_loc, pmax=pmax,perm=perm)})
	#columns <- c("REGION", "GENE","SNP_AP","SNP_IM","CHR","BP","P")
	yank.df <- thing[columns]
	write.table(yank.df, file = yank_loc, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)

#thing <- apply(regions, 1, function(x){do.stuff(x, table_loc=table_loc,out_folder=out_folder, yank_loc=yank_loc)})
#print(thing)
}

yank_all_cells <- function(region_loc="/home/jkb4y/work/data/Region_Lists/hg19/achilleas_all_09202012.txt",perm=FALSE, filter=1e-4, exon=FALSE, cis=TRUE, expansion=1){
	if (cis){ pmax = 0.05/1125}
	else {pmax = 0.05/55336}
	for (cell in c("B","CD4","CD8","MONO","NK")){
		region_yank(cell.type=cell,region_loc=region_loc,perm=perm, pmax=pmax, exon=exon, filter=filter,cis=cis, expansion=expansion)
	
	}
	
}
