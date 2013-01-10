

filter_p <-function(df, pmax, perm){
	if (perm){
		df$P = df$EMP1
	}
	sub <- as.numeric(df$P) <= as.numeric(pmax)
	df.filter <-df[sub,]
	if (perm){
		df.filter <- df.filter[ , -which(names(df.filter) %in% c("P"))]
	}
	return(df.filter)

}
write_out <- function(df, out_loc, sep="\t"){
	write.table(df, file = out_loc, append = FALSE, quote = FALSE, sep = sep,
            	eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)
}




write_filter <- function(cell, perm=FALSE, exon=TRUE, pmax=1e-4){
	#basefolder = '/home/jkb4y/cphgdesk_share/Achilleas/cis-eQTLs_Oct2012/data'
	basefolder = '/home/jkb4y/ubs/work/data/Achilleas/cis-eQTLs_Oct2012/'
	permflag = ''
	
	if (perm){permflag = '_perm'}
	if (exon){
		filename <- paste0(cell, '_exon_cis_eqtl',permflag,'_results_noNA.txt')
		outname <- paste0(cell, '_exon_cis_eqtl', permflag, '_results_filtered_1e-4.txt')
		}
	else{
		filename <- paste0(cell,'_cis_transcript_eqtls',permflag,'_all_new_noNA.txt')
		outname <- paste0(cell,'_cis_transcript_eqtls',permflag,'_filtered_1e-4.txt')
		}
	file_loc <- file.path(basefolder, filename)
	out_loc <- file.path(basefolder, outname)
	results <- read.table(file_loc,T)
	results <- filter_p(results, pmax, perm)
	write_out(results, out_loc, sep=' ')



}

filter_out_region <- function(df, region, expansion=0){
	MBstart = as.numeric(region$region_start) - expansion
	MBend = as.numeric(region$region_end) + expansion
	chrom =as.numeric(region$chr)
	BPstart <- MBstart  * 1e6
	BPend <- MBend  * 1e6
	sub <- as.numeric(df$CHR) == chrom &
	 		as.numeric(df$BP) >= BPstart &
	  		as.numeric(df$BP) <= BPend;
	sub_df <- df[!sub, ]
	#print(head(region_results))
	return(sub_df)

}


inter_region_check <- function(results, regions.df, pmax, perm, expansion){
	results <- filter_p(results, pmax, perm)
	for (i in 1:nrow(regions.df))  {
		region <- regions.df[i,]
		print(paste("Now removing region", as.character(region$region_id)))
		results <- filter_out_region(results, region, expansion)
	}
	return(results)
}

table.locate <- function(cell, perm=FALSE, exon=FALSE){
	tablefolder = '/home/jkb4y/ubs/work/data/Achilleas/cis-eQTLs_Oct2012/'
	permflag = ''
	if (perm){permflag = '_perm'}
	if (exon){
		tablename <- paste0(cell, '_exon_cis_eqtl',permflag,'_results_noNA_JB.txt')
		}
	else{
		tablename <- paste0(cell,'_cis_transcript_eqtls',permflag,'_all_new_noNA_JB.txt')
		}
	table_loc = file.path(tablefolder, tablename)
}

write_inter_region <- function(cell,perm=FALSE, exon=FALSE, pmax=1e-4, expansion=0){

	tablefolder = '/home/jkb4y/ubs/work/data/Achilleas/cis-eQTLs_Oct2012/'
	outbase = '/home/jkb4y/ubs/work/results/Achilleas/cis-eQTLs_Oct2012/'
	region_loc = "/home/jkb4y/work/data/Region_Lists/hg19/achilleas_all_09202012.txt"
	regions.df <- read.table(region_loc,T,strip.white = TRUE)
	permflag = ''
	endtag = '_interregion.tbl'
	if (perm){permflag = '_perm'}
	transflag = 'transcript'
	if (exon){
		transflag = 'exon'
		tablename <- paste0(cell, '_exon_cis_eqtl',permflag,'_results_noNA.txt')
		outname <- paste0(cell, '_exon_cis_eqtl', permflag, endtag)
		}
	else{
		transflag = 'transcript'
		tablename <- paste0(cell,'_cis_transcript_eqtls',permflag,'_all_new_noNA.txt')
		outname <- paste0(cell,'_cis_transcript_eqtls',permflag, endtag)
		}
	table_loc = file.path(tablefolder, tablename)
	outfolder = file.path(outbase, transflag, paste0(cell,permflag),'RegionYank')
	out_loc = file.path(outfolder, outname)
	results <- read.table(table_loc, T)
	results <- inter_region_check(results, regions.df, pmax, perm, expansion)
	write_out(results, out_loc)

}

inter_region_all_cells <- function(perm=FALSE, exon=FALSE, pmax=1e-4, expansion=1){
	for (cell in c("B","CD4","CD8","MONO","NK")){
	write_inter_region(cell, perm, exon, pmax, expansion)
	}
}



filter_for_snp <- function(snp, df){
	sub <- as.character(df$SNP_LZ_JB) == snp
	df_sub <- df[sub, ]
	df_sub <- df_sub[!is.na(df_sub$CHR),] 
	#print(head(df_sub))
	#exit()
	return(df_sub)
}

make.df <- function(snp_vector, eqtl_df){
	print(snp_vector)
	snp <- as.character(snp_vector$conditional_SNP)
	print(snp)
	out_df <- filter_for_snp(snp, eqtl_df)
	#na.omit(filter_for_snp(snp, eqtl_df))
	#print(out_df)
	if (!nrow(out_df) == 0){
	out_df$conditional_SNP <- snp
	out_df$Region_ID <- as.character(snp_vector$Region_ID)
	out_df$SNP. <- as.character(snp_vector$SNP.)
	}
	return(out_df)

}

filter_for_snps <- function(cell, perm=FALSE, exon=FALSE){
	library("plyr")
	permflag = ''
	transflag = 'transcript'
	endtag = '_intersect_SI.tbl'
	if(perm){permflag = '_perm'}
	outname <- paste0(cell,'_cis_transcript_eqtls',permflag, endtag)
	if(exon){
		transflag = 'exon'
		outname <- paste0(cell, '_exon_cis_eqtl',permflag,endtag)
		}
		
	snp_list_loc <- '/home/jkb4y/cphgdesk_share/Achilleas/cis-eQTLs_Oct2012/data/SI_SNP_LIST_IC_NOMHC_noDup_20121005_Query2.txt'
	out_folder <- file.path('/home/jkb4y/cphgdesk_share/Achilleas/cis-eQTLs_Oct2012/', transflag, paste0(cell,permflag))
	out_loc <- file.path(out_folder, outname)
	eqtl_loc <- table.locate(cell, perm, exon)
	eqtl_table <- read.table(eqtl_loc, T)
	eqtl_table <- eqtl_table[ ,-which(names(eqtl_table) %in% c('NMISS','TEST','EXTREME_LZ','annotation'))]
	print(head(eqtl_table))
	snp_table <- read.table(snp_list_loc, T,strip.white = TRUE)
	snp_table <- snp_table[ ,c('conditional_SNP','Region_ID','SNP.')]
	#print(head(snp_table))
	thing <- adply(snp_table,1,function(x){make.df(x,eqtl_table)})
	names(thing)[which(names(thing)=="SNP_AP")] = "SNP_eQTL"
	names(thing)[which(names(thing)=="SNP.")] = "SNP_TOP"
	names(thing)[which(names(thing)=="conditional_SNP")] = "SI_SNP"
	thing <- thing[ ,-which(names(thing) %in% c('SNP_LZ_JB'))]
	
	write_out(thing, out_loc)
	

}
