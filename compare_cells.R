make_manageable_table <- function(file_loc, cell_type, wanted_cols=c('P','A1','BETA'), merge_cols=c('SNP','GENE'),filter='1e-4', cis=TRUE, exon=FALSE){
	#cell_table <- read.table(file_loc, T)
	#source('gene_lists.R')
	cell_table <- read_and_filter(file_loc,as.numeric(filter),perm=FALSE)
	cell_table$SNP <- cell_table$SNP_AP
	names(cell_table)[names(cell_table)=="SNP_LZ_JB"] <- "SNP_JB"
	#cell_table$SNP_JB <- cell_table$SNP_LZ_JB
	print(colnames(cell_table))
	cell_table <- cell_table[ , c(merge_cols,wanted_cols)]
	for (name in wanted_cols){
		new_name <- paste0(cell_type, '.',name)
		names(cell_table)[which(names(cell_table)==name)] <- new_name
	}
	print(head(cell_table))
	return(cell_table)
}



locate_table <- function(cell, exon, cis, tablefolder=NA){
	#tablefolder = '/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Feb2013_pcaCorrected'
	filter_tail = paste0('_filtered_','1e-04','.txt')
	if (exon){
		exonflag = 'exon'
		}
	else{ exonflag = 'transcript'}
	if (cis){ cisflag = 'Cis'}
	else {cisflag = 'Trans'}
	
	#if (cis & ! exon){filter_tail = '.txt'}
	filename = paste0(exonflag,cisflag,cell,filter_tail)
		
	#else if (! cis){
	#filename = paste0(cell,'_all_data_JB.txt')
	#}
	#else {
	#filename = paste0(cell,'_cis_transcript_eqtls',permflag,filter_tail)
	#}
	file_loc = file.path(tablefolder, filename)
	return(file_loc)

}
locate_merge_out <- function(basefolder= NA,perm, exon, cis, filter, type='merge', outfolder=NA){
	out_loc = NULL
if (is.na(outfolder)){
	#basefolder = '/home/jkb4y/ubs/work/results/Achilleas/eQTLs_Oct2012'
	#basefolder= '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected'
	if (!cis){cisfolder = 'trans'; cis_flag = 'trans_'}
	else{cisfolder = 'cis'; cis_flag = 'cis_'}
	if (! cis){ cis_flag = 'trans_'}

	if (exon){transfolder = "exon"}
	else{transfolder = "transcript"}
	if (perm){permflag <- "_perm"}
	else{permflag <- ""}
	outfolder <- file.path(basefolder, transfolder, cisfolder)
	}
	if (type=='merge'){
	out_loc <- file.path(outfolder, paste0(cis_flag,'merge_',filter,'.tbl'))
	}
	else if (type == 'noNA'){out_loc <- file.path(outfolder, paste0(cis_flag, 'NAabove1e-4_merge_',filter,'.tbl'))
	}
	else if (type=='thin'){out_loc <- file.path(outfolder, paste0(cis_flag,'thin_merge_',filter,'.tbl'))
	}
	else if (type =='multigene'){out_loc <- file.path(outfolder, paste0(cis_flag,'multigene_region_merge_',filter,'.tbl'))
	}
	else if (type =='countsum'){out_loc <- file.path(outfolder, paste0(cis_flag,'count_summary_',filter,'.tbl'))
	}
	else if (type =='countcross'){out_loc <- file.path(outfolder, paste0(cis_flag,'count_crosstab_',filter,'.tbl'))
	}
	else if (type =='si'){out_loc <- file.path(outfolder, paste0(cis_flag,'SI_NAabove1e-4_merge.tbl'))
	}
	else if (type =='folder'){out_loc <- outfolder
	}
	print(out_loc)
	
	return(out_loc)

}


write_out <-function(df, out_loc){
	#out_loc <- locate_merge_out(perm=perm,exon=exon,cis=cis,filter=filter,type=type, outfolder=outfolder)
	write.table(df, file = out_loc, append = FALSE, quote = FALSE, sep = "\t",
            	eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)


}

read_and_filter <- function(table_loc, filter, perm){
	df <- read.table(table_loc, header=T, stringsAsFactors=FALSE, sep='\t')
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




add_SI_value <- function(df, si_loc){
	si_table <- read.table(si_loc, T)
	si_list <- si_table$conditional_SNP
	df$SI <- ifelse(df$SNP_JB %in% si_list, TRUE, FALSE)
	print(head(df))
return(df)
}

add_literature_info <- function(df, gwas_loc){
	gwas <- read.csv(gwas_loc,T,sep='\t',quote="")
	gwas$gwas_pos <- paste0('chr',as.character(gwas$Chr_id),":",as.character(gwas$Chr_pos))
	#gwas$disease_list <- gwas$Disease.Trait[gwas$gwas_pos == gwas$gwas_pos,]
	gwas$disease_list <- sapply(gwas$gwas_pos, FUN=function(x){thing <- unique(gwas[gwas$gwas_pos == x, c('Disease.Trait')]);out <- paste(thing, collapse = ';'); return(out)})

	gwas$link_list <- sapply(gwas$gwas_pos, FUN=function(x){thing <- unique(gwas[gwas$gwas_pos == x, c('Link')]);out <- paste(thing, collapse = ' ; '); return(out)})

	df$chr_pos <- paste0('chr',as.character(df$CHR),":",as.character(df$BP))
	df$GWAS <- ifelse(df$chr_pos  %in% gwas$gwas_pos, TRUE, FALSE)
	#df$GWAS <- ifelse(paste0(as.character(df$CHR),":",as.character(df$BP))  %in% gwas$gwas_pos, TRUE, FALSE)
	#df$chr_pos <- paste0('chr',as.character(df$CHR),":",as.character(df$BP))
	df$GWAS_DISEASES <- NA
	df$GWAS_LINKS <- NA
	counter = 1
	countermax = length(unique(gwas$gwas_pos))
	for (snp in unique(gwas$gwas_pos)){
		print(paste0(counter, "/",countermax))
		df$GWAS_DISEASES[df$chr_pos == snp] <- gwas[gwas$gwas_pos == snp,c('disease_list')][[1]]
		df$GWAS_LINKS[df$chr_pos == snp] <- gwas[gwas$gwas_pos == snp,c('link_list')][[1]]
		counter = counter + 1
		}
	df <- df[,-which(colnames(df) %in% c('chr_pos'))]
	gwas <- NULL
	#df$GWAS_DISEASES <- ifelse(df$chr_pos  %in% gwas$gwas_pos, gwas[gwas$gwas_pos == df$chr_pos [[1]], c('disease_list')], NA)
	print(head(df))
	return(df)
	}


add_r2_info <- function(df, r2_loc,yank_loc){
	r2 <- read.table(r2_loc,T)
	r2 <- r2[,c('SNP_B','SNP_A','R2')]
	#print(head(r2))
	yank <- read.table(yank_loc,T)
	yank <- yank[yank$p.value <= 3.5e-7, c('snp_name','p.value')]
	#print(head(yank))
	yank_snps <- as.character(yank$snp_name)
	print(nrow(r2))
	r2 <- r2[which(r2$SNP_A %in% yank_snps), ]
	yank_snps <- NULL
	yank <- NULL
	print(nrow(r2))
	sharedsnps = unique(as.character(df$SNP_JB[df$SNP_JB %in% r2$SNP_B]))
	r2 <- r2[which(r2$SNP_B %in% sharedsnps),]
	print(nrow(r2))
	print(head(sharedsnps))
	df$R2 <- NA
	df$csnp <- NA
	#df[,c('R2','csnp')] <- sapply(df$SNP_JB, FUN=function(x){ if (x %in% sharedsnps){r2val <- as.character(r2$R2[r2$SNP_B == x]); csnpval <- as.character(r2$SNP_A[r2$SNP_B == x]); print(x)} else{ r2val = NA; csnpval = NA}; return(c(r2val,csnpval))})
	counter = 1
	countermax = length(sharedsnps)
	for (snp in sharedsnps){
		print(paste0(counter, "/",countermax))
		print(snp)
		df$R2[df$SNP_JB == snp] <- as.character(r2$R2[r2$SNP_B == snp])
		df$csnp[df$SNP_JB == snp] <- as.character(r2$SNP_A[r2$SNP_B == snp])
		#print(df$csnp[df$SNP_JB == snp])
		counter = counter + 1
	}
	
	#overlap<- which(df$SNP_JB %in% r2$SNP_B & df$REGION != '17q21.31_z')
	#snp_overlap <- df$SNP_JB[overlap]
	#df$R2 <- NA
	#df$csnp <- NA
	#counter = 1
	#countermax = length(overlap)
	#for (index in overlap){
	#	print(paste0(counter, "/",countermax))
		#print(as.character(df$SNP_JB[[index]]))
	#	df$R2[index] <- r2$R2[r2$SNP_B == as.character(df$SNP_JB[[index]])]
	#	df$csnp[index] <- as.character(r2$SNP_A[r2$SNP_B == as.character(df$SNP_JB[[index]])])
	#	counter = counter + 1
	#}
	r2 <- NULL
	#snp_overlap <- NULL
	#df$R2 <- ifelse(df$SNP_JB %in% r2$SNP_B, r2[r2$SNP_B == df$SNP_JB, c('R2') ][[1]], NA)
	print(head(df))
	return(df)
	


}

match_region <- function(region, df){
	#print(region)
	sub <- as.numeric(df$CHR) == as.numeric(region$chr) & as.numeric(df$BP) > as.numeric(region$region_start)*1e6 & as.numeric(df$BP) < as.numeric(region$region_end)*1e6
	df[sub,'REGION'] <- as.character(region$region_id)
	return(df)
}

add_regions <- function(df, region_loc){
	regions <- read.table(region_loc, T)
	rcount = nrow(regions)
	#print(rcount)
	
	#print(regions)
	for (i in 1:rcount){
		print(regions[i,])
		df <- match_region(regions[i,],df)	
	}
	print(head(df))
	return(df)
}

determine_merge_cols <- function(exon){
	
	if (! exon){
		merge_cols <- c('SNP','GENE','CHR','BP', 'SNP_JB','SNP_IM')
		}
	else{
		merge_cols <- c('SNP','GENE','EXON','CHR','BP', 'SNP_JB','SNP_IM')
	
		}
	return(merge_cols)
}

create_empty_df <- function(col_names){
	df <- data.frame(t(rep(NA,length(col_names))))
	names(df) <- col_names
	df <- df[-1,]
	return(df)
}



merge_tables <- function(outfolder=NA, perm=FALSE, exon=FALSE, filter='1e-4',region_loc='/home/jkb4y/work/data/Region_Lists/hg19/achilleas_all_09202012.txt',cis=TRUE, si_loc='/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected/data/si_SNP_NoMHC_20121024_Query.txt', gwas_loc="/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected/data/gwascatalog_20130111.txt",r2_loc='/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/eurmeta_LD/all_regions_r2_0.ld',yank_loc='/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/RegionYank/eurmeta_yank.tbl', basefolder=NA, tablefolder = NA, ...){
#merge_tables <- function(outfolder=NA, perm=FALSE, exon=FALSE, filter='1e-4',region_loc='/home/jkb4y/work/data/Region_Lists/hg19/achilleas_all_09202012.txt',cis=TRUE, si_loc='/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected/data/si_SNP_NoMHC_20121024_Query.txt', gwas_loc="/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected/data/gwascatalog_20130111.txt",r2_loc='/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/eurmeta_LD/all_regions_r2_0.ld',yank_loc='/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/RegionYank/eurmeta_yank.tbl', basefolder='/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected', tablefolder = '/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Feb2013_pcaCorrected', ...){
	library(plyr)
	cell_types <- c('B','CD4','CD8','MONO','NK')
	#merge_cols <- c('SNP','GENE','CHR','BP', 'SNP_JB')
	#if (exon){ merge_cols <- c(merge_cols,'EXON')}
	wanted_cols <- c('P','BETA','A1')
	merge_cols <- determine_merge_cols(exon)
	empty_df <- create_empty_df(merge_cols)
	if (! exon){
		merge_cols <- c('SNP','GENE','CHR','BP', 'SNP_JB','SNP_IM')
		mamatable <- data.frame(SNP = character(0), GENE = character(0),
			CHR = character(0), BP = character(0),SNP_JB = character(0),SNP_IM = character(0))
		}
	else{
		merge_cols <- c('SNP','GENE','EXON','CHR','BP', 'SNP_JB','SNP_IM')
		mamatable <- data.frame(SNP = character(0), GENE = character(0),
			EXON = character(0), CHR = character(0), BP = character(0),SNP_JB = character(0),SNP_IM = character(0))
	
		}
	for (cell_type in cell_types){
		table_loc <- locate_table(cell_type, exon, cis, tablefolder)
		print(table_loc)
		cell_table <- make_manageable_table(table_loc, cell_type, wanted_cols=wanted_cols,merge_cols=merge_cols, filter, cis, exon)
		mamatable <-  merge(mamatable, cell_table, by=merge_cols, all=TRUE, all.x=TRUE, all.y=TRUE)
		cell_table <- NULL
	}
	
	
	mamatable <- add_regions(mamatable, region_loc)
	mamatable <- add_SI_value(mamatable, si_loc)
	mamatable <- add_r2_info(mamatable, r2_loc, yank_loc)
	
	mamanona <- mamatable[,c(merge_cols,'REGION','SI','R2','csnp')]
	
	names(mamatable)[names(mamatable)=="BP"] <- "bp"
	#names(mamanona)[names(mamanona)=="BP"] <- "bp"
	#mamatable <- add_literature_info(mamatable, gwas_loc)
	
	
	sumtable <- summarize_mama(mamatable)
	sum_loc <- locate_merge_out(basefolder=basefolder,perm=perm,exon=exon,cis=cis,filter=filter,type='countsum', outfolder=outfolder)
	print(sum_loc)
	write.table(sumtable,file=sum_loc, row.names=FALSE, col.names=TRUE, sep='\t',quote=FALSE)
	sumtable <- NULL
	
	crosstab <- crosstab_mama(mamatable)
	cross_loc <- locate_merge_out(basefolder=basefolder,perm=perm,exon=exon,cis=cis,filter=filter,type='countcross', outfolder=outfolder)
	print(cross_loc)
	write.table(crosstab,file=cross_loc, row.names=FALSE, col.names=TRUE, sep='\t',quote=FALSE)
	crosstab <- NULL
	
	print(head(mamatable))
	mama_loc <- locate_merge_out(basefolder=basefolder,perm=perm,exon=exon,cis=cis,filter=filter,type='merge', outfolder=outfolder)
	write_out(mamatable, mama_loc)
	multigenes <- find_multiregs(mamatable)
	
	#multigenes <- add_regions(regions, multigenes)
	multi_loc <- locate_merge_out(basefolder=basefolder,perm=perm,exon=exon,cis=cis,filter=filter,type='multigene', outfolder=outfolder)
	write_out(multigenes, multi_loc)
	multigenes <- NULL
	
	minimama <- thin_out(mamatable, cell_types)
	thin_loc <- locate_merge_out(basefolder=basefolder,perm=perm,exon=exon,cis=cis,filter=filter,type='thin', outfolder=outfolder)
	write_out(minimama, out_loc=thin_loc)
	minimama <- NULL
	#return(mamatable)
	
	#make merges without a p-value threshold
	mamatable <-NULL
	#create si_list
	si_table = read.table(si_loc, T)
	si_list <- si_table$conditional_SNP
	si_merge <- empty_df
	names(si_merge)[names(si_merge)=="bp"] <- "BP"
	
	for (cell_type in cell_types){
		table_loc <- locate_table(cell_type, exon, cis, tablefolder)
		print(table_loc)
		cell_table <- make_manageable_table(table_loc, cell_type, wanted_cols=wanted_cols,merge_cols=merge_cols, filter=1, cis, exon)
		mamanona <-  merge(mamanona, cell_table, by=merge_cols, all.x=TRUE, all.y=FALSE)
		cell_table <- cell_table[cell_table$SNP_JB %in% si_list, ]
		si_merge <- merge(si_merge,cell_table, by=merge_cols, all=TRUE)
		print(head(si_merge))
		cell_table <- NULL
	}
	nona_loc <- locate_merge_out(basefolder=basefolder,perm=perm,exon=exon,cis=cis,filter=filter,type='noNA', outfolder=outfolder)
	
	names(mamanona)[names(mamanona)=="BP"] <- "bp"
	#extra_snp_name_cols = c(which(names(mamanona)=="SNP_IM"),which(names(mamanona)=="SNP_JB"))
	extra_snp_names = c("SNP_IM","SNP_JB")
	col_order = c(names(mamanona)[-which(names(mamanona) %in% extra_snp_names)],extra_snp_names)
	mamanona <- mamanona[ , col_order]
	print(colnames(mamanona))
	write_out(mamanona, out_loc=nona_loc)
	mamanona <- NULL
	
	si_merge_loc <- locate_merge_out(basefolder=basefolder,perm=perm,exon=exon,cis=cis,filter=filter,type='si', outfolder=outfolder)
	names(si_merge)[names(si_merge)=="BP"] <- "bp"
	extra_snp_names = c("SNP_IM","SNP_JB")
	col_order = c(names(si_merge)[-which(names(si_merge) %in% extra_snp_names)],extra_snp_names)
	si_merge <- si_merge[ , col_order]
	print(colnames(si_merge))
	write_out(si_merge, out_loc=si_merge_loc)
	si_merge <- NULL
}

find_multiregs <- function(df){
	multigene_snps <- df[duplicated(df$SNP),]
	#snps_sub <- df$SNP %in% multigenes$SNP
	multigenes <-df[df$SNP %in% unique(multigene_snps$SNP),]
	#print(head(multigenes))
	return(multigenes)

}

thin_out <- function(total_table, cell_types){
#print("Made it.")
 col.count <- ncol(total_table)
 countNAs <- apply(total_table, 1, function(x) sum(is.na(x)))
 #print(head(countNAs))
 total_table$NAs <- countNAs
 #print(head(total_table))
 sub <- total_table$NAs < 12
 total_table <- total_table[sub, 1:col.count]
 #print(head(total_table))
return(total_table)
}


locate_yank <- function(basefolder, cell, perm, exon, cis){
	cisfolder = 'cis'
	transfolder = 'transcript'
	if (perm){
		cell = paste0(cell,'_perm')
	}
	if (exon){
	transfolder = 'exon'
	}
	if (!cis){
	cisfolder = 'trans'
	basefolder <- file.path(basefolder, transfolder, cisfolder, cell, 'RegionYank')
	filename <- paste0(cell, '_R_yank.tbl')
	file_loc <- file.path(basefolder, filename)
	return(file_loc)
	}

}



double_summary_helper <- function(cell_type, df, cell_list){
	sub <- df[!is.na(df[ ,paste0(cell_type,'.P')]), ]
	sumline <- data.frame(CELL_TYPE = cell_type,
							B.SNP_COUNT=character(1),B.GENE_COUNT=character(1),B.UNIQUE_SNP_COUNT=character(1),
							CD4.SNP_COUNT=character(1),CD4.GENE_COUNT=character(1),CD4.UNIQUE_SNP_COUNT=character(1),
							CD8.SNP_COUNT=character(1),CD8.GENE_COUNT=character(1),CD8.UNIQUE_SNP_COUNT=character(1),
							MONO.SNP_COUNT=character(1),MONO.GENE_COUNT=character(1),MONO.UNIQUE_SNP_COUNT=character(1),
							NK.SNP_COUNT=character(1),NK.GENE_COUNT=character(1),NK.UNIQUE_SNP_COUNT=character(1))
	for (cell in cell_list){
		sub_sub <- sub[!is.na(sub[ ,paste0(cell, '.P')]), ]
		sumline[ ,paste0(cell, '.SNP_COUNT')] <- nrow(sub_sub)
		sumline[ ,paste0(cell, '.GENE_COUNT')] <- length(unique(sub_sub$GENE))
		sumline[ ,paste0(cell, '.UNIQUE_SNP_COUNT')] <- length(unique(sub_sub$SNP))	
	}
	print(sumline)
	return(sumline)
}

summary_helper <- function(cell, df){
sub <- df[!is.na(df[ ,paste0(cell,'.P')]), ]

sumline <- data.frame(CELL_TYPE = cell, SNP_COUNT = nrow(sub), GENE_COUNT = length(unique(sub$GENE)), UNIQUE_SNP_COUNT = length(unique(sub$SNP)))
print(sumline)
return(sumline)
}

summarize_mama <- function(df){
	cells <-c('B','CD4','CD8','MONO','NK')
	sumtable <- ldply(cells,.fun=function(x){summary_helper(x, df)})
	return(sumtable)
	
}
crosstab_mama <- function(df){
	cells <-c('B','CD4','CD8','MONO','NK')
	crosstab <- ldply(cells,.fun=function(x){double_summary_helper(x, df, cells)})
	return(crosstab)
}

compare_merges <- function(exon=FALSE, basefolder='/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected/'){
	if (exon){ exonflag = 'exon'}
	else {exonflag = 'transcript'}
	outfolder = file.path(basefolder, exonflag)
	filename = paste0(exonflag,'_cis_trans_merge.tbl')
	cistable_loc = file.path(outfolder, 'cis','cis_merge_4.45e-5.tbl')
	transtable_loc = file.path(outfolder, 'trans','trans_merge_9.05e-7.tbl')
	cistable = read.table(cistable_loc, T, sep='\t', stringsAsFactors=F)
	transtable = read.table(transtable_loc, T, sep='\t',stringsAsFactors=F)

	#RENAME the columns that aren't being merged
	trans_indices <- c(grep(".P", names(transtable), fixed=TRUE),grep(".A1", names(transtable), fixed=TRUE), grep(".BETA", names(cistable), fixed=TRUE),grep("GENE", names(transtable)))
	cis_indices <- c(grep(".P", names(cistable), fixed=TRUE),grep(".A1", names(cistable), fixed=TRUE), grep(".BETA", names(cistable), fixed=TRUE), grep("GENE",names(cistable)))
	colnames(transtable)[trans_indices] <- paste0(names(transtable)[trans_indices],'.TRANS')
	colnames(cistable)[cis_indices] <- paste0(names(cistable)[cis_indices],'.CIS')
	
	#Discover the names of the ones that are to be merged!
	trans_merge_cols <- colnames(transtable)[-trans_indices]
	cis_merge_cols <- colnames(cistable)[-cis_indices]
	mergemerge <-  merge(cistable, transtable, by=trans_merge_cols, all=TRUE)
	print(head(mergemerge))
	return(mergemerge)

}

merge_overlap <- function(exon=FALSE, basefolder='/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected/'){
	dict_loc = '/home/jkb4y/work/data/annot_hg19_extended.dict'
	annot_dict <- read.table(dict_loc, stringsAsFactors=F)
	colnames(annot_dict) <- c('imchip', 'lz', 'rs', 'hg19_chr', 'hg19_pos','hg18_chr','hg18_pos','chr_band', 'FunctionGVS', 'weimin', 'FunctionDBSNP', 'GeneList', 'DistanceToSplice', 'AminoAcids')
	annot_dict <- annot_dict[,which(names(annot_dict) %in% c('imchip', 'FunctionGVS', 'GeneList', 'AminoAcids', 'FunctionDBSNP', 'DistanceToSplice'))]
	
	
	if (exon){ exonflag = 'exon'}
	else {exonflag = 'transcript'}
	outfolder = file.path(basefolder, exonflag)
	#filename = paste0(exonflag,'_cis_trans_merge.tbl')
	cistable_loc = file.path(outfolder, 'cis','cis_merge_4.45e-5.tbl')
	transtable_loc = file.path(outfolder, 'trans','trans_merge_9.05e-7.tbl')
	cistable = read.table(cistable_loc, T, sep='\t', stringsAsFactors=F)
	transtable = read.table(transtable_loc, T, sep='\t',stringsAsFactors=F)
	
	#FIND UNIQUE SNPS in CIS
	cis_snps <- unique(cistable$SNP)
	sub_transtable <- 	transtable[transtable$SNP %in% cis_snps,]
	
	
	#FIND UNIQUE SNPs IN TRANS
	trans_snps <- unique(transtable$SNP)
	sub_cistable <- cistable[cistable$SNP %in% trans_snps,]
	
	#ADD DICT INFO
	sub_cis_merge <- merge(sub_cistable, annot_dict, by.x='SNP_IM', by.y='imchip', all.x=TRUE, all.y=FALSE)
	sub_trans_merge <- merge(sub_transtable, annot_dict, by.x='SNP_IM', by.y='imchip', all.x=TRUE, all.y=FALSE)
	
	#WRITE TRANS ALSO IN CIS
	transoutloc = file.path(outfolder, paste0(exonflag,'_TRANS_overlap.tbl'))
	cisoutloc = file.path(outfolder, paste0(exonflag,'_CIS_overlap.tbl'))
	write_out(sub_trans_merge, transoutloc)
	write_out(sub_cis_merge, cisoutloc)
	
	

}



compare_pca_and_non <- function(basebasefolder='/home/jkb4y/cphgdesk_share/Achilleas', filter='', exon=FALSE, cis=TRUE, perm=FALSE){
	
	
	non_pcaCorrected_base <- file.path(basebasefolder, 'eQTLs_Feb2013')
	pcaCorrected_base <- paste0(non_pcaCorrected_base, '_pcaCorrected')
	
	
	non_pcaCorrected_loc <- locate_merge_out(type='merge',basefolder=non_pcaCorrected_base,exon=exon,cis=cis, filter=filter, perm=perm)
	pcaCorrected_loc <- locate_merge_out(type='merge',basefolder=pcaCorrected_base,exon=exon,cis=cis, filter=filter, perm=perm)
	pcacorr = substr(pcaCorrected_loc, 1, nchar(pcaCorrected_loc) - 4)
	nonpcacorr = substr(pcaCorrected_loc, 1, nchar(non_pcaCorrected_loc) - 4)
	out_loc = paste0(pcacorr , '_compared_to_nonPCAcorrected.txt')
	outfolder <- locate_merge_out(type='folder',basefolder=pcaCorrected_base, exon=exon,cis=cis,filter=filter,perm=perm)
	
	non_corrected <- read.table(non_pcaCorrected_loc, T, stringsAsFactors=F)
	corrected <- read.table(pcaCorrected_loc, T, stringsAsFactors=F)
	
	if (!exon){
		matched_cols <- c('SNP','GENE')
		#FIND SNP/GENE PAIRS in NON CORRECTED
		non_pairs <- with(non_corrected, paste0(SNP, GENE))
		#FIND SNP/GENE PAIRS in PCA CORRECTED
		pairs <- with(corrected, paste0(SNP, GENE))
		}
	else{
		matched_cols <- c('SNP','GENE','EXON')
		#FIND SNP/GENE PAIRS in NON CORRECTED
		non_pairs <- with(non_corrected, paste0(SNP, GENE, EXON))
		#FIND SNP/GENE PAIRS in PCA CORRECTED
		pairs <- with(corrected, paste0(SNP, GENE, EXON))
		}
	non_corrected$pairs <- non_pairs
	corrected$pairs <- pairs
	
	non_genes <- unique(non_corrected[,'GENE'])
	genes <- unique(corrected[,'GENE'])
	
	non_snps <- unique(non_corrected[,'SNP'])
	snps <- unique(corrected[,'SNP'])
	
	sink(out_loc)
	print('TRUE is the number of GENE/SNP pairs that can be found in both the pcaCorrected and the non_pcaCorrected merges, FALSE is the number of GENE/SNP pairs in the pcaCorrected merge that are NOT found in the non_pcaCorrected merge:')
	summary <- summary(pairs %in% non_pairs)
	names(summary)[names(summary)=="TRUE"] <- "GENE/SNP PAIRS IN BOTH"
	names(summary)[names(summary)=="FALSE"] <- "GENE/SNP PAIRS IN PCA CORRECTED ONLY"
	print(summary)
	

	print('TRUE is the number of GENE/SNP pairs that can be found in both the pcaCorrected and the non_pcaCorrected merges, FALSE is the number of GENE/SNP pairs in the non_pcaCorrected merge that are NOT found in the pcaCorrected merge:')
	
	summary <- summary(non_pairs %in% pairs)
	names(summary)[names(summary)=="TRUE"] <- "GENE/SNP PAIRS IN BOTH"
	names(summary)[names(summary)=="FALSE"] <- "GENE/SNP PAIRS IN NON CORRECTED ONLY"
	print(summary)
	
	print('TRUE is the number of GENES that can be found in both the pcaCorrected and the non_pcaCorrected merges, FALSE is the number of GENES in the pcaCorrected merge that are NOT found in the non_pcaCorrected merge:')
	print(summary(genes %in% non_genes))

	
	print('TRUE is the number of GENES that can be found in both the pcaCorrected and the non_pcaCorrected merges, FALSE is the number of GENES in the non_pcaCorrected merge that are NOT found in the pcaCorrected merge:')
	print(summary(non_genes %in% genes))

	print('TRUE is the number of individual SNPS that can be found in both the pcaCorrected and the non_pcaCorrected merges, FALSE is the number of individual SNPs in the pcaCorrected merge that are NOT found in the non_pcaCorrected merge:')
	print(summary(snps %in% non_snps))

	print('TRUE is the number of individual SNPS that can be found in both the pcaCorrected and the non_pcaCorrected merges, FALSE is the number of individual SNPs in the non_pcaCorrected merge that are NOT found in the pcaCorrected merge:')
	print(summary(non_snps %in% snps))	
	
	
	sink()
	
	#BUILD ANNOTATION DICTIONARY
	dict_loc = '/home/jkb4y/work/data/annot_hg19_extended.dict'
	annot_dict <- read.table(dict_loc, stringsAsFactors=F)
	colnames(annot_dict) <- c('imchip', 'lz', 'rs', 'hg19_chr', 'hg19_pos','hg18_chr','hg18_pos','chr_band', 'FunctionGVS', 'weimin', 'FunctionDBSNP', 'GeneList', 'DistanceToSplice', 'AminoAcids')
	annot_dict <- annot_dict[,which(names(annot_dict) %in% c('imchip', 'FunctionGVS', 'GeneList', 'AminoAcids', 'FunctionDBSNP', 'DistanceToSplice'))]
	
	
	#SUBSET THE TABLES WITH THE INTERSECTIONS
	sub_non_corrected <- non_corrected[!(non_corrected$pairs %in% pairs),]
	sub_corrected <- corrected[!(corrected$pairs %in% non_pairs),]
	
	#ADD DICT INFO
	sub_non_corrected <- merge(sub_non_corrected, annot_dict, by.x='SNP_IM', by.y='imchip', all.x=TRUE, all.y=FALSE)
	sub_corrected <- merge(sub_corrected, annot_dict, by.x='SNP_IM', by.y='imchip', all.x=TRUE, all.y=FALSE)
	
	#WRITE TRANS ALSO IN CIS
	correctedoutloc = file.path(outfolder, 'PCA_merge_unique_pairs.tbl')
	noncorrectedoutloc = file.path(outfolder,'nonPCA_merge_unique_pairs.tbl')
	write_out(sub_corrected, correctedoutloc)
	write_out(sub_non_corrected, noncorrectedoutloc)

}

do_all_the_things <- function(basefolder='/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected/', tablefolder='/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Feb2013_pcaCorrected', exon=FALSE, outfolder=NA, perm=FALSE, gwas_loc="/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected/data/gwascatalog_20130111.txt",r2_loc='/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/eurmeta_LD/all_regions_r2_0.ld',yank_loc='/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/RegionYank/eurmeta_yank.tbl'){
	#run merge tables for cis and trans
	cislist = c(TRUE, FALSE)
	for (cis in cislist){
		if (cis){ filter = '4.45e-5'}
		else {filter = '9.05e-7'}
		merge_tables(filter=filter, cis=cis, exon=exon, tablefolder=tablefolder,basefolder=basefolder,outfolder=outfolder, gwas_loc=gwas_loc,r2_loc=r2_loc,yank_loc=yank_loc,perm=perm)
		compare_pca_and_non(basebasefolder='/home/jkb4y/cphgdesk_share/Achilleas/', exon=exon, cis=cis, filter=filter)
	}
	compare_merges(exon=exon, basefolder=basefolder)
	merge_overlap(exon=exon, basefolder=basefolder)



}
