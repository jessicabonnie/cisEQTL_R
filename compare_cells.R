make_manageable_table <- function(file_loc, cell_type, wanted_cols=c('P','A1','BETA'), merge_cols=c('SNP','GENE'),filter='1e-4', cis=TRUE){
	#cell_table <- read.table(file_loc, T)
	source('gene_lists.R')
	cell_table <- read_and_filter(file_loc,as.numeric(filter),perm=FALSE, cis)
	cell_table$SNP <- cell_table$SNP_AP
	cell_table$SNP_JB <- cell_table$SNP_LZ_JB
	print(colnames(cell_table))
	cell_table <- cell_table[ , c(merge_cols,wanted_cols)]
	for (name in wanted_cols){
		new_name <- paste0(cell_type, '.',name)
		names(cell_table)[which(names(cell_table)==name)] <- new_name
	}
	print(head(cell_table))
	return(cell_table)
}

big_fish_eats_little <- function(total_table, cell_table){
	merge_by <- c('SNP','GENE')
	total_table <- merge(total_table, cell_table, by=merge_by, all=TRUE)




}

locate_table <- function(cell, perm, exon, cis){
	tablefolder = '/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Oct2012'
	filter_tail = paste0('_filtered_','1e-4','_JB','.txt')
	permflag = ''
	if (perm){ permflag = '_perm'}
	if (exon){
		filename = paste0(cell,'_exon_cis_eqtl',permflag,'_results',filter_tail)
	}
	else if (! cis){
	filename = paste0(cell,'_all_data_JB.txt')
	}
	else {
	filename = paste0(cell,'_cis_transcript_eqtls',permflag,filter_tail)
	}
	file_loc = file.path(tablefolder, filename)
	return(file_loc)

}
locate_merge_out <- function(perm, exon, cis, filter, type='merge', outfolder=NA){
if (is.na(outfolder)){
	#basefolder = '/home/jkb4y/ubs/work/results/Achilleas/eQTLs_Oct2012'
	basefolder= '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Oct2012'
	if (!cis){cisfolder = 'trans'; cis_flag = 'trans_'}
	else{cisfolder = 'cis'; cis_flag = 'cis_'}
	if (! cis){ cis_flag = 'trans_'}

	if (exon){transfolder = "exon"}
	else{transfolder = "transcript"}
	if (perm){permflag <- "_perm"}
	else{permflag <- ""}
	outfolder <- file.path(basefolder, cisfolder, transfolder)
	}
	if (type=='merge'){
	out_loc <- file.path(outfolder, paste0(cis_flag,'merge_',filter,'.tbl'))
	}
	else if (type=='thin'){out_loc <- file.path(outfolder, paste0(cis_flag,'thin_merge_',filter,'.tbl'))
	}
	else if (type =='multigene'){out_loc <- file.path(outfolder, paste0(cis_flag,'multigene_region_merge_',filter,'.tbl'))
	}
	else if (type =='countsum'){out_loc <- file.path(outfolder, paste0(cis_flag,'count_summary_',filter,'.tbl'))
	}
	else if (type =='countcross'){out_loc <- file.path(outfolder, paste0(cis_flag,'count_crosstab_',filter,'.tbl'))
	}
	print(out_loc)
	return(out_loc)

}


write_out <-function(df, out_loc){
	#out_loc <- locate_merge_out(perm=perm,exon=exon,cis=cis,filter=filter,type=type, outfolder=outfolder)
	write.table(df, file = out_loc, append = FALSE, quote = FALSE, sep = "\t",
            	eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)


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

merge_tables <- function(outfolder=NA, perm=FALSE, exon=FALSE, filter='1e-4',region_loc='/home/jkb4y/work/data/Region_Lists/hg19/achilleas_all_09202012.txt',cis=TRUE, si_loc='/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Oct2012/data/si_SNP_NoMHC_20121024_Query.txt', gwas_loc="/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Oct2012/data/gwascatalog_11262012.txt",r2_loc='/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/eurmeta_LD/all_regions_r2_0.ld',yank_loc='/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/RegionYank/eurmeta_yank.tbl'){
	library(plyr)
	cell_types <- c('B','CD4','CD8','MONO','NK')
	#merge_cols <- c('SNP','GENE','CHR','BP', 'SNP_JB')
	#if (exon){ merge_cols <- c(merge_cols,'EXON')}
	wanted_cols <- c('P','BETA','A1')
	if (! exon){
		merge_cols <- c('SNP','GENE','CHR','BP', 'SNP_JB')
		mamatable <- data.frame(SNP = character(0), GENE = character(0),
			CHR = character(0), BP = character(0),SNP_JB = character(0))
		}
	else{
		merge_cols <- c('SNP','GENE','EXON','CHR','BP', 'SNP_JB')
		mamatable <- data.frame(SNP = character(0), GENE = character(0),
			EXON = character(0), CHR = character(0), BP = character(0),SNP_JB = character(0))
	
		}
	for (cell_type in cell_types){
		table_loc <- locate_table(cell_type, perm, exon, cis)
		print(table_loc)
		cell_table <- make_manageable_table(table_loc, cell_type, wanted_cols=wanted_cols,merge_cols=merge_cols, filter, cis)
		mamatable <-  merge(mamatable, cell_table, by=merge_cols, all=TRUE, all.x=TRUE, all.y=TRUE)
		cell_table <- NULL
	}
	
	
	mamatable <- add_regions(mamatable, region_loc)
	mamatable <- add_SI_value(mamatable, si_loc)
	mamatable <- add_r2_info(mamatable, r2_loc, yank_loc)
	#mamatable <- add_literature_info(mamatable, gwas_loc)
	
	
	sumtable <- summarize_mama(mamatable)
	sum_loc <- locate_merge_out(perm=perm,exon=exon,cis=cis,filter=filter,type='countsum', outfolder=outfolder)
	print(sum_loc)
	write.table(sumtable,file=sum_loc, row.names=FALSE, col.names=TRUE, sep='\t',quote=FALSE)
	sumtable <- NULL
	
	crosstab <- crosstab_mama(mamatable)
	cross_loc <- locate_merge_out(perm=perm,exon=exon,cis=cis,filter=filter,type='countcross', outfolder=outfolder)
	print(cross_loc)
	write.table(crosstab,file=cross_loc, row.names=FALSE, col.names=TRUE, sep='\t',quote=FALSE)
	crosstab <- NULL
	
	print(head(mamatable))
	mama_loc <- locate_merge_out(perm=perm,exon=exon,cis=cis,filter=filter,type='merge', outfolder=outfolder)
	write_out(mamatable, mama_loc)
	multigenes <- find_multiregs(mamatable)
	
	#multigenes <- add_regions(regions, multigenes)
	multi_loc <- locate_merge_out(perm=perm,exon=exon,cis=cis,filter=filter,type='multigene', outfolder=outfolder)
	write_out(multigenes, multi_loc)
	multigenes <- NULL
	
	minimama <- thin_out(mamatable, cell_types)
	thin_loc <- locate_merge_out(perm=perm,exon=exon,cis=cis,filter=filter,type='thin', outfolder=outfolder)
	write_out(minimama, out_loc=thin_loc)
	minimama <- NULL
	return(mamatable)
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
	basefolder <- file.path(basefolder, cisfolder, transfolder, cell, 'RegionYank')
	filename <- paste0(cell, '_R_yank.tbl')
	file_loc <- file.path(basefolder, filename)
	return(file_loc)
	}

}

compare <- function(outfolder, basefolder, perm=FALSE, exon=FALSE){
	cells <-c('B','CD4','CD8','MONO','NK')
	subcells <- cells
	for (cell in cells){
		cell_loc <- locate_yank(basefolder, cell, perm, exon)
		cell_table <- read.table(cell_loc, T)
		subcells <- subcells[cells != cell]
		for (subcell in subcells){
			subcell_loc <- locate_yank(basefolder, subcell, perm, exon)
			subcell_table <- read.table(subcell_loc, T)
			
			compare_pair(cell_yank, subcell_yank)
		
		}
	
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
