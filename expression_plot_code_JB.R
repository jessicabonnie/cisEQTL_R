
read_genotype_table <- function(genotype_loc, graph_table){
	
	genotypes = read.table(genotype_loc, T)
	genotypes <- genotypes[,colnames(genotypes) %in% colnames(genotypes[1:4]) | colnames(genotypes) %in% graph_table$SNP]
	print(colnames(genotypes))
	return(genotypes)
	}
main <- function(outfolder="/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013",cis=TRUE,exon=FALSE, genotypes=NULL, summary.only=TRUE){
	basefolder <- "/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013"
	datafolder <- file.path(basefolder, 'data')
	alleles_loc <- file.path(datafolder, "MmAlleles.txt")
	
	pat_table_loc <- file.path(datafolder, "Table3a_IC_V14_20121118.aq_JB.txt")
	pat_table <- read.table(pat_table_loc, T, sep='\t', stringsAsFactors=FALSE)
	print(paste("Pat's table is ",nrow(pat_table),"rows."))
	suna_table_loc <- file.path(datafolder, "trans_eQTL_SI-SNPs.txt")
	suna_table <- read.table(suna_table_loc, T, sep='\t', stringsAsFactors=FALSE)

	print(paste("Suna's first table is ",nrow(suna_table),"rows."))
	genotype_loc <- file.path(datafolder,"talk_genotypes.txt")
	suna_table2_loc <- file.path(datafolder, "cis_eQTL_SNP_list_expression_plot_group2.txt")
	suna_table2 <- read.table(suna_table2_loc, T, sep='\t', stringsAsFactors=FALSE)
	print(paste("Suna's second table is ",nrow(suna_table2),"rows."))
	
	if (exon){
		exonflag = 'exon'
		#pmax = "9.05e-7"
		thinflag = "_thin"
		#thinflag = ""
		}
	else{ 
		exonflag = 'transcript'
		#pmax = "4.45e-5"
		thinflag = ""}
	basefolder <- file.path(basefolder, exonflag)
	if (cis){
		cisflag <- "cis"
		pmax = "4.45e-5"
		}
	else {
		cisflag <- "trans"
		pmax = "9.05e-7"
		}
	basefolder <- file.path(basefolder, cisflag)
	merge_loc <- file.path(basefolder,paste0(cisflag,thinflag,"_merge_",pmax,"_lit.tbl"))
	print(merge_loc)
	outfolder <- file.path(basefolder, 'expression_plots')
	folder_exists(outfolder)
	merge_table <- read.delim(merge_loc, T, sep='\t', stringsAsFactors=FALSE)
	print("read merge_table")

	merge_table <- merge_table[
					!is.na(merge_table$GWAS_DISEASES) | 
					!is.na(merge_table$EXTRA_DISEASES) | 
					as.character(merge_table$SI) == 'TRUE'|as.character(merge_table$SNP) %in% as.character(pat_table$BestSNP) |
					as.character(merge_table$SNP) %in% as.character(suna_table$SNP) |as.character(merge_table$SNP) %in% as.character(suna_table2$SNP), which(colnames(merge_table) %in% c('SNP','GENE','REGION','B.P','CD4.P','CD8.P','MONO.P','NK.P'))]
	#merge_table <- suna_table[ which(colnames(suna_table) %in% c('SNP','GENE','REGION','B.P','CD4.P','CD8.P','MONO.P','NK.P'))]
	colnames(merge_table) <- sub(".P","",colnames(merge_table), fixed = TRUE)
	#print(colnames(merge_table))
	print("These are the SNPs from Pat's list that are in the table:")
	print(merge_table[as.character(merge_table$SNP) %in% as.character(pat_table$BestSNP),])
	print("These are the SNPs from Suna's lists that are in the table:")
	print(merge_table[as.character(merge_table$SNP) %in% as.character(suna_table$SNP),])
	print(merge_table[as.character(merge_table$SNP) %in% as.character(suna_table2$SNP),])
	
	
	if (is.null(genotypes)){
		genotypes = read.table(genotype_loc, T)
		#print(colnames(genotypes))
	}
	
	genotype_sub <- genotypes[,colnames(genotypes) %in% colnames(genotypes[1:4]) | colnames(genotypes) %in% merge_table$SNP]
	#print(genotype_sub)
	allele_table <-read.table(alleles_loc, header=T, sep=' ')
	print("read allele_table")
	allele_table <- allele_table[allele_table$SNP %in% merge_table$SNP,]
	
	do_cells(genotypes=genotype_sub, merge_table=merge_table, allele_table=allele_table, outfolder=outfolder)
	do(genotypes=genotype_sub, graph_table=merge_table, allele_table=allele_table,outfolder=outfolder, summary.only=summary.only)
	
	
	}

do_cells <-function(genotypes, merge_table, allele_table, outfolder){
	pdf.options(width=11, height=8.5)
	#graph_helper <- "/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Jan2013/data/eQTL_graph_helper.txt"
	#outfolder <- "/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Jan2013/transcript/expression_plots/"
	summaryfolder <-file.path(outfolder, 'summary')
	folder_exists(summaryfolder)
	
	#graph_table <- read.table(graph_helper, T)
	
	#alleles_loc <- file.path(outfolder, "MmAlleles.txt")
	#allele_table = read.table(alleles_loc,T)
	#allele_table <- allele_table[allele_table$SNP %in% graph_table$SNP,]
	
	cells <- list('B','CD4','CD8','MONO','NK')
	snps <- unique(as.character(merge_table$SNP))
	
	lapply(snps, per_snp_across_cells, summaryfolder=summaryfolder, graph_table=merge_table, genotypes=genotypes, allele_table=allele_table,cells=cells)
	
	}


determine_par <-  function(sub_graph,cells){
	columns = nrow(sub_graph)
	rows = 0
	for (cell in cells){
		if (length(na.omit(unique(sub_graph[[cell]]))) > 0){
			rows = rows + 1
			}
		}
	print(columns)
	print(rows)
	if (columns > 4){
		columns = 4
		}
	par=c(rows,columns)
	return(par)
}

per_snp_across_cells <- function(snp, summaryfolder,graph_table, genotypes, allele_table, cells){
	sub_graph <- graph_table[graph_table$SNP == snp,]
	region = as.character(sub_graph$REGION[1])
	print(region)
	region_folder = file.path(summaryfolder, region)
	folder_exists(region_folder)
	pdf.name = file.path(region_folder, paste0(snp, '_AllCells.pdf'))
	print(pdf.name)
	pdf(pdf.name)
	#par(mfrow=c(2,3))
	par(mfrow=determine_par(sub_graph,cells))
	for (cell in cells){
		load(paste0("/m/cphg-expr1/cphg-expr1/temp_data/data/BRI_eqtl_data/raw/", cell,".RData"))
		cell.data <- get(cell)
		cell.order = match(genotype.ids, genotypes$IID)
		draw_per_snp(snp,graph_table, cell, cell.order,cell.data, genotypes, allele_table, across_cells=TRUE)
		}
	dev.off()
	
	
	
	}


do <- function(genotypes, graph_table,allele_table,outfolder, summary.only=TRUE){
	pdf.options(width=11, height=8.5)
	library('plyr')
#graph_helper <- "/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Jan2013/data/eQTL_graph_helper.txt"
#outfolder <- "/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Jan2013/transcript/expression_plots/"
#graph_table <- read.table(graph_helper, T)

#alleles_loc <- file.path(outfolder, "MmAlleles.txt")
#allele_table = read.table(alleles_loc,T)
#allele_table <- allele_table[allele_table$SNP %in% graph_table$SNP,]

	cells <- list('B','CD4','CD8','MONO','NK')
	
	for (cell in cells){
		load(paste0("/m/cphg-expr1/cphg-expr1/temp_data/data/BRI_eqtl_data/raw/", cell,".RData"))
		cellfolder = file.path(outfolder, cell)
		if (! file.exists(cellfolder)){
			dir.create(cellfolder)
			}
		cell.data <- get(cell)
		cell.order = match(genotype.ids, genotypes$IID)
		if (!summary.only){
			adply(graph_table, 1, draw_individual_pdf, cell=cell, cell.order=cell.order, cellfolder=cellfolder, cell.data=cell.data, genotypes=genotypes, allele_table)
			}
		lapply(unique(as.character(graph_table$SNP)), draw_per_snp_pdf, graph_table, cell=cell, cell.order=cell.order, outfolder=outfolder, cell.data=cell.data, genotypes=genotypes, allele_table)
		print(paste("Finished cell: ",cell))
	
		}


	}


find_alleles <- function(allele_table, snp){
	allele_line <- allele_table[allele_table$SNP == snp,]
	print(allele_line)

	alleles <- list(as.character(allele_line$A1[1]),as.character(allele_line$A2)[1])
	return(alleles)
	}

draw_per_snp_pdf <-function(snp,graph_table, cell, cell.order,outfolder,cell.data, genotypes, allele_table){
	library('plyr')
	summaryfolder <- file.path(outfolder, 'summary')
	cellfolder <- file.path(outfolder, cell)
	folder_exists(cellfolder)
	folder_exists(summaryfolder)
	sub_graph <- graph_table[graph_table$SNP == snp,]
	if (length(na.omit(unique(sub_graph[[cell]]))) > 0){
		region = as.character(sub_graph$REGION[1])
		#region_folder = file.path(summaryfolder, region)
		region_folder = file.path(cellfolder, region)
		folder_exists(region_folder)
		pdf.name = file.path(region_folder, paste0(snp,'_',cell,'cells.pdf'))
		pdf(pdf.name)
		count = length(na.omit(sub_graph[[cell]]))
		if (count > 3){
			par(mfrow=c(2,2))
			}
		else{
			par(mfrow=c(1,count))
			}
		adply(sub_graph, 1, draw_individual_plot, cell=cell, cell.order=cell.order, cell.data=cell.data, genotypes=genotypes, allele_table)
		dev.off()
		}
	

}
draw_per_snp <-function(snp,graph_table, cell, cell.order,cell.data, genotypes, allele_table, across_cells=FALSE){
library('plyr')
	sub_graph <- graph_table[graph_table$SNP == snp,]
	if (length(na.omit(unique(sub_graph[[cell]]))) > 0){
	adply(sub_graph, 1, draw_individual_plot, cell=cell, cell.order=cell.order, cell.data=cell.data, genotypes=genotypes, allele_table, across_cells=across_cells)
	}

}

draw_gene_vs_gene_single <- function(gene1, gene2,cell){
	load(paste0("/m/cphg-expr1/cphg-expr1/temp_data/data/BRI_eqtl_data/raw/", cell,".RData"))
	cell.data <- get(cell)
	print(gene1)
	print(gene2)
	if (! (gene1 %in% row.names(cell.data) & gene2 %in% row.names(cell.data))){
		print(paste("Cannot draw ",gene1,"vs",gene2,"for",cell,"cells, because one or both is missing in this cell type."))
		}
	else{
		print(length(cell.data[gene1,]))
		print(length(cell.data[gene2,]))
		
		plot( cell.data[gene1,] ~ cell.data[gene2,],ylab = paste(gene1,"Expression"), xlab = paste(gene2,"Expression"), main = paste(cell,"Cells:",gene1,"vs",gene2,"Expression"), pch=19, col="purple")
		#points( cell.data[gene,] ~ jitter(as.factor(genotypes[cell.order,c(snp)])), pch=19, col="darkcyan")
		#points( cell.data[gene1,] ~ cell.data[gene2,], pch=19, col="purple")
		}

	}
draw_gene_vs_gene_multicell <- function(gene1, gene2){
	print(gene1)
	print(gene2)
	cells = c('B','CD4','CD8','MONO','NK')
	colors = c('navy','green','darkcyan','purple','red')
	graph <- data.frame(gene1=numeric(0),gene2=numeric(0),cell=(character(0)))
	counter = 0
	for (cell in cells){
	load(paste0("/m/cphg-expr1/cphg-expr1/temp_data/data/BRI_eqtl_data/raw/", cell,".RData"))
	cell.data <- get(cell)
	print(cell)
		if (! (gene1 %in% row.names(cell.data) & gene2 %in% row.names(cell.data))){
		print(paste("Cannot draw ",gene1,"vs",gene2,"for",cell,"cells, because one or both is missing in this cell type."))
		}
	else{
		df <-data.frame(gene1=cell.data[gene1,],gene2=cell.data[gene2,], cell=cell)
		print(head(df))
		graph <-rbind(graph,df)
		print(head(graph))
		}
		pch=19
		color.factor = as.factor(graph$cell)
		plot( graph$gene1 ~ graph$gene2,ylab = paste(gene1,"Expression"), xlab = paste(gene2,"Expression"), main = paste("All Cells:",gene1,"vs",gene2,"Expression"), pch=pch, col=colors[color.factor])
		legend(x="topleft", c(levels(graph$cell)), col=colors, 
       inset=0.05, bty="n", pch=pch)

	}
	}	
draw_gene_vs_gene_jpegs <- function(){
	library('plyr')
	outfolder = '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013/gene_expression_plots/'
	interest_table <- read.table('/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013/gene_expression_plots/TEST_cis-trans_correlation.txt',T, stringsAsFactors=F)
	do_things <-function(line, outfolder){
		gene1 = as.character(line$Cis_GENE)
		gene2 = as.character(line$Trans_GENE)
		cell = as.character(line$Cell_Type)
		#DRAW ONE CELL TYPE PLOT
		onecell_loc = file.path(outfolder, paste0(cell,'_',gene1,'_',gene2,'.jpeg'))
		jpeg(onecell_loc, width=650, height=425, quality=100)
		draw_gene_vs_gene_single(gene1=gene1,gene2=gene2,cell=cell)
		dev.off()
		#DRAW ALL CELLTYPES FOR THIS
		allcell_loc = file.path(outfolder, paste0(gene1,'_',gene2,'.jpeg'))
		jpeg(allcell_loc, width=650, height=425, quality=100)
		draw_gene_vs_gene_multicell(gene1=gene1,gene2=gene2)
		dev.off()
	
		}
	adply(interest_table, 1, do_things, outfolder)


	}
	
draw_plot <- function(snp, gene,cell,cell.order,cell.data, genotypes, allele_table){
	print(gene)
	print(snp)
	if (! gene %in% row.names(cell.data)){
		print("HERHERHEHREHEHREHRHERHE")

		}
	else{
		print(length(cell.data[gene,]))
		print(length(genotypes[cell.order,c(snp)]))
		
		alleles = find_alleles(allele_table,snp)
		allele = paste(alleles,collapse='/')
		
		plot( cell.data[gene,] ~ as.factor(genotypes[cell.order,c(snp)]),ylab = paste(gene,"expression"), xlab =paste0(snp), main = paste(cell,"Cells: Expression vs Genotype"), xaxt="n", outline=F)
		axis(1, at=1:3, lab=c(paste0(alleles[1],alleles[1]),paste0(alleles[1],alleles[2]),paste0(alleles[2],alleles[2])))
		#points( cell.data[gene,] ~ jitter(as.factor(genotypes[cell.order,c(snp)])), pch=19, col="darkcyan")
		points( cell.data[gene,] ~ jitter(as.numeric(as.factor(genotypes[cell.order,c(snp)])), .75), pch=19, col="darkcyan")
		}

	}

draw_individual_plot <- function(line,cell,cell.order,cell.data, genotypes, allele_table, across_cells=FALSE){
	#print(line)
	if (!is.na(line[[cell]]) | across_cells){
		gene = as.character(line$GENE)
		snp = as.character(line$SNP)
		draw_plot(snp, gene,cell,cell.order,cell.data, genotypes, allele_table)
		}

	}

draw_individual_pdf <- function(line,cell,cell.order,cellfolder,cell.data, genotypes, allele_table){
	if (!is.na(line[[cell]])){
		gene = as.character(line$GENE)
		snp = as.character(line$SNP)
		region = as.character(line$REGION[1])
		region_folder = file.path(cellfolder, region)
		folder_exists(region_folder)
		pdf.name = file.path(region_folder,paste0(snp,"_",gene,"_",cell,"cells.pdf"))
		
		pdf(pdf.name)
		draw_plot(snp, gene,cell,cell.order,cell.data, genotypes, allele_table)
		dev.off()
		}

	}

folder_exists <- function(folder){
 	if (! file.exists(folder)){
		dir.create(folder)
		}
	}
	
