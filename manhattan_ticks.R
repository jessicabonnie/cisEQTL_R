
# Original manhattan code forked from:
# Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Last updated: Tuesday, April 19, 2011
# R code for making manhattan plots and QQ plots from plink output files.
# manhattan() with GWAS data this can take a lot of memory, recommended for use on 64bit machines only, for now.
# Altnernatively, use bmanhattan() , i.e., base manhattan. uses base graphics. way faster.


## This is for testing purposes.
# set.seed(42)
# nchr=23
# nsnps=1000
# d=data.frame(
# SNP=sapply(1:(nchr*nsnps), function(x) paste("rs",x,sep='')),
# CHR=rep(1:nchr,each=nsnps),
# BP=rep(1:nsnps,nchr),
# P=runif(nchr*nsnps)
# )
# annotatesnps <- d$SNP[7550:7750]


#plot all of Achilleas cis-eQTL data files
achilleas_all <- function(flag='', file_type='jpeg', ymin=3, onepage=FALSE, filetype='jpeg', ymax=15, genomewideline=T, pca=F, exon=NULL, cis=NULL, ...){
	pcaflag = ''
	if (pca){
		pcaflag = '_pcaCorrected'
		}
 	datafolder = paste0("/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Feb2013", pcaflag)
 	outfolder = paste0("/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013", pcaflag)
	draw = TRUE
	both = c(TRUE,FALSE)
	if (is.null(cis)){
		cislist = both}
	else {cislist = c(cis)}
	if (is.null(exon)){
		exonlist = both}
	else {exonlist = c(exon)}
	
	#if (onepage){
	#	draw=FALSE
		#base_out_loc =file.path("/home/jkb4y/ubs/work/results/Achilleas","cis_trans_Manhattans")
	#	base_out_loc =file.path("/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013/","Manhattans","all_Manhattans")
	#	if (filetype== 'jpeg'){
	#		jpg_out_loc = paste0(base_out_loc,".jpeg")
	#		jpeg(jpg_out_loc, width=1100, height=850, quality=100, title="Manhattan Plots")
	#		}
	#	else if (filetype == 'ps'){
	#		ps_out_loc = paste0(base_out_loc,".ps")
	#		postscript(file=ps_out_loc, paper="letter", horizontal=F, title="Manhattan Plots")
	#		}
	#	par(mfcol=c(5,2),mar=c(3.1,4.1,2.1,2.1))
	#}
cisp = -log10(.05/1125)
transp = -log10(.05/55336)
	for (exon in exonlist){
		if (onepage){
		if (exon){title_flag = 'Exon'}
		else {title_flag = 'Transcript'}
		draw=FALSE
		#base_out_loc =file.path("/home/jkb4y/ubs/work/results/Achilleas","cis_trans_Manhattans")
		base_out_loc =file.path(outfolder,"Manhattans","all_Manhattans")
		if (filetype== 'jpeg'){
			jpg_out_loc = paste0(base_out_loc,title_flag,".jpeg")
			jpeg(jpg_out_loc, width=1100, height=850, quality=100, title=paste(title_flag,"Manhattan Plots"))
			}
		else if (filetype == 'ps'){
			ps_out_loc = paste0(base_out_loc,title_flag,".ps")
			postscript(file=ps_out_loc, paper="letter", horizontal=F, title=paste(title_flag,"Manhattan Plots"))
			}
		par(mfcol=c(5,2),mar=c(3.1,4.1,2.1,2.1))
	}
		for (cell_type in c("B","CD4","CD8","MONO","NK")){
		for (cis in cislist){
			if (! exon){ ymin = 2}
			else {ymin = 3}
				if (genomewideline!=F){
					if (cis){
					genomewideline=cisp}
					else {genomewideline=transp}
					}
				print(ymin)
				achilleas_eqtl_draw(cell_type=cell_type,exon=exon, cis=cis, flag=flag, ymax=ymax, draw=draw, ymin=ymin, datafolder=datafolder, genomewideline=genomewideline, outfolder=outfolder,...)
	}
	}
	
	if (onepage){
		dev.off()
		}
	}

}

# plot Achilleas' EQTL data, taking a perm flag if it is a perm file
achilleas_eqtl_draw <-function(cell_type = NULL, exon = FALSE, cis = TRUE, flag = NULL, file_type = 'jpeg',colors=c("navy","darkorange"), draw=TRUE, ymin=0, datafolder= NA, genomewideline=T,outfolder=NA,...){

filter_label = ''
if (exon){
	exonflag = 'exon'
	exonlabel = 'Exon'
	filter_label = '_filtered_1e-04'
	}
else {
	exonflag = 'transcript'
	exonlabel = 'Transcript'
	filter_label = '_filtered_1e-02'
	}
if (cis){
	cislabel = 'Cis'
	cisflag = 'cis'
	#filter_label = ''
	}
else{
	cislabel = 'Trans'
	cisflag = 'trans'
	}


data_name = paste0(exonflag,cislabel,cell_type, filter_label,'_manhattan.txt')

#if (exon==TRUE){
#	arrayflag = 'exon'
#	data_name = paste0(cell_type,'_exon_cis_eqtl',permflag,'_results_noNA_JB.txt')
#	arraytitle = 'Exon'
#	}
#data_name = paste0(cell_type,'_cis_transcript_eqtls',permflag,'_all_new_JB.txt')
data_loc = file.path(datafolder, data_name)
print(data_loc)

#out_folder =file.path("/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Jan2013",exonflag,cisflag,cell_type,"Manhattans")
out_folder = file.path(outfolder,"Manhattans_test")
if (! file.exists(out_folder)){
		dir.create(out_folder)
		}
out_name_base = paste0(cell_type,"_",exonflag,"_",cisflag,"_eqtls")
pdf_out_name = paste0(out_name_base,'.pdf')
pdf_out_loc = file.path(out_folder, pdf_out_name)
png_out_loc = file.path(out_folder, paste0(out_name_base,'.png'))
ps_out_loc = file.path(out_folder, paste0(out_name_base,'.ps'))
jpg_out_loc = file.path(out_folder, paste0(out_name_base,'.jpeg'))

#out_loc = file.path(out_folder, out_name)
	print('Draw is')
	print(draw)
if (draw){	print("Draw is")
	print(draw)
	}
title = paste(cell_type,"Cells:",exonlabel,cislabel,"eQTLs")




results <- read.table(data_loc, header=T, stringsAsFactors=F)
#if (! cis & ! exon){
#if (genomewideline == F){
#threshold = 10 ^ (-1*(ymin +1))
#}
#else {threshold = 10 ^ (-1*(genomewideline -1))}

#results <- alter_table(results, threshold=threshold)
#print(head(results))
#}
#if (! exon){
#results <- read.table(data_loc, T)
#}
#else{results <- read.table(data_loc, T)}
#print(tail(results))
#postscript(file=ps_out_loc, paper="letter", horizontal=T, title=title)
#pdf(out_loc, title=title, width=11, height=8.5) #,paper="letter", horizontal=T)
#png(png_out_loc, quality=100)
if (file_type == 'ps' & draw){
	#print(tail(results))
	postscript(file=ps_out_loc, paper="letter", horizontal=T, title=title)
	}
else if (file_type == 'pdf' & draw){
	pdf(pdf_out_loc, title=title, width=11, height=8.5) #,paper="letter", horizontal=T)
	}
else if (file_type == 'png' & draw) {
	png(png_out_loc, quality=100)
	}
else if (file_type == 'jpeg' & draw){
	jpeg(jpg_out_loc, width=1100, height=850, quality=100)
	}
	print(title)
manhattan(results, colors=colors, main=title, suggestiveline=F, genomewideline=genomewideline, pch=16, limitchromosomes=1:23,ymin=ymin,...)
if (draw){
	dev.off()
	}

}


# plot Achilleas' trans-EQTL data
achilleas_trans_acting <-function(cell_type = NULL, flag = NULL, file_type = 'jpeg',colors=c("navy","darkorange"), genomewideline=F,ymin=0,draw=TRUE, data_folder,...){
#data_folder = "/home/jkb4y/ubs/work/data/Achilleas/transcript_eqtl_results_Oct12"

data_name = paste0(cell_type,'')

data_loc = file.path(data_folder, data_name)
print(data_loc)

out_folder =file.path("/home/jkb4y/ubs/work/results/Achilleas/transcript_eqtl_results_Oct12",cell_type,"Manhattans")
if (! file.exists(out_folder)){
		dir.create(out_folder)
		}
if (is.null(flag)) { addflag=''}
else{ addflag = paste0('_',flag)}
pdf_out_name = paste0(cell_type, "_trans_eqtls",addflag,".pdf")
pdf_out_loc = file.path(out_folder, pdf_out_name)
png_out_loc = file.path(out_folder, paste0(cell_type, "_trans_eqtls",addflag,".png"))
ps_out_loc = file.path(out_folder, paste0(cell_type, "_trans_eqtls",addflag,".ps"))
jpg_out_loc = file.path(out_folder, paste0(cell_type, "_trans_eqtls",addflag,".jpeg"))


title = paste(cell_type,"Cells: Trans eQTLs")


results <- read.table(data_loc, T)
results <- results[results$TRANS.CIS == "trans", ]
print(head(results))
if (file_type == 'ps' & draw){
#print(tail(results))
postscript(file=ps_out_loc, paper="letter", horizontal=T, title=title)
}
else if (file_type == 'pdf' & draw){
pdf(pdf_out_loc, title=title, width=11, height=8.5) #,paper="letter", horizontal=T)
}
else if (file_type == 'png' & draw) {
png(png_out_loc, quality=100)
}
else if (file_type == 'jpeg' & draw){
jpeg(jpg_out_loc, width=1100, height=850, quality=100)
}
manhattan(results, colors=colors, main=title, suggestiveline=F, genomewideline=genomewideline, pch=16, limitchromosomes=1:23,ymin=ymin, ...)
if (draw){
	dev.off()
	}

}

#function to alter the table to be drawable within the human lifespan
alter_table <- function(df, threshold){
library('plyr')
	snp.list <- unique(df$SNP)
	subtable <- ldply(snp.list,.fun=function(x){alter_helper(as.character(x), df, threshold)})
	return(subtable)

}

#helper function
alter_helper <-function(snp, df, threshold){
	
	sub = df[as.character(df$SNP)==as.character(snp), ]
	#print(nrow(sub))
	top = as.numeric(sub$P[order(sub$P)[1]])
	bottom = as.numeric(sub$P[order(sub$P)[nrow(sub)]])
	sub = sub[as.numeric(sub$P) < as.numeric(threshold) | as.numeric(sub$P) == top,]
	#print(sub)
	#print(top)
	#print(nrow(sub))
	#print(threshold)
	return(sub)

}


# manhattan plot using base graphics
manhattan <- function(dataframe, colors=c("navy","darkorange"), ymax="max", limitchromosomes=1:23,
 suggestiveline=F, genomewideline=F, annotate=NULL, ymin=0,...) {
#suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, ...) {
    d=dataframe
    #print(names(d))
    if (!("CHR" %in% names(d) & "BP" %in% names(d) )) stop("Make sure your data frame contains columns CHR and BP")
    if (!("P" %in% names(d) | "EMP1" %in% names(d))) stop("Make sure your data frame contains EITHER column P or column EMP1")
    if("EMP1" %in% names(d)) d$P = d$EMP1;
    
    if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
    if (ymin == 0){
    pmax = 1}
    else{
    pmax = 10^-(ymin)
    }
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=pmax)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    ticks=NULL
    lastbase=0
    colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
    if (ymax=="max") ymax<-ceiling(max(d$logp))
    if (ymax<8) ymax<-8
    
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
    } else {
        for (i in unique(d$CHR)) {
          if (i==1) {
     d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
     } else {
     
     lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
     #print(lastbase)
     d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
     }
     ordered = order(d[d$CHR==i,]$BP)
     len = length(ordered)
	 first <- as.numeric(d[d$CHR==i,]$pos[ ordered[1] ]);
     last <- as.numeric(d[d$CHR==i,]$pos[ ordered[len]] );
     mid = (last - first)/2
     middle = first + mid
     #print(c(first, last, middle))
     ordered = NULL
     ticks=c(ticks, middle)
     #ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
     #print(ticks)
     }
    }
    
    if (numchroms==1) {
        #with(d, plot(pos, logp, ylim=c(ymin,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
        with(d, plot(pos, logp, ylim=c(ymin,ymax), ylab=expression(-log[10](p)), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    }	else {
        #with(d, plot(pos, logp, ylim=c(ymin,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
        with(d, plot(pos, logp, ylim=c(ymin,ymax), ylab=expression(-log[10](p)), xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=unique(d$CHR), ...)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
            icol=icol+1
     }
    }
    
    if (!is.null(annotate)) {
        d.annotate=d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, logp, col="green3", ...))
    }
    
    if (suggestiveline) abline(h=suggestiveline, col="blue")
    if (genomewideline) abline(h=genomewideline, col="red")
}


## Make a pretty QQ plot of p-values
qq = function(pvector, ...) {
if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
o = -log10(sort(pvector,decreasing=F))
e = -log10( ppoints(length(pvector) ))
plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
abline(0,1,col="red")
}

## Make a pretty QQ plot of p-values
qq_jb = function(pvector, ...) {
if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
o = -log10(sort(pvector,decreasing=F))
e = -log10( ppoints(length(pvector) ))
lim_max = ceiling(max(c(e,o)))
plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](p)), ylab=expression(Observed~~-log[10](p)), xlim=c(0,lim_max), ylim=c(0,lim_max),xaxs='i',yaxs='i', ...)
abline(0,1,col="red")
}
### OLD GGPLOT2 CODE ###

# manhattan plot using ggplot2
gg.manhattan = function(dataframe, title=NULL, max.y="max", suggestiveline=0, genomewideline=-log10(5e-8), size.x.labels=9, size.y.labels=10, annotate=F, SNPlist=NULL) {
library(ggplot2)
    if (annotate & is.null(SNPlist)) stop("You requested annotation but provided no SNPlist!")
d=dataframe
#limit to only chrs 1-23?
d=d[d$CHR %in% 1:23, ]
if ("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d) ) {
d=na.omit(d)
d=d[d$P>0 & d$P<=1, ]
d$logp = -log10(d$P)
d$pos=NA
ticks=NULL
lastbase=0
#new 2010-05-10
numchroms=length(unique(d$CHR))
if (numchroms==1) {
d$pos=d$BP
} else {

for (i in unique(d$CHR)) {
if (i==1) {
d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
}	else {
lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
}
ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
}
ticklim=c(min(d$pos),max(d$pos))

}
mycols=rep(c("gray10","gray60"),max(d$CHR))
if (max.y=="max") maxy=ceiling(max(d$logp)) else maxy=max.y
if (maxy<8) maxy=8
if (annotate) d.annotate=d[as.numeric(substr(d$SNP,3,100)) %in% SNPlist, ]
if (numchroms==1) {
plot=qplot(pos,logp,data=d,ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"))
}	else {
plot=qplot(pos,logp,data=d, ylab=expression(-log[10](italic(p))) , colour=factor(CHR))
plot=plot+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(d$CHR)))
plot=plot+scale_y_continuous(limits=c(0,maxy), breaks=1:maxy, labels=1:maxy)
plot=plot+scale_colour_manual(value=mycols)
}
if (annotate) plot=plot + geom_point(data=d.annotate, colour=I("green3"))
plot=plot + opts(legend.position = "none")
plot=plot + opts(title=title)
plot=plot+opts(
panel.background=theme_blank(),
panel.grid.minor=theme_blank(),
axis.text.x=theme_text(size=size.x.labels, colour="grey50"),
axis.text.y=theme_text(size=size.y.labels, colour="grey50"),
axis.ticks=theme_segment(colour=NA)
)
if (suggestiveline) plot=plot+geom_hline(yintercept=suggestiveline,colour="blue", alpha=I(1/3))
if (genomewideline) plot=plot+geom_hline(yintercept=genomewideline,colour="red")
plot
}	else {
stop("Make sure your data frame contains columns CHR, BP, and P")
}
}

gg.qq = function(pvector, title=NULL, spartan=F) {
library(ggplot2)
o = -log10(sort(pvector,decreasing=F))
#e = -log10( 1:length(o)/length(o) )
e = -log10( ppoints(length(pvector) ))
plot=qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + stat_abline(intercept=0,slope=1, col="red")
plot=plot+opts(title=title)
plot=plot+scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
plot=plot+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
if (spartan) plot=plot+opts(panel.background=theme_rect(col="grey50"), panel.grid.minor=theme_blank())
plot
}

gg.qqman = function(data="plinkresults") {
myqqplot = ggqq(data$P)
mymanplot = ggmanhattan(data)
ggsave(file="qqplot.png",myqqplot,w=5,h=5,dpi=100)
ggsave(file="manhattan.png",mymanplot,width=12,height=9,dpi=100)
}

gg.qqmanall= function(command="ls *assoc") {
filelist=system(command,intern=T)
datalist=NULL
for (i in filelist) {datalist[[i]]=read.table(i,T)}
highestneglogp=ceiling(max(sapply(datalist, function(df) max(na.omit(-log10(df$P))))))
print(paste("Highest -log10(P) = ",highestneglogp),quote=F)
start=Sys.time()
for (i in names(datalist)) {
myqqplot=ggqq(datalist[[i]]$P, title=i)
ggsave(file=paste("qqplot-", i, ".png", sep=""),myqqplot, width=5, height=5,dpi=100)
mymanplot=ggmanhattan(datalist[[i]], title=i, max.y=highestneglogp)
ggsave(file=paste("manhattan-", i, ".png", sep=""),mymanplot,width=12,height=9,dpi=100)
}
end=Sys.time()
print(elapsed<-end-start)
}
