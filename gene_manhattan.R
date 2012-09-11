# Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Last updated: Tuesday, April19, 2011
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




read_assoc <-function(address){

results <- read.table(address,T)
return(results)
}

achilleas <-function(data_loc,out_loc, title='Manhattan Plot',colors=c("navy","orange"),shapes=c(1,2,3,4,5,6,7,8,9,19,20,21,22,23,24), genomewideline=-log10(3.5e-7),...){
results <- read_assoc(data_loc)
print(tail(results))
postscript(file=out_loc, paper="letter", horizontal=T)
manhattan(results, colors=colors, shapes=shapes, main=title, suggestiveline=F, genomewideline=genomewideline)
dev.off()

}

# manhattan plot using base graphics
manhattan <- function(dataframe, colors=c("navy","orange","chartreuse4"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=-log10(3.5e-7),shapes=c(1,2,3,4,5,6,7,8,9,19,20,21,22,23,24), annotate=NULL, ...) {
	library(plyr)
    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "GENE" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, GENE and SNP")
    if (!("P" %in% names(d) | "EMP1" %in% names(d))) stop("Make sure your data frame contains EITHER column P or column EMP1")
    if("EMP1" %in% names(d)) d$P = d$EMP1; print("DONE")
    if("SNP_AP" %in% names(d)) d$SNP = d$SNP_AP; print("DONE")
    
    if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d=data.frame(indices=1:nrow(d),d)
    d$logp = -log10(d$P)
    print(head(d))
    d$pos=NA
    ticks=NULL
    lastbase=0
    colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
    shapes <-rep(shapes,length(unique(d$GENE)))[1:length(unique(d$GENE))]
    glabels <- ifelse(d$logp > 6,paste(d$GENE,'\n',d$SNP),'')
    


    shape.factor = as.factor(d$GENE)
    color.factor = as.factor(d$CHR)
    #ymax = ifelse(ymax=="max",ceiling(max(d$logp)),ymax)
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
     d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
     }
     ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
     }
    }
    
    if (numchroms==1) {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    } else {
    	#MaxP <- tapply(d$logp, d$GENE, max)
    	#MaxP <- tapply(seq_len(nrow(d)),d$SNP,function(i)d$GENE[i][which.max(d$logp[i])])
    	#print(MaxP)
		#names(MaxP) <- c("SNP","logp")
    	#MaxGene <- merge(d, MaxP)
    	blah1 = function(a){if (a$logp < 6){return(NA)}
    	else{ return(a$indices[which.max(a$logp)])}
    	}
		blah.func <- function(df){
    	df$V1 <- ifelse(df$logp < 5,NA,df$indices[which.max(df$logp)])
    	#print("working...")
    	return(df)}    	
    	blah = function(a){ifelse(a$logp < 6,return(NA),return(a$indices[which.max(a$logp)]))}
    	
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=unique(d$CHR), ...)
        with(d, points(pos,logp,col=colors[color.factor],pch=shapes[shape.factor],...))
        print("I'm in 1")
        index.frame <- ddply(d,.(GENE),blah.func)
        print("I'm in 2")
        index.frame = na.omit(index.frame)
        print(head(index.frame))
        print("I'm in 3")
        #print(d[index.frame$V1,])
        print("I'm in 4")
        d.genes = d[index.frame$V1,]
        print(head(d.genes))
        print("I'm in 5")
        with(d.genes, text(x=pos, y = logp+.5,labels = paste(SNP,'\n',GENE), adj = NULL, col="purple",vfont = NULL,font = NULL))
      #  g <-subset(d,with(d,tapply(seq_len(nrow(d)),SNP,function(i)GENE[i][which.max(logp[i])])))
        #with(d,tapply(seq_len(nrow(d)),SNP,function(i)GENE[i][which.max(logp[i])]),text(x=pos, y = logp+.5,labels = paste(SNP,'\n',GENE), adj = NULL, col="red",vfont = NULL,font = NULL))
        #with(d,	tapply(seq_len(nrow(d)),SNP,function(i)GENE[i][which.max(logp[i])]))
#        icol=1
#        for (i in unique(d$CHR)) {
#            with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
#            icol=icol+1
#        }
    }
    
    if (!is.null(annotate)) {
        d.annotate=d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, logp, col="green3", ...))
    }
    
    if (suggestiveline) abline(h=suggestiveline, col="blue")
    if (genomewideline) abline(h=genomewideline, col="red")
}


