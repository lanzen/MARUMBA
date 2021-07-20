setwd("~/projects/MARUMBA/")
require(vegan)
source("~/kode/R/filtering.R")
source("~/kode/R/diversity.r")
source("~/kode/R/correlationTests.r")
source("~/kode/R/mergeOTUTable.R")
source("~/kode/R/taxaplot.R")

md = read.csv("metadata.csv",header=T,row.names=1,
              colClasses=c("character","factor","Date",
                           "character","factor","numeric","factor"))
summary(md)
# Site         Date             Parameters             Rep      VolFiltered     Fraction 
# Antartic :20   Min.   :2019-01-07   Length:149         Rep 1  :39   Min.   : 0.735   0.1 :72  
# AZTI     :14   1st Qu.:2019-09-18   Class :character   -      :22   1st Qu.: 5.000   0.22:77  
# Blanes   :13   Median :2019-10-22   Mode  :character   Rep 2  :12   Median : 7.500            
# KAEC     :96   Mean   :2019-10-03                      Rep 3  : 9   Mean   : 6.795            
# Xixon    : 6   3rd Qu.:2019-12-16                      Setup 1: 6   3rd Qu.:10.000            
# Max.   :2020-02-10                      Setup 2: 6   Max.   :10.000            
# NA's   :6                               (Other):55   NA's   :1                               

otus.t =  read.table("SWARM_20210716/CREST_LULU/SWARM_table_curated.tsv", row.names=1, header=T, 
                     check.names = F,sep="\t")
names(otus.t) = gsub("X","",names(otus.t))

tax.16S=data.frame(row.names=row.names(otus.t), 
                   classification = otus.t$classification)
bestTx = array(dim=dim(tax.16S)[1])
for (i in 1:length(bestTx)) bestTx[i] <- tail(unlist(strsplit(as.character(tax.16S[i,]), 
                                                              split=";", 
                                                              fixed=TRUE)), 1)
tax.16S$bestTx = bestTx

table(row.names(md) %in% names(otus.t))  # 149

otus = as.data.frame(t(otus.t[,row.names(md)]))

dim(otus) # 149 samples, 17,005 OTUs (5000 new)
summary(rowSums(otus)) # 64 -- 535k reads per dataset
sum(otus) # 35,504,300 reads


write.csv(file = "Unfiltered_div.csv",data.frame(reads=rowSums(otus),
                                                     richness=specnumber(otus)),
          quote=F)

# ---- Filter probable cross-contamination reads ----

otus = filterCrossContaminants2(otus)
sum(otus) # 35,501,800 (2500 removed)

# ---- Filtering all unclassified reads and separate into euk. and prokaryotes  -----

# Remove unclassified reads
toDel = row.names(tax.16S)[grep("No hits",tax.16S$classification)]
otus = otus[,!(names(otus) %in% toDel)]
dim(otus) # 16,563 (432 OTUs removed)
sum(otus) # 35,429,783 -> with 70k reads

# for (taxonDel in c("No hits","Eukaryota")){
#   toDel = row.names(tax.16S)[grep(taxonDel,tax.16S$classification)]
#   if (length(toDel)>0){
#     otus = otus[,!(names(otus) %in% toDel)]
#     print(paste("Removing ",taxonDel))
#     print(dim(otus)[2])
#     print(sum(otus))
#   }
# }

eukSWARMs = row.names(tax.16S)[grep("Eukaryota",tax.16S$classification)]
otus.euk = otus[,names(otus) %in% eukSWARMs]
dim(otus.euk) # 3026 euk 
sum(otus.euk) # 1,881,249

otus.prok = otus[,!names(otus) %in% eukSWARMs]
dim(otus.prok) # 13,537 euk 
sum(otus.prok) # 33,548,534

eukShare = rowSums(otus.euk) / rowSums(otus)
boxplot(eukShare*100~md$Site, ylab = "share of eukaryotic reads (%)")

eukShareDiv = specnumber(otus.euk) / specnumber(otus)
boxplot(eukShareDiv*100~md$Site, ylab = "share of eukaryotic OTUs (%)")


# --- Alpha diversity ----

require(ggplot2)

writeDivStats("Diversity_clean.csv", otus)
writeDivStats("Diversity_prokaryotic.csv", otus.prok)
writeDivStats("Diversity_eukaryotic.csv", otus.euk)

div.prok = read.csv("Diversity_prokaryotic.csv",header=T,row.names=1)
div.prok = div.prok[,c(1,3,5:7,9:10)]
table(row.names(div.prok) == row.names(md))
div.prok$Site = md$Site
div.prok$Fraction = md$Fraction

pdf("img/RRichness_bacterial_site.pdf",width=4,height=5)
ggplot(div.prok, aes(x=Site, y=Rarefied.richness)) + geom_boxplot(notch=T, outlier.colour="red", outlier.shape=16,
                                                 outlier.size=3) + geom_jitter(shape=16, position=position_jitter(0.2))

#boxplot(div.prok$Rarefied.richness~md$Site, main="Rarefied richness per site")
dev.off()

pdf("img/RRichness_fraction.pdf",width=2.5,height=5)
ggplot(div.prok, aes(x=Fraction, y=Rarefied.richness)) + geom_boxplot(notch=T, outlier.colour="red", outlier.shape=16,
                                                                   outlier.size=3) + geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()

pdf("img/RRichness_Site_and_fraciton",width=8.5,height=5)
div.prok$SiteFraction = paste(div.prok$Site, div.prok$Fraction,sep="/")
md$SiteFraction = paste(md$Site, md$Fraction,sep="/")
ggplot(div.prok, aes(x=SiteFraction, y=Rarefied.richness)) + geom_boxplot(notch=T, outlier.colour="red", outlier.shape=16,
                                                                      outlier.size=3) + geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()

pdf("img/Shannon_Site_and_fraciton",width=8.5,height=5)
ggplot(div.prok, aes(x=SiteFraction, y=H)) + geom_boxplot(notch=T, outlier.colour="red", outlier.shape=16,
                                                                          outlier.size=3) + geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()

rvf = lm(div.prok$Rarefied.richness~md$VolFiltered)
summary(rvf) #R2=.08 p=4E-4
plot(div.prok$Rarefied.richness~md$VolFiltered, xlab="Volume filtered (L)",
     ylab="Rarefied richness",col=as.numeric(md$Fraction),pch=as.numeric(md$Fraction))
abline(rvf,col="grey")
legend("topleft",pch=c(1,2),col=c(1,2),legend=c("0.1 uM", "0.22 uM"))

printANOVA(md[,c(1,5:7)], div.prok[,c(1:7)], a=.01)

## Rarefaction curves

pdf("img/Rarefaction_ultramicro_fraction.pdf",width = 9, height = 9)
otus.prok.micro = otus.prok[md$Fraction=="0.1",]
md.micro = md[md$Fraction=="0.1",]
rarecurve(otus.prok.micro, step = 1000, col=as.numeric(md.micro$Site)+1,cex=.6)
legend("bottomright",col=c(2:6),cex=.6,ncol=5,lwd=2,
       legend=sort(unique(md$Site)))
dev.off()


pdf("img/Rarefaction_ultramicro_fraction_logX.pdf",width = 9, height = 9)
rarecurve(otus.prok.micro, step = 1000, col=as.numeric(md.micro$Site)+1,cex=.6,
          log="x",xlim=range(2E3,5E5))
legend("bottomright",col=c(2:6),cex=.6,ncol=5,lwd=2,
       legend=sort(unique(md$Site)))
dev.off()


# ---- Filter rare OTUs by max abundance and export new OTU tables for CREST ----------

summary(rowSums(otus.prok))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 21788  107710  210379  225158  309930  530870 

# (earlier used 1E-4 due to 60k reads being sequened but most were euk. 
# "Still this is the real sequencing depth, and not 22k", not sure anymore)

otus.prok.ab = dropRareByMaxAbundance(otus.prok, 2E-4)
dim(otus.prok.ab) # 4015 OTUs of  9587
write.table(as.data.frame(t(otus.prok.ab)), file = "SWARM_20210716/SWARM_table_prok_abFiltered.tsv",
          sep="\t",quote=F)

# ---- Taxa barplots -------


gi = data.frame(row.names=paste(row.names(md),md$Site,md$Date),md$SiteFraction)

otus.prok.ab.t = as.data.frame(t(otus.prok.ab))
taxa.ab = tax.16S[row.names(otus.prok.ab.t),]
assignments.prok.ab = mergeOTUTable(otus.prok.ab.t, metadata = taxa.ab, by="bestTx")
assignments.prok.ab.ra = decostand(as.data.frame(t(assignments.prok.ab)),method="total")

pdf(file = "img/Assignments_Prok.pdf", width=25,height = 9)
taxaplot(30,gi,assignments.prok.ab.ra)
dev.off()

# ---- Taxa barplot for Euk ---
otus.euk.ab = dropRareByMaxAbundance(otus.euk,2E-4)
dim(otus.euk.ab) #2370 OTUs (of 2430)
otus.euk.ab.t = as.data.frame(t(otus.euk.ab))
taxa.ab = tax.16S[row.names(otus.euk.ab.t),]
assignments.euk.ab = mergeOTUTable(otus.euk.ab.t, metadata = taxa.ab, by="bestTx")
assignments.euk.ab.ra = decostand(as.data.frame(t(assignments.euk.ab)),method="total")

pdf(file = "img/Assignments_euk.pdf", width=25,height = 9)
taxaplot(30,gi,assignments.euk.ab.ra)
dev.off()

# ---- Intermediate ranks ------

ra.prok = read.table("SWARM_20210716/CREST_abFiltered/Relative_Abundance.tsv",sep="\t",
                     header=T,row.names=3, check.names = F)


ranks = data.frame(rank=c("domain","kingdom","phylum","class","order","family","genus","species"),
                   levels=c(2,11,21,21,21,21,21,3))


for (i in c(1:8)){
        r=as.character(ranks$rank[i])
        leveltaxa = as.data.frame(t(ra.prok[ra.prok$Rank==r,-c(1:2)]))
        print(paste(mean(rowSums(leveltaxa))*100,"% classified at rank",r))
        pdf(paste("img/",r,".pdf",sep=""),height=9,width=25)
        taxaplot(ranks$levels[i],gi,leveltaxa)
        dev.off()
}





# ---- NMDS ---------

otus.prok.ab.ra = decostand(otus.prok.ab,method="total")
nmds=metaMDS(otus.prok.ab.ra)
# Run 6 stress 0.1418382 
# ... Procrustes: rmse 1.271145e-05  max resid 0.0001010861 
# ... Similar to previous best

pdf("img/NMDS.pdf")
ordiplot(nmds,type="none")#,xlim=range(-4,2),ylim=range(-1,1))
ordihull(nmds,md$SiteFraction,label=T,col="grey", draw="line", cex=.6)
points(nmds,pch=as.numeric(md$Fraction),col=as.numeric(md$Site)+1)
legend("topright",pch=c(1,2,0,rep(1,5)),
       col=c(1,1,0,2:6),cex=.8,
       legend=c("0.1 uM", "0.2 uM", "",as.character(sort(unique(md$Site)))))
dev.off()

