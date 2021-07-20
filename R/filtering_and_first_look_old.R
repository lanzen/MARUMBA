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
# Site            Date             Parameters           Rep      VolFiltered        Fraction     
# Antartic :20   2019-10-22:38   BL1C          :  2   -          :22   Min.   : 0.735   Min.   :0.1000  
# AZTI     :14   2019-09-18:14   BL1D          :  2   Rep 1      :12   1st Qu.: 3.190   1st Qu.:0.1000  
# Blanes   :13   2019-09-17:13   BL1E          :  2   Rep 2      :12   Median : 9.500   Median :0.2200  
# KAEC     :60   2020-02-10:12   100m 10L 0.1 1:  1   Setup 1    : 6   Mean   : 6.836   Mean   :0.1627  
# Xixon    : 6   2019-09-24:10   100m 10L 0.1 2:  1   Setup 2    : 6   3rd Qu.:10.000   3rd Qu.:0.2200  
# (Other)   :20   100m 10L 0.2 1:  1   Replicate A: 5   Max.   :10.000   Max.   :0.2200  
# NA's      : 6   (Other)       :104   (Other)    :50   NA's   :1                        

otus.t =  read.table("SWARM_20210322/CREST_LULU/SWARM_table_curated.tsv", row.names=1, header=T, 
                     check.names = F,sep="\t")
names(otus.t) = gsub("X","",names(otus.t))

tax.16S=data.frame(row.names=row.names(otus.t), 
                   classification = otus.t$classification)
bestTx = array(dim=dim(tax.16S)[1])
for (i in 1:length(bestTx)) bestTx[i] <- tail(unlist(strsplit(as.character(tax.16S[i,]), 
                                                              split=";", 
                                                              fixed=TRUE)), 1)
tax.16S$bestTx = bestTx

table(row.names(md) %in% names(otus.t))  # Yes


otus = as.data.frame(t(otus.t[,row.names(md)]))

dim(otus) # 113 samples, 12,483 OTUs
summary(rowSums(otus)) # 64 -- 328k reads per dataset
sum(otus) # 19,893,442 reads


write.csv(file = "Unfiltered_div.csv",data.frame(reads=rowSums(otus),
                                                     richness=specnumber(otus)),
          quote=F)

# ---- Filter probable cross-contamination reads ----

otus = filterCrossContaminants2(otus)
sum(otus) # 19,891,554 (1900 removed)

# ---- Filtering all unclassified reads and separate into euk. and prokaryotes  -----

# Remove unclassified reads
toDel = row.names(tax.16S)[grep("No hits",tax.16S$classification)]
otus = otus[,!(names(otus) %in% toDel)]
dim(otus) # 12,017 483 (466 OTUs removed)
sum(otus) # with 70k reads

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
dim(otus.euk) # 2430 euk 
sum(otus.euk) # 1,491,669

otus.prok = otus[,!names(otus) %in% eukSWARMs]
dim(otus.prok) # 9587 euk 
sum(otus.prok) # 18,333,044

eukShare = rowSums(otus.euk) / rowSums(otus)
boxplot(eukShare*100~md$Site, ylab = "share of eukaryotic reads (%)")

eukShareDiv = specnumber(otus.euk) / specnumber(otus)
boxplot(eukShareDiv*100~md$Site, ylab = "share of eukaryotic OTUs (%)")


# --- Alpha diversity ----

writeDivStats("Diversity_clean.csv", otus)
writeDivStats("Diversity_prokaryotic.csv", otus.prok)
writeDivStats("Diversity_eukaryotic.csv", otus.euk)

div.prok = read.csv("Diversity_prokaryotic.csv",header=T,row.names=1)
div.prok = div.prok[,c(1,3,5:7,9:10)]

boxplot(div.prok$Rarefied.richness~md$Site, main="Rarefied richness per site")

boxplot(div.prok$Rarefied.richness~md$Fraction, main="Rarefied richness per fraction",
        notch=T)

md$SiteFraction = paste(md$Site, md$Fraction)
boxplot(div.prok$Rarefied.richness~md$SiteFraction, 
        main="Rarefied richness per site and fraction",las=2)

rvf = lm(div.prok$Rarefied.richness~md$VolFiltered)
summary(rvf) #R2=.18 p=2E-6
plot(div.prok$Rarefied.richness~md$VolFiltered, xlab="Volume filtered (L)",
     ylab="Rarefied richness",col=as.numeric(md$Fraction),pch=as.numeric(md$Fraction))
abline(rvf,col="grey")
legend("topleft",pch=c(1,2),col=c(1,2),legend=c("0.1 uM", "0.22 uM"))

printANOVA(md[,c(1,6:7)], div.prok, a=.01)

## Rarefaction curves

otus.prok.micro = otus.prok[md$Fraction=="0.1",]
md.micro = md[md$Fraction=="0.1",]
rarecurve(otus.prok.micro, step = 1000, col=as.numeric(md$Site)+1,cex=.6)
legend("bottomright",col=c(2:6),cex=.6,ncol=5,lwd=2,
       legend=sort(unique(md$Site)))

# ---- Filter rare OTUs by max abundance and export new OTU tables for CREST ----------

summary(rowSums(otus.prok))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 22319   98699  142409  162239  228745  315930 
# Very big difference to total because where the lowest had 63k reads! 
# Still this is the real sequencing depth, and not 22k

# We require at least 0.01% abundance in at least one sample. 
# This corresponds to 6 reads in the smallest dataset

otus.prok.ab = dropRareByMaxAbundance(otus.prok, 1E-4)
dim(otus.prok.ab) # 4773 OTUs of  9587
write.table(as.data.frame(t(otus.prok.ab)), file = "SWARM_20210322/SWARM_table_prok_abFiltered.tsv",
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
otus.euk.ab = dropRareByMaxAbundance(otus.euk,1E-4)
dim(otus.euk.ab) #2370 OTUs (of 2430)
otus.euk.ab.t = as.data.frame(t(otus.euk.ab))
taxa.ab = tax.16S[row.names(otus.euk.ab.t),]
assignments.euk.ab = mergeOTUTable(otus.euk.ab.t, metadata = taxa.ab, by="bestTx")
assignments.euk.ab.ra = decostand(as.data.frame(t(assignments.euk.ab)),method="total")

pdf(file = "img/Assignments_euk.pdf", width=25,height = 9)
taxaplot(30,gi,assignments.euk.ab.ra)
dev.off()

# ---- NMDS ---------

otus.prok.ab.ra = decostand(otus.prok.ab,method="total")
nmds=metaMDS(otus.prok.ab.ra)
# Run 6 stress 0.1448751 
# ... New best solution
# ... procrustes: rmse 8.991444e-05  max resid 0.0002549636 
# *** Solution reached

ordiplot(nmds,type="none")#,xlim=range(-4,2),ylim=range(-1,1))
ordihull(nmds,md$SiteFraction,label=T,col="grey", draw="line", cex=.6)
points(nmds,pch=as.numeric(md$Fraction),col=as.numeric(md$Site)+1)
legend("bottomright",pch=c(1,2,0,rep(1,5)),
       col=c(1,1,0,2:6),cex=.8,
       legend=c("0.1 uM", "0.2 uM", "",as.character(sort(unique(md$Site)))))

