# Script to replicate Figure 1 of the manuscript 

# libraries & functions ---------------------------------------------------
library(phyloseq)
library(picante)
library(UpSetR)
library(ggpubr)
library(RColorBrewer)
set.seed(333)

# colors for the plot
mancols <- brewer.pal(n = 8, name = 'Dark2')
pvalTransform <- list(cutpoints = c(0, 0.0001, 0.001,0.01, 0.05, 1), symbols = c("****", "***", "**", "*",  NA))
# FDR (ns: p > 0.05;*: p <= 0.05;**: p <= 0.01;***: p < = 0.001;****: p <= 0.0001)

# input data --------------------------------------------------------------
# load anonymized phyloseq object with 16S rRNA information
tab_16S <- readRDS("ASV_PhyloseqObject.rds")


# main --------------------------------------------------------------------

## obtain number of ASV
# calculate pairwise comparisons between groups and correct for multiple testing
tab_nASV <- plot_richness(tab_16S, x="Population",measures = "Observed")$data
tab_nASV$Population <- factor(x = tab_nASV$Population, levels = c("Agta", "Palanan", "BaYaka"),ordered = T)

pairwisePvals <- pairwise.t.test(x = tab_nASV$value,g = tab_nASV$Population,p.adjust.method = "bonferroni")
adjustedPvals <- as.vector(na.omit(as.vector(pairwisePvals$p.value)))
mycomparisons <- combn(x = levels(tab_nASV$Population),m = 2,simplify = F)

fig1A <- ggplot(tab_nASV, aes(x = Population, y = value)) +
  geom_violin(aes(fill=Population,alpha = 0.7),trim = F,show.legend = F) + geom_boxplot(color="white",width=0.1, aes(fill=Population,alpha = 0.7), show.legend = F) + stat_compare_means(comparisons = mycomparisons, label = 'p.signif', method = 't.test', tip.length = 0.0001,p.adjust.method='fdr')+
  theme_pubr() + scale_fill_manual(values = mancols[c(1,6,4)]) + ylab("Number of ASVs") + xlab("") +
  theme(plot.title = element_text(face = "bold"),legend.position = "none",strip.text = element_text(size = 18),axis.text =  element_text(size=16, colour = "black"),axis.title.y = element_text(size=18),axis.ticks = element_line(colour = "black"))  
## ggpubr does not apply multiple test correction, replace values
fig1A$layers[[3]]$stat_params$annotations <- as.vector(cut(x = adjustedPvals,breaks = pvalTransform$cutpoints,labels = pvalTransform$symbols))
fig1A

## calculate Faith phylogenetic diversity
# calculate pairwise comparisons between groups and correct for multiple testing
tab_pd <- picante::pd(t(as.data.frame(tab_16S@otu_table)), tab_16S@phy_tree,include.root=T) 

## add population information
tab_pd$id <- rownames(tab_pd)
tab_pd <- merge(tab_pd,tab_nASV[c("id", "Population")], by='id')

pairwisePvals <- pairwise.t.test(x = tab_pd$PD,g = tab_pd$Population,p.adjust.method = "bonferroni")
adjustedPvals <- as.vector(na.omit(as.vector(pairwisePvals$p.value)))
mycomparisons <- combn(x = levels(tab_pd$Population),m = 2,simplify = F)


fig1B <- ggplot(tab_pd, aes(x = Population, y = PD)) +
  geom_violin(aes(fill=Population,alpha = 0.7),trim = F,show.legend = F) + geom_boxplot(color="white",width=0.1, aes(fill=Population,alpha = 0.7), show.legend = F) + stat_compare_means(comparisons = mycomparisons, label = 'p.signif', method = 't.test', tip.length = 0.0001,p.adjust.method='fdr')+
  theme_pubr() + scale_fill_manual(values = mancols[c(1,6,4)]) + ylab("Faith's Phylogenetic Diversity") + xlab("") +
  theme(plot.title = element_text(face = "bold"),legend.position = "none",strip.text = element_text(size = 18),axis.text =  element_text(size=16, colour = "black"),axis.title.y = element_text(size=18),axis.ticks = element_line(colour = "black"))  
fig1B$layers[[3]]$stat_params$annotations <- as.vector(cut(x = adjustedPvals,breaks = pvalTransform$cutpoints,labels = pvalTransform$symbols))
fig1B


## Calculate number of shared and exclusive ASVs between populations
sampleSize <- 10 # sample size of the subsets
npermut <- 100 # number of permutations

# store the number of shared ASV between pairs of populations
resAB <- numeric(npermut)
resAC <- numeric(npermut)
resCB <- numeric(npermut)

# store the number of private ASVs in each population
pA <- numeric(npermut)
pB <- numeric(npermut)
pC <- numeric(npermut)


for (i in seq(npermut)){
  ## grab ASV present only in that subset of Agta
  table_pop <- subset_samples(tab_16S, id %in% sample(tab_16S@sam_data$id[which(tab_16S@sam_data$Population == "Agta")],size = sampleSize,replace = FALSE))
  table_pop <- prune_taxa(taxa_sums(table_pop) > 0 , table_pop)
  Agta <- rownames(otu_table(table_pop, taxa_are_rows = T))
  
  ## grab ASV present only in that subset of BaYaka
  table_pop <- subset_samples(tab_16S, id %in% sample(tab_16S@sam_data$id[which(tab_16S@sam_data$Population == "BaYaka")],size = sampleSize,replace = FALSE))
  table_pop <- prune_taxa(taxa_sums(table_pop) > 0 , table_pop)
  BaYaka <- rownames(otu_table(table_pop, taxa_are_rows = T))
  
  ## grab ASV present only in that subset of Palanan
  table_pop <- subset_samples(tab_16S, id %in% sample(tab_16S@sam_data$id[which(tab_16S@sam_data$Population == "Palanan")],size = sampleSize,replace = FALSE))
  table_pop <- prune_taxa(taxa_sums(table_pop) > 0 , table_pop)
  Palanan <- rownames(otu_table(table_pop, taxa_are_rows = T))
  
  ## intersect each pair of populations
  resAB[i] <- length(intersect(Agta,BaYaka))
  resAC[i] <- length(intersect(Agta,Palanan))
  resCB[i] <- length(intersect(Palanan,BaYaka))
  
  ## number of private bacteria
  tab_permut <- upset(fromList(list(Agta=Agta,BaYaka=BaYaka,Palanan=Palanan)))$New_data
  tab_permut$bacteria <- unique(unlist(list(Agta=Agta,BaYaka=BaYaka,Palanan=Palanan)))
  
  pA[i] <- length(tab_permut$bacteria[which(tab_permut$Agta == 1 & tab_permut$BaYaka == 0 & tab_permut$Palanan == 0 )])
  pB[i] <- length(tab_permut$bacteria[which(tab_permut$Agta == 0 & tab_permut$BaYaka == 1 & tab_permut$Palanan == 0 )])  
  pC[i] <- length(tab_permut$bacteria[which(tab_permut$Agta == 0 & tab_permut$BaYaka == 0 & tab_permut$Palanan == 1 )])  
  
}

## calculate number of shared ASVs
tab_overlap <- data.frame( Overlap = c(resAB, resAC, resCB),
                   Intersect = c(rep("Agta\nBaYaka",npermut), 
                                 rep("Agta\nPalanan",npermut), 
                                 rep("BaYaka\nPalanan", npermut) ))

pairwisePvals <- pairwise.t.test(x = tab_overlap$Overlap,g = tab_overlap$Intersect,p.adjust.method = "bonferroni")
adjustedPvals <- as.vector(na.omit(as.vector(pairwisePvals$p.value)))
mycomparisons <- combn(x = levels(tab_overlap$Intersect),m = 2,simplify = F)


fig1C <- ggplot(tab_overlap, aes(x = Intersect, y = Overlap)) + 
  geom_violin(aes(fill=Intersect,alpha = 0.7),trim = F,show.legend = F) + geom_boxplot(color="white",width=0.1, aes(fill=Intersect,alpha = 0.7), show.legend = F) +
  stat_compare_means(comparisons = mycomparisons, label = 'p.signif', method = 't.test', tip.length = 0.0001,p.adjust.method='fdr',size=5)+
  theme_pubr() + ylab("Number of shared ASVs") + xlab("") + 
  scale_fill_manual(values = mancols[c(2,2,2)]) +
  theme(plot.title = element_text(face = "bold"),legend.position = "none",strip.text = element_text(size = 18),axis.text =  element_text(size=16, colour = "black"),axis.title.y = element_text(size=18),axis.ticks = element_line(colour = "black")) 
fig1C$layers[[3]]$stat_params$annotations <- as.vector(cut(x = adjustedPvals,breaks = pvalTransform$cutpoints,labels = pvalTransform$symbols))
fig1C

## calculate number of private ASVs
datPriv <- data.frame( Private = c(pA, pB, pC),
                       Pop = c(rep("Agta",npermut), 
                               rep("BaYaka",npermut), 
                               rep("Palanan", npermut) ))

datPriv$Pop <- factor(x = datPriv$Pop, levels = c("Agta", "Palanan", "BaYaka"), ordered = T)

pairwisePvals <- pairwise.t.test(x = datPriv$Private,g = datPriv$Pop,p.adjust.method = "bonferroni")
adjustedPvals <- as.vector(na.omit(as.vector(pairwisePvals$p.value)))
mycomparisons <- combn(x = levels(datPriv$Pop),m = 2,simplify = F)

fig1D <- ggplot(datPriv, aes(x = Pop, y = Private)) +  
  geom_violin(aes(fill=Pop,alpha = 0.7),trim = F,show.legend = F) + geom_boxplot(color="white",width=0.08, aes(fill=Pop,alpha = 0.7), show.legend = F) + 
  stat_compare_means(comparisons = mycomparisons, label = 'p.signif', method = 't.test', tip.length = 0.0001,p.adjust.method='fdr',size=5)+
  theme_pubr() + ylab("Number of private ASVs") + xlab("") +
  scale_fill_manual(values = mancols[c(1,6,4)]) +
  theme(plot.title = element_text(face = "bold"),legend.position = "none",strip.text = element_text(size = 18),axis.text =  element_text(size=16, colour = "black"),axis.title.y = element_text(size=18),axis.ticks = element_line(colour = "black")) 
fig1D$layers[[3]]$stat_params$annotations <- as.vector(cut(x = adjustedPvals,breaks = pvalTransform$cutpoints,labels = pvalTransform$symbols))
fig1D

## final figure
fig1 <- ggarrange(plotlist = list(fig1A,fig1B,fig1C,fig1D),nrow = 2,ncol = 2,labels = c("a", "b", "c","d"),font.label = list(face="bold",family="arial",size=16))
ggsave(plot = fig1,filename = "Figure1.png" ,width = 10,height = 10)


