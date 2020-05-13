


########## code for manuscript: 
# "Ecological insights from gut microbiome networks: linking taxa distributions to co-occurrence networks across vertebrates"

# Risely et al 2020

library(phyloseq)
library(ggplot2)
library(vegan)
library(directlabels)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(ape)
library(gridExtra)
library(ade4)
library(plyr)
library(tidyr)
library(data.table)
library(stringr)
library(expss)
library(ranacapa)


library("igraph")
library("qgraph")
library("vegan")
library("MCL")
library("SpiecEasi")
library("ape")
library("reshape2")
library("expss")
library("microbiome")



########################### UPLOAD DATA ##################
########################### UPLOAD DATA ##################
########################### UPLOAD DATA ##################
########################### UPLOAD DATA ##################
########################### UPLOAD DATA ##################


merged_7species_unrarefied<-readRDS("data_7species_unrarefied.rds") #phyloseq object containing data for 7 species
shorebird_unrarefied<-readRDS("data_shorebird_unrarefied.rds") #phyloseq object containing data for Red-necked stint


### merge for plotting (without tree, because stint were sequenced using different primers)
## this means for some part of the processing stint have to be processed seperately from the other 7 species


table<-otu_table(merged_7species_unrarefied)
map<-sample_data(merged_7species_unrarefied)
taxonomy<-tax_table(merged_7species_unrarefied)

merged_8species <-merge_phyloseq(table, map, taxonomy)


table<-otu_table(shorebird_unrarefied)
map<-sample_data(shorebird_unrarefied)
taxonomy<-tax_table(shorebird_unrarefied)

stint <-merge_phyloseq(table, map, taxonomy)

merged_8species<-merge_phyloseq(merged_8species, stint) #phylo object 8 species with no tree


################## SUPPLEMENTARY FIGURE 2
################## SUPPLEMENTARY FIGURE 2
################## SUPPLEMENTARY FIGURE 2
################## SUPPLEMENTARY FIGURE 2

############################ rarefaction curves ########################
############################ rarefaction curves ########################
############################ rarefaction curves ########################
############################ rarefaction curves ########################



sample_data(merged_8species)$Species<- factor(sample_data(merged_8species)$Species, levels = c("human", "meerkat", "deer", "carollia","spinyrat","mouselemur","flamingo","Red-necked stint"))


p <- ranacapa::ggrare(merged_8species, step = 500,   se = FALSE)+ xlim(c(0,50000))

p + xlim(c(0,50000))

p1<- p + facet_wrap(~Species, scales="free_y", ncol=4)+geom_vline(xintercept=10000)+theme_bw()+theme(legend.position = "none")

rarefaction_fig<-p1


rarefaction_fig1<- rarefaction_fig+geom_line(aes(col=Species))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
 # theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=14))+ylab("Number of ASVs")+xlab("Number of reads")+
  scale_color_brewer(palette = "Dark2")

rarefaction_fig1



############################### rarefy to 10000 ###################
############################### rarefy to 10000 ###################
############################### rarefy to 10000 ###################
############################### rarefy to 10000 ###################

merged_7species<- rarefy_even_depth(merged_7species_unrarefied, sample.size = 10000, rngsee = 100, replace = TRUE, trimOTUs=TRUE,verbose=TRUE)
merged_7species<-prune_taxa(taxa_sums(merged_7species)>0, merged_7species)

###############


shorebird<- rarefy_even_depth(shorebird_unrarefied, sample.size = 10000, rngsee = 100, replace = TRUE, trimOTUs=TRUE,verbose=TRUE)
shorebird<-prune_taxa(taxa_sums(shorebird)>0, shorebird)

############# generate phyloseq object for all species by merging the two rarefied objects (rarefying seperately changes number of ASVs)

table<-otu_table(merged_7species)
map<-sample_data(merged_7species)
taxonomy<-tax_table(merged_7species)

merged_8species_rare <-merge_phyloseq(table, map, taxonomy)


table<-otu_table(shorebird)
map<-sample_data(shorebird)
taxonomy<-tax_table(shorebird)

stint <-merge_phyloseq(table, map, taxonomy)

merged_8species_rare<-merge_phyloseq(merged_8species_rare, stint) #rarefied phylo object 8 species no tree



###############################  species accumulation curve ######
###############################  species accumulation curve ######
###############################  species accumulation curve ######
###############################  species accumulation curve ######
###############################  species accumulation curve ######

SAClist<-list()

uniq <- unique(sample_data(merged_8species)$Species)



for (i in 1:length(uniq)){
  
  
data_1<-subset_samples(merged_8species, Species == uniq[i])

data_1<-prune_taxa(taxa_sums(data_1)>0, data_1)
data_1_matrix<-data.frame(t(data.frame(otu_table(data_1))))
data_1_specaccum<-specaccum(data_1_matrix, method="random", permutations = 500)


sac_df<- data_1_specaccum$sites
sac_df<-data.frame(sac_df)
names(sac_df)[1]<-"Site"
sac_df$Richness <-  data_1_specaccum$richness
sac_df$SD <-  data_1_specaccum$sd


sac_total_estimated<-specpool(data_1_matrix)
sac_df$Total <- sac_total_estimated$boot

sac_df$Species <- as.character(sample_data(data_1)$Species[1])
head(sac_df)
SAClist[[i]]<-sac_df

}


sac_df_all<-do.call(rbind, SAClist)

head(sac_df_all)


######### species accumulation figure
######### species accumulation figure
######### species accumulation figure



ggplot(sac_df_all, aes(x = Site, y = Richness, group = Species))+
  geom_line(alpha=0.7, linetype = "dashed")+
  geom_point( aes(col=Species), size = 2)+
  theme_bw()+
  xlab("Number of individuals sampled")+
  ylab("Number of ASVs")+
  theme(text=element_text(size=14))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(legend.position = "none")+ 
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Species accumulation curves")+
  ylim(c(0,6000))+xlim(c(0,50))+
  geom_dl(aes(label = Species), method = list("last.points", cex = 0.8, vjust = -0.2, hjust = 0.5)) 




######################################### PREVALENCE/OCCUPANCY FREQUENCY ##################
######################################### PREVALENCE/OCCUPANCY FREQUENCY ##################


## function for prevalence

prevalence <- function(physeq, add_tax = TRUE){
  
  ## Check if taxa are rows
  trows <- taxa_are_rows(physeq)
  
  ## Extract OTU table
  otutab <- as.data.frame(otu_table(physeq))
  
  ## Transpose OTU table (species should be arranged by rows)
  if(trows == FALSE){
    otutab <- t(otutab)
  }
  
  ## Estimate prevalence (number of samples with OTU present)
  prevdf <- apply(X = otutab,
    #MARGIN = ifelse(trows, yes = 1, no = 2),  # for a non-transposed data
    MARGIN = 1,
    FUN = function(x){sum(x > 0)})
  
  ## Add total and average read counts per OTU
  prevdf <- data.frame(Prevalence = prevdf,
    TotalAbundance = taxa_sums(physeq),
    MeanAbundance = rowMeans(otutab),
    MedianAbundance = apply(otutab, 1, median))
  
  ## Add taxonomy table
  if(add_tax == TRUE && !is.null(tax_table(physeq, errorIfNULL = F))){
    prevdf <- cbind(prevdf, tax_table(physeq))
  }
  return(prevdf)
}

#################### loop to calculate prevalence and abundance of every ASV per host species

Prevlist<-list()

uniq <- unique(sample_data(merged_8species_rare)$Species)

for (i in 1:length(uniq)){
  
  data_1<-subset_samples(merged_8species_rare, Species == uniq[i])
  data_1<-prune_taxa(taxa_sums(data_1)>0, data_1)
  
  occupancy_abundance<-prevalence(data_1)
  occupancy_abundance$host_species <- as.character(sample_data(data_1)$Species[1])
  occupancy_abundance$RelAbundance<- (occupancy_abundance$TotalAbundance/sum(occupancy_abundance$TotalAbundance))
  occupancy_abundance$RelPrev<-(occupancy_abundance$Prevalence / length(sample_data(data_1)$feature.id))
  occupancy_abundance$RelAbundanceInd<-(occupancy_abundance$TotalAbundance/ occupancy_abundance$Prevalence/10000)
  occupancy_abundance$ASV<-taxa_names(data_1) #this is important! Don't use 'row.names' as it sneakily adds a 1 at the end if it occurs twice
  occupancy_abundance$Sample_size<-length(unique(sample_data(data_1)$feature.id))
  
  Prevlist[[i]]<-occupancy_abundance
  
}

#combine

occupancy_abundance_df<-do.call(rbind, Prevlist)


#occupancy_abundance2 <- read.csv("C:/Users/risel/Dropbox/Sommer postdoc/Core microbiome project/Analysis/MERGED/occupancy_abundance_final.csv", sep = ",")[,-1]


######################################## CO-OCCURRENCE NETWORKS #####################
######################################## CO-OCCURRENCE NETWORKS #####################
######################################## CO-OCCURRENCE NETWORKS #####################
######################################## CO-OCCURRENCE NETWORKS #####################
######################################## CO-OCCURRENCE NETWORKS #####################
######################################## CO-OCCURRENCE NETWORKS #####################


##################### generate loop to generate co-occurrence networks and extract node and vertex info for each species

#memory.limit(size = 30000)

############### NOT RUN ################ (takes many days with 500 iterations)
############### NOT RUN ################
############### NOT RUN ################
############### NOT RUN ################
############### NOT RUN ################

Net.list<-list()

##########

set.seed(100)

uniq <- unique(sample_data(merged_8species)$Species) #unrarefied data


for (i in 1:length(uniq)){
  
  data<-subset_samples(merged_8species, Species==uniq[i])
  data<-prune_taxa(taxa_sums(data)>0, data)
  
  data_df<-subset(occupancy_abundance_df, host_species==uniq[i])
  
  
  data10.df<-subset(data_df, RelPrev >= 0.10)
  asv10<-data10.df$ASV
  
  
  data10<-prune_taxa(asv10, data) #subset only taxa we want from the rarefied dataset
  data10<-prune_taxa(taxa_sums(data10)>0, data10)
  
  #### start network analysis
  
  data10.df<-data.frame(otu_table(data10))
  
  
  ##need to transpose dataset so it is in the right format for network analysis
  
  data10.df<-t(data10.df)
  
  
  
  ##spiec easi MB method
  
  
  se.mb.data10 <- spiec.easi(data10.df, 
    method='mb', 
    lambda.min.ratio=1e-2,
    nlambda=20, 
    pulsar.params=list(rep.num=500))
  
  
  ##get weights
  
  sebeta <- symBeta(getOptBeta(se.mb.data10 ), mode='maxabs')
  diag(sebeta) <- 0
  elist.mb     <- summary(sebeta)
  #str(elist.mb) #x are weights
  # hist(elist.mb[,3], main='', xlab='edge weights')
  weights<-elist.mb$x
  
  ###concert to igraph object
  
  net.data10  <- adj2igraph(Matrix::drop0(getRefit(se.mb.data10)),
    edge.attr=list(weight=weights),
    rmEmptyNodes = FALSE, 
    diag = FALSE,
    vertex.attr = list(name=colnames(data10.df)))
  
  
  
  ###ok lets start building up our network object with useful information
  
  E(net.data10)$weight1<-weights
  E(net.data10)$dir<-ifelse(E(net.data10)$weight1 > 0 , "pos", "neg") #1 = positive 0 = negative
  E(net.data10)$weight<-abs(E(net.data10)$weight1)
  
  ##add taxonomy
  
  taxonomy<-data.frame(tax_table(data10))
  
  
  family<-as.character(taxonomy$Family)
  order<-as.character(taxonomy$Order)
  genus<-as.character(taxonomy$Genus)
  
  
  V(net.data10)$family<-family
  V(net.data10)$order<-order
  V(net.data10)$genus<-genus
  
  
  #Calculate network indices for each ASV (hub score, betweenness, degree, and closeness centrality)
  
  
  net.cn <- closeness(net.data10) #closeness centrality
  net.bn <- betweenness(net.data10) #betweenness centrality
  net.deg<-degree(net.data10) #degree
  net.hs <- hub_score(net.data10)$vector
  
  
  V(net.data10)$closeness<-net.cn
  V(net.data10)$betweenness<-net.bn
  V(net.data10)$degree<-net.deg
  V(net.data10)$hubbiness<-net.hs
  
  
  # add weighted degree
  
  net.weighted_deg <- strength(net.data10)
  V(net.data10)$weighted_deg<-net.weighted_deg
  
  V(net.data10)$host_species<-as.character(sample_data(data)$Species[1])
  
  Net.list[[i]]<-net.data10
  
}

names(Net.list)<-uniq  


####################### END RUN #####################
####################### END RUN #####################
####################### END RUN #####################
####################### END RUN #####################

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

# instead of running the above, import list of networks

#saveRDS(Net.list, "C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Core microbiome project\\Analysis\\Network objects\\all_networks_list.RDS")
Net.list<-readRDS("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Core microbiome project\\Analysis\\Network objects\\all_networks_list.RDS")


########################################################## generate vertex and edge info - loop ################
########################################################## generate vertex and edge info - loop ################
########################################################## generate vertex and edge info - loop ################
########################################################## generate vertex and edge info - loop ################
########################################################## generate vertex and edge info - loop ################
########################################################## generate vertex and edge info - loop ################


Net.list7<-Net.list[-8] #remove shorebird


Edge.list<-list()
Vertex.list<-list()

##########


uniq <- unique(sample_data(merged_7species)$Species)


for (i in 1:length(uniq)){
  
  net_1<-Net.list7[[i]]

  edges_df<-igraph::as_long_data_frame(net_1)
  
  head(edges_df)
  

  ## add relatedness. Need trees
  
  data<-subset_samples(merged_7species, Species==uniq[i])
  data<-prune_taxa(taxa_sums(data)>0, data)
  
  data_df<-subset(occupancy_abundance_df, host_species==uniq[i])
  
  
  data10.df<-subset(data_df, RelPrev >= 0.10)
  asv10<-data10.df$ASV
  
  
  data10<-prune_taxa(asv10, data) #subset only taxa we want from the rarefied dataset
  data10<-prune_taxa(taxa_sums(data10)>0, data10)
  
  
  otu_table<-t(data.frame(otu_table(data10)))
  tree<-phy_tree(data10)
  PatristicDistMatrix<-cophenetic.phylo(tree)
  
  relatedness<-reshape2::melt(PatristicDistMatrix)
  head(relatedness)
  
  colnames(relatedness)<-c("From","To","Relatedness")
  
  
  ########################################
  
  edges_df$Edge<-paste(edges_df$from_name, "-", edges_df$to_name)
  
  relatedness$Edge<-paste(relatedness$From, "-", relatedness$To)
  
  edges_df$Relatedness<-vlookup(edges_df$Edge, relatedness, lookup_column = "Edge", result_column = "Relatedness")
  
  head(edges_df)
  
  
  #######################################################################################################################
  
  
  
  vertices_df<-igraph::as_data_frame(net_1, what = "vertices")
  head(vertices_df)
  
  edges<-edges_df
  
  
  edges_long <-  edges[,c("dir","from_name","to_name")]
  edges_long <- gather(edges_long, To_From, ASV, c(from_name, to_name), factor_key=TRUE)
  
  
  positive_edges<-subset(edges_long, dir == "pos")
  negative_edges<-subset(edges_long, dir == "neg")
  
  
  sum_pos<-ddply(positive_edges, .(ASV), summarize, sum_pos= length(dir))
  sum_neg<-ddply(negative_edges, .(ASV), summarize, sum_neg= length(dir))
  
  merged<-merge(sum_pos, sum_neg, by = c("ASV"), all = T)
  head(merged)
  
  
  vertices_df1<-merge(vertices_df, merged, by.x = c("name"), by.y = c("ASV"), all = T)
  
  head(vertices_df1)
  
  vertices_df1$sum_pos[is.na(vertices_df1$sum_pos)]<-0 #make NAs into zeros
  vertices_df1$sum_neg[is.na(vertices_df1$sum_neg)]<-0
  
  
  ########## calculate freq won by looking at which ASV has the highest local abundance
  
  str(data10.df)
  
  traits<-data10.df[,c("ASV", "RelAbundance", "RelAbundanceInd", "RelPrev")]
  
  edges_df1<-merge(edges_df, traits, by.x=c("from_name"), by.y = c("ASV"))
  
  head(edges_df1)
  
  colnames(edges_df1)[26]<-"From_abundance"
  colnames(edges_df1)[27]<-"From_abundanceInd"
  colnames(edges_df1)[28]<-"From_prevalence"
  
  edges_df1<-merge(edges_df1, traits, by.x=c("to_name"), by.y = c("ASV"))
  
  colnames(edges_df1)[29]<-"To_abundance"
  colnames(edges_df1)[30]<-"To_abundanceInd"
  colnames(edges_df1)[31]<-"To_prevalence"
  
  head(edges_df1)
  
  #######################################
  
  negative<-subset(edges_df1, dir =="neg")
  
  negative$from_name<-as.character(negative$from_name)
  negative$to_name<-as.character(negative$to_name)
  
  negative$winner<-ifelse(negative$From_abundanceInd > negative$To_abundanceInd, negative$from_name, negative$to_name)
  negative$loser<-ifelse(negative$From_abundanceInd < negative$To_abundanceInd, negative$from_name, negative$to_name)
  
  
  winners<-data.frame(table(negative$winner))
  losers<-data.frame(table(negative$loser))
  
  
  winners<-winners[order(-winners$Freq),]
  losers<-losers[order(-losers$Freq),]

  colnames(winners)<-c("ASV","Freq_won")
  colnames(losers)<-c("ASV","Freq_lost")
  
  vertices_df1$Freq_won<-vlookup(vertices_df1$name, winners, lookup_column = "ASV", result_column = "Freq_won")
  vertices_df1$Freq_lost<-vlookup(vertices_df1$name, losers, lookup_column = "ASV", result_column = "Freq_lost")
  

  vertices_df1$Freq_won[is.na(vertices_df1$Freq_won)]<-0
  vertices_df1$Freq_lost[is.na(vertices_df1$Freq_lost)]<-0
  
  vertices_df1$Percent_won<-vertices_df1$Freq_won/(vertices_df1$Freq_won+vertices_df1$Freq_lost)
  
  ####################
  ####################
  ####################
  
  vertices_df1$host_species<-uniq[i]
  edges_df1$host_species<-uniq[i]
  
  Vertex.list[[i]]<-vertices_df1
  Edge.list[[i]]<-edges_df1
  
}


###### combine

vertex_attributes_df<-do.call(rbind, Vertex.list)
edge_attributes_df<-do.call(rbind, Edge.list)


############################ repeat for shorebird ##############################################
############################ repeat for shorebird ##############################################
############################ repeat for shorebird ##############################################
############################ repeat for shorebird ##############################################

net_1<-Net.list$`Red-necked stint`

edges_df<-igraph::as_long_data_frame(net_1)

head(edges_df)

########## add relatedness

data<-shorebird
data<-prune_taxa(taxa_sums(data)>0, data)

data_df<-subset(occupancy_abundance_df, host_species=="Red-necked stint")


data10.df<-subset(data_df, RelPrev >= 0.10)
asv10<-data10.df$ASV


data10<-prune_taxa(asv10, data) #subset only taxa we want from the rarefied dataset
data10<-prune_taxa(taxa_sums(data10)>0, data10)


otu_table<-t(data.frame(otu_table(data10)))
tree<-phy_tree(data10)
PatristicDistMatrix<-cophenetic.phylo(tree)

relatedness<-reshape2::melt(PatristicDistMatrix)
head(relatedness)

colnames(relatedness)<-c("From","To","Relatedness")


########################################

edges_df$Edge<-paste(edges_df$from_name, "-", edges_df$to_name)

relatedness$Edge<-paste(relatedness$From, "-", relatedness$To)

edges_df$Relatedness<-vlookup(edges_df$Edge, relatedness, lookup_column = "Edge", result_column = "Relatedness")

head(edges_df)


#######################################################################################################################



vertices_df<-igraph::as_data_frame(net_1, what = "vertices")
head(vertices_df)

edges<-edges_df


edges_long <-  edges[,c("dir","from_name","to_name")]
edges_long <- gather(edges_long, To_From, ASV, c(from_name, to_name), factor_key=TRUE)


positive_edges<-subset(edges_long, dir == "pos")
negative_edges<-subset(edges_long, dir == "neg")


sum_pos<-ddply(positive_edges, .(ASV), summarize, sum_pos= length(dir))
sum_neg<-ddply(negative_edges, .(ASV), summarize, sum_neg= length(dir))

merged<-merge(sum_pos, sum_neg, by = c("ASV"), all = T)
head(merged)


vertices_df1<-merge(vertices_df, merged, by.x = c("name"), by.y = c("ASV"), all = T)

head(vertices_df1)

vertices_df1$sum_pos[is.na(vertices_df1$sum_pos)]<-0
vertices_df1$sum_neg[is.na(vertices_df1$sum_neg)]<-0


########## freq won
str(data10.df)

traits<-data10.df[,c("ASV", "RelAbundance", "RelAbundanceInd", "RelPrev")]

edges_df1<-merge(edges_df, traits, by.x=c("from_name"), by.y = c("ASV"))

head(edges_df1)

colnames(edges_df1)[26]<-"From_abundance"
colnames(edges_df1)[27]<-"From_abundanceInd"
colnames(edges_df1)[28]<-"From_prevalence"

edges_df1<-merge(edges_df1, traits, by.x=c("to_name"), by.y = c("ASV"))

colnames(edges_df1)[29]<-"To_abundance"
colnames(edges_df1)[30]<-"To_abundanceInd"
colnames(edges_df1)[31]<-"To_prevalence"

head(edges_df1)

#######################################

negative<-subset(edges_df1, dir =="neg")

negative$from_name<-as.character(negative$from_name)
negative$to_name<-as.character(negative$to_name)

negative$winner<-ifelse(negative$From_abundanceInd > negative$To_abundanceInd, negative$from_name, negative$to_name)
negative$loser<-ifelse(negative$From_abundanceInd < negative$To_abundanceInd, negative$from_name, negative$to_name)


winners<-data.frame(table(negative$winner))
losers<-data.frame(table(negative$loser))


winners<-winners[order(-winners$Freq),]
losers<-losers[order(-losers$Freq),]


colnames(winners)<-c("ASV","Freq_won")
colnames(losers)<-c("ASV","Freq_lost")

vertices_df1$Freq_won<-vlookup(vertices_df1$name, winners, lookup_column = "ASV", result_column = "Freq_won")
vertices_df1$Freq_lost<-vlookup(vertices_df1$name, losers, lookup_column = "ASV", result_column = "Freq_lost")

head(vertices_df1)


vertices_df1$Freq_won[is.na(vertices_df1$Freq_won)]<-0
vertices_df1$Freq_lost[is.na(vertices_df1$Freq_lost)]<-0

vertices_df1$Percent_won<-vertices_df1$Freq_won/(vertices_df1$Freq_won+vertices_df1$Freq_lost)

####################
####################
####################

vertices_df1$host_species<-"Red-necked stint"
edges_df1$host_species<-"Red-necked stint"


####combine

vertex_attributes_df<-rbind(vertex_attributes_df, vertices_df1)
edge_attributes_df<-rbind(edge_attributes_df, edges_df1)

##change name of first column to 'ASV'

names(vertex_attributes_df)[1]<-"ASV"

str(vertex_attributes_df)
str(edge_attributes_df)


str(occupancy_abundance_df)


######################################## add abundance and prevalence data from occupancy_abundance_df

variables_to_add<-occupancy_abundance_df[,c(1:4,12:17)]

occupancy_abundance_2<-merge(vertex_attributes_df, variables_to_add, by = c("ASV", "host_species"), all.x = T, all.y = F)


str(occupancy_abundance_2)

###########################

################################################# STATS AND ANALYSIS #########################################
################################################# STATS AND ANALYSIS #########################################
################################################# STATS AND ANALYSIS #########################################
################################################# STATS AND ANALYSIS #########################################
################################################# STATS AND ANALYSIS #########################################
################################################# STATS AND ANALYSIS #########################################
################################################# STATS AND ANALYSIS #########################################
################################################# STATS AND ANALYSIS #########################################




library(effects)
library(MuMIn)
library(lme4)
library(lmerTest)
library(standardize)
library(ggplot2)
library(gam)
library(tidyverse)
library(expss)
library(plyr)
library(broom)
library(performance)



#order levels

str(occupancy_abundance_2)
table(occupancy_abundance_2$host_species)


occupancy_abundance_2$host_species<-factor(occupancy_abundance_2$host_species, levels = c("human", "meerkat", "deer", "carollia","spinyrat","mouselemur","flamingo","Red-necked stint"))

#rename

levels(occupancy_abundance_2$host_species)[levels(occupancy_abundance_2$host_species)=="human"] <- "Human"
levels(occupancy_abundance_2$host_species)[levels(occupancy_abundance_2$host_species)=="meerkat"] <- "Meerkat"
levels(occupancy_abundance_2$host_species)[levels(occupancy_abundance_2$host_species)=="deer"] <- "Deer"
levels(occupancy_abundance_2$host_species)[levels(occupancy_abundance_2$host_species)=="carollia"] <- "Carollia"
levels(occupancy_abundance_2$host_species)[levels(occupancy_abundance_2$host_species)=="spinyrat"] <- "Spinyrat"
levels(occupancy_abundance_2$host_species)[levels(occupancy_abundance_2$host_species)=="mouselemur"] <- "Mouselemur"
levels(occupancy_abundance_2$host_species)[levels(occupancy_abundance_2$host_species)=="flamingo"] <- "Flamingo"
levels(occupancy_abundance_2$host_species)[levels(occupancy_abundance_2$host_species)=="Red-necked stint"] <- "Stint"


###replace NA's with zeros


occupancy_abundance_2$closeness[is.na(occupancy_abundance_2$closeness)] <- 0
occupancy_abundance_2$betweenness[is.na(occupancy_abundance_2$betweenness)] <- 0
occupancy_abundance_2$degree[is.na(occupancy_abundance_2$degree)] <- 0
occupancy_abundance_2$hubbiness[is.na(occupancy_abundance_2$hubbiness)] <- 0
occupancy_abundance_2$weighted_deg[is.na(occupancy_abundance_2$weighted_deg)] <- 0
occupancy_abundance_2$sum_neg[is.na(occupancy_abundance_2$sum_neg)] <- 0
occupancy_abundance_2$sum_pos[is.na(occupancy_abundance_2$sum_pos)] <- 0


head(occupancy_abundance_2)

#generate data without deer, to test whether the deer dataset is driving any trends (due to much higher number of positive and negative edges)

occupancy_abundance_reduced<-subset(occupancy_abundance_2, host_species != "Deer")




########## add column showing whether ASV has high betweenness (is a linker)

betweenness_SD<-ddply(occupancy_abundance_2, .(host_species), summarize, mean=mean(betweenness), std_dev=sd(betweenness))
betweenness_SD$std_dev_upper<-betweenness_SD$mean+betweenness_SD$std_dev
str(betweenness_SD)


occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Human" & occupancy_abundance_2$betweenness > betweenness_SD[1,4], TRUE, NA)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Meerkat" & occupancy_abundance_2$betweenness > betweenness_SD[2,4], TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Deer" & occupancy_abundance_2$betweenness > betweenness_SD[3,4], TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Carollia" & occupancy_abundance_2$betweenness > betweenness_SD[4,4], TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Spinyrat" & occupancy_abundance_2$betweenness > betweenness_SD[5,4], TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Mouselemur" & occupancy_abundance_2$betweenness > betweenness_SD[6,4], TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Flamingo" & occupancy_abundance_2$betweenness > betweenness_SD[7,4], TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Stint" & occupancy_abundance_2$betweenness > betweenness_SD[8,4], TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(is.na(occupancy_abundance_2$linker), FALSE,occupancy_abundance_2$linker )


head(occupancy_abundance_2)

##############################################

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

## look at ditributions of variables

ggplot(occupancy_abundance_2, aes(x = RelPrev))+
  geom_histogram(alpha = 0.6, col = "black",fill = "grey")+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 14))+
  xlab("Occupancy frequency")+
  ylab("Count")+
  facet_wrap(~host_species, ncol = 4, scales = "free_y")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())


ggplot(occupancy_abundance_2, aes(x = log(RelAbundance)))+
  geom_histogram(alpha = 0.6, col = "black",fill = "grey")+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 14))+
  xlab("Log abundance")+
  ylab("Count")+
  facet_wrap(~host_species, ncol = 4, scales = "free_y")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())




ggplot(occupancy_abundance_2, aes(x = sum_pos))+
  geom_histogram(alpha = 0.6, col = "black",fill = "grey")+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 14))+
  xlab("No. positive associations")+
  ylab("Count")+
  facet_wrap(~host_species, ncol = 4, scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

ggplot(occupancy_abundance_2, aes(x = sum_neg))+
  geom_histogram(alpha = 0.6, col = "black",fill = "grey")+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 14))+
  xlab("No. negative associations")+
  ylab("Count")+
  facet_wrap(~host_species, ncol = 4, scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

ggplot(occupancy_abundance_2, aes(x = Freq_won))+
  geom_histogram(alpha = 0.6, col = "black",fill = "grey")+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 14))+
  xlab("No. negative associations won")+
  ylab("Count")+
  facet_wrap(~host_species, ncol = 4, scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

### supplementary figure 1

ggplot(occupancy_abundance_2, aes(x = betweenness))+
  geom_histogram(alpha = 0.6, col = "black",aes(fill = linker))+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 18))+
  xlab("Betweenness")+
  ylab("Count")+
  facet_wrap(~host_species, ncol = 4, scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  scale_fill_manual(values = c('grey','red'))+theme(legend.position = 'none')


###########################################################################################################################################
###########################################################################################################################################
##################################### STATISTICS ##########################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################


####### PREDICTION ONE - OCCUPANCY FREQUENCY

############# global model

global_occupancy<-glm(cbind(Prevalence, Sample_size - Prevalence) ~
    
    host_species*log(RelAbundance)+
    host_species*sum_neg+
    host_species*sum_pos+
    log(MeanAbundance)*sum_neg+
    log(MeanAbundance)*sum_pos,
  
  family = binomial("logit"), 
  data = occupancy_abundance_2, 
  na.action = "na.fail")



########## check psudo R2

R2logit<- function(y,model){
  R2<- 1-(model$deviance/model$null.deviance)
  return(R2)
}

R2logit(Prevalence,global_occupancy)


##################### AIV comparison ### supplementary table S1

AIC_results_occupancy<- dredge(global_occupancy)

########## supplemtary table S1

head(AIC_results_occupancy, 20)

#write.csv(AIC_results[1:20], "AIC_results_occupancy.csv")


########################### top model is global model


model_occupancy2<-glm(cbind(Prevalence, Sample_size - Prevalence) ~
    
    host_species*log(RelAbundance)+
    host_species*sum_neg+
    host_species*sum_pos+
    log(RelAbundance)*sum_neg+
    log(RelAbundance)*sum_pos,
  
  family = binomial("logit"), 
  data = occupancy_abundance_2, 
  na.action = "na.fail")



########## check psudo R2


R2logit(Prevalence,model_occupancy2)


#############################################

summary(model_occupancy2)

#convert to odds ratio

model.df <- tidy(model_occupancy2)  # Convert model to dataframe for easy manipulation
model.df

model_estimates_occupancy<-model.df %>% 
  mutate(or = exp(estimate),  # Odds ratio/gradient
    var.diag = diag(vcov(model_occupancy2)),  # Variance of each coefficient
    or.se = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 

########### supplementary table S6

model_estimates_occupancy
print(tbl_df(model_estimates_occupancy), n=40)

#write.csv(model_estimates_occupancy, "model_results_occupancy_editabundance.csv")

############################### model fit


par(mfrow=c(2,2))
plot(model_occupancy2)
dev.off()


# visualise model R2

occupancy_abundance_2$predicted_values_restricted<-as.numeric(fitted(model_occupancy2))

ggplot(occupancy_abundance_2, aes(x = RelPrev, y = predicted_values_restricted))+
  geom_point( alpha = 0.5)+
  facet_wrap(~host_species,  ncol = 4)+
  geom_smooth(method = "lm")+theme_bw()+
  scale_size(range = c(1,5))+
  ylab("Predicted occupancy frequency")+
  xlab("Observed occupancy frequency")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

summary(lm(occupancy_abundance_2$RelPrev~occupancy_abundance_2$predicted_values_restricted))$adj.r.squared #62



################################################################### plot



### generate predictions for abundance = 0.03%
occupancy_abundance_2$predicted_values_0.0003<-predict(model_occupancy2, newdata=base::transform(occupancy_abundance_2,  RelAbundance=runif(4141, 0.00025,0.00035)))

### generate predictions for abundance = 1%

occupancy_abundance_2$predicted_values_0.01<-predict(model_occupancy2,  newdata=base::transform(occupancy_abundance_2, RelAbundance=runif(4141,0.005, 0.015)))


######## convert for logit function to predicted occupancy 

occupancy_abundance_2$predicted_values_0.0003<-exp(occupancy_abundance_2$predicted_values_0.0003)/(1+exp(occupancy_abundance_2$predicted_values_0.0003))
occupancy_abundance_2$predicted_values_0.01<-exp(occupancy_abundance_2$predicted_values_0.01)/(1+exp(occupancy_abundance_2$predicted_values_0.01))


######## long format

prediction_plots_df <- gather(occupancy_abundance_2, Abundance_control, Predicted_occupancy, c(predicted_values_0.0003:predicted_values_0.01, predicted_values_0.0003), factor_key=TRUE)
head(prediction_plots_df)

#unique(prediction_plots_df$Abundance_control)
prediction_plots_df$Abundance_control<-factor(prediction_plots_df$Abundance_control, levels = c("predicted_values_0.01", "predicted_values_0.0003"))


ggplot(prediction_plots_df, aes(y = Predicted_occupancy, x = sum_pos))+
  
  geom_jitter(alpha = 0.7, size = 0.7, width = 0.1)+
  geom_smooth(method = "lm") +
  facet_wrap(~Abundance_control+host_species, ncol = 8, scales = "free")+
  #ggtitle("Phylogenetic uniquness ~ Predicted occupancy when controlling for abundance")+
  theme_classic()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=12))+
  #scale_y_log10(labels = point) +
  theme(axis.text = element_text(size = 8))+
  theme(axis.title = element_text(size = 16))+
  xlab("No. positive associations")+
  ylab("Predicted occupancy frequency")+
  theme(legend.position = "none")

ggplot(prediction_plots_df, aes(y = Predicted_occupancy, x = sum_neg))+
  
  geom_jitter(alpha = 0.7, size = 0.7, width = 0.1)+
  geom_smooth(method = "lm") +
  facet_wrap(~Abundance_control+host_species, ncol = 8, scales = "free")+
  #ggtitle("Phylogenetic uniquness ~ Predicted occupancy when controlling for abundance")+
  theme_classic()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=12))+
  #scale_y_log10(labels = point) +
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 16))+
  xlab("No. negative associations")+
  ylab("Predicted occupancy frequency")+
  theme(legend.position = "none")



ggplot(prediction_plots_df, aes(y = Predicted_occupancy, x = linker))+
  
  geom_boxplot()+
  geom_smooth(method = "lm") +
  facet_wrap(~Abundance_control+host_species, ncol = 8, scales = "free")+
  #ggtitle("Phylogenetic uniquness ~ Predicted occupancy when controlling for abundance")+
  theme_classic()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=12))+
  #scale_y_log10(labels = point) +
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 16))+
  xlab("Probability of having high betweenness centrality")+
  ylab("Predicted occupancy frequency")+
  theme(legend.position = "none")

################################################################################################


######################## FIGURE 1A AND B

NotFancy <- function(l) {
  l <- format(l, scientific = FALSE)
  parse(text=l)
}

ggplot(occupancy_abundance_2, aes(x = RelAbundance, y = RelPrev))+
  geom_point(alpha = 0.7, width = 0.1)+
  scale_x_log10(labels = NotFancy)+theme_bw()+
  xlab("Relative abundance")+
  ylab("Occupancy frequency")+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 14))




ggplot(occupancy_abundance_2, aes(x = RelPrev))+
  geom_histogram(alpha = 0.6, col = "black",fill = "grey")+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 14))+
  xlab("Occupancy frequency")+
  ylab("Count")+
  facet_wrap(~host_species, ncol = 4, scales = "free_y")+
  theme_bw()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

## Supplemtary figure S3

ggplot(occupancy_abundance_2, aes(x = RelAbundance, y = RelPrev))+
  geom_point(alpha = 0.7)+
  scale_x_log10(labels = NotFancy)+theme_bw()+
  xlab("Relative abundance")+
  ylab("Occupancy frequency")+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 14))+
  facet_wrap(~host_species, ncol = 4)+
  geom_vline(xintercept = 0.01, linetype = "dashed", size = 1.5, col = "orange")+
  geom_vline(xintercept = 0.0003,  linetype = "dashed", size = 1.5, col = "orange")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())




#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

## repeat but exclude red deer, because this data set has exponentally larger numbers of positive and negative edges and therefore influence model


model_restrict<-glm(cbind(Prevalence, Sample_size - Prevalence) ~
    
    host_species*log(RelAbundance)+
    host_species*sum_neg+
    host_species*sum_pos+
    log(RelAbundance)*sum_neg+
    log(RelAbundance)*sum_pos,
  
  family = binomial("logit"), 
  data = occupancy_abundance_reduced, 
  na.action = "na.fail")



########################### top model statistics (top model is global model)


summary(model_restrict)

#convert to odds ratio

model.df <- tidy(model_restrict)  # Convert model to dataframe for easy manipulation
model.df

model_estimates_occupancy_NoDeer<-model.df %>% 
  mutate(or = exp(estimate),  # Odds ratio/gradient
    var.diag = diag(vcov(model_restrict)),  # Variance of each coefficient
    or.se = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 

########### supplementary table S7

model_estimates_occupancy_NoDeer

#write.csv(model_estimates_occupancy_NoDeer, "model_results_occupancy_NoDeer.csv")





##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################


####### test hyporthsis that specialists should have lower prevalence

specialist_analysis<-subset(occupancy_abundance_2, host_species=="Deer" | host_species == "Meerkat")

ggplot(specialist_analysis, aes(x = Freq_won))+geom_histogram()+facet_wrap(~host_species)

zeros<-subset(specialist_analysis, Freq_won == 0)

table(zeros$host_species)
table(specialist_analysis$host_species)
hist(specialist_analysis$Freq_won)

#########################################################################

### fit model 

global_specialist_model<-glm(Freq_won ~ 
    RelPrev +
    log(RelAbundance)+
    sum_neg +
    host_species+
    sum_pos+
    host_species*log(RelAbundance)+
    host_species*RelPrev+
    host_species*linker,
  data = specialist_analysis, 
  family = poisson(link = "log"),
  na.action = "na.fail")



######## AIC comparison

summary(global_specialist_model)

### supplementary table S2

AIC_comparison<-dredge(global_specialist_model)

head(AIC_comparison, 20)


#write.csv(AIC_comparison[1:20], "AIC_results_Specialists.csv")

############################## TOP MODEL (TOP MODEL doen't have betweenness/linker variables)


specialist_model<-glm(Freq_won ~ 
    RelPrev +
    log(RelAbundance)+
    
    sum_neg +
    host_species+
    sum_pos+
    host_species*RelPrev+
    host_species*log(RelAbundance),
  
  data = specialist_analysis, 
  family = poisson(link = "log"),
  na.action = "na.fail")


summary(specialist_model)

#check over dispersion
performance::check_overdispersion(specialist_model)


AIC_comparison<-dredge(specialist_model)

head(AIC_comparison, 20)

#convert to odds ratio

model.df <- tidy(specialist_model)  # Convert model to dataframe for easy manipulation
model.df

model_estimates_highcompetitors<-model.df %>% 
  mutate(or = exp(estimate),  # Odds ratio/gradient
    var.diag = diag(vcov(specialist_model)),  # Variance of each coefficient
    or.se = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 

######### supplementary table S8

model_estimates_highcompetitors


#write.csv(model_estimates_highcompetitors, "model_results_specialists.csv")

#################################################################

############################### model fit


par(mfrow=c(2,2))
plot(specialist_model)
dev.off()

########## check psudo R2


R2logit(Freq_won,specialist_model)


# visualise model R2
# visualise model R2

specialist_analysis$predicted_values_highcompetitors<-as.numeric(fitted(specialist_model))

ggplot(specialist_analysis, aes(x = Freq_won, y = predicted_values_highcompetitors))+
  geom_point( alpha = 0.5)+
  facet_wrap(~host_species,  ncol = 4, scales = "free")+
  geom_smooth(method = "lm")+theme_bw()+
  scale_size(range = c(1,5))+
  ylab("Predicted no. directed negative associations")+
  xlab("Observed no. directed negative associations")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

### R2

summary(lm(specialist_analysis$predicted_values_highcompetitors~specialist_analysis$Freq_won))$adj.r.squared #62


################################################# PLOT

############# FIGURE 3


ggplot(specialist_analysis, aes(x = RelAbundance, y = predicted_values_highcompetitors))+
  geom_point( alpha = 0.5)+
  
  facet_wrap(~host_species,  ncol = 4, scales = "free")+
  geom_smooth(method = "lm")+theme_bw()+
  scale_x_log10(labels = NotFancy)+
  ylab("Predicted no. directed negative associations")+
  xlab("Relative abundance")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())


mean(specialist_analysis$RelAbundance)

## HOLD ABUNDANCE AT 1%


specialist_analysis$predicted_values_specialist_0.01<-predict(specialist_model, newdata=base::transform(specialist_analysis, RelAbundance=0.01))


ggplot(specialist_analysis, aes(x = RelPrev, y = predicted_values_specialist_0.01))+
  geom_point( alpha = 0.5)+
  facet_wrap(~host_species,  ncol = 4, scales = 'free')+
  geom_smooth(method = "lm")+theme_bw()+
  ylab("Predicted no. directed negative associations")+
  xlab("Occupancy frequency")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())





#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################



################# what about taxa characterised by many positive associations??

generalist_model<-glm(sum_pos ~ 
    RelPrev+
    log(RelAbundance)+
    sum_neg+
    host_species*log(RelAbundance)+
    host_species*RelPrev+
    host_species*linker,
  
  data = occupancy_abundance_2, 
  family = poisson(link = "log"),
  na.action = "na.fail")



######## AIC comparison

summary(generalist_model)

#plot(allEffects(generalist_model))

### supplementary table S3

AIC_comparison<-dredge(generalist_model)

head(AIC_comparison, 20)

#write.csv(AIC_comparison[1:20], "AIC_results_generalists.csv")

############################## TOP MODEL (TOP MODEL IS GLOBAL MODEL)

########## check psudo R2


R2logit(sum_pos,generalist_model)


#convert to odds ratio

model.df <- tidy(generalist_model)  # Convert model to dataframe for easy manipulation
model.df

generalist_model_estimates<-model.df %>% 
  mutate(or = exp(estimate),  # Odds ratio/gradient
    var.diag = diag(vcov(generalist_model)),  # Variance of each coefficient
    or.se = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 

######### supplementary table S9

print(tbl_df(generalist_model_estimates), n=40)

#write.csv(generalist_model_estimates, "model_results_generalists.csv")


############################### model fit

########## check psudo R2


R2logit(sum_pos,generalist_model)


##################

par(mfrow=c(2,2))
plot(generalist_model)
dev.off()


# visualise model R2

occupancy_abundance_2$predicted_values_generalists<-as.numeric(fitted(generalist_model))

ggplot(occupancy_abundance_2, aes(x = sum_pos, y = predicted_values_generalists))+
  geom_point( alpha = 0.5, size = 1)+
  # facet_wrap(~host_species,  ncol = 4, scales = "free")+
  geom_smooth(method = "lm")+theme_bw()+
  scale_size(range = c(1,5))+
  ylab("Predicted no. positive associations")+
  xlab("Observed no. positive associations")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

### R2

summary(lm(occupancy_abundance_2$predicted_values_generalists~occupancy_abundance_2$sum_pos))$adj.r.squared #62

##########################################

############# PLOTS


## HOLD ABUNDANCE AT 1% and linker = TRUE

newdata=transform(occupancy_abundance_2, linker=TRUE)
newdata1=transform(newdata, RelAbundance=runif(4141,0.005,0.015))
newdata2=transform(newdata, RelAbundance=runif(4141,0.0003,0.0004))

newdata3=transform(occupancy_abundance_2, RelAbundance=runif(4141,0.005,0.015))
newdata4=transform(occupancy_abundance_2, RelAbundance=runif(4141,0.0003,0.0004))

occupancy_abundance_2$predicted_values_generalist_linker<-predict(generalist_model, newdata=newdata)
occupancy_abundance_2$predicted_values_generalist_0.01<-predict(generalist_model, newdata=newdata1)
occupancy_abundance_2$predicted_values_generalist_0.0003<-predict(generalist_model, newdata=newdata2)

occupancy_abundance_2$predicted_values_generalist_0.01X<-predict(generalist_model, newdata=newdata3)
occupancy_abundance_2$predicted_values_generalist_0.0003X<-predict(generalist_model, newdata=newdata4)


####### supplementary figure S5




ggplot(occupancy_abundance_2, aes(x = RelPrev, y = predicted_values_generalist_linker))+
  geom_point( alpha = 0.5, size = 1)+
  facet_wrap(~host_species,  ncol =4, scales = "free")+
  geom_smooth(method = "lm")+
  theme_bw()+
  ylab("Predicted no. positive associations")+
  xlab("Occupancy frequency")+
  scale_x_continuous(labels =NotFancy)+
  theme(axis.text = element_text(size = 8))+
  theme(axis.title = element_text(size = 14))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())


ggplot(occupancy_abundance_2, aes(x = RelAbundance, y = predicted_values_generalist_linker))+
  geom_point( alpha = 0.5, size = 1)+
  facet_wrap(~host_species,  ncol =4, scales = "free")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_log10(labels = NotFancy)+
  ylab("Predicted no. positive associations")+
  xlab("Relative abundance")+
  theme(axis.text = element_text(size = 8))+
  theme(axis.title = element_text(size = 14))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

###### Figure 4

ggplot(occupancy_abundance_2, aes(x = linker, y = predicted_values_generalists))+
  geom_jitter(alpha = 0.3, width = 0.2, size = 1)+
  geom_boxplot(fill = 'lightgrey', alpha = 0.5)+
  facet_wrap(~host_species,  ncol =4, scales = 'free_y')+
  geom_smooth(method = "lm")+
  theme_bw()+
  ylab("Predicted no. positive associations")+
  xlab("Betweenness centrality")+
  theme(axis.text = element_text(size = 8))+
  theme(axis.title = element_text(size = 14))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

#### testing competiton related-hypothesis

edge_properties1<- read.csv("C:/Users/Alice_R/Dropbox/Sommer postdoc/Core microbiome project/Analysis/MERGED/edge_data_risely.csv", sep=";")

####################################################################

edge_positive<-subset(edge_properties1, dir == "pos")
edge_negative<-subset(edge_properties1, dir == "neg")

table(edge_negative$host_species)
table(edge_properties1$host_species)


model_edge_positive<- lm(weight~ 
    Relatedness+
    host_species+
    RelAbundance_diff+
    RelAbundance_sum+
    RelPrev_diff+
    RelPrev_sum+
    degree_diff+
    degree_sum,
  data = edge_positive,
  na.action = "na.fail")



######## AIC comparison


#plot(allEffects(model_edge_positive))

######## supplmentary table S4

AIC_comparison<-dredge(model_edge_positive)

model.avg(AIC_comparison, delta<3)

#write.csv(AIC_comparison[1:20], "AIC_results_Relatedness.csv")

##############################

## top model

model_edge_positive<- lm(weight~ 
    Relatedness+
    host_species+
    RelAbundance_sum ,
  data = edge_positive,
  na.action = "na.fail")



summary(model_edge_positive)
plot(allEffects(model_edge_positive))

#convert to odds ratio

model.df <- tidy(model_edge_positive)  # Convert model to dataframe for easy manipulation
model.df

model_estimates_positive<-model.df %>% 
  mutate(or = exp(estimate),  # Odds ratio/gradient
    var.diag = diag(vcov(model_edge_positive)),  # Variance of each coefficient
    or.se = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 


model_estimates_positive


#write.csv(odds_ratio, "model_results_relatedness_pos.csv")

######## R2


edge_positive$predicted_values_positive<-as.numeric(fitted(model_edge_positive))

summary(lm(edge_positive$weight~edge_positive$predicted_values_positive))$adj.r.squared 

##### model explanatory power pretty much nada




#################################################################


model_edge_negative<- lm(weight~ 
    Relatedness+
    host_species+
    RelAbundance_diff+
    RelAbundance_sum+
    RelPrev_diff+
    RelPrev_sum+
    degree_diff+
    degree_sum,
  data = edge_negative,
  na.action = "na.fail")



######## AIC comparison


#plot(allEffects(model_edge_positive))

### supplemtary table S5

AIC_comparison<-dredge(model_edge_negative)

model.avg(AIC_comparison, delta<3)


### no good predictors.

#write.csv(AIC_comparison[1:20], "AIC_results_Relatedness_neg.csv")

##############################



