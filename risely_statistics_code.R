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


occupancy_abundance_2<- read.csv("C:/Users/Alice_R/Dropbox/Sommer postdoc/Core microbiome project/Analysis/MERGED/occupancy_data_risely.csv", sep=";")

#order levels

occupancy_abundance_2$host_species<-factor(occupancy_abundance_2$host_species, levels = c("Human", "Meerkat", "Deer", "Carollia","Spinyrat","Mouselemur","Flamingo","Stint"))

head(occupancy_abundance_2)

#generate data without deer, to test whether the deer dataset is driving any trends (due to much higher number of positive and negative edges)

occupancy_abundance_reduced<-subset(occupancy_abundance_2, host_species != "Deer")


########## add column showing whether ASV has high betweenness (is a linker)

betweenness_SD<-ddply(occupancy_abundance_2, .(host_species), summarize, mean=mean(betweenness), std_dev=sd(betweenness))
betweenness_SD$std_dev_upper<-betweenness_SD$mean+betweenness_SD$std_dev



occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Human" & occupancy_abundance_2$betweenness > 3862, TRUE, NA)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Meerkat" & occupancy_abundance_2$betweenness > 4113, TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Deer" & occupancy_abundance_2$betweenness > 14635, TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Carollia" & occupancy_abundance_2$betweenness > 722, TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Spinyrat" & occupancy_abundance_2$betweenness > 6153, TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Mouselemur" & occupancy_abundance_2$betweenness > 1061, TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Flamingo" & occupancy_abundance_2$betweenness > 785, TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(occupancy_abundance_2$host_species=="Stint" & occupancy_abundance_2$betweenness > 3070, TRUE, occupancy_abundance_2$linker)
occupancy_abundance_2$linker<-ifelse(is.na(occupancy_abundance_2$linker), FALSE,occupancy_abundance_2$linker )




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


ggplot(occupancy_abundance_2, aes(x = log(MeanAbundance)))+
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
                        
                      host_species*log(MeanAbundance)+
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
                      
                      host_species*log(MeanAbundance)+
                      host_species*sum_neg+
                      host_species*sum_pos+
                      log(MeanAbundance)*sum_neg+
                      log(MeanAbundance)*sum_pos,
                    
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
occupancy_abundance_2$predicted_values_0.0003<-predict(model_occupancy2, newdata=transform(occupancy_abundance_2, 
                                                          MeanAbundance=runif(4141, 0.00025,0.00035)))

### generate predictions for abundance = 1%

occupancy_abundance_2$predicted_values_0.01<-predict(model_occupancy2, 
                                                     newdata=transform(occupancy_abundance_2, MeanAbundance=runif(4141,0.005, 0.015)))


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

ggplot(occupancy_abundance_2, aes(x = MeanAbundance, y = RelPrev))+
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

ggplot(occupancy_abundance_2, aes(x = MeanAbundance, y = RelPrev))+
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
                      
                      host_species*log(MeanAbundance)+
                      host_species*sum_neg+
                      host_species*sum_pos+
                      log(MeanAbundance)*sum_neg+
                      log(MeanAbundance)*sum_pos,
                    
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
max(specialist_analysis$Freq_won)

#########################################################################

### fit model 

global_specialist_model<-glm(Freq_won ~ 
                        RelPrev +
                        log(MeanAbundance)+
                        sum_neg +
                        host_species+
                        sum_pos+
                        host_species*log(MeanAbundance)+
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
                        log(MeanAbundance)+
                        sum_neg +
                        host_species+
                        sum_pos+
                        host_species*log(MeanAbundance)+
                        host_species*RelPrev,
                      data = specialist_analysis, 
                      family = poisson(link = "log"),
                      na.action = "na.fail")


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


ggplot(specialist_analysis, aes(x = MeanAbundance, y = predicted_values_highcompetitors))+
  geom_point( alpha = 0.5)+
  
  facet_wrap(~host_species,  ncol = 4, scales = "free")+
  geom_smooth(method = "lm")+theme_bw()+
  scale_x_log10(labels = NotFancy)+
  ylab("Predicted no. directed negative associations")+
  xlab("Relative abundance")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())


mean(specialist_analysis$MeanAbundance)

## HOLD ABUNDANCE AT 1%


specialist_analysis$predicted_values_specialist_0.01<-predict(specialist_model, newdata=transform(specialist_analysis, MeanAbundance=0.01))


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
                        log(MeanAbundance)+
                        sum_neg+
                        host_species*log(MeanAbundance)+
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
newdata1=transform(newdata, MeanAbundance=runif(4141,0.005,0.015))
newdata2=transform(newdata, MeanAbundance=runif(4141,0.0003,0.0004))

newdata3=transform(occupancy_abundance_2, MeanAbundance=runif(4141,0.005,0.015))
newdata4=transform(occupancy_abundance_2, MeanAbundance=runif(4141,0.0003,0.0004))

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


ggplot(occupancy_abundance_2, aes(x = MeanAbundance, y = predicted_values_generalist_linker))+
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
                           MeanAbundance_diff+
                           MeanAbundance_sum+
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
                           MeanAbundance_sum ,
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
                           MeanAbundance_diff+
                           MeanAbundance_sum+
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



