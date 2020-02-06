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

occupancy_abundance_2$host_species<-factor(occupancy_abundance_2$host_species, levels = c("Human", "Meerkat", "Deer", "Carollia","Spinyrat","Mouselemur","Flamingo","Stint"))

head(occupancy_abundance_2)


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################


############# global model


model_restrict<-glm(cbind(Prevalence, Sample_size - Prevalence) ~
                      
                      host_species*log(MeanAbundance)+
                      host_species*sum_neg+
                      host_species*sum_pos+
                      log(MeanAbundance)*sum_neg+
                      log(MeanAbundance)*sum_pos,
                    
                    family = binomial("logit"), 
                    data = occupancy_abundance_2, 
                    na.action = "na.fail")


##################### AIV comparison ### supplementary table S1

AIC_results_occupancy<- dredge(model_restrict)


########## supplemtary table S2

head(AIC_results_occupancy, 20)

#write.csv(AIC_results[1:20], "AIC_results_occupancy.csv")


########################### top model statistics (top model is global model)


summary(model_restrict)

#convert to odds ratio

model.df <- tidy(model_restrict)  # Convert model to dataframe for easy manipulation
model.df

model_estimates_occupancy<-model.df %>% 
  mutate(or = exp(estimate),  # Odds ratio/gradient
         var.diag = diag(vcov(model_restrict)),  # Variance of each coefficient
         or.se = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 

########### supplementary table S5

model_estimates_occupancy

#write.csv(odds_ratio, "model_results_occupancy.csv")



############################### model fit


par(mfrow=c(2,2))
plot(model_restrict)

dev.off()


# visualise model R2

occupancy_abundance_2$predicted_values_restricted<-as.numeric(fitted(model_restrict))

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
occupancy_abundance_2$predicted_values_0.0003<-predict(model_restrict, newdata=transform(occupancy_abundance_2, MeanAbundance=0.0003))

### generate predictions for abundance = 1%

occupancy_abundance_2$predicted_values_0.01<-predict(model_restrict, newdata=transform(occupancy_abundance_2, MeanAbundance=0.01))


######## convert for logit function to predicted occupancy 

occupancy_abundance_2$predicted_values_0.0003<-exp(occupancy_abundance_2$predicted_values_0.0003)/(1+exp(occupancy_abundance_2$predicted_values_0.0003))
occupancy_abundance_2$predicted_values_0.01<-exp(occupancy_abundance_2$predicted_values_0.01)/(1+exp(occupancy_abundance_2$predicted_values_0.01))


######## long format

prediction_plots_df <- gather(occupancy_abundance_2, Abundance_control, Predicted_occupancy, c(predicted_values_0.0003:predicted_values_0.01, predicted_values_0.0003), factor_key=TRUE)
head(prediction_plots_df)

unique(prediction_plots_df$Abundance_control)


ggplot(prediction_plots_df, aes(y = Predicted_occupancy, x = sum_pos))+
  
  geom_jitter(alpha = 0.4)+
  geom_smooth(method = "lm") +
  facet_wrap(~Abundance_control+host_species, ncol = 8, scales = "free")+
  #ggtitle("Phylogenetic uniquness ~ Predicted occupancy when controlling for abundance")+
  theme_classic()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=12))+
  #scale_y_log10(labels = point) +
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 16))+
  xlab("No. positive associations")+
  ylab("Predicted occupancy frequency")+
  theme(legend.position = "none")

ggplot(prediction_plots_df, aes(y = Predicted_occupancy, x = sum_neg))+
  
  geom_jitter(alpha = 0.4)+
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



######################## FIGURE 1A AND B

NotFancy <- function(l) {
  l <- format(l, scientific = FALSE)
  parse(text=l)
}

ggplot(occupancy_abundance_2, aes(x = MeanAbundance, y = RelPrev))+
  geom_point(alpha = 0.7)+
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



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################


####### test hyporthsis that specialists should have lower prevalence

specialist_analysis<-subset(occupancy_abundance_2, host_species=="Deer" | host_species == "Meerkat")


### fit model 

specialist_model<-glm(Freq_won ~ 
                        RelPrev +
                        log(MeanAbundance)+
                        sum_neg+
                        host_species+
                        sum_pos,
                      data = specialist_analysis, 
                      family = poisson(link = "log"),
                      na.action = "na.fail")




#plot(specialist_model)

######## AIC comparison

summary(specialist_model)

plot(allEffects(specialist_model))

### supplementary table S2

AIC_comparison<-dredge(specialist_model)

#write.csv(AIC_comparison[1:20], "AIC_results_Freq_won.csv")

############################## TOP MODEL (TOP MODEL IS GLOBAL MODEL)


#convert to odds ratio

model.df <- tidy(specialist_model)  # Convert model to dataframe for easy manipulation
model.df

model_estimates_highcompetitors<-model.df %>% 
  mutate(or = exp(estimate),  # Odds ratio/gradient
         var.diag = diag(vcov(specialist_model)),  # Variance of each coefficient
         or.se = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 

######### supplementary table S6

model_estimates_highcompetitors


#write.csv(odds_ratio, "model_results_freq_won.csv")

#################################################################

############################### model fit


par(mfrow=c(2,2))
plot(specialist_model)
dev.off()


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

############# FIGURE 4


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
  facet_wrap(~host_species,  ncol = 4, scales = "free")+
  geom_smooth(method = "lm")+theme_bw()+
  ylab("Predicted no. directed negative associations")+
  xlab("Occupancy frequency")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
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

######## supplmentary table S3

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

### supplemtary table S4

AIC_comparison<-dredge(model_edge_negative)

model.avg(AIC_comparison, delta<3)


### no good predictors.

#write.csv(AIC_comparison[1:20], "AIC_results_Relatedness_neg.csv")

##############################

