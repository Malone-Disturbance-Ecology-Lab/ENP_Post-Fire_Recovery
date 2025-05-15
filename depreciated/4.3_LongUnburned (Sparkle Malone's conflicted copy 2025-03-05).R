
library(tidyverse)
library(randomForest)
library(ggplot2)
library(ggpubr)


setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
setwd("./Baseline")
load( file="BL_train_test.RDATA")

setwd('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline')
load(file="NDVI_rf.RDATA")

NDVI_rf

# Add NDVI Model to the training and testing dataframes. 
train <- train %>% mutate( NDVI.Model = predict(NDVI_rf,  cur_data()))
test <- test %>% mutate( NDVI.Model = predict(NDVI_rf,  cur_data()))

# Extract the long unburned (TF < 2 and TSF >= 15)
# Based on the training data:
# TF.2 <- train %>% filter( TotalFires < 2  )
#TF.2$ptID %>% length/ train$ptID %>% length *100
#long.unburned <- TF.2 %>% filter( TSF >= 15  ) 
#long.unburned$ptID %>% length/ train$ptID %>% length *100

# Based on the test data:
TF.2 <- test %>% filter( TotalFires < 2  )
TF.2$ptID %>% length/ test$ptID %>% length *100
long.unburned <- TF.2 %>% filter( TSF > 15  ) 
long.unburned$ptID %>% length/ test$ptID %>% length *100

ideal <- train %>% filter( TotalFires == 9) %>% reframe(NIR.SWIR1 = mean(NIR.SWIR1), 
                                                        SWIR1.SWIR2 = mean(SWIR1.SWIR2))
# Set the target for the long unburned:
long.unburned.target <-long.unburned %>% filter( TotalFires < 2, TSF > 15) %>% mutate(
  TotalFires = 9,
  NIR.SWIR1 = ideal$NIR.SWIR1,
  SWIR1.SWIR2= ideal$SWIR1.SWIR2,
  NDVI.Model = predict(NDVI_rf,  cur_data()))


long.unburned
lm( NDVI.Model ~ NDVI, data =train) %>% summary
lm( NDVI.Model ~ NDVI, data =long.unburned) %>% summary
lm( NDVI.Model ~ NDVI, data =long.unburned.target) %>% summary

p.ndvi.1 <- long.unburned %>% ggplot(aes( x= NDVI, y= NDVI.Model)) + 
  geom_point( aes( x= NDVI, y= NDVI.Model), alpha=0.1, ) + 
  stat_smooth(method = "lm", se=FALSE, color="goldenrod1", formula = y ~ x) + 
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(..p.label..), sep = "~`,`~")), # adds R^2 and p-value
           r.accuracy = 0.01,
           p.accuracy = 0.001,
           label.x = 0.2, label.y = 0.6, size = 6) +
  stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
                        label.x = 0.2, label.y = 0.5, size =6) +
  geom_abline(intercept = 0, slope = 1, col = 'red',linetype="dashed") +
  xlab("Observed NDVI") + ylab('NDVI Model')  + theme_bw()+
  theme(text=element_text(size=20))

p.ndvi.2 <-long.unburned.target %>% ggplot(aes( x= NDVI, y= NDVI.Model)) + 
  geom_point( aes( x= NDVI, y= NDVI.Model), alpha=0.1) + 
  stat_smooth(method = "lm", se=FALSE, color="goldenrod1", formula = y ~ x) + 
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(..p.label..), sep = "~`,`~")), # adds R^2 and p-value
           r.accuracy = 0.01,
           p.accuracy = 0.001,
           label.x = 0.2, label.y = 0.6, size = 6) +
  stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
                        label.x = 0.2, label.y = 0.5, size =6) +
  geom_abline(intercept = 0, slope = 1, col = 'red',linetype="dashed") +
  xlab("Observed NDVI") + ylab(expression('NDVI'[Enhanced])) + theme_bw() +
  theme(text=element_text(size=20)) 

setwd('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/Figures')

png('LongUnburned_NDVI.png',
    width = 600, height = 300)
ggarrange(p.ndvi.1, p.ndvi.2, labels=c("A", "B"),
          nrow=1, ncol=2 )
dev.off()

# Explore how the long unburned impacts the recovery dataframe: ####

setwd('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod')
# Load data
load("./Recovery/Recov_Rates.RDATA")

summary(Recov_Rate)

long.unburned.rec <- Recov_Rate %>% filter(Prev.Int > 15, 
                                                     TotalFires <2)
# How much of the recovery data is long unburned?
long.unburned.rec$ptID %>% length /Recov_Rate$ptID %>% length*100


long.unburned.rec.high <- long.unburned.rec %>% filter(PreNDVI >= model.NDVI)

# How much of the long unburned had a higher pre-NDVI than target?
long.unburned.rec.high$ptID %>% length/long.unburned.rec$ptID %>% length*100

# How much of the long unburned had a higher pre-NDVI than target?
long.unburned.rec.high$ptID %>% length/Recov_Rate$ptID %>% length *100

mean(long.unburned.rec.high$Rec100_Yrs, na.rm=T)
mean(long.unburned.rec.high$Rec80_Yrs, na.rm=T)
mean(long.unburned.rec.high$Rec70_Yrs, na.rm=T)

mean(Recov_Rate$Rec100_Yrs, na.rm=T)
mean(Recov_Rate$Rec80_Yrs, na.rm=T)
mean(Recov_Rate$Rec70_Yrs, na.rm=T)



df.1 <- long.unburned.rec.high %>% select(ptID, Rec70_Yrs) %>% mutate(Recovery= "70%" , Rec = Rec70_Yrs) %>% select(ptID, Rec, Recovery)%>% mutate(Fire = 'Long Unburned')
df.2 <- long.unburned.rec.high %>% select(ptID, Rec80_Yrs) %>% mutate(Recovery= "80%" , Rec = Rec80_Yrs) %>% select(ptID, Rec, Recovery)%>% mutate(Fire = 'Long Unburned')
df.3 <- long.unburned.rec.high %>% select(ptID, Rec90_Yrs) %>% mutate(Recovery= "90%" , Rec = Rec90_Yrs) %>% select(ptID, Rec, Recovery)%>% mutate(Fire = 'Long Unburned')
df.4 <- long.unburned.rec.high %>% select(ptID, Rec100_Yrs) %>% mutate(Recovery= "100%" , Rec = Rec100_Yrs) %>% select(ptID, Rec, Recovery)%>% mutate(Fire = 'Long Unburned')


df.1f <- Recov_Rate %>% filter(Prev.Int < 15, 
                               TotalFires >2) %>%  select(ptID, Rec70_Yrs) %>% mutate(Recovery= "70%" , Rec = Rec70_Yrs) %>% select(ptID, Rec, Recovery) %>% mutate(Fire = 'Maintained')
df.2f <- Recov_Rate %>% filter(Prev.Int < 15, 
                              TotalFires >2) %>% select(ptID, Rec80_Yrs) %>% mutate(Recovery= "80%" , Rec = Rec80_Yrs) %>% select(ptID, Rec, Recovery)%>% mutate(Fire = 'Maintained')
df.3f <- Recov_Rate %>% filter(Prev.Int < 15, 
                               TotalFires >2) %>% select(ptID, Rec90_Yrs) %>% mutate(Recovery= "90%" , Rec = Rec90_Yrs) %>% select(ptID, Rec, Recovery)%>% mutate(Fire = 'Maintained')
df.4f <- Recov_Rate %>% filter(Prev.Int < 15, 
                               TotalFires >2)%>% select(ptID, Rec100_Yrs) %>% mutate(Recovery= "100%" , Rec = Rec100_Yrs) %>% select(ptID, Rec, Recovery)%>% mutate(Fire = 'Maintained')


longUn <- rbind( df.1, df.2, df.3, df.4, df.1f, df.2f, df.3f, df.4f)
longUn$Recovery <- factor(longUn$Recovery , levels=
                            c("70%", "80%", "90%", "100%"))

longUn$Recovery %>% levels
plot.LUNB1 <- longUn %>%  ggplot( aes(x=Recovery, y=Rec, col = Fire)) + 
  geom_violin(trim=FALSE) + stat_summary(fun.y=mean, geom="point", shape=1, size=1) +
  ylab("Recovery Time (Years)") + theme_bw() + xlab("Recovery Threshold") +
  scale_color_manual(values=c("black", "goldenrod1" )) +
  theme(text=element_text(size=20))

# Resulted in recovery rates much faster than what is observed on the landscape. 

setwd('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/Figures')

png('LongUnburned_RecoveryTime.png',
    width = 400, height = 300)
plot.LUNB1
dev.off()
