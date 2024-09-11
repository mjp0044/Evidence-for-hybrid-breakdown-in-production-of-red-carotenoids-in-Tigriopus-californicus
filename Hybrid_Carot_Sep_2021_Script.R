#### Tigriopus parental and hybrid carotenoid mito analysis #### 
#Updated last on 9-2-21

library(tidyverse)
library(lme4)
library(cowplot)
library(lmerTest)
library(agricolae)
library(MASS)
library(emmeans)
library(ggpubr)
library(MuMIn)
library(car)
library(reshape2)
library(plotrix)
library(ggjoy)

#Model analysis using lmer class objects

dat <- read.csv("astax.ril_R.csv")

dat %>% 
  mutate(sex = fct_recode(sex, "Females" = "F",
                          "Males" = "M")) -> dat

#Subset out experimental groups of interest
dat.algae = subset(dat, food == "Algae")

#Set up data for plotting crosses with BR, SD, and BUF
dat %>% 
  filter(food == "Algae",
         cross %in% c( "BUFxSD",
                       "BUFxBUF",
                       "BRxSD",
                       "BRxBR",
                       "SDxSD"))-> all.dat


all.dat$ril.id<- factor(all.dat$ril.id,
                        levels = c("BUF", "BUFSD19", "BUFSD24","BUFSD4", "SD", 
                                   "BRSD45", "BRSD50", "BRSD56","BR"))
all.dat$group <-as.factor(all.dat$group)


#Set up data for plotting crosses with PES, AB, and CAT
dat %>% 
  filter(food == "Algae",
         cross %in% c( "ABxAB",
                       "PESxPES",
                       "CATxCAT",
                       "PESxAB",
                       "CATxAB", 
                       "ABxCAT"))-> all.dat2

all.dat2$ril.id<- factor(all.dat2$ril.id,
                        levels = c("CAT", "CATAB27", "ABCAT11","AB", "PESAB20", 
                                   "PES"))
all.dat2$group <-as.factor(all.dat2$group)


#Filter data set to set up groups of individual crosses
#Set up group 1 (BUF x SD)
dat %>% 
  filter(food == "Algae",
         cross %in% c( "BUFxSD",
                       "BUFxBUF",
                        "SDxSD"))-> group.1

group.1$ril.id<- factor(group.1$ril.id,
                        levels = c("BUF", "BUFSD19", "BUFSD24","BUFSD4", "SD"))

#Set up group 2 (BR x SD)
dat %>% 
  filter(food == "Algae",
         cross %in% c( "BRxSD",
                       "BRxBR",
                       "SDxSD"))-> group.2

group.2$ril.id<- factor(group.2$ril.id, 
                        levels = c("BR", "BRSD45", "BRSD50","BRSD56", "SD"))

#Set up group 3 (CAT x AB)
dat %>% 
  filter(food == "Algae",
         cross %in% c( "ABxAB",
                       "ABxCAT",
                       "CATxAB", "CATxCAT"))-> group.3

group.3$ril.id<- factor(group.3$ril.id, 
                        levels = c("AB", "ABCAT11", "CATAB27", "CAT"))

#Set up group 4 (PES x AB) with only males
dat %>% 
  filter(food == "Algae", 
         cross %in% c( "ABxAB",
                       "PESxAB", "PESxPES"))-> group.4

group.4$ril.id<- factor(group.4$ril.id, 
                        levels = c("AB", "PESAB20", "PES"))

#### Comparing Parental vs Hybrid Astaxanthin cross by cross ####

### put the crosses in a specific order
group.1$cross<- factor(group.1$cross,
                        levels = c("BUFxBUF", "BUFxSD", "SDxSD"))

###Comparing BUF, SD, and hybrids averaged
# set the contrasts we want to test, from a priori prediction that hybrids will be lower than parentals

m.1 <- rbind(c(1,-1,0), # This compares the first cross (BUFxBUF) to the sencond cross (BUFxSD))
           c(0,-1,1)) # This compares the third cross (SDxSD) to the second cross (BUFxSD)

invm.1 <-ginv(m.1) # create a generalize inverse matrix 

mod.1 <- lm(astax.mg ~ cross , data = subset(group.1, sex == "Males"),
              contrasts = list(cross = invm.1))
summary(mod.1)

### BUF, SD, and their individual Hybrids
mod.2 <-lm(astax.mg ~ ril.id*sex, data = group.1)

### Making multiple comparisons and correcting for type 1 error
emmeans(mod.2, pairwise ~ ril.id | sex)

comp.g1<- emmeans(mod.2, pairwise ~ ril.id | sex)
write.csv(comp.g1$emmeans, file = "group1.plot.means.csv", row.names = FALSE)

#### Comparing BR, SD, and their hybrids averaged
### put the crosses in a specific order
group.2$cross<- factor(group.2$cross,
                       levels = c("BRxBR", "BRxSD", "SDxSD"))

# set the contrasts we want to test, from a priori prediction that hybrids will be lower than parentals
m.3 <- rbind(c(1,-1,0), # This compares the first cross (BRxBR) to the sencond cross (BRxSD))
             c(0,-1,1)) # This compares the third cross (SDxSD) to the second cross (BRxSD)

invm.3 <-ginv(m.3) # create a generalize inverse matrix 

mod.3 <- lm(astax.mg ~ cross , data = subset(group.2, sex =="Males"),
              contrasts = list(cross = invm.3))
summary(mod.3)


#### Comparing BR and SD Parentals to each RIL
mod.4 <-lm(astax.mg ~ ril.id * sex, data = group.2)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.4, pairwise ~ ril.id | sex)

comp.g2<- emmeans(mod.4, pairwise ~ ril.id | sex)
write.csv(comp.g2$emmeans, file = "group2.plot.means.csv", row.names = FALSE)


###Comparing CAT, AB, and hybrids averaged (not enough female data to include, so just male data. No sex effect coded)
#Fit model
mod.5 <- lmer(astax.mg ~ cross + (1|ril.id), data = group.3)
summary(mod.5)
#Relevel reference group to either CATxCAT or ABxAB
group.3$cross = relevel(group.3$cross, ref = "ABxAB")

### CAT, AB, and their individual Hybrids
mod.6 <-lm(astax.mg ~ ril.id, data = group.3)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.6, pairwise ~ ril.id )
confint(emmeans(mod.6, pairwise ~ ril.id ))

comp.g3<- emmeans(mod.6, pairwise ~ ril.id)
write.csv(comp.g3$emmeans, file = "group3.plot.means.csv", row.names = FALSE)

###Comparing PES, AB, and hybrids averaged
#Fit model
mod.7 <- lmer(astax.mg ~ cross + (1|ril.id), data = group.4)
summary(mod.7)
#Relevel reference group to either CATxCAT or ABxAB
group.4$cross = relevel(group.4$cross, ref = "ABxAB")

### PES, AB, and their individual Hybrids
mod.8 <-lm(astax.mg ~ ril.id * sex, data = group.4)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.8, pairwise ~ ril.id | sex)
confint(emmeans(mod.8, pairwise ~ ril.id | sex))

comp.g4<- emmeans(mod.8, pairwise ~ ril.id | sex)
write.csv(comp.g4$emmeans, file = "group4.plot.means.csv", row.names = FALSE)

####Make figure for BUF, BR and SD

#Cross.plot.means.csv made from combining group plot means from above into single summary excel file. 
plot.dat <- read.csv("cross.plot.means.asta.csv")

plot.dat$ril.id<- factor(plot.dat$ril.id,
                        levels = c("BUF", "BUFSD19", "BUFSD24","BUFSD4", "SD", 
                                   "BRSD45", "BRSD50", "BRSD56","BR"))
plot.dat %>% 
  mutate(sex = fct_recode(sex, "Females" = "F",
                          "Males" = "M")) -> plot.dat

#Make plot for comparing BUF, SD, and BR with individual RILs

pbbsmales <- subset(all.dat, sex == "Males") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = subset(plot.dat, sex == "Males"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat, sex =="Males"), mapping=aes(xmin=1, xmax=4.8, ymin=0.361, ymax=0.446), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat, sex =="Males"), mapping=aes(xmin=1.2, xmax=8.8, ymin=0.303, ymax=0.389), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat, sex =="Males"), mapping=aes(xmin=5.2, xmax=9.0, ymin=0.147, ymax=0.245), fill="slateblue1", alpha=0.03, linetype = "blank") +
  
  geom_linerange(data= subset(plot.dat, sex == "Males"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = astax.mg, col = group), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(0,0.75)+
  scale_color_manual(values = c("chartreuse4","orangered1", "lightsteelblue4","slateblue1","goldenrod3" ))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")

pbbsfemales <- subset(all.dat, sex == "Females") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = subset(plot.dat, sex == "Females"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Female confidence bars
  geom_rect(data= subset(plot.dat, sex =="Females"), mapping=aes(xmin=1, xmax=4.8, ymin=0.162, ymax=0.248), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat, sex =="Females"), mapping=aes(xmin=1.2, xmax=8.8, ymin=0.1045, ymax=0.226), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat, sex =="Females"), mapping=aes(xmin=5.2, xmax=9.0, ymin=0.142, ymax=0.240), fill="slateblue1", alpha=0.03, linetype = "blank") +
  
  geom_linerange(data= subset(plot.dat, sex == "Females"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = astax.mg, col = group), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  scale_color_manual(values = c("chartreuse4","orangered1", "lightsteelblue4","slateblue1","goldenrod3" ))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")

####Make figure for AB, CAT, and PES

#Cross.plot.means2.csv made from combining group plot means from above into single summary excel file. 
plot.dat2 <- read.csv("cross.plot.means2.asta.csv")

plot.dat2$ril.id<- factor(plot.dat2$ril.id,
                         levels = c("CAT", "CATAB27", "ABCAT11","AB", "PESAB20", 
                                    "PES"))

#Make plot for comparing AB, CAT, and PES with individual RILs
all.dat2 %>%
  filter(sex == "Males")-> all.dat2.males


pacpmales <- subset(all.dat2, sex =="Males") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data= subset(plot.dat2, sex == "Males"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat2, sex == "Males"), mapping=aes(xmin=1, xmax=3.8, ymin=0.098, ymax=0.365), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat2, sex == "Males"), mapping=aes(xmin=1.2, xmax=5.8, ymin=0.211, ymax=0.344), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat2, sex == "Males"), mapping=aes(xmin=4.2, xmax=6, ymin=0.342, ymax=0.466), fill="slateblue1", alpha=0.03, linetype = "blank") +
  geom_linerange(data= subset(plot.dat2, sex == "Males"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = astax.mg, col = ril.id), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(0,0.75)+
  scale_color_manual(values = c("chartreuse4","orangered1", "purple" , "lightsteelblue4","goldenrod3","slateblue1"))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("Cross") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")

pacpfemales <- subset(all.dat2, sex =="Females") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data= subset(plot.dat2, sex == "Females"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat2, sex == "Females"), mapping=aes(xmin=1, xmax=2, ymin=0.219, ymax=0.336), fill="slateblue1", alpha=0.07, linetype = "blank") +
  geom_linerange(data= subset(plot.dat2, sex == "Females"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = astax.mg, col = ril.id), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(0.1,0.4)+
  scale_color_manual(values = c("goldenrod3","slateblue1"))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("Cross") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")

figuremales <- ggarrange(pbbsmales, pacpmales, ncol=1, nrow=2)

jpeg(file = "Cross Carot Males for Figure 2.jpg", units = "in", width = 8.5, height = 6, res = 500)
figuremales
dev.off()

figurefemales <- ggarrange(pbbsfemales, pacpfemales, labels = c("A", "B"), ncol=1, nrow=2)

tiff(file = "Fig S4.tif", units = "in", width = 8.5, height = 6, res = 300)
figurefemales
dev.off()

#Overall hybrid vs parental
mod.overall <- lmer(astax.cop~gen + (1|ril.id), data = subset(dat, sex == "Males"))
summary(mod.overall)
confint(mod.overall)

gen <- c("Hybrid","Parental")
mean <- c(0.24138, 0.29641)
lower.cl <- c(0.16368346,0.20435962)
upper.cl <- c(0.3189984,0.3883551)

plot.dat.overall <- data.frame(gen, mean, lower.cl, upper.cl)

jpeg(file = "Overall.jpg", units = "in", width = 5, height = 5, res = 500)
subset(dat, sex =="Males")  %>% 
  ggplot(aes(x = gen))+
  geom_point(data= plot.dat.overall, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.overall, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = astax.mg, col = gen), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  theme_cowplot()+
  xlab("All Lines") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")+
  scale_color_manual(values=c("red","blue"))+
  annotate("text", 1.5, 0.6, label = expression(italic(beta)*"=0.055,"~italic(p)*"=0.393"), size = 5)
dev.off()

#BR, BUF, SD hybrid vs parental
mod.overall.BBS <- lmer(astax.mg~gen + (1|ril.id), data = subset(all.dat, sex == "Males"))
summary(mod.overall.BBS)
confint(mod.overall.BBS)

gen <- c("Hybrid","Parental")
mean <- c(0.16056, 0.31506)
lower.cl <- c(0.09086198,0.22016042)
upper.cl <- c(0.22874362,0.40996139)

plot.dat.overall.BBS <- data.frame(gen, mean, lower.cl, upper.cl)

jpeg(file = "Overall BR BUF SD MALES.jpg", units = "in", width = 5, height = 5, res = 500)
subset(all.dat, sex =="Males") %>% 
  ggplot(aes(x = gen))+
  geom_point(data= plot.dat.overall.BBS, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.overall.BBS, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = astax.mg, col = gen), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  ylim(0,0.75)+
  theme_cowplot()+
  xlab("BR, BUF, SD") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.position = "")+
  scale_color_manual(values=c("red","blue"))+
  annotate("text", 1.5, 0.6, label = expression(italic(beta)*"=0.155,"~italic(p)*"=0.041"), size = 6)
dev.off()

#AB, CAT, PES hybrid vs parental
mod.overall.ACP <- lmer(astax.mg~gen + (1|ril.id), data = subset(all.dat2, sex == "Males"))
summary(mod.overall.ACP)
confint(mod.overall.ACP)

gen <- c("Hybrid","Parental")
mean <- c(0.41564, 0.30630)
lower.cl <- c(0.32459726,0.22103449)
upper.cl <- c(0.5095524,0.3906706)

plot.dat.overall.ACP <- data.frame(gen, mean, lower.cl, upper.cl)

jpeg(file = "Overall AB-CAT-PES.jpg", units = "in", width = 5, height = 5, res = 500)
subset(all.dat2, sex =="Males") %>% 
  ggplot(aes(x = gen))+
  geom_point(data= plot.dat.overall.ACP, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.overall.ACP, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = astax.mg, col = gen), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  ylim(0,0.75)+
  theme_cowplot()+
  xlab("AB, CAT, PES") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.position = "")+
  scale_color_manual(values=c("red","blue"))+
  annotate("text", 1.5, 0.6, label = expression(italic(beta)*"=-0.109,"~italic(p)*"=0.169"), size = 6)
dev.off()








####Repeat comparison across lines using dietary carotenoids b-carotene and hydroxyechineneon/echinenone as point of comparison####

### BUF, SD, and their individual Hybrids
mod.2.dietary <-lm(total.dietary ~ ril.id*sex, data = group.1)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.2.dietary, pairwise ~ ril.id | sex)
confint(emmeans(mod.2.dietary, pairwise ~ ril.id | sex))

#### Comparing BR and SD Parentals to each RIL
mod.4.dietary <-lm(total.dietary ~ ril.id * sex, data = group.2)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.4.dietary, pairwise ~ ril.id | sex)
confint(emmeans(mod.4.dietary, pairwise ~ ril.id | sex))

### CAT, AB, and their individual Hybrids
mod.6.dietary <-lm(total.dietary ~ ril.id, data = group.3)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.6.dietary, pairwise ~ ril.id )
confint(emmeans(mod.6.dietary, pairwise ~ ril.id ))

### PES, AB, and their individual Hybrids
mod.8.dietary <-lm(total.dietary ~ ril.id * sex, data = group.4)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.8.dietary, pairwise ~ ril.id | sex)
confint(emmeans(mod.8.dietary, pairwise ~ ril.id | sex))

#Take means and CI's from above tables and combine into data sheet for plotting
#Read plotting data sheet in to R
plot.dat.dietary <- read.csv("cross.plot.means.dietary.csv")

plot.dat.dietary$ril.id<- factor(plot.dat.dietary$ril.id,
                         levels = c("BUF", "BUFSD19", "BUFSD24","BUFSD4", "SD", 
                                    "BRSD45", "BRSD50", "BRSD56","BR"))

plot.dat.dietary %>% 
  mutate(sex = fct_recode(sex, "Females" = "F",
                          "Males" = "M")) -> plot.dat.dietary

#Make plot for comparing BUF, SD, and BR with individual RILs
pbbsmales.dietary <- subset(all.dat, sex == "Males") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = subset(plot.dat.dietary, sex == "Males"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat.dietary, sex =="Males"), mapping=aes(xmin=1, xmax=4.8, ymin=0.05562, ymax=0.0921), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.dietary, sex =="Males"), mapping=aes(xmin=1.2, xmax=8.8, ymin=0.02451, ymax=0.0609), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.dietary, sex =="Males"), mapping=aes(xmin=5.2, xmax=9.0, ymin=0.06493, ymax=0.0938), fill="slateblue1", alpha=0.03, linetype = "blank") +
  
  geom_linerange(data= subset(plot.dat.dietary, sex == "Males"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = total.dietary, col = group), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(0,0.2)+
  scale_color_manual(values = c("chartreuse4","orangered1", "lightsteelblue4","slateblue1","goldenrod3" ))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

pbbsfemales.dietary <- subset(all.dat, sex == "Females") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = subset(plot.dat.dietary, sex == "Females"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Female confidence bars
  geom_rect(data= subset(plot.dat.dietary, sex =="Females"), mapping=aes(xmin=1, xmax=4.8, ymin=0.02005, ymax=0.0565), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.dietary, sex =="Females"), mapping=aes(xmin=1.2, xmax=8.8, ymin=0.00428, ymax=0.0558), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.dietary, sex =="Females"), mapping=aes(xmin=5.2, xmax=9.0, ymin=0.0184, ymax=0.0472), fill="slateblue1", alpha=0.03, linetype = "blank") +
  
  geom_linerange(data= subset(plot.dat.dietary, sex == "Females"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = total.dietary, col = group), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  scale_color_manual(values = c("chartreuse4","orangered1", "lightsteelblue4","slateblue1","goldenrod3" ))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

####Make figure for AB, CAT, and PES

#Cross.plot.means2.csv made from combining group plot means from above into single summary excel file. 
plot.dat2.dietary <- read.csv("cross.plot.means2.dietary.csv")

plot.dat2.dietary$ril.id<- factor(plot.dat2.dietary$ril.id,
                          levels = c("CAT", "CATAB27", "ABCAT11","AB", "PESAB20", 
                                     "PES"))

#Make plot for comparing AB, CAT, and PES with individual RILs
all.dat2 %>%
  filter(sex == "Males")-> all.dat2.males


pacpmales.dietary <- subset(all.dat2, sex =="Males") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data= subset(plot.dat2.dietary, sex == "Males"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat2.dietary, sex == "Males"), mapping=aes(xmin=1, xmax=3.8, ymin=0.0613, ymax=0.1104), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat2.dietary, sex == "Males"), mapping=aes(xmin=1.2, xmax=5.8, ymin=0.1091, ymax=0.15), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat2.dietary, sex == "Males"), mapping=aes(xmin=4.2, xmax=6, ymin=0.0817, ymax=0.1226), fill="slateblue1", alpha=0.03, linetype = "blank") +
  geom_linerange(data= subset(plot.dat2.dietary, sex == "Males"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = total.dietary, col = ril.id), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(0,0.2)+
  scale_color_manual(values = c("chartreuse4","orangered1", "purple" , "lightsteelblue4","goldenrod3","slateblue1"))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("Cross") +
  ylab(bquote('Dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

pacpfemales.dietary <- subset(all.dat2, sex =="Females") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data= subset(plot.dat2.dietary, sex == "Females"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat2.dietary, sex == "Females"), mapping=aes(xmin=1, xmax=2, ymin=0.0211, ymax=0.0668), fill="slateblue1", alpha=0.07, linetype = "blank") +
  geom_linerange(data= subset(plot.dat2.dietary, sex == "Females"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = total.dietary, col = ril.id), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(0.0,0.2)+
  scale_color_manual(values = c("goldenrod3","slateblue1"))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("Cross") +
  ylab(bquote('Dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

figuremales.dietary <- ggarrange(pbbsmales.dietary, pacpmales.dietary, ncol=1, nrow=2)

jpeg(file = "Cross Carot Dietary Males.jpg", units = "in", width = 8.5, height = 6, res = 500)
figuremales.dietary
dev.off()

figurefemales.dietary <- ggarrange(pbbsfemales.dietary, pacpfemales.dietary, labels = c("A", "B"), ncol=1, nrow=2)

tiff(file = "Cross Carot Dietary Females.tif", units = "in", width = 8.5, height = 6.5, res = 300)
figurefemales.dietary
dev.off()



#Overall hybrid vs parental in dietary carotenoids

#BR, BUF, SD hybrid vs parental
mod.overall.BBS.dietary <- lmer(total.dietary~gen + (1|ril.id), data = subset(all.dat, sex == "Males"))
summary(mod.overall.BBS.dietary)
confint(mod.overall.BBS.dietary)

gen <- c("Hybrid","Parental")
mean <- c(0.07928, 0.06530)
lower.cl <- c(0.05528392,0.03341783)
upper.cl <- c(0.10252494,0.09719133)

plot.dat.overall.BBS.dietary <- data.frame(gen, mean, lower.cl, upper.cl)

jpeg(file = "Overall BR BUF SD MALES dietary.jpg", units = "in", width = 5, height = 5, res = 500)
subset(all.dat, sex =="Males") %>% 
  ggplot(aes(x = gen))+
  geom_point(data= plot.dat.overall.BBS.dietary, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.overall.BBS.dietary, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = total.dietary, col = gen), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  ylim(0,0.2)+
  theme_cowplot()+
  xlab("BR, BUF, SD") +
  ylab(bquote('Dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.position = "")+
  scale_color_manual(values=c("red","blue"))+
  annotate("text", 1.5, 0.18, label = expression(italic(beta)*"=-0.014,"~italic(p)*"=0.521"), size = 6)
dev.off()

#AB, CAT, PES hybrid vs parental
mod.overall.ACP.dietary <- lmer(total.dietary~gen + (1|ril.id), data = subset(all.dat2, sex == "Males"))
summary(mod.overall.ACP.dietary)
confint(mod.overall.ACP.dietary)

gen <- c("Hybrid","Parental")
mean <- c(0.08799, 0.10621)
lower.cl <- c(0.05954632,0.07989336)
upper.cl <- c(0.11557163,0.13235752)

plot.dat.overall.ACP.dietary <- data.frame(gen, mean, lower.cl, upper.cl)

jpeg(file = "Overall AB-CAT-PES dietary.jpg", units = "in", width = 5, height = 5, res = 500)
subset(all.dat2, sex =="Males") %>% 
  ggplot(aes(x = gen))+
  geom_point(data= plot.dat.overall.ACP.dietary, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.overall.ACP.dietary, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = total.dietary, col = gen), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  ylim(0,0.2)+
  theme_cowplot()+
  xlab("AB, CAT, PES") +
  ylab(bquote('Dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.position = "")+
  scale_color_manual(values=c("red","blue"))+
  annotate("text", 1.5, 0.18, label = expression(italic(beta)*"=0.182,"~italic(p)*"=0.423"), size = 6)
dev.off()



####Comparison across lines in ratio of asta to dietary carotenoids####

### BUF, SD, and their individual Hybrids
mod.2.ratio <-lm(asta.to.dietary ~ ril.id*sex, data = group.1)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.2.ratio, pairwise ~ ril.id | sex)
confint(emmeans(mod.2.ratio, pairwise ~ ril.id | sex))

#### Comparing BR and SD Parentals to each RIL
mod.4.ratio <-lm(asta.to.dietary ~ ril.id * sex, data = group.2)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.4.ratio, pairwise ~ ril.id | sex)
confint(emmeans(mod.4.ratio, pairwise ~ ril.id | sex))

### CAT, AB, and their individual Hybrids
mod.6.ratio <-lm(asta.to.dietary ~ ril.id, data = group.3)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.6.ratio, pairwise ~ ril.id )
confint(emmeans(mod.6.ratio, pairwise ~ ril.id ))

### PES, AB, and their individual Hybrids
mod.8.ratio <-lm(asta.to.dietary ~ ril.id * sex, data = group.4)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.8.ratio, pairwise ~ ril.id | sex)
confint(emmeans(mod.8.ratio, pairwise ~ ril.id | sex))


#Take means and CI's from above tables and combine into data sheet for plotting
#Read plotting data sheet in to R
plot.dat.ratio <- read.csv("cross.plot.means.ratio.csv")

plot.dat.ratio$ril.id<- factor(plot.dat.ratio$ril.id,
                                 levels = c("BUF", "BUFSD19", "BUFSD24","BUFSD4", "SD", 
                                            "BRSD45", "BRSD50", "BRSD56","BR"))

plot.dat.ratio %>% 
  mutate(sex = fct_recode(sex, "Females" = "F",
                          "Males" = "M")) -> plot.dat.ratio

#Make plot for comparing BUF, SD, and BR with individual RILs
pbbsmales.ratio <- subset(all.dat, sex == "Males") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = subset(plot.dat.ratio, sex == "Males"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat.ratio, sex =="Males"), mapping=aes(xmin=1, xmax=4.8, ymin=4.532, ymax=6.74), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.ratio, sex =="Males"), mapping=aes(xmin=1.2, xmax=8.8, ymin=7.2225, ymax=9.43), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.ratio, sex =="Males"), mapping=aes(xmin=5.2, xmax=9.0, ymin=1.413, ymax=3.67), fill="slateblue1", alpha=0.03, linetype = "blank") +
  
  geom_linerange(data= subset(plot.dat.ratio, sex == "Males"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = asta.to.dietary, col = group), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(-2.2,22)+
  scale_color_manual(values = c("chartreuse4","orangered1", "lightsteelblue4","slateblue1","goldenrod3" ))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Astaxanthin:dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

pbbsfemales.ratio <- subset(all.dat, sex == "Females") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = subset(plot.dat.ratio, sex == "Females"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Female confidence bars
  geom_rect(data= subset(plot.dat.ratio, sex =="Females"), mapping=aes(xmin=1, xmax=4.8, ymin=4.5308, ymax=6.55), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.ratio, sex =="Females"), mapping=aes(xmin=1.2, xmax=8.8, ymin=4.0024, ymax=7.12), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.ratio, sex =="Females"), mapping=aes(xmin=5.2, xmax=9.0, ymin=5.098, ymax=7.53), fill="slateblue1", alpha=0.03, linetype = "blank") +
  
  geom_linerange(data= subset(plot.dat.ratio, sex == "Females"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = asta.to.dietary, col = group), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(0,10)+
  scale_color_manual(values = c("chartreuse4","orangered1", "lightsteelblue4","slateblue1","goldenrod3" ))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Astaxanthin:dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

####Make figure for AB, CAT, and PES

#Cross.plot.means2.csv made from combining group plot means from above into single summary excel file. 
plot.dat2.ratio <- read.csv("cross.plot.means2.ratio.csv")

plot.dat2.ratio$ril.id<- factor(plot.dat2.ratio$ril.id,
                                  levels = c("CAT", "CATAB27", "ABCAT11","AB", "PESAB20", 
                                             "PES"))

#Make plot for comparing AB, CAT, and PES with individual RILs
all.dat2 %>%
  filter(sex == "Males")-> all.dat2.males


pacpmales.ratio <- subset(all.dat2, sex =="Males") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data= subset(plot.dat2.ratio, sex == "Males"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat2.ratio, sex == "Males"), mapping=aes(xmin=1, xmax=3.8, ymin=-2.17, ymax=7.57), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat2.ratio, sex == "Males"), mapping=aes(xmin=1.2, xmax=5.8, ymin=1.21, ymax=3.23), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat2.ratio, sex == "Males"), mapping=aes(xmin=4.2, xmax=6, ymin=3.12, ymax=5.14), fill="slateblue1", alpha=0.03, linetype = "blank") +
  geom_linerange(data= subset(plot.dat2.ratio, sex == "Males"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = asta.to.dietary, col = ril.id), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(-2.2,22)+
  scale_color_manual(values = c("chartreuse4","orangered1", "purple" , "lightsteelblue4","goldenrod3","slateblue1"))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("Cross") +
  ylab(bquote('Astaxanthin:dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

pacpfemales.ratio <- subset(all.dat2, sex =="Females") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data= subset(plot.dat2.ratio, sex == "Females"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat2.ratio, sex == "Females"), mapping=aes(xmin=1, xmax=2, ymin=5.5, ymax=7.76), fill="slateblue1", alpha=0.07, linetype = "blank") +
  geom_linerange(data= subset(plot.dat2.ratio, sex == "Females"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = asta.to.dietary, col = ril.id), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(0.0,10)+
  scale_color_manual(values = c("goldenrod3","slateblue1"))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("Cross") +
  ylab(bquote('Astaxanthin:dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

figuremales.ratio <- ggarrange(pbbsmales.ratio, pacpmales.ratio, ncol=1, nrow=2)

jpeg(file = "Cross Carot Ratio Males.jpg", units = "in", width = 8.5, height = 6, res = 500)
figuremales.ratio
dev.off()

figurefemales.ratio <- ggarrange(pbbsfemales.ratio, pacpfemales.ratio, labels = c("A", "B"), ncol=1, nrow=2)

tiff(file = "Cross Carot Ratio Females.tif", units = "in", width = 8.5, height = 7, res = 300)
figurefemales.ratio
dev.off()



#Overall hybrid vs parental in carotenoid ratio

#BR, BUF, SD hybrid vs parental
mod.overall.BBS.ratio <- lmer(asta.to.dietary~gen + (1|ril.id), data = subset(all.dat, sex == "Males"))
summary(mod.overall.BBS.ratio)
confint(mod.overall.BBS.ratio)

#Run as block
gen <- c("Hybrid","Parental")
mean <- c(2.3197, 5.4995)
lower.cl <- c(0.7878819,3.3958620)
upper.cl <- c(3.838796,7.603170)
plot.dat.overall.BBS.ratio <- data.frame(gen, mean, lower.cl, upper.cl)

jpeg(file = "Overall BR BUF SD MALES ratio.jpg", units = "in", width = 5, height = 5, res = 500)
subset(all.dat, sex =="Males") %>% 
  ggplot(aes(x = gen))+
  geom_point(data= plot.dat.overall.BBS.ratio, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.overall.BBS.ratio, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = asta.to.dietary, col = gen), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  ylim(-2.2,22)+
  theme_cowplot()+
  xlab("BR, BUF, SD") +
  ylab(bquote('Astaxanthin:dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 16),
        legend.position = "")+
  scale_color_manual(values=c("red","blue"))+
  annotate("text", 1.5, 20, label = expression(italic(beta)*"=3.18,"~italic(p)*"=0.051"), size = 6)
dev.off()

#AB, CAT, PES hybrid vs parental
mod.overall.ACP.ratio <- lmer(asta.to.dietary~gen + (1|ril.id), data = subset(all.dat2, sex == "Males"))
summary(mod.overall.ACP.ratio)
confint(mod.overall.ACP.ratio)

#Run as block
gen <- c("Hybrid","Parental")
mean <- c(6.858, 3.020)
lower.cl <- c(3.1275762,-0.5725172)
upper.cl <- c(10.713982,6.609991)
plot.dat.overall.ACP.ratio <- data.frame(gen, mean, lower.cl, upper.cl)

jpeg(file = "Overall AB-CAT-PES ratio.jpg", units = "in", width = 5, height = 5, res = 500)
subset(all.dat2, sex =="Males") %>% 
  ggplot(aes(x = gen))+
  geom_point(data= plot.dat.overall.ACP.ratio, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.overall.ACP.ratio, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = asta.to.dietary, col = gen), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  ylim(-2.2,22)+
  theme_cowplot()+
  xlab("AB, CAT, PES") +
  ylab(bquote('Astaxanthin:dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 16),
        legend.position = "")+
  scale_color_manual(values=c("red","blue"))+
  annotate("text", 1.5, 20, label = expression(italic(beta)*"=-3.84,"~italic(p)*"=0.254"), size = 6)
dev.off()

####Comparison across lines in ratio of asta to hydroxyechinenone ####

### BUF, SD, and their individual Hybrids
mod.2.hydroxy <-lm(asta.to.hydroxy ~ ril.id*sex, data = group.1)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.2.hydroxy, pairwise ~ ril.id | sex)

#### Comparing BR and SD Parentals to each RIL
mod.4.hydroxy <-lm(asta.to.hydroxy ~ ril.id * sex, data = group.2)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.4.hydroxy, pairwise ~ ril.id | sex)

### CAT, AB, and their individual Hybrids
mod.6.hydroxy <-lm(asta.to.hydroxy ~ ril.id, data = group.3)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.6.hydroxy, pairwise ~ ril.id )

### PES, AB, and their individual Hybrids
mod.8.hydroxy <-lm(asta.to.hydroxy ~ ril.id * sex, data = group.4)
### Making multiple comparisons and correcting for type 1 error
emmeans(mod.8.hydroxy, pairwise ~ ril.id | sex)


#Take means and CI's from above tables and combine into data sheet for plotting
#Read plotting data sheet in to R
plot.dat.hydroxy <- read.csv("cross.plot.means.hydroxy.csv")

plot.dat.hydroxy$ril.id<- factor(plot.dat.hydroxy$ril.id,
                               levels = c("BUF", "BUFSD19", "BUFSD24","BUFSD4", "SD", 
                                          "BRSD45", "BRSD50", "BRSD56","BR"))

plot.dat.hydroxy %>% 
  mutate(sex = fct_recode(sex, "Females" = "F",
                          "Males" = "M")) -> plot.dat.hydroxy

#Make plot for comparing BUF, SD, and BR with individual RILs
pbbsmales.hydroxy <- subset(all.dat, sex == "Males") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = subset(plot.dat.hydroxy, sex == "Males"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat.hydroxy, sex =="Males"), mapping=aes(xmin=1, xmax=4.8, ymin=10.6293, ymax=13.38), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.hydroxy, sex =="Males"), mapping=aes(xmin=1.2, xmax=8.8, ymin=8.5582, ymax=11.31), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.hydroxy, sex =="Males"), mapping=aes(xmin=5.2, xmax=9.0, ymin=4.03, ymax=6.56), fill="slateblue1", alpha=0.03, linetype = "blank") +
  geom_linerange(data= subset(plot.dat.hydroxy, sex == "Males"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = asta.to.hydroxy, col = group), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(-2.2,22)+
  scale_color_manual(values = c("chartreuse4","orangered1", "lightsteelblue4","slateblue1","goldenrod3" ))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Astaxanthin:Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

pbbsfemales.hydroxy <- subset(all.dat, sex == "Females") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = subset(plot.dat.hydroxy, sex == "Females"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Female confidence bars
  geom_rect(data= subset(plot.dat.hydroxy, sex =="Females"), mapping=aes(xmin=1, xmax=4.8, ymin=7.2745, ymax=10.03), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.hydroxy, sex =="Females"), mapping=aes(xmin=1.2, xmax=8.8, ymin=4.6345, ymax=8.53), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat.hydroxy, sex =="Females"), mapping=aes(xmin=5.2, xmax=9.0, ymin=5.58, ymax=8.11), fill="slateblue1", alpha=0.03, linetype = "blank") +
  
  geom_linerange(data= subset(plot.dat.hydroxy, sex == "Females"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = asta.to.hydroxy, col = group), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(0,12)+
  scale_color_manual(values = c("chartreuse4","orangered1", "lightsteelblue4","slateblue1","goldenrod3" ))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Astaxanthin:Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

####Make figure for AB, CAT, and PES

#Cross.plot.means2.csv made from combining group plot means from above into single summary excel file. 
plot.dat2.hydroxy <- read.csv("cross.plot.means2.hydroxy.csv")

plot.dat2.hydroxy$ril.id<- factor(plot.dat2.hydroxy$ril.id,
                                levels = c("CAT", "CATAB27", "ABCAT11","AB", "PESAB20", 
                                           "PES"))

#Make plot for comparing AB, CAT, and PES with individual RILs
all.dat2 %>%
  filter(sex == "Males")-> all.dat2.males


pacpmales.hydroxy <- subset(all.dat2, sex =="Males") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data= subset(plot.dat2.hydroxy, sex == "Males"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat2.hydroxy, sex == "Males"), mapping=aes(xmin=1, xmax=3.8, ymin=-1.615, ymax=8.41), fill="chartreuse4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat2.hydroxy, sex == "Males"), mapping=aes(xmin=1.2, xmax=5.8, ymin=2.22, ymax=4.77), fill="lightsteelblue4", alpha=0.03, linetype = "blank") +
  geom_rect(data= subset(plot.dat2.hydroxy, sex == "Males"), mapping=aes(xmin=4.2, xmax=6, ymin=6.64, ymax=9.19), fill="slateblue1", alpha=0.03, linetype = "blank") +
  geom_linerange(data= subset(plot.dat2.hydroxy, sex == "Males"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = asta.to.hydroxy, col = ril.id), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(-2.2,22)+
  scale_color_manual(values = c("chartreuse4","orangered1", "purple" , "lightsteelblue4","goldenrod3","slateblue1"))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("Cross") +
  ylab(bquote('Astaxanthin:Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

pacpfemales.hydroxy <- subset(all.dat2, sex =="Females") %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data= subset(plot.dat2.hydroxy, sex == "Females"), aes(y = emmean, shape = gen),size = 7, alpha = 0.7)+
  #Male confidence bars
  geom_rect(data= subset(plot.dat2.hydroxy, sex == "Females"), mapping=aes(xmin=1, xmax=2, ymin=6.02, ymax=8.87), fill="slateblue1", alpha=0.07, linetype = "blank") +
  geom_linerange(data= subset(plot.dat2.hydroxy, sex == "Females"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = asta.to.hydroxy, col = ril.id), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  ylim(0.0,12)+
  scale_color_manual(values = c("goldenrod3","slateblue1"))+
  scale_shape_manual(values = c(15, 19),
                     labels = c("Parental", "Hybrid"),
                     name = "")+
  theme_cowplot()+
  xlab("Cross") +
  ylab(bquote('Astaxanthin:Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "")

figuremales.hydroxy <- ggarrange(pbbsmales.hydroxy, pacpmales.hydroxy, ncol=1, nrow=2)

jpeg(file = "Cross Carot Hydroxy Males.jpg", units = "in", width = 8.5, height = 6.2, res = 500)
figuremales.hydroxy
dev.off()

figurefemales.hydroxy <- ggarrange(pbbsfemales.hydroxy, pacpfemales.hydroxy, labels = c("A", "B"), ncol=1, nrow=2)

tiff(file = "Cross Carot Hydroxy Females.tif", units = "in", width = 8.5, height = 7, res = 300)
figurefemales.hydroxy
dev.off()



#Overall hybrid vs parental in carotenoid ratio

#BR, BUF, SD hybrid vs parental
mod.overall.BBS.hydroxy <- lmer(asta.to.hydroxy~gen + (1|ril.id), data = subset(all.dat, sex == "Males"))
summary(mod.overall.BBS.hydroxy)
confint(mod.overall.BBS.hydroxy)

#Run as block
gen <- c("Hybrid","Parental")
mean <- c(3.8345, 9.0796)
lower.cl <- c(1.944928,6.518019)
upper.cl <- c(5.691162,11.641279)
plot.dat.overall.BBS.hydroxy <- data.frame(gen, mean, lower.cl, upper.cl)

jpeg(file = "Overall BR BUF SD MALES hydroxy.jpg", units = "in", width = 5, height = 5, res = 500)
subset(all.dat, sex =="Males") %>% 
  ggplot(aes(x = gen))+
  geom_point(data= plot.dat.overall.BBS.hydroxy, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.overall.BBS.hydroxy, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = asta.to.hydroxy, col = gen), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  ylim(-2.2,22)+
  theme_cowplot()+
  xlab("BR, BUF, SD") +
  ylab(bquote('Astaxanthin:Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 16),
        legend.position = "")+
  scale_color_manual(values=c("red","blue"))+
  annotate("text", 1.5, 20, label = expression(italic(beta)*"=5.24,"~italic(p)*"=0.016"), size = 6)
dev.off()

#AB, CAT, PES hybrid vs parental
mod.overall.ACP.hydroxy <- lmer(asta.to.hydroxy~gen + (1|ril.id), data = subset(all.dat2, sex == "Males"))
summary(mod.overall.ACP.hydroxy)
confint(mod.overall.ACP.hydroxy)

#Run as block
gen <- c("Hybrid","Parental")
mean <- c(9.232, 4.970)
lower.cl <- c(6.003307,1.924931)
upper.cl <- c(12.579581,7.997852)
plot.dat.overall.ACP.hydroxy <- data.frame(gen, mean, lower.cl, upper.cl)

jpeg(file = "Overall AB-CAT-PES hydroxy.jpg", units = "in", width = 5, height = 5, res = 500)
subset(all.dat2, sex =="Males") %>% 
  ggplot(aes(x = gen))+
  geom_point(data= plot.dat.overall.ACP.hydroxy, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.overall.ACP.hydroxy, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = asta.to.hydroxy, col = gen), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  ylim(-2.2,22)+
  theme_cowplot()+
  xlab("AB, CAT, PES") +
  ylab(bquote('Astaxanthin:Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 16),
        legend.position = "")+
  scale_color_manual(values=c("red","blue"))+
  annotate("text", 1.5, 20, label = expression(italic(beta)*"=-4.262,"~italic(p)*"=0.148"), size = 6)
dev.off()


###Compare beta carotene to hydroxyechinenone accumulated
#Subset desired samples from data frame
dat.algae_sub = dat.algae[,c(6,13,15)]
##Then rearrange your data frame
dd.algae = melt(dat.algae_sub, id=c("ril.id"))

mod.beta.vs.hydroxy <- lmer(value~variable + (1|ril.id), data = dd.algae)
summary(mod.beta.vs.hydroxy)
emmeans(mod.beta.vs.hydroxy, pairwise ~ variable)
comp.beta.vs.hydroxy<- emmeans(mod.beta.vs.hydroxy, pairwise ~ variable)
plot.beta.vs.hydroxy <- as.data.frame(comp.beta.vs.hydroxy$emmeans)

jpeg(file = "Beta vs hydroxy dotplot.jpg", units = "in", width = 5, height = 5, res = 500)
dd.algae %>% 
  ggplot(aes(x = variable))+
  geom_point(aes(y = value, colour = variable), size = 4, shape = 19, alpha = 0.5, 
             position = position_dodge2(width = .4))+
  scale_color_manual(values=c("goldenrod3","orangered3"))+
  scale_x_discrete(labels = c('Beta-carotene','Hydroxyechinenone'))+
  geom_point(data= plot.beta.vs.hydroxy, aes(y = emmean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.beta.vs.hydroxy, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  ylim(-0.01,0.10)+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Concentration ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 16),
        legend.position = "")+
  annotate("text", 1.5, -0.01, label = expression(italic(beta)*"=0.040,"~italic(p)*"<0.001"), size = 6)
dev.off()

#Ridgeline plot
jpeg(file = "Beta vs hydroxy ridgeline.jpg", units = "in", width = 6, height = 4, res = 500)
ggplot(dd.algae, aes(x=value, y=variable, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.5, quantile_lines = TRUE, size =1, vline_color = "black", color = 'white')+
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="yellow", high="red")+
  xlab(bquote('Concentration ('~mu*'g'~'mg'^-1*') across all samples')) +ylab("")+
  scale_y_discrete(labels = c('Beta-carotene','Hydroxyechinenone'))+
  theme(legend.position = "non")+
  annotate("text", x = 0.1, y = 2.5, label = expression(italic(n)*"=224,"~italic(p)*"<0.001"), size = 4)

#### Testing for Sex differences by line in BR, SD, and BUF ####
####
#Overall male vs female
mod.sex <- lmer(astax.mg~sex + (1|ril.id), data = dat)
summary(mod.sex)
confint(mod.sex)

sex <- c("Females","Males")
mean <- c(0.20042, 0.26316)
lower.cl <- c(0.1713827,0.2294541)
upper.cl <- c(0.2390727,0.2872513)

plot.dat.sex <- data.frame(sex, mean, lower.cl, upper.cl)

jpeg(file = "sex.jpg", units = "in", width = 5, height = 5, res = 500)
dat %>% 
  ggplot(aes(x = sex))+
  geom_point(data= plot.dat.sex, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.sex, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = astax.mg, col = sex), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  theme_cowplot()+
  xlab("All Lines") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")+
  annotate("text", 1.5, 0.6, label = expression(italic(beta)*"=0.05692,"~italic(p)*"<0.001"), size = 5)
dev.off()

#BR, BUF, SD male vs female
mod.sex.BBS <- lmer(astax.mg~sex + (1|mito.type), data = all.dat)
summary(mod.sex.BBS)
r.squaredGLMM(mod.sex.BBS)
confint(mod.sex.BBS)

sex <- c("Females","Males")
mean <- c(0.19079, 0.23682)
lower.cl <- c(0.1602567,0.2099811)
upper.cl <- c(0.2213306,0.2636680)

plot.dat.sex.BBS <- data.frame(sex, mean, lower.cl, upper.cl)

jpeg(file = "sex BR BUF SD.jpg", units = "in", width = 5, height = 5, res = 500)
all.dat %>% 
  ggplot(aes(x = sex))+
  geom_point(data= plot.dat.sex.BBS, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.sex.BBS, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = astax.mg, col = sex), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  theme_cowplot()+
  xlab("BR, BUF, SD") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")
dev.off()

#AB, CAT, PES hybrid vs parental
mod.sex.ACP <- lmer(astax.mg~sex + (1|mito.type), data = all.dat2)
summary(mod.sex.ACP)
r.squaredGLMM(mod.sex.ACP)
confint(mod.sex.ACP)

sex <- c("Females","Males")
mean <- c(0.29969, 0.34810)
lower.cl <- c(0.2341061,0.3057647)
upper.cl <- c(0.3652813,0.3904379)

plot.dat.sex.ACP <- data.frame(sex, mean, lower.cl, upper.cl)

jpeg(file = "sex AB-CAT-PES.jpg", units = "in", width = 5, height = 5, res = 500)
all.dat2 %>% 
  ggplot(aes(x = sex))+
  geom_point(data= plot.dat.sex.ACP, aes(y = mean),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat.sex.ACP, aes(ymin=lower.cl, ymax = upper.cl), lwd = 1)+
  geom_point(aes(y = astax.mg, col = sex), size = 4, shape = 19, alpha = 0.3, 
             position = position_dodge2(width = .4))+
  theme_cowplot()+
  xlab("AB, CAT, PES") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")
dev.off()



####Comparing carotenoid production to ATP production among generations####
all.males = subset(dat.algae, sex == "Males")

summary <- all.males[,c(6,11,16,17,18)] %>% group_by(ril.id) %>% summarise_each(funs(mean,std.error))

#Take averages and standard errors and paste into "astax-complex-dev-comp-all-norm-males.csv" file
write.csv(summary, file = "carotenoid averages.csv", row.names = F)

#enter C1 ATP, C2 ATP, and CS values to each line and add generation info (in excel). These values come from emmeans summary comparisons across lines below. 
#Renamed astax-complex-comp-all.csv to "astax-complex-dev-comp-all-norm-males.csv"
line.dat <- read.csv("astax-complex-dev-comp-all-norm-males.csv")
#Standard errors in this data sheet come from pairwise comparisons among lines, corrected for multiple comparisons and included a random effect of mito type

#Transform and scale/center data normalized ATP data
line.dat$atpc1.tf <- log(line.dat$atpc1.new)
line.dat$atpc2.tf <- log(line.dat$atpc2.new)
line.dat$atpc1 <- scale(line.dat$atpc1.tf, scale = TRUE, center = TRUE)
line.dat$atpc2 <- scale(line.dat$atpc2.tf, scale = TRUE, center = TRUE)

#Transform and scale/center data non-normalized ATP data
line.dat$atpc1.nonnorm.tf <- log(line.dat$atpc1.nonnorm)
line.dat$atpc2.nonnorm.tf <- log(line.dat$atpc2.nonnorm)
line.dat$atpc1.nonnorm.scale <- scale(line.dat$atpc1.nonnorm.tf, scale = TRUE, center = TRUE)
line.dat$atpc2.nonnorm.scale <- scale(line.dat$atpc2.nonnorm.tf, scale = TRUE, center = TRUE)

#Transform and scale/center data normalized ATP data
line.dat$ratio.tf <- log(line.dat$ratio)
line.dat$ratio.scale <- scale(line.dat$ratio.tf, scale = TRUE, center = TRUE)

#Transform and scale/center data normalized ATP data
line.dat$hydroxy.ratio.tf <- log(line.dat$hydroxy.ratio)
line.dat$hydroxy.ratio.scale <- scale(line.dat$hydroxy.ratio.tf, scale = TRUE, center = TRUE)

# Comparing carotennoid production to atp production 
#ASTAXANTHIN
mod.10 <- lmer(avg ~ atpc1 + gen+ (1|cross), data= line.dat)
summary(mod.10)
r.squaredGLMM(mod.10)
#Removed outlier from row 6, BRSD56 close to cooks distance limit
mod.11 <- lmer(avg ~ atpc2 + gen+ (1|cross), data= line.dat[-c(6),])
summary(mod.11)
r.squaredGLMM(mod.11)

#ASTAXANTHIN : DIETARY CAROTENOIDS
mod.10.ratio <- lmer(ratio.scale ~ atpc1 + gen +(1|cross), data= line.dat)
summary(mod.10.ratio)
r.squaredGLMM(mod.10.ratio)
#Removed outlier from row 6, BRSD56 close to cooks distance limit
mod.11.ratio <- lmer(ratio.scale ~ atpc2 + gen+ (1|cross), data= line.dat[-c(6),])
summary(mod.11.ratio)
r.squaredGLMM(mod.11.ratio)

#ASTAXANTHIN : DIETARY CAROTENOIDS
mod.10.hyd.ratio <- lmer(hydroxy.ratio.scale ~ atpc1 + gen +(1|cross), data= line.dat)
summary(mod.10.hyd.ratio)
r.squaredGLMM(mod.10.hyd.ratio)
#Removed outlier from row 6, BRSD56 close to cooks distance limit
mod.11.hyd.ratio <- lmer(hydroxy.ratio.scale ~ atpc2 + gen+ (1|cross), data= line.dat[-c(6),])
summary(mod.11.hyd.ratio)
r.squaredGLMM(mod.11.hyd.ratio)



#comparing ATP production in Parentals vs Hybrids
mod.12 <- lm(atpc1.new ~gen, data= line.dat)
summary(mod.12)

mod.13 <- lm (atpc2.new ~ gen, data= line.dat)
summary(mod.13)

#plot comparing ATP production in Parentals vs Hybrids
mod.12 <- lm(atpc1.new ~gen-1, data= line.dat)
summary(mod.12)
confint(mod.12)

mod.13 <- lm (atpc2.new ~ gen-1, data= line.dat)
summary(mod.13)
confint(mod.13)

gen <- c("Hybrid","Parental")
average.atp.c1 <- c(565.91, 487.18)
average.atp.c2 <- c(943.9,670.3)
lower.cl.atp.c1 <- c(404.4033,312.7324)
upper.cl.atp.c1 <- c(727.4189,661.6291)
lower.cl.atp.c2 <- c(676.6697,381.6450)
upper.cl.atp.c2 <- c(1211.171,958.972)
plot.dat4 <- data.frame(gen, average.atp.c1, average.atp.c2, lower.cl.atp.c1, upper.cl.atp.c1, lower.cl.atp.c2, upper.cl.atp.c2)

p1 = line.dat %>% 
  ggplot(aes(x = gen)) +
  geom_point(data = plot.dat4, aes(y = average.atp.c1),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat4, aes(ymin=lower.cl.atp.c1, ymax = upper.cl.atp.c1), lwd = 1)+
  geom_point(aes(y = atpc1.new, col = gen),size = 5, alpha =0.7, position = position_dodge2(width = .3))+
  xlab("ATP complex 1") +
  ylab(bquote('ATP log nmol('~'min'^-1*')'))  +
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))+
  theme(legend.position = "none")+
  annotate("text", 1.5, 800, label = expression(italic(beta)*"=-78.73,"~italic(p)*"=0.481"), size = 5)

p2 = line.dat %>% 
  ggplot(aes(x = gen)) +
  geom_point(data = plot.dat4, aes(y = average.atp.c2),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat4, aes(ymin=lower.cl.atp.c2, ymax = upper.cl.atp.c2), lwd = 1)+
  geom_point(aes(y = atpc2.new, color=gen),size = 5, alpha =0.7, position = position_dodge2(width = .3))+
  xlab("ATP complex 2") +
  ylab(bquote('ATP log nmol('~'min'^-1*')'))  +
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))+
  theme(legend.position = "none")+
  annotate("text", 1.5, 1500, label = expression(italic(beta)*"=-273.6,"~italic(p)*"=0.154"), size = 5)

p3 = line.dat %>% 
  ggplot(aes(x = atpc1.new, y=avg)) +
  geom_point(aes(y = avg, col = gen),size = 5, alpha =0.5)+
  geom_smooth(method='lm')+
  geom_errorbarh(xmax= line.dat$atpc1.new + line.dat$atpc1.new.se, xmin = line.dat$atpc1.new - line.dat$atpc1.new.se)+
  geom_errorbar(ymax=line.dat$avg + (line.dat$avg.se), ymin = line.dat$avg - (line.dat$avg.se))+
  xlim(200,900)+
  xlab(bquote('CI ATP Synthesis (nmol'~'min'^-1*')')) +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))+
  guides(color = guide_legend(title="Group"))+
  annotate("text", 750, 0.45, label =  bquote('p=0.018, R'^2*'=0.448'), size=5)

p4 = line.dat[-c(6),] %>% 
  ggplot(aes(x = atpc2.new, y=avg)) +
  geom_point(aes(y = avg,  col = gen),size = 5, alpha =0.5)+
  geom_smooth(method='lm')+
  #Removed outlier from row 6, BRSD56 close to cooks distance limit
  geom_errorbarh(xmax= line.dat[-c(6),]$atpc2.new + line.dat[-c(6),]$atpc2.new.se, xmin = line.dat[-c(6),]$atpc2.new-line.dat[-c(6),]$atpc2.new.se)+
  geom_errorbar(ymax=line.dat[-c(6),]$avg + (line.dat[-c(6),]$avg.se), ymin = line.dat[-c(6),]$avg - (line.dat[-c(6),]$avg.se))+
  xlab(bquote('CII ATP Synthesis (nmol'~'min'^-1*')')) +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))+
  guides(color = guide_legend(title="Group"))+
  annotate("text", 900, 0.45, label = bquote('p=0.412, R'^2*'=0.440'), size = 5)

p3.ratio = line.dat %>% 
  ggplot(aes(x = atpc1.new, y=ratio)) +
  geom_point(aes(y = ratio, col = gen),size = 5, alpha =0.5)+
  geom_smooth(method='lm')+
  geom_errorbarh(xmax= line.dat$atpc1.new + line.dat$atpc1.new.se, xmin = line.dat$atpc1.new - line.dat$atpc1.new.se)+
  geom_errorbar(ymax=line.dat$ratio + (line.dat$ratio.se), ymin = line.dat$ratio - (line.dat$ratio.se))+
  xlim(200,900)+
  xlab(bquote('CI ATP Synthesis (nmol'~'min'^-1*')')) +
  ylab(bquote('Astaxanthin:Dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 10))+
  guides(color = guide_legend(title="Group"))+
  annotate("text", 400, 10, label =  bquote('p=0.178, R'^2*'=0.294'), size=4)

p4.ratio = line.dat[-c(6),] %>% 
  ggplot(aes(x = atpc2.new, y=ratio)) +
  geom_point(aes(y = ratio,  col = gen),size = 5, alpha =0.5)+
  geom_smooth(method='lm')+
  #Removed outlier from row 6, BRSD56 close to cooks distance limit
  geom_errorbarh(xmax= line.dat[-c(6),]$atpc2.new + line.dat[-c(6),]$atpc2.new.se, xmin = line.dat[-c(6),]$atpc2.new-line.dat[-c(6),]$atpc2.new.se)+
  geom_errorbar(ymax=line.dat[-c(6),]$ratio + (line.dat[-c(6),]$ratio.se), ymin = line.dat[-c(6),]$ratio - (line.dat[-c(6),]$ratio.se))+
  xlab(bquote('CII ATP Synthesis (nmol'~'min'^-1*')')) +
  ylab(bquote('Astaxanthin:Dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 10))+
  guides(color = guide_legend(title="Group"))+
  annotate("text", 550, 10, label = bquote('p=0.778, R'^2*'=0.117'), size = 4)

p3.hydroxy.ratio = line.dat %>% 
  ggplot(aes(x = atpc1.new, y=hydroxy.ratio)) +
  geom_point(aes(y = hydroxy.ratio, col = gen),size = 5, alpha =0.5)+
  geom_smooth(method='lm')+
  geom_errorbarh(xmax= line.dat$atpc1.new + line.dat$atpc1.new.se, xmin = line.dat$atpc1.new - line.dat$atpc1.new.se)+
  geom_errorbar(ymax=line.dat$hydroxy.ratio + (line.dat$hydroxy.ratio.se), ymin = line.dat$hydroxy.ratio - (line.dat$hydroxy.ratio.se))+
  xlim(200,900)+
  xlab(bquote('CI ATP Synthesis (nmol'~'min'^-1*')')) +
  ylab(bquote('Astaxanthin:Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 10))+
  guides(color = guide_legend(title="Group"))+
  annotate("text", 450, 0, label =  bquote('p=0.075, R'^2*'=0.340'), size=4)

p4.hydroxy.ratio = line.dat[-c(6),] %>% 
  ggplot(aes(x = atpc2.new, y=hydroxy.ratio)) +
  geom_point(aes(y = hydroxy.ratio,  col = gen),size = 5, alpha =0.5)+
  geom_smooth(method='lm')+
  #Removed outlier from row 6, BRSD56 close to cooks distance limit
  geom_errorbarh(xmax= line.dat[-c(6),]$atpc2.new + line.dat[-c(6),]$atpc2.new.se, xmin = line.dat[-c(6),]$atpc2.new-line.dat[-c(6),]$atpc2.new.se)+
  geom_errorbar(ymax=line.dat[-c(6),]$hydroxy.ratio + (line.dat[-c(6),]$hydroxy.ratio.se), ymin = line.dat[-c(6),]$hydroxy.ratio - (line.dat[-c(6),]$hydroxy.ratio.se))+
  xlab(bquote('CII ATP Synthesis (nmol'~'min'^-1*')')) +
  ylab(bquote('Astaxanthin:Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 10))+
  guides(color = guide_legend(title="Group"))+
  annotate("text", 700, 0, label = bquote('p=0.165, R'^2*'=0.233'), size = 4)

#Combine plots and export
figure.hyb.parental <- ggarrange(p1, p2, labels = c("A", "B"), ncol=1, nrow=2, common.legend = TRUE, legend = "right")

jpeg(file="ATP_hybrid_parental.jpg", units="in", width=6, height=10, res=500)
figure.hyb.parental
dev.off()

####Compare ATP by line####
atp.dat <- read.csv(file="ATP.dat.norm.csv")
#Remove top line with messed up ATP assay for AB
atp.dat <- atp.dat[-c(1),]

#Set up group 1 (BR, BUF, and SD)
atp.dat %>% 
  filter(cross %in% c( "BUFSD", "BRSD", "BR",
                       "BUF",
                       "SD"))-> group.1.atp

group.1.atp$ril.id<- factor(group.1.atp$ril.id,
                        levels = c("BUF", "BUFSD19", "BUFSD24","BUFSD4", "SD", "BRSD50", "BRSD56", "BR"))

#Set up group 2 (AB, CAT, and PES)
atp.dat %>% 
  filter(cross %in% c( "CAT", "CATAB", "ABCAT",
                       "AB", "PESAB",
                       "PES"))-> group.2.atp

group.2.atp$ril.id<- factor(group.2.atp$ril.id, 
                        levels = c("CAT", "CATAB27", "AB ","PESAB20", "PES"))

#Complex 1 BR, BUF, and SD
mod.atpc1.bbs <-lm(atpc1 ~ ril.id , data = group.1.atp)
### Making multiple comparisons and correcting for type 1 error
comp.atp.c1.bbs <- emmeans(mod.atpc1.bbs, pairwise ~ ril.id)
atp.c1.bbs <- as.data.frame(comp.atp.c1.bbs$emmeans)
atp.c1.bbs$line = c("Parental", "Hybrid","Hybrid", "Hybrid", "Parental", "Hybrid", 
                    "Hybrid", "Parental")

#Complex 1 AB, CAT and PES
mod.atpc1.acp <-lm(atpc1 ~ ril.id , data = group.2.atp)
### Making multiple comparisons and correcting for type 1 error
comp.atp.c1.acp <- emmeans(mod.atpc1.acp, pairwise ~ ril.id)
atp.c1.acp <- as.data.frame(comp.atp.c1.acp$emmeans)
atp.c1.acp$line = c("Parental", "Hybrid","Parental", "Hybrid", "Parental")

#Complex2 BR, BUF, SD
mod.atpc2.bbs <-lm(atpc2 ~ ril.id , data = group.1.atp)
### Making multiple comparisons and correcting for type 1 error
comp.atp.c2.bbs <- emmeans(mod.atpc2.bbs, pairwise ~ ril.id)
atp.c2.bbs <- as.data.frame(comp.atp.c2.bbs$emmeans)
atp.c2.bbs$line = c("Parental", "Hybrid","Hybrid", "Hybrid", "Parental", "Hybrid", 
                    "Hybrid", "Parental")

#Complex 2 AB, CAT and PES
mod.atpc2.acp <-lm(atpc2 ~ ril.id , data = group.2.atp)
### Making multiple comparisons and correcting for type 1 error
comp.atp.c2.acp <- emmeans(mod.atpc2.acp, pairwise ~ ril.id)
atp.c2.acp <- as.data.frame(comp.atp.c2.acp$emmeans)
atp.c2.acp$line = c("Parental", "Hybrid","Parental", "Hybrid", "Parental")


#Make plot for comparing complex 1 in BBS
ATPC1BBS_groups <- group.1.atp %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = atp.c1.bbs, aes(y = emmean),size = 7, alpha = 0.7)+
  geom_linerange(data= atp.c1.bbs, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = atpc1, col = cross), size = 4, shape = 19, alpha = 0.5, 
             position = position_dodge2(width = .2))+
  theme_cowplot()+
  xlab("Complex 1") +
  ylab(bquote('ATP nmol/min'))  +
  scale_color_manual(values = c("slateblue1","goldenrod3", "chartreuse4","orangered1","lightsteelblue4" ))+
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.position = "")+
  guides(color = guide_legend(title="Group:"))


#Make plot for comparing complex 1 in BBS
ATPC1ACP_groups <- group.2.atp %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = atp.c1.acp, aes(y = emmean),size = 7, alpha = 0.8)+
  geom_linerange(data= atp.c1.acp, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = atpc1, col = cross), size = 4, shape = 19, alpha = 0.5, 
             position = position_dodge2(width = .4))+
  theme_cowplot()+
  xlab("Complex 1") +
  ylab(bquote('ATP nmol/min'))  +
  scale_color_manual(values = c("lightsteelblue4","chartreuse4", "orangered1","slateblue1","goldenrod3" ))+
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.position = "")+
  guides(color = guide_legend(title="Group:"))


#Make plot for comparing complex 2 in BBS
ATPC2BBS_groups <- group.1.atp %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = atp.c2.bbs, aes(y = emmean),size = 7, alpha = 0.8)+
  geom_linerange(data= atp.c2.bbs, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = atpc2, col = cross), size = 4, shape = 19, alpha = 0.4, 
             position = position_dodge2(width = .4))+
  theme_cowplot()+
  xlab("Complex 2") +
  ylab(bquote('ATP nmol/min'))  +
  scale_color_manual(values = c("slateblue1","goldenrod3", "chartreuse4","orangered1","lightsteelblue4" ))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")+
guides(color = guide_legend(title="Group:"))

#Make plot for comparing complex 2 in ACP
ATPC2ACP_groups<- group.2.atp %>% 
  ggplot(aes(x = ril.id))+
  geom_point(data = atp.c2.acp, aes(y = emmean),size = 7, alpha = 0.8)+
  geom_linerange(data= atp.c2.acp, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = atpc2, col = cross), size = 4, shape = 19, alpha = 0.4, 
             position = position_dodge2(width = .4))+
  theme_cowplot()+
  xlab("Complex 2") +
  ylab(bquote('ATP nmol/min'))  +
  scale_color_manual(values = c("lightsteelblue4","chartreuse4", "orangered1","slateblue1","goldenrod3" ))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")+
  guides(color = guide_legend(title="Group:"))


#### Compare ATP to dev time ####
#log transform dev time data and then center and scale
line.dat$devtime.tf = log(line.dat$devtime)
line.dat$devtimescale = scale(line.dat$devtime.tf, scale = TRUE, center = TRUE)

mod.atp1.dev <- lmer(devtimescale ~ atpc1 + gen + (1|cross), data= line.dat)
summary(mod.atp1.dev)
r.squaredGLMM(mod.atp1.dev)

#Removed outlier from row 6, BRSD56 close to cooks distance limit
mod.atp2.dev <- lmer(devtimescale ~ atpc2 + gen + (1|cross), data= line.dat[-c(6),])
summary(mod.atp2.dev)
r.squaredGLMM(mod.atp2.dev)

p5 = line.dat %>% 
  ggplot(aes(x = atpc1.new, y=devtime)) +
  geom_point(aes(y = devtime, col = gen),size = 5, alpha =0.5)+
  geom_smooth(method='lm')+
  geom_errorbarh(xmax= line.dat$atpc1.new + line.dat$atpc1.new.se, xmin = line.dat$atpc1.new - line.dat$atpc1.new.se)+
  geom_errorbar(ymax=line.dat$devtime + (line.dat$devtime.se), ymin = line.dat$devtime - (line.dat$devtime.se))+
  xlab(bquote('CI ATP Synthesis (nmol'~'min'^-1*')')) +
  ylab("Development Time (days)")  +
  ylim(6,9)+
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))+
  guides(color = guide_legend(title="Group"))+
  annotate("text", 500, 6.1, label = bquote('p=0.118, R'^2*'=0.323'), size=5)

p6 = line.dat[-c(6),] %>% 
  ggplot(aes(x = atpc2.new, y=devtime)) +
  geom_point(aes(y = devtime,  col = gen),size = 5, alpha =0.5)+
  geom_smooth(method='lm')+
  #Removed outlier from row 6, BRSD56 close to cooks distance limit
  geom_errorbarh(xmax= line.dat[-c(6),]$atpc2.new + line.dat[-c(6),]$atpc2.new.se, xmin = line.dat[-c(6),]$atpc2.new - line.dat[-c(6),]$atpc2.new.se)+
  geom_errorbar(ymax=line.dat[-c(6),]$devtime + (line.dat[-c(6),]$devtime.se), ymin = line.dat[-c(6),]$devtime - (line.dat[-c(6),]$devtime.se))+
  xlab(bquote('CII ATP Synthesis (nmol'~'min'^-1*')')) +
  ylab("Development Time (days)")  +
  ylim(5.7,9)+
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))+
  guides(color = guide_legend(title="Group"))+
  annotate("text", 700, 6.1, label = bquote('p=0.445, R'^2*'=0.206'), size = 5)


####Compare astax content by dev time and gen ####
# Comparing carotenoid production to dev time by generation
#ASTAXANTHIN
mod.dev.astax <- lmer(avg ~ devtimescale + gen + (1|cross), data= line.dat)
summary(mod.dev.astax)
r.squaredGLMM(mod.dev.astax)

#CAROTENOID RATIO
mod.dev.ratio <- lmer(ratio.scale ~ devtimescale + gen + (1|cross), data= line.dat)
summary(mod.dev.ratio)
r.squaredGLMM(mod.dev.ratio)

#CAROTENOID RATIO JUST HYDROXY
mod.dev.hydroxy.ratio <- lmer(hydroxy.ratio.scale ~ devtimescale + gen + (1|cross), data= line.dat)
summary(mod.dev.hydroxy.ratio)
r.squaredGLMM(mod.dev.hydroxy.ratio)

#comparing dev time in Parentals vs Hybrids
mod.dev <- lm(devtime ~gen, data= line.dat)
summary(mod.dev)

#plot comparing dev time in Parentals vs Hybrids
mod.dev <- lm(devtime ~gen-1, data= line.dat)
summary(mod.dev)
confint(mod.dev)

gen <- c("Hybrid","Parental")
average.dev <- c(7.6367, 7.2075)
lower.cl.dev <- c(6.988690,6.413894)
upper.cl.dev <- c(8.284644,8.001106)
plot.dat4.2 <- data.frame(gen, average.dev, lower.cl.dev, upper.cl.dev)

#output to file just to look at later if needed
jpeg(file = "Cross by dev time.jpg", units = "in", width = 5.5, height = 5.5, res = 500)
line.dat %>% 
  ggplot(aes(x = gen)) +
  geom_point(data = plot.dat4.2, aes(y = average.dev),size = 7, alpha = 0.7)+
  geom_linerange(data= plot.dat4.2, aes(ymin=lower.cl.dev, ymax = upper.cl.dev), lwd = 1)+
  geom_point(aes(y = devtime),size = 4, alpha =0.3, position = position_dodge2(width = .3))+
  xlab("") +
  ylab(bquote('Dev. time (days)'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))+
  theme_cowplot()+
  annotate("text", 1.5, 10, label = expression(italic(beta)*"=-0.429,"~italic(p)*"=0.362"), size = 5)
dev.off()

#Line graph of astaxanthin vs devtime
astax.dev <- line.dat %>% 
  ggplot(aes(x = devtime, y = avg)) +
  geom_point(aes(y = avg, col = gen),size = 5, alpha =0.6)+
  geom_smooth(method='lm')+
  geom_errorbarh(xmax= line.dat$devtime + line.dat$devtime.se, xmin = line.dat$devtime-line.dat$devtime.se)+
  geom_errorbar(ymax=line.dat$avg + (line.dat$avg.se), ymin = line.dat$avg - (line.dat$avg.se))+
  xlim(6.9,9.0)+
  theme_cowplot()+
  xlab("Dev. Time (Days)") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))+
  guides(color = guide_legend(title="Group"))+
  
  theme(legend.position = "top", legend.direction = "horizontal")+
  annotate("text", 7.5, 0.6, label = bquote('p=0.178, R'^2*'=0.728'), size = 5)

#Line graph of astaxanthin : dietary carotenoids ratio vs devtime
ratio.dev <- line.dat %>% 
  ggplot(aes(x = devtime, y = ratio)) +
  geom_point(aes(y = ratio, col = gen),size = 5, alpha =0.6)+
  geom_smooth(method='lm')+
  geom_errorbarh(xmax= line.dat$devtime + line.dat$devtime.se, xmin = line.dat$devtime-line.dat$devtime.se)+
  geom_errorbar(ymax=line.dat$ratio + (line.dat$ratio.se), ymin = line.dat$ratio - (line.dat$ratio.se))+
  theme_cowplot()+
  xlab("Dev. Time (Days)") +
  ylab(bquote('Astaxanthin:Dietary carotenoids ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 10))+
  guides(color = guide_legend(title="Group"))+
  theme(legend.position = "top", legend.direction = "horizontal")+
  annotate("text", 9, 14, label = bquote('p=0.006, R'^2*'=0.806'), size = 4)

#Line graph of astaxanthin : hydroxy ratio vs devtime
hydroxy.ratio.dev <- line.dat %>% 
  ggplot(aes(x = devtime, y = hydroxy.ratio)) +
  geom_point(aes(y = hydroxy.ratio, col = gen),size = 5, alpha =0.6)+
  geom_smooth(method='lm')+
  geom_errorbarh(xmax= line.dat$devtime + line.dat$devtime.se, xmin = line.dat$devtime-line.dat$devtime.se)+
  geom_errorbar(ymax=line.dat$hydroxy.ratio + (line.dat$hydroxy.ratio.se), ymin = line.dat$hydroxy.ratio - (line.dat$hydroxy.ratio.se))+
  theme_cowplot()+
  xlab("Dev. Time (Days)") +
  ylab(bquote('Astaxanthin:hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 10))+
  guides(color = guide_legend(title="Group"))+
  theme(legend.position = "top", legend.direction = "horizontal")+
  annotate("text", 9, 15, label = bquote('p=0.023, R'^2*'=0.510'), size = 4)


####Compare devtime across all lines and make figure ####
#Read in data
devtime.dat = read.csv(file = "devtime.dat.csv")

#Set up group 1 (BR, BUF, and SD)
devtime.dat %>% 
  filter(line %in% c( "BR", "BRSD50", "BRSD56", "BUFSD24", "BUFSD4",
                       "BUF",
                       "SD"))-> group.1.devtime

group.1.devtime$line<- factor(group.1.devtime$line,
                            levels = c("BUF", "BUFSD24","BUFSD4", "SD", "BRSD50", "BRSD56", "BR"))

#Set up group 2 (AB, CAT, and PES)
devtime.dat %>% 
  filter(line %in% c("CATAB27", "ABCAT11", "PESAB20",
                       "PES"))-> group.2.devtime

group.2.devtime$line<- factor(group.2.devtime$line, 
                            levels = c("ABCAT11", "CATAB27","PESAB20", "PES"))

#Dev time BR, BUF, and SD
mod.devtime.bbs <-lm(dph ~ line , data = group.1.devtime)
### Making multiple comparisons and correcting for type 1 error
comp.devtime.bbs <- emmeans(mod.devtime.bbs, pairwise ~ line)
devtime.bbs <- as.data.frame(comp.devtime.bbs$emmeans)
devtime.bbs$group = c("Parental", "Hybrid","Hybrid", "Parental", "Hybrid", 
                    "Hybrid", "Parental")

#Devtime AB, CAT and PES
mod.devtime.acp <-lm(dph ~ line , data = group.2.devtime)
### Making multiple comparisons and correcting for type 1 error
comp.devtime.acp <- emmeans(mod.devtime.acp, pairwise ~ line)
devtime.acp <- as.data.frame(comp.devtime.acp$emmeans)
devtime.acp$group = c("Hybrid","Hybrid", "Hybrid", "Parental")


#Make comparison plot for BBS
devtime_across_lines_BBS <- group.1.devtime %>% 
  ggplot(aes(x = line))+
  geom_point(data = devtime.bbs, aes(y = emmean),size = 7, alpha = 0.8)+
  geom_linerange(data= devtime.bbs, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = dph, col = cross), size = 4, shape = 19, alpha = 0.5, 
             position = position_dodge2(width = .2))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Development time (days)'))  +
  scale_color_manual(values = c("slateblue1","goldenrod3", "chartreuse4","orangered1","lightsteelblue4" ))+
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.position = "")+
  guides(color = guide_legend(title="Group:"))


#Make comparison plot for ACP
devtime_across_lines_ACP <- group.2.devtime %>% 
  ggplot(aes(x = line))+
  geom_point(data = devtime.acp, aes(y = emmean),size = 7, alpha = 0.8)+
  geom_linerange(data= devtime.acp, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = dph, col = cross), size = 4, shape = 19, alpha = 0.5, 
             position = position_dodge2(width = .2))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Development time (days)'))  +
  scale_color_manual(values = c("purple1","orangered1", "slateblue1","goldenrod3" ))+
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.position = "")+
  guides(color = guide_legend(title="Group:"))


####Look at relationship between CS and ATP####
### Comparing CS to ATP
#Transform and scale
line.dat$CS.tf = log(line.dat$CS)
line.dat$CS.scale = scale(line.dat$CS.tf, scale = TRUE, center = TRUE)


#Comparing CS to ATP
#Complex 1
mod.atp1.cs <- lmer(atpc1.nonnorm.scale~CS.scale + gen + (1|cross), data = line.dat)
summary(mod.atp1.cs)
r.squaredGLMM(mod.atp1.cs)

p7 <- line.dat %>%
  ggplot(aes(x = atpc1.nonnorm, y=CS)) +
  geom_point(aes(y = CS, col = gen),size = 5, alpha =0.5)+
  geom_smooth(method='glm')+
  geom_errorbarh(xmax= line.dat$atpc1.nonnorm + line.dat$atpc1.nonnorm.se, xmin = line.dat$atpc1.nonnorm-line.dat$atpc1.nonnorm.se)+
  geom_errorbar(ymax=line.dat$CS + (line.dat$CS.se), ymin = line.dat$CS - (line.dat$CS.se))+
  ylab(bquote('Citrate (nmol'~'min'^-1*')')) +
  xlab(bquote('Uncorrected CI ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size =18))+
  guides(color = guide_legend(title="Group"))+
  annotate("text", 109, 0.16, label = bquote('p=0.963, R'^2*'=0.298'), size=5)

#ATP complex 2
mod.atp2.cs <- lmer(atpc2.nonnorm.scale~CS.scale + gen+ (1|cross), data = line.dat)
summary(mod.atp2.cs)
r.squaredGLMM(mod.atp2.cs)


p8<- line.dat %>%
  ggplot(aes(x = atpc2.nonnorm, y=CS)) +
  geom_point(aes(y = CS, col = gen),size = 5, alpha =0.5)+
  geom_smooth(method='glm')+
  #Removed outlier from row 6, BRSD56 close to cooks distance limit
  geom_errorbarh(xmax= line.dat$atpc2.nonnorm + line.dat$atpc2.nonnorm.se, xmin = line.dat$atpc2.nonnorm-line.dat$atpc2.nonnorm.se)+
  geom_errorbar(ymax=line.dat$CS + (line.dat$CS.se), ymin = line.dat$CS - (line.dat$CS.se))+
  ylab(bquote('Citrate (nmol'~'min'^-1*')')) +
  xlab(bquote('Uncorrected CII ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size =18))+
  guides(color = guide_legend(title="Group"))+
  annotate("text", 300, 0.14, label = bquote('p=0.979, R'^2*'=0.456'), size=5)


#### SD x SCN cross ####

#Read in astaxanthin data
dat.sdscn <- read.csv("astax.hyb.alt.SDxSCN.csv")
dat.sdscn %>% 
  mutate(sex = fct_recode(sex, "Males" = "M")) -> dat.sdscn

#Read in devtime+astaxanthin experiment data 
devtime.dat.SDSCN = read.csv(file = "devtime.astax.SDxSCN.csv")

#Read in summary data across lines/generations
dtime.astax.atp.sdscn <- read.csv("devtime.hyb.astax.atp.sdscn.csv")
#Subset data frame
dtime.astax.atp.sdscn.subset = dtime.astax.atp.sdscn[c(3,6:8),]


#Transform and scale/center devtime data
dtime.astax.atp.sdscn$dph.tf <- log(dtime.astax.atp.sdscn$dph)
dtime.astax.atp.sdscn$dph.tfs <- scale(dtime.astax.atp.sdscn$dph.tf, scale = TRUE, center = TRUE)
#Transform and scale/center CS data
dtime.astax.atp.sdscn$cs.tf <- log(dtime.astax.atp.sdscn$cs)
dtime.astax.atp.sdscn$cs.tfs <- scale(dtime.astax.atp.sdscn$cs.tf, scale = TRUE, center = TRUE)

#Transform and scale/center ATP data
dtime.astax.atp.sdscn$atpc1.tf <- log(dtime.astax.atp.sdscn$c1)
dtime.astax.atp.sdscn$atpc1 <- scale(dtime.astax.atp.sdscn$atpc1.tf, scale = TRUE, center = TRUE)
dtime.astax.atp.sdscn$atpc2.tf <- log(dtime.astax.atp.sdscn$c2)
dtime.astax.atp.sdscn$atpc2 <- scale(dtime.astax.atp.sdscn$atpc2.tf, scale = TRUE, center = TRUE)

dtime.astax.atp.sdscn$atpc1.nonnorm.tf <- log(dtime.astax.atp.sdscn$c1.nonnorm)
dtime.astax.atp.sdscn$atpc1.nonnorm <- scale(dtime.astax.atp.sdscn$atpc1.nonnorm.tf, scale = TRUE, center = TRUE)
dtime.astax.atp.sdscn$atpc2.nonnorm.tf <- log(dtime.astax.atp.sdscn$c2.nonnorm)
dtime.astax.atp.sdscn$atpc2.nonnorm <- scale(dtime.astax.atp.sdscn$atpc2.nonnorm.tf, scale = TRUE, center = TRUE)



####Just astaxanthin comparison across lines

atp <- read.csv("atp.sdscn.csv")
atp %>% 
  mutate(Population = fct_recode(Population, "SDxSCN F3" = "W F3")) -> atp
atp %>% 
  mutate(Population = fct_recode(Population, "SCNxSD F3" = "X F3")) -> atp

#Set up data for plotting hybrids vs. parentals
#make new dat with only SD, SCN, and hybrids
dat.sdscn %>% 
  filter(cross %in% c( "SDxSD",
                       "SDxSCN",
                       "SCNxSCN",
                       "SCNxSD"))-> group.1.sdscn

#put lines in order
group.1.sdscn$hyb.id<- factor(group.1.sdscn$hyb.id,
                        levels = c("SD", "SCN", "SDxSCN F1","SDxSCN F2", "SDxSCN F3", "SCNxSD F1", 
                                   "SCNxSD F2", "SCNxSD F3"))
#put crosses in order
group.1.sdscn$cross<- factor(group.1.sdscn$cross,
                       levels = c("SCNxSCN", "SDxSD", "SCNxSD", "SDxSCN"))

#put generations in order
group.1.sdscn$gen<- factor(group.1.sdscn$gen,
                             levels = c("Parental", "F1", "F2", "F3"))

#make linear model of asta by line
mod.1.sdscn <-lm(astax.mg ~ hyb.id, data = group.1.sdscn)
emmeans(mod.1.sdscn, pairwise~hyb.id)

#take means of model and turn into csv
comp.g1.sdscn <- emmeans(mod.1.sdscn, pairwise ~ hyb.id)
plot.dat.sdscn <- as.data.frame(comp.g1.sdscn$emmeans) 

#make plot of asta by LINE
tiff(file = "Fig S7.tif", units = "in", width = 6.5, height = 4.5, res = 300)
group.1.sdscn %>% 
  ggplot(aes(x = hyb.id))+
  geom_point(data= plot.dat.sdscn, aes(y = emmean),size = 7, alpha = 0.7, shape=c(15, 15, 19, 19, 19, 19, 19, 19))+
  geom_point(aes(y = astax.mg, col = hyb.id), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  geom_linerange(data= plot.dat.sdscn, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_rect(mapping=aes(xmin=.8, xmax=4.5, ymin=0.108, ymax=0.161), fill="mediumblue", alpha= 0.002, linetype = "blank") +
  geom_rect(mapping=aes(xmin=4.9, xmax=8.5, ymin=0.122, ymax=0.199), fill="mediumblue", alpha= 0.002, linetype = "blank") +
  scale_color_manual(values = c("mediumblue", "mediumblue", "firebrick3","firebrick3","firebrick3", "firebrick3","firebrick3", "firebrick3"))+
  scale_shape_manual(values = c(15))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=45, hjust =1))+
  xlab("") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")
dev.off()

#make linear model of asta by GENERATION while controlling for mitotype
mod.1.sdscn.gen <-lmer(astax.mg ~ gen + (1|mito.type), data = group.1.sdscn)
emmeans(mod.1.sdscn.gen, pairwise~gen)
confint(emmeans(mod.1.sdscn.gen, pairwise~gen))
#take means of model and turn into csv
comp.g1.sdscn.gen <- emmeans(mod.1.sdscn.gen, pairwise ~ gen)
plot.dat.sdscn.gen <- as.data.frame(comp.g1.sdscn.gen$emmeans)

#put generations in order in plot data
plot.dat.sdscn.gen$gen<- factor(plot.dat.sdscn.gen$gen,
                           levels = c("Parental", "F1", "F2", "F3"))

#make plot of asta by GENERATION
tiff(file = "Fig 4.tif", units = "in", width = 8.5, height = 5.5, res = 300)
group.1.sdscn %>% 
  ggplot(aes(x = gen))+
  geom_point(data= plot.dat.sdscn.gen, aes(y = emmean),size = 7, alpha = 0.7, shape=c(15, 19, 19, 19))+
  geom_point(aes(y = astax.mg, col = gen), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  geom_linerange(data= plot.dat.sdscn.gen, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_rect(mapping=aes(xmin=.8, xmax=4.5, ymin=0.10842, ymax=0.177), fill="mediumblue", alpha= 0.002, linetype = "blank") +
  scale_color_manual(values = c("mediumblue", "firebrick3","firebrick3","firebrick3"))+
  scale_shape_manual(values = c(15))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "")

dev.off()

#### Devtime analysis

#Compare across all lines and make figure

#Run model for lines
mod.devtime.SDSCN <-lm(dph ~ line, data = devtime.dat.SDSCN)
### Making multiple comparisons and correcting for type 1 error
comp.devtime.SDSCN <- emmeans(mod.devtime.SDSCN, pairwise ~ line)
devtime.comparison.SDSCN <- as.data.frame(comp.devtime.SDSCN$emmeans)
devtime.comparison.SDSCN$group = c("Parental", "Hybrid", "Hybrid","Hybrid", "Parental", "Hybrid", "Hybrid", 
                                   "Hybrid")

#Make comparison plot for all lines
devtime_across_lines_SDSCN <-devtime.dat.SDSCN %>% 
  ggplot(aes(x = line))+
  geom_point(data = devtime.comparison.SDSCN, aes(y = emmean, col = group),size = 7, alpha = 0.8)+
  geom_rect(mapping=aes(xmin=1, xmax=4.1, ymin=4.89, ymax=7.89), fill="mediumblue", alpha= 0.002, linetype = "blank") +
  geom_rect(mapping=aes(xmin=4.98, xmax=8.1, ymin=4.48, ymax=7.83), fill="mediumblue", alpha= 0.002, linetype = "blank") +
  geom_linerange(data= devtime.comparison.SDSCN, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = dph), size = 4, shape = 19, alpha = 0.5, 
             position = position_dodge2(width = .2))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Development time (days)'))  +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.position = "top")+
  guides(color = guide_legend(title="Group"))

#Run model for GENERATIONS
mod.devtime.SDSCN.gen <-lm(dph ~ gen, data = devtime.dat.SDSCN)
### Making multiple comparisons and correcting for type 1 error
comp.devtime.SDSCN.gen <- emmeans(mod.devtime.SDSCN.gen, pairwise ~ gen)
confint(emmeans(mod.devtime.SDSCN.gen, pairwise ~ gen))
devtime.comparison.SDSCN.gen <- as.data.frame(comp.devtime.SDSCN.gen$emmeans)
devtime.comparison.SDSCN.gen$group = c("Parental", "Hybrid","Hybrid", "Hybrid")

#Make comparison plot for GENERATIONS
#put generations in order in plot data
devtime.dat.SDSCN$gen<- factor(devtime.dat.SDSCN$gen,
                                levels = c("Parental", "F1", "F2", "F3"))

devtime.comparison.SDSCN.gen$gen<- factor(devtime.comparison.SDSCN.gen$gen,
                               levels = c("Parental", "F1", "F2", "F3"))


devtime_across_gen_SDSCN <- devtime.dat.SDSCN %>% 
  ggplot(aes(x = gen))+
  geom_point(data = devtime.comparison.SDSCN.gen, aes(y = emmean),size = 7, alpha = 0.8, shape = c(15, 19, 19, 19))+
  #geom_rect(mapping=aes(xmin=1, xmax=4.1, ymin=5.15, ymax=7.40), fill="mediumblue", alpha= 0.002, linetype = "blank") +
  geom_linerange(data= devtime.comparison.SDSCN.gen, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  geom_point(aes(y = dph, col = gen), size = 4, shape = 19, alpha = 0.5, 
             position = position_dodge2(width = .2))+
  scale_color_manual(values = c("firebrick3","firebrick3", "firebrick3", "mediumblue"))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Development time (days)'))  +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.position = "top")+
  guides(color = guide_legend(title="Group"))


##Line graphs for averages of each line and generation through F3

#Devtime vs astaxanthin
SDSCNdevtimeast <- ggplot(dtime.astax.atp.sdscn, aes(x=dph, y=asta))+
  geom_point(aes(y=asta, col= gen), size =6)+
  geom_smooth(method="lm", level=.95)+
  geom_errorbarh(xmax= dtime.astax.atp.sdscn$dph + dtime.astax.atp.sdscn$dph.sd, xmin = dtime.astax.atp.sdscn$dph-dtime.astax.atp.sdscn$dph.sd)+
  geom_errorbar(ymax=dtime.astax.atp.sdscn$asta + dtime.astax.atp.sdscn$asta.se, ymin = dtime.astax.atp.sdscn$asta - dtime.astax.atp.sdscn$asta.se)+
  xlab("Development Time (days)")+
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))+
  theme_cowplot()+
  theme(legend.position = "top")+
  guides(color = guide_legend(title="Generation:"))+
  xlim(5,10)+
  annotate("text", 8.5, 0.20, label = bquote('p=0.016, R'^2*'=0.59'), size = 4)



#Astaxanthin vs devtime across all generations
mod.dtimeavg.hyb <- lm(asta~dph.tfs, data=dtime.astax.atp.sdscn)
summary(mod.dtimeavg.hyb)


#Test for relationship between astax and development time subsetting each generation
##Line graph for all replicates within each line
#Parental
mod.devtime.asta.SDSCN.parental = lm(astax~dph, data = subset(devtime.dat.SDSCN, gen == "Parental"))
summary(mod.devtime.asta.SDSCN.parental)

#F1
mod.devtime.asta.SDSCN.f1 = lm(astax~dph, data = subset(devtime.dat.SDSCN, gen == "F1"))
summary(mod.devtime.asta.SDSCN.f1)

#F2
mod.devtime.asta.SDSCN.f2 = lm(astax~dph, data = subset(devtime.dat.SDSCN, gen == "F2"))
summary(mod.devtime.asta.SDSCN.f2)

#F3
mod.devtime.asta.SDSCN.f3 = lm(astax~dph, data = subset(devtime.dat.SDSCN, gen == "F3"))
summary(mod.devtime.asta.SDSCN.f3)

jpeg(file = "ALL devtime vs asta SD x SCN.jpg", units = "in", width = 6.5, height = 5.5, res = 500)
ggplot(devtime.dat.SDSCN, aes(x=dph, y=astax, col = gen))+
  geom_point(aes(y=astax))+
  geom_smooth(method="lm")+
  geom_errorbarh(xmax= devtime.dat.SDSCN$dph + devtime.dat.SDSCN$dph.sd, xmin = devtime.dat.SDSCN$dph-devtime.dat.SDSCN$dph.sd)+
  scale_color_manual(name = "", values = c("blue", "goldenrod", "black", "red"))+
  xlab("Development time (days)")+
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))+
  theme_cowplot()+
  annotate("text", 12, 0.31, label = bquote('p=0.346, R'^2*'=0.002'), size = 4, col = "blue")+
  annotate("text", 12, 0.29, label = bquote('p=0.025, R'^2*'=0.26'), size = 4, col = "goldenrod")+
  annotate("text", 12, 0.27, label = bquote('p=0.909, R'^2*'=-0.12'), size = 4, col = "black")+
  annotate("text", 12, 0.25, label = bquote('p=0.345, R'^2*'=-0.002'), size = 4, col = "red")+
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.position = "top")+
  guides(color = guide_legend(title="Group"))
dev.off()


####ATP analysis

### Hybrid ATP
atp %>% 
  filter(Food == "Tetraselmis",
         Population %in% c( "SD",                            
                            "SCN",
                            "SDxSCN F3",
                            "SCNxSD F3")) -> hyb.atp
hyb.atp$Population <- factor(hyb.atp$Population, levels = c("SD", "SCN", "SDxSCN F3", "SCNxSD F3"))

#Complex 1 emmeans
mod.hyb.atptetra <- lm(c1~Population, data=hyb.atp) 
comp.hyb.atptetra <- emmeans(mod.hyb.atptetra, pairwise~Population)
plot.dat.hyb.atptetra <- as.data.frame(comp.hyb.atptetra$emmeans)
plot.dat.hyb.atptetra$Population <- factor(plot.dat.hyb.atptetra$Population, levels = c("SD", "SCN", "SDxSCN F3", "SCNxSD F3"))

#Complex 2 emmeans
mod.hyb.atptetra2 <- lm(c2~Population, data=hyb.atp) 
comp.hyb.atptetra2 <- emmeans(mod.hyb.atptetra2, pairwise~Population)
plot.dat.hyb.atptetra2 <- as.data.frame(comp.hyb.atptetra2$emmeans)
plot.dat.hyb.atptetra2$Population <- factor(plot.dat.hyb.atptetra$Population, levels = c("SD", "SCN", "SDxSCN F3", "SCNxSD F3"))

#Comparison graph among lines
#Complex 1
f3atp1 <- 
  hyb.atp %>%
  ggplot(aes(x= Population))+
  geom_point(data=plot.dat.hyb.atptetra, aes(y=emmean), size = 7, alpha = .7, shape = c(15, 15, 19, 19))+
  geom_point(aes(y = c1, col = Population), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  geom_linerange(data= plot.dat.hyb.atptetra, aes(ymin=emmean-SE, ymax = emmean+SE), lwd = 1)+
  scale_color_manual(values = c("mediumblue","mediumblue", "firebrick3", "firebrick3"))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('CI ATP Synthesis (nmol '~' min'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")

#Complex 2
f3atp2 <- 
  hyb.atp %>%
  ggplot(aes(x= Population))+
  geom_point(data=plot.dat.hyb.atptetra2, aes(y=emmean), size = 7, alpha = .7, shape = c(15, 15, 19, 19))+
  geom_point(aes(y = c2, col = Population), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  geom_linerange(data= plot.dat.hyb.atptetra2, aes(ymin=emmean-SE, ymax = emmean+SE), lwd = 1)+
  scale_color_manual(values = c("mediumblue","mediumblue", "firebrick3", "firebrick3"))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('CII ATP Synthesis (nmol '~' min'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")



#Atp1 vs dev time
atpdevtimeSDSCN <- ggplot(dtime.astax.atp.sdscn.subset, aes(x=c1, y=dph))+
  geom_point(aes(y=dph, col= gen), size=7)+
  geom_smooth(method="lm", level=.95)+
  geom_errorbar(ymax= dtime.astax.atp.sdscn.subset$dph + dtime.astax.atp.sdscn.subset$dph.sd, ymin = dtime.astax.atp.sdscn.subset$dph-dtime.astax.atp.sdscn.subset$dph.sd)+
  geom_errorbarh(xmax=dtime.astax.atp.sdscn.subset$c1 + (dtime.astax.atp.sdscn.subset$c1.se), xmin = dtime.astax.atp.sdscn.subset$c1 - (dtime.astax.atp.sdscn.subset$c1.se))+
  ylab("Development Time (days)")+
  xlab(bquote('CI ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme_cowplot()+
  theme(legend.position = "top")+
  guides(color = guide_legend(title="Generation:"))+
  annotate("text", 455, 11, label = bquote('p=0.087, R'^2*'=0.751'), size = 4)

mod.dtime.atp <- lm(atpc1~dph.tfs, data=dtime.astax.atp.sdscn)
summary(mod.dtime.atp)

#Atp2 vs dev time
atpdevtimeSDSCN2 <- ggplot(dtime.astax.atp.sdscn.subset, aes(x=c2, y=dph))+
  geom_point(aes(y=dph, col= gen), size=7)+
  geom_smooth(method="lm", level=.95)+
  geom_errorbar(ymax= dtime.astax.atp.sdscn.subset$dph + dtime.astax.atp.sdscn.subset$dph.sd, ymin = dtime.astax.atp.sdscn.subset$dph-dtime.astax.atp.sdscn.subset$dph.sd)+
  geom_errorbarh(xmax=dtime.astax.atp.sdscn.subset$c2 + (dtime.astax.atp.sdscn.subset$c2.se), xmin = dtime.astax.atp.sdscn.subset$c2 - (dtime.astax.atp.sdscn.subset$c2.se))+
  ylab("Development Time (days)")+
  xlab(bquote('CII ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme_cowplot()+
  theme(legend.position = "top")+
  guides(color = guide_legend(title="Generation:"))+
  annotate("text", 1030, 13, label = bquote('p=0.837, R'^2*'=-0.46'), size = 4)

mod.dtime.atp2 <- lm(atpc2~dph.tfs, data=dtime.astax.atp.sdscn)
summary(mod.dtime.atp2)

#ATP vs astaxanthin
f3atpast <- 
  ggplot(dtime.astax.atp.sdscn.subset, aes(x=c1, y=asta))+
  geom_point(aes(y=asta, col= gen), size = 7)+
  geom_smooth(method="lm", level = .95)+
  geom_errorbarh(xmax= dtime.astax.atp.sdscn.subset$c1 + (dtime.astax.atp.sdscn.subset$c1.se), xmin = dtime.astax.atp.sdscn.subset$c1- (dtime.astax.atp.sdscn.subset$c1.se))+
  geom_errorbar(ymax=dtime.astax.atp.sdscn.subset$asta + (dtime.astax.atp.sdscn.subset$asta.se), ymin = dtime.astax.atp.sdscn.subset$asta - (dtime.astax.atp.sdscn.subset$asta.se))+
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))+
  xlab(bquote('CI ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme_cowplot()+
  theme(legend.position = "top")+
  guides(color = guide_legend(title="Generation:"))+
  ylim(-0.03, .25)+
  xlim(350,550)+
  annotate("text", 400, 0.23, label = bquote('p=0.058, R'^2*'=0.84'), size = 4)


#ATP2 vs astaxanthin
f3atpast2 <- 
  ggplot(dtime.astax.atp.sdscn.subset, aes(x=c2, y=asta))+
  geom_point(aes(y=asta, col= gen), size = 7)+
  geom_smooth(method="lm", level = .95)+
  geom_errorbarh(xmax= dtime.astax.atp.sdscn.subset$c2 + (dtime.astax.atp.sdscn.subset$c2.se), xmin = dtime.astax.atp.sdscn.subset$c2- (dtime.astax.atp.sdscn.subset$c2.se))+
  geom_errorbar(ymax=dtime.astax.atp.sdscn.subset$asta + (dtime.astax.atp.sdscn.subset$asta.se), ymin = dtime.astax.atp.sdscn.subset$asta - (dtime.astax.atp.sdscn.subset$asta.se))+
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))+
  xlab(bquote('CII ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme_cowplot()+
  theme(legend.position = "top")+
  guides(color = guide_legend(title="Generation:"))+
  xlim(950,1230)+
  annotate("text", 1150, 0.3, label = bquote('p=0.388, R'^2*'=0.06'), size = 4)

mod.astax.atp <- lm(asta~atpc1, data=dtime.astax.atp.sdscn)
summary(mod.astax.atp)

mod.astax.atp2 <- lm(asta~atpc2, data=dtime.astax.atp.sdscn)
summary(mod.astax.atp2)


### Citrate synthase vs ATP

#CS vs atp1
f3atpcs <- 
  ggplot(dtime.astax.atp.sdscn.subset, aes(x=c1.nonnorm, y=cs))+
  geom_point(aes(y=cs, col= gen), size = 7)+
  geom_smooth(method="lm", level = .95)+
  geom_errorbarh(xmax= dtime.astax.atp.sdscn.subset$c1.nonnorm + (dtime.astax.atp.sdscn.subset$c1.nonnorm.se), xmin = dtime.astax.atp.sdscn.subset$c1.nonnorm- (dtime.astax.atp.sdscn.subset$c1.nonnorm.se))+
  geom_errorbar(ymax=dtime.astax.atp.sdscn.subset$cs + (dtime.astax.atp.sdscn.subset$cs.se), ymin = dtime.astax.atp.sdscn.subset$cs - (dtime.astax.atp.sdscn.subset$cs.se))+
  ylab(bquote('Citrate (nmol'~'min'^-1*')'))+
  xlab(bquote('Unocorrected CI ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme_cowplot()+
  theme(legend.position = "top")+
  guides(color = guide_legend(title="Generation:"))+
  annotate("text", 200, 0.35, label = bquote('p=0.023, R'^2*'=0.933'), size = 4)

#Complex 1 vs ATP model
mod.atp.cs <- lm(cs.tfs~atpc1.nonnorm, data=dtime.astax.atp.sdscn)
summary(mod.atp.cs)

#CS vs atp2
f3atpcs2 <- 
  ggplot(dtime.astax.atp.sdscn.subset, aes(x=c2.nonnorm, y=cs))+
  geom_point(aes(y=cs, col= gen), size = 7)+
  geom_smooth(method="lm", level = .95)+
  geom_errorbarh(xmax= dtime.astax.atp.sdscn.subset$c2.nonnorm + (dtime.astax.atp.sdscn.subset$c2.nonnorm.se), xmin = dtime.astax.atp.sdscn.subset$c2.nonnorm- (dtime.astax.atp.sdscn.subset$c2.nonnorm.se))+
  geom_errorbar(ymax=dtime.astax.atp.sdscn.subset$cs + (dtime.astax.atp.sdscn.subset$cs.se), ymin = dtime.astax.atp.sdscn.subset$cs - (dtime.astax.atp.sdscn.subset$cs.se))+
  ylab(bquote('Citrate (nmol'~'min'^-1*')'))+
  xlab(bquote('Uncorrected CII ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme_cowplot()+
  theme(legend.position = "top")+
  guides(color = guide_legend(title="Generation:"))+
  annotate("text", 425, 0.63, label = bquote('p=0.371, R'^2*'=0.093'), size = 4)

#Run model
mod.atp.cs2 <- lm(cs.tfs~atpc2.nonnorm, data=dtime.astax.atp.sdscn)
summary(mod.atp.cs2)


####Effect of diet analysis in BUFSD19 and BUF ####
#ATP Both
buf19both.1 <- read.csv("diet.effect.csv")

mod.buf19atpboth <- lm(c1~Population, data=buf19both.1) 
comp.buf19atpboth <- emmeans(mod.buf19atpboth, pairwise~Population)
plot.dat.buf19atpboth <- as.data.frame(comp.buf19atpboth$emmeans)

plot.dat.buf19atpboth$Population <- factor(plot.dat.buf19atpboth$Population, levels = c("BUF Tetra.", "BUFSD19 Tetra.", "BUF Yeast", "BUFSD19 Yeast"))
plot.dat.buf19atpboth$gen <- c("parental", "parental", "hybrid", "hybrid")

##Complex 1
#Parental line
b19atp_BUF <- subset(buf19both.1, gen == "parental") %>%
  ggplot(aes(x= Population))+
  geom_point(data=subset(plot.dat.buf19atpboth, gen == "parental"), aes(y=emmean), size = 7, alpha = .7, shape = c(15, 19))+
  geom_point(aes(y = c1, col = Population), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  geom_linerange(data= subset(plot.dat.buf19atpboth, gen == "parental"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  scale_color_manual(values = c("mediumblue", "firebrick3"))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('CI ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")+
  annotate("text", 2, 650, label = bquote('t=2.134, p=0.198'), size = 5)

#Hybrid line
b19atp_BUFSD19 <- subset(buf19both.1, gen == "hybrid") %>%
  ggplot(aes(x= Population))+
  geom_point(data=subset(plot.dat.buf19atpboth, gen == "hybrid"), aes(y=emmean), size = 7, alpha = .7, shape = c(15, 19))+
  geom_point(aes(y = c1, col = Population), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  geom_linerange(data= subset(plot.dat.buf19atpboth, gen == "hybrid"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  scale_color_manual(values = c("mediumblue", "firebrick3"))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('CI ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")+
  annotate("text", 2, 850, label = bquote('t=3.454, p=0.022'), size = 5)


##Complex 2##
mod.buf19atpboth.2 <- lm(c2~Population, data=buf19both.1) 
comp.buf19atpboth.2 <- emmeans(mod.buf19atpboth.2, pairwise~Population)
plot.dat.buf19atpboth.2 <- as.data.frame(comp.buf19atpboth.2$emmeans)

plot.dat.buf19atpboth.2$Population <- factor(plot.dat.buf19atpboth.2$Population, levels = c("BUF Tetra.", "BUFSD19 Tetra.", "BUF Yeast", "BUFSD19 Yeast"))
plot.dat.buf19atpboth.2$gen <- c("parental", "parental", "hybrid", "hybrid")

#Parental line
b19atp_BUF.2 <- subset(buf19both.1, gen == "parental") %>%
  ggplot(aes(x= Population))+
  geom_point(data=subset(plot.dat.buf19atpboth.2, gen == "parental"), aes(y=emmean), size = 7, alpha = .7, shape = c(15, 19))+
  geom_point(aes(y = c2, col = Population), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  geom_linerange(data= subset(plot.dat.buf19atpboth.2, gen == "parental"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  scale_color_manual(values = c("mediumblue", "firebrick3"))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('CII ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")+
  annotate("text", 2, 1000, label = bquote('t=5.121, p=0.001'), size = 5)

#Hybrid line
b19atp_BUFSD19.2 <- subset(buf19both.1, gen == "hybrid") %>%
  ggplot(aes(x= Population))+
  geom_point(data=subset(plot.dat.buf19atpboth.2, gen == "hybrid"), aes(y=emmean), size = 7, alpha = .7, shape = c(15, 19))+
  geom_point(aes(y = c2, col = Population), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  geom_linerange(data= subset(plot.dat.buf19atpboth.2, gen == "hybrid"), aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  scale_color_manual(values = c("mediumblue", "firebrick3"))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('CII ATP Synthesis (nmol'~'min'^-1*')'))  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "")+
  annotate("text", 2, 1400, label = bquote('t=4.180, p=0.006'), size = 5)


## Parental comparison between diets ##

#Devtime
diet.comparison.dat.1 <- read.csv("diet_comparison_data.csv")

diet.mod.1 <- lmer(dph~diet + (1|line), data = diet.comparison.dat.1)
summary(diet.mod.1)

mod.diet.dev <- lm(dph~diet, data=diet.comparison.dat.1) 
comp.diet.dev <- emmeans(mod.diet.dev, pairwise~diet)
plot.dat.diet.dev <- as.data.frame(comp.diet.dev$emmeans)

diet.dev <- diet.comparison.dat.1 %>%
  ggplot(aes(x= diet))+
  geom_point(data=plot.dat.diet.dev, aes(y=emmean), size = 7, alpha = .7, shape = c(15, 19))+
  geom_point(aes(y = dph, col = diet), size = 4, shape = 19, alpha = 0.6, 
             position = position_dodge2(width = .4))+
  geom_linerange(data= plot.dat.diet.dev, aes(ymin=lower.CL, ymax = upper.CL), lwd = 1)+
  scale_color_manual(values = c("mediumblue", "firebrick3"))+
  theme_cowplot()+
  xlab("") +
  ylab(bquote('Days post hatch to copepodid (C1) stage'))  +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 12),
        legend.position = "")+
  annotate("text", 1, 11, label = bquote('t=2.164, p=0.036'), size = 5)


#Summary figure for diet effect on ATP and development
summary_diet_effect <- ggarrange(b19atp_BUF, b19atp_BUFSD19, b19atp_BUF.2, b19atp_BUFSD19.2, diet.dev,  
                                 labels = c("A", "B", "C", "D", "E"), ncol=2, nrow=3, common.legend = , legend = "none")

tiff(file="Fig S1.tif", units="in", width=8, height=12, res=300)
summary_diet_effect
dev.off()


#### Make summary figures for astaxanthin and fitness correlations ####

#Astaxanthin relationship summary Figure 3
fig3 <- ggarrange(p3, p4, astax.dev,  labels = c("A", "B", "C"), ncol=1, nrow=3, common.legend = TRUE, legend = "top")

jpeg(file="Fig3.jpg", units="in", width=5, height=12, res=500)
fig3
dev.off()

#Astaxanthin:dietary carotenoids and astaxanthin:hydroxy relationship summary Figure 
ratio.summary <- ggarrange(p3.ratio, p4.ratio, ratio.dev, 
                                     labels = c("A", "B", "C"), ncol=1, nrow=3, common.legend = TRUE, legend = "top")

jpeg(file="Ratio summary figure.jpg", units="in", width=4, height=11.5, res=300)
ratio.summary
dev.off()


#Fitness variables relationships summary Figure S6
summaryfigureSI <- ggarrange(p5, p6, p7, p8,  labels = c("A", "B", "C", "D"), ncol=2, nrow=2, common.legend = TRUE, legend = "top")

tiff(file="Fig S6.tif", units="in", width=9, height=9, res=300)
summaryfigureSI
dev.off()

#Make summary figure for ATPC1 and ATPC2 across lines Figure S2
summary_ATP_across_lines <- ggarrange(ATPC1BBS_groups, ATPC1ACP_groups, ATPC2BBS_groups,ATPC2ACP_groups,  
                                      labels = c("A", "B", "C", "D"), ncol=1, nrow=4, common.legend = FALSE, legend = "right")

tiff(file="Fig S2.tif", units="in", width=9, height=12, res=300)
summary_ATP_across_lines
dev.off()

#Summary figure for devtime across lines Figure S3
summary_devtime_across_lines <- ggarrange(devtime_across_lines_BBS, devtime_across_lines_ACP,  
                                          labels = c("A", "B"), ncol=1, nrow=2, common.legend = FALSE, legend = "right")
tiff(file="Fig S3.tif", units="in", width=7, height=7, res=300)
summary_devtime_across_lines
dev.off()



#Astaxanthiin relationships summary for SDxSCN cross
summaryfigureSDSCN_asta <- ggarrange(f3atpast, f3atpast2, SDSCNdevtimeast,  labels = c("A", "B", "C"), ncol=1, nrow=3, common.legend = FALSE, legend = "right")

tiff(file="Fig5.tif", units="in", width=6, height=11, res=300)
summaryfigureSDSCN_asta
dev.off()

#Fitness variable relationship summary for SD x SCN cross
summaryfigureSDSCN_fitness <- ggarrange(f3atpcs, f3atpcs2, atpdevtimeSDSCN, atpdevtimeSDSCN2, labels = c("A", "B", "C", "D"), ncol=2, nrow=2, common.legend = TRUE, legend = "top")

tiff(file="Fig S8.tif", units="in", width=9, height=9, res=300)
summaryfigureSDSCN_fitness
dev.off()

#Summary figure for ATP and devtime across SD x SCN 
summary_devtime_ATP_SDSCN <- ggarrange(f3atp1, f3atp2, devtime_across_gen_SDSCN,  
                                             labels = c("A", "B", "C"), ncol=1, nrow=3, common.legend = FALSE, legend = "none")

tiff(file="Fig S9.tif", units="in", width=5, height=11.5, res=300)
summary_devtime_ATP_SDSCN
dev.off()


#### Phylogeny for summary figure ####
library(phytools) #To build trees and make branches ultrametric

tree1='((((BUF, SD, BR,),AB),CAT),(SCN),PES));'
tree1=read.newick(text=tree1)
plot(tree1)

#Convert tree to be ultrametric for plotting purposes
tree1=compute.brlen(tree1,method="Grafen")
is.ultrametric(tree1)

jpeg("Tigriopus tree for figure 1.jpg", units = "in", height = 3, width = 1, res=300)
par(mar=c(0,0,0,0))
plot(tree1)
dev.off()






##### Hybrid Breakdown Astaxanthin Analysis for AB x SD RILs #########

require(ggplot2)
require(nlme)
require (car)
require(stats)
require(lawstat)
library(tidyverse)
library(cowplot)
transdat= read.csv("ABxSD_May_evo_transf.csv")
transdat$gen = as.factor(transdat$gen) ### Added because of R v 4.0 update
transdat$gen = relevel(transdat$gen, ref="Parental")

##################################################################################
######################## Mixed-effects model with equal variance not assumed ############
######################## Overall Hybrid vs Parental ######################
########################################################################################

log.memod= lme(log.astax~ gen ,
               data= transdat, random = ~1|line, weights=varIdent(form=~1|gen), na.action=na.omit)
summary(log.memod)
#Get confidence limits
-0.1690282-(1.96*0.07703210)
-0.1690282+(1.96*0.07703210)


VarCorr(log.memod)
ranef(log.memod)
intervals(log.memod)
plot(log.memod)
hist(resid(log.memod))

log.dens <- density(resid(log.memod))
plot(log.dens)

##Contrasts between parental lines and hybrid lines
transdat$analysisID = relevel(transdat$analysisID, ref = "Parental")
mod.pc = lm(log.astax~analysisID, data = transdat)
summary(mod.pc)
confint(mod.pc)

#############################################################################
################### Mixed-effects model with equal variance assumed ########
############################################################################

log.memod2= lme(log.astax~ gen  ,
                data= transdat, random = ~1|line, na.action=na.omit)
summary(log.memod2)


#################          Comparing the two models        #####################
#################        with equal vs unequal variances   #####################
################################################################################

model.test <-anova(log.memod2, log.memod )

print(model.test)


#################################################################################
##################    Variance test    ##########################################
#################################################################################


plot( log.memod2, resid(., type = "p") ~ fitted(.) | gen,
      id = 0.05, adj = -0.3 )

bartlett.test(log.astax~gen, data=transdat)



############# PLot figure S5 ####################################
#Make smaller dataset with values condensed 
plotdat= read.csv("5-3-2017-estimates-logtransf.csv")
plotdat2 = plotdat[!(plotdat$line %in% c("P-1", "P-2", "P-3", "P-4",
                                         "P-5", "P-6", "P-7", "P-8", "P-9")),]

plotdat2 %>%
  group_by(line, gen) %>%
  summarize(est = mean(est), se=mean(se)) -> plotdat3


avgplot <- ggplot(plotdat2, aes(x=line))+
  geom_point(data = plotdat3, aes(y=est, fill = "black"),size =8,shape =16, alpha = 0.7) +
  geom_point(aes(y = log.astax), col = "red", size = 4, alpha = 0.7) +
  scale_shape(solid=TRUE) +
  scale_size_continuous(range= c(4,20)) +
  geom_errorbar(aes(ymin = est-se, ymax = est+se), width = 0.1) +
  theme_cowplot() +
  scale_x_discrete(labels = c("1a" = "Parental",
                              "1b" = "Hybrid", "A1" = "A", "A5" = "B", "B2" = "C", "B3" = "D", "D2" = "E", "D3" = "F", "D6" = "G", "E2" = "H")) +
  geom_vline(xintercept= 2.5,linetype = "solid", lwd=1) +
  geom_hline(aes(yintercept = plotdat2$est[1]),linetype = "dashed", lwd=1)+
  xlab(" Line") + ylab("log Astaxanthin concentration (ng/mg)")+
  theme(strip.text = element_text(size = 16),
        axis.text = element_text(size=14), axis.title = element_text(size=17))+
  theme(legend.position = "")
avgplot

tiff(file = "Fig S5.tif", units = "in", width = 8, height = 6, res = 300)
avgplot
dev.off()

#### Extra figure showing median and IQRs because of skewed residual distribution ####

transdat$ng.astax = transdat$astax.mass*1000

jpeg(file = "Median Hybrid Parental Carot.jpg", units = "in", width = 5.5, height = 5.5, res = 300)

transdat %>% 
  group_by(gen) %>% 
  summarise(med.astax = median(ng.astax),
            liqr = quantile(ng.astax, prob = 0.25),
            hiqr = quantile(ng.astax, prob = 0.75)) %>% 
  ggplot(aes(x = gen)) +
  geom_point(data = transdat, aes(y = ng.astax, col = gen),size = 6, 
             alpha =0.7)+
  geom_point( aes(y = med.astax, fill ="black"), size = 12, shape =15, alpha = 0.7,
              position = position_dodge(width = 0.4))+
  geom_errorbar( aes(ymin = liqr, ymax = hiqr),
                 position = position_dodge(width = 0.4), width = 0.1) +
  theme_cowplot()+
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Parental", "Hybrid"),
                     name = "") +
  scale_fill_manual(values = c("blue", "red"),
                    labels = c("Parental", "Hybrid"),
                    name = "") +
  xlab("") +
  ylab("Astaxanthin (ng/mg)") +
  #scale_y_continuous(breaks = seq(0,350,50), limits = c(0,350))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.title.y = element_text(size = 20, margin = margin(0,20,0,0)))  +
  theme(legend.position = "none")
dev.off()

kruskal.test(ng.astax ~ gen, data = transdat)



####Begin plotting of HPLC chromatograms####
datumhplc = read.csv(file = "HPLC_chromatograms.csv")

#Trim first 12 minutes off to remove solvent front peak and reduce white space of graph
#trim off last 6 minutes to remove equilibration period and reduce white space of graph
datumhplctrim = subset(datumhplc, Time >=13 & Time <=26)

#Subset desired samples from data frame
dd_sub = datumhplctrim[,c(1,2,3,4,5)]
dd_subtetra = datumhplctrim[,c(1,6)]

#You will notice numbers after the Y variable call for the standard mix. This is just an adjustment to 
#all of the values of the intensity column for that standard to correct for the natural baseline drift
#of the HPLC system

#Adjust y values for baseline drift during hplc
dd_sub$Intensitystandardmix = dd_sub$Intensitystandardmix+200

##Then rearrange your data frame
dd = melt(dd_sub, id=c("Time"))
ddtetra = melt(dd_subtetra, id=c("Time"))

#Make plot
chromatograms <- ggplot(dd) + geom_line(aes(x=Time, y=value, colour=variable)) +
  scale_colour_manual(name="",values=c("black","goldenrod3", "red3", "blue2"),
                      labels=c("Tigriopus Sample 131 (low ratio example)","Tigriopus Sample 142 (high ratio example)", "Hydroxyechinenone standard", "Standard Mix"))+
  theme_cowplot()+
  scale_y_continuous(breaks=seq(0,4000,1000))+
  scale_x_continuous(breaks=seq(15,27,3))+
  ylab("Intensity (mv)")+ xlab("Time (min)")+
  annotate("text", x = 16, y = 2800, label = "1",size = 4, color ="black")+
  annotate("text", x = 17, y = 3300, label = "2", size = 4, color ="black")+
  annotate("text", x = 18.5, y = 3300, label = "3", size = 4, color ="black")+
  annotate("text", x = 20.5, y = 1800, label = "4", size = 4, color ="black")+
  annotate("text", x = 23.5, y = 2500, label = "5",  size = 4, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 2, ncol = 2))+
  theme(legend.position = "top", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
#View plot
chromatograms


#Make plot #2
chromatogramtetra <- ggplot(ddtetra) + geom_line(aes(x=Time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("darkgreen"),
                      labels=c("Tetraselmis algae"))+
  theme_cowplot()+
  scale_y_continuous(breaks=seq(0,40000,5000))+
  scale_x_continuous(breaks=seq(15,27,3))+
  ylab("Intensity (mv)")+ xlab("Time (min)")+
  annotate("text", x = 13.7, y = 5000, label = "1",size = 4, color ="black")+
  annotate("text", x = 15.5, y = 5000, label = "2", size = 4, color ="black")+
  annotate("text", x = 17, y = 33500, label = "3", size = 4, color ="black")+
  annotate("text", x = 20, y = 34000, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 13000, label = "5",  size = 4, color ="black")+
  annotate("text", x = 24, y = 24000, label = "6",  size = 4, color ="black")+
  annotate("text", x = 14, y = 35000, label = "Tetraselmis algae",  size = 4, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 4, ncol = 1))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
#View plot
chromatogramtetra



figurehplcchroma <- ggarrange(chromatograms, chromatogramtetra,
                     labels = c("A", "B"),
                     font.label = list(size=12), hjust = -0.1,vjust = 1.5, 
                     ncol=1, nrow=2)

jpeg(file="HPLC_chromatograms.jpg", units="in", width=8, height=10, res=500)
figurehplcchroma
dev.off()