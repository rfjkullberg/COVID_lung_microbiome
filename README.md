# Lung microbiota in COVID-19-related ARDS

This is the code used for the main analyses in "Lung microbiota of critically ill COVID-19 patients are associated with non-resolving acute respiratory distress syndrome" (submitted). For questions: Bob Kullberg at r.f.j.kullberg@amsterdamumc.nl

## Step 1 - Load libraries
```
library(tidyverse)
library(yingtools2)
library(phyloseq)
library(vegan)
library(microbiome)
library(RColorBrewer) 
library(cowplot)
library(ggpubr) 
library(ggrepel)
library(scales)
library(DESeq2)
library(randomForest)
library(readxl)
library(survival)
library(survminer)
library(cmprsk)
```

## Step 2 - Load data (phyloseq file)
```
ps <- readRDS("~/Documents/COVID ARDS/Data/phyloseq-covid-ards.RDS") # Main phyloseq file
metadata <- read_excel("~/Documents/COVID ARDS/Data/metadata-covid-ards.xlsx") # Metadata file
microbialburden <- read_excel("~/Documents/COVID ARDS/Data/microbialburden-covid-ards.xlsx") # Data file with the microbial (bacterial/fungal) burden for all samples (patients and negative controls)
```

Samples that did not pass quality checks (details provided in the manuscript) were already filtered out using the following code:
```
metadata <- metadata %>%
  filter(reads_mapped > 199)
sample_data(ps) <- set.samp(metadata)
```

```
table(metadata$neg_control) # 16S rRNA gene sequencing data of 163 samples from COVID-19 patients and 32 negative experimental controls. 
table(microbialburden$neg_control) # qPCR data of 163 samples from COVID-19 patients and 46 negative experimental controls.  
ps # 2592 taxa and 195 samples
```

## Step 3 - Lung microbiota of BAL samples are distinct from negative experimental controls

First, we compared the microbial burden of our BAL samples with negative experimental control samples. 
```
comparisons <- list(c("y", "n"))
colours <- c("#1b849e","#9e351b")

# Bacterial
ggplot(microbialburden, aes(x=neg_control, y=`16s_quant`, fill=neg_control))+ 
  geom_boxplot(outlier.shape = NA, alpha = .20) + 
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+ 
  theme_cowplot(12)+  
  scale_fill_manual(values = colours) + 
  scale_color_manual(values = colours) + 
  xlab("")+ 
  ylab("Bacterial DNA Burden")+ 
  scale_y_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+ 
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(),
        legend.title = element_blank())+
  labs(fill = "")

# Fungal
microbialburden %>%
  mutate(pseudocount_18s_quant = (`18s_quant`+(0.5*5.80510960001517E-07))) %>% # add pseudocount (50% of lowest value) to enable log transformations
  ggplot(aes(x=neg_control, y=`pseudocount_18s_quant`, fill=neg_control))+ 
  geom_boxplot(outlier.shape = NA, alpha = .20) + 
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+ 
  theme_cowplot(12)+  
  scale_fill_manual(values = colours) + 
  scale_color_manual(values = colours) + 
  xlab("")+ 
  ylab("Fungal DNA Burden")+ 
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+ 
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(),
        legend.title = element_blank())+
  labs(fill = "")
rm(microbialburden)
```
```
ggplot(metadata, aes(x=neg_control, y=reads_mapped, fill=neg_control))+ 
  geom_boxplot(outlier.shape = NA, alpha = .20) + 
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+ 
  theme_cowplot(12)+  
  scale_fill_manual(values = colours) + 
  scale_color_manual(values = colours) + 
  xlab("")+ 
  ylab("High-quality reads")+
  scale_y_continuous(trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(),
        legend.title = element_blank())+
labs(fill = "")
```

Next, we compared the community composition of BAL samples and negative controls. 
```
set.seed(711)
wunifrac <- phyloseq::distance(ps, method = "wunifrac")
ord.wuni <- ordinate(ps, method = 'PCoA', distance = wunifrac) 
wuni_df <- as.data.frame(as.matrix(ord.wuni$vectors))
wuni_df$sample <- row.names(wuni_df) 
wuni_df <- wuni_df %>%
  left_join(metadata) %>%
  select(Axis.1, Axis.2, sample, neg_control)
centroids <- aggregate(cbind(Axis.1,Axis.2)~ neg_control,data= wuni_df, mean)
wuni_df <- merge(wuni_df,centroids,by="neg_control",suffixes=c("",".centroid"))
ggplot(wuni_df, aes(Axis.1, Axis.2, color = neg_control)) +
  geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, color= neg_control), alpha = 0.3)+ 
  geom_point(data=wuni_df,aes(color=neg_control),size=3.5,alpha=1.0) + 
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, label=c("Patients", "Negative Controls")), size=6, fill = ggplot2::alpha(c("white"),0.76)) +
  theme_bw() +
  xlab("Axis 1") +
  ylab("Axis 2") + 
  scale_color_manual(values=colours) +
  theme(legend.position = "none")
adonis2((wunifrac ~ neg_control), by = 'margin', data = metadata, permutations =9999)
rm(wunifrac, wuni_df, ord.wuni, centroids)

set.seed(711)
unifrac <- phyloseq::distance(ps, method = "unifrac")
ord.uni <- ordinate(ps, method = 'PCoA', distance = unifrac) 
uni_df <- as.data.frame(as.matrix(ord.uni$vectors))
uni_df$sample <- row.names(uni_df) 
uni_df <- uni_df %>%
  left_join(metadata) %>%
  select(Axis.1, Axis.2, sample, neg_control)
centroids <- aggregate(cbind(Axis.1,Axis.2)~ neg_control,data= uni_df, mean)
uni_df <- merge(uni_df,centroids,by="neg_control",suffixes=c("",".centroid"))
ggplot(uni_df, aes(Axis.1, Axis.2, color = neg_control)) +
  geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, color= neg_control), alpha = 0.3)+ 
  geom_point(data=uni_df,aes(color=neg_control),size=3.5,alpha=1.0) + 
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, label=c("Patients", "Negative Controls")), size=6, fill = ggplot2::alpha(c("white"),0.76)) +
  theme_bw() +
  xlab("Axis 1") +
  ylab("Axis 2") + 
  scale_color_manual(values=colours) +
  theme(legend.position = "none")
adonis2((unifrac ~ neg_control), by = 'margin', data = metadata, permutations =9999)
rm(unifrac, uni_df, ord.uni, centroids)
rm(comparisons, colours)
```

```
tax <- get.tax(ps)
ps_gen <- phyloseq::tax_glom(ps, "Genus", NArm = TRUE)
taxgen <- get.tax(ps_gen) %>%
  select(-Species)
for (i in 1:6){taxgen[,i] <- as.character(taxgen[,i])} 
taxgen[is.na(taxgen)] <- "unspecified" 
tax_table(ps_gen)<-set.tax(taxgen) 
topN <- 15 
most_abundant_gen <- as.data.frame(sort(taxa_sums(ps_gen), TRUE)[1:topN])
most_abundant_gen$otu <- row.names(most_abundant_gen)
composition <- get.otu.melt(ps_gen, filter.zero =F) %>%
  mutate(Genus = if_else(otu %in% most_abundant_gen$otu, Genus, "Others")) %>%
  mutate(Genus = fct_relevel(Genus, "Massilia", "Staphylococcus", "Ralstonia", "Streptococcus",
                             "Prevotella", "Enterococcus", "Pseudarthrobacter", "Methylobacterium", 
                             "Alloprevotella",  "Klebsiella","Bifidobacterium","Lactobacillus",
                             "Peptostreptococcus","Prevotella_7","Porphyromonas", "Others")) %>%
  mutate(Phylum = if_else(Phylum == "Firmicutes", "Firmicutes", 
                          if_else(Phylum == "Bacteroidetes", "Bacteroidetes", 
                                  if_else(Phylum == "Proteobacteria", "Proteobacteria", 
                                          if_else(Phylum == "Actinobacteria", "Actinobacteria", "Other"))))) %>%
  mutate(Phylum = if_else(Genus == "Others", "Other", Phylum)) %>%
  mutate(Phylum = fct_relevel(Phylum, "Proteobacteria", "Firmicutes", "Bacteroidetes", 
                              "Actinobacteria", "Other")) 
comp_neg <- composition %>%
  group_by(Phylum, Genus, neg_control) %>%
  dplyr::summarize(mean = mean(pctseqs), sd = sd(pctseqs))
getPalette <- c("#849e1b",  "#9e351b", "#1b849e", "#bea93b", "#808080")

ggplot(comp_neg, aes(x=Genus, y=mean, fill=Phylum))+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(~neg_control) +
  scale_fill_manual(values = getPalette, na.value="grey50")+
  theme_bw(base_size=15)+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Mean relative abundance")+ 
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
```

We used a Wilcoxon rank sum test (with Benjamini-Hochberg correction for multiple comparisons) to identify differences in the relative abundance of genera between groups.
```
p.genus <- composition %>%
  group_by(Genus) %>%
  summarize(p = wilcox.test(pctseqs ~ neg_control)$p.value)

p.genus$BH = p.adjust(p.genus$p, method = "BH")
rm(comp_neg, composition, most_abundant_gen, p.genus, ps_gen, tax, taxgen, topN, getPalette)
```

## Step 4 - Lung bacterial and fungal burden are associated with clinical outcomes of COVID-19-related ARDS
Create metadata and phyloseq file with one sample per patient and without negative experimental controls
```
m <- metadata %>%
  filter(bal_nr == "1")
P <- ps
sample_data(P) <- set.samp(m)
```

We compared patients who reached the primary outcome with patients who died within 60 days or were still mechanically ventilated at day 60. 
```
alpha <- estimate_richness(P)
alpha$sample <- row.names(alpha)
alpha$sample <- gsub("X","",as.character(alpha$sample))

df <- left_join(alpha, m) %>%
  mutate(vfd60 = as.factor(if_else(ventfreedays60_intub == "0", "Deceased or Intubated", "Extubated"))) %>%
  mutate(event = if_else(mortality_d60_intub == "1", "Deceased",  # if a patients is deceased at day 60, the event is 'deceased'
                         if_else(ventfreedays60_intub == "0", "Intubated", "Extubated"))) %>% # if a patients is not deceased after 60d, they are either still intubated (=zero ventilator-free days) or extubated
  mutate(event = as.factor(event)) %>%
  mutate(event = fct_relevel(event, "Intubated", "Extubated", "Deceased")) %>%
  mutate(time = if_else(mechvent_length > 60, 60, mechvent_length)) %>%
  mutate(time = as.numeric(time)) %>%
  mutate(log16s = log10(as.numeric(`16s`))) %>%
  mutate(`18s_pseudocount` = `18s`+0.009018145) %>% # add pseudocount (50% of lowest value) to enable log transformations
  mutate(log18s = log10(as.numeric(`18s_pseudocount`))) %>%
  mutate(tertiles = ntile(`log16s`, 3)) %>%
  mutate(tertiles = if_else(tertiles == 1, 'Low burden', if_else(tertiles == 2, 'Intermediate burden', 'High burden'))) %>%
  mutate(tertiles = as.factor(tertiles)) %>%
  mutate(tertiles_fung = ntile(`log18s`, 3)) %>%
  mutate(tertiles_fung = if_else(tertiles_fung == 1, 'Low fungal burden', if_else(tertiles_fung == 2, 'Intermediate fungal burden', 'High fungal burden'))) %>%
  mutate(tertiles_fung = as.factor(tertiles_fung)) %>%
  mutate(tertiles_alpha = ntile(Shannon, 3)) %>%
  mutate(tertiles_alpha = if_else(tertiles_alpha == 1, 'Low diversity', if_else(tertiles_alpha == 2, 'Intermediate diversity', 'High diversity'))) %>%
  mutate(tertiles_alpha = as.factor(tertiles_alpha))
df$tertiles <- relevel(df$tertiles, ref = "Low burden")
df$tertiles_fung <- relevel(df$tertiles_fung, ref = "Low fungal burden")
df$tertiles_alpha <- relevel(df$tertiles_alpha, ref = "Low diversity")
df$event <- relevel(df$event, ref = "Intubated")
df$bal_bact <- relevel(df$bal_bact, ref = "0")
df$vfd60 <- relevel(df$vfd60, ref = "Extubated")

rm(alpha)
comparisons <- list(c("Extubated", "Deceased or Intubated"))

df %>%
  ggplot(aes(x = vfd60, y = `16s`, fill = vfd60)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+
  scale_y_continuous(trans=log_epsilon_trans(epsilon=0.5),limits=c(0,750))+ 
  geom_hline(linetype = 2, yintercept = 1.32) +
  geom_hline(linetype = 2, yintercept = 0.293) +
  scale_colour_manual(values = c("#297568","#752936"))+
  scale_fill_manual(values = c("#297568","#752936"))+
  theme_bw()+
  xlab("")+
  ylab("Bacterial burden")+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+
  theme(legend.position = "none")

df %>%
  ggplot(aes(x = vfd60, y = `18s_pseudocount`, fill = vfd60)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+
  scale_y_continuous(trans=log_epsilon_trans(epsilon=0.1),limits=c(0,5000))+  
  geom_hline(linetype = 2, yintercept = 0.63) +
  geom_hline(linetype = 2, yintercept = 0.1625) +
  scale_colour_manual(values = c("#297568","#752936"))+
  scale_fill_manual(values = c("#297568","#752936"))+
  theme_bw()+
  xlab("")+
  ylab("Fungal burden")+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+
  theme(legend.position = "none")

df %>%
  ggplot(aes(x = vfd60, y = Shannon, fill = vfd60)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+
  geom_hline(linetype = 2, yintercept = 2.167) +
  geom_hline(linetype = 2, yintercept = 1.5873) +
  scale_colour_manual(values = c("#297568","#752936"))+
  scale_fill_manual(values = c("#297568","#752936"))+
  theme_bw()+
  xlab("")+
  ylab("Shannon diversity")+
  theme(legend.position = "none")
```

For time-to-event analyses, we caculated SHRs using competing risk regression models (compete package). 
```
extdat <- finegray(Surv(time, event) ~ ., data=df, etype="Extubated")

summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ log16s, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ tertiles, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ log18s, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ tertiles_fung, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ Shannon, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ tertiles_alpha, data=extdat,weight= fgwt))

df <- df %>%
  mutate(tertiles_n = as.numeric(if_else(tertiles == "Low burden", "1", if_else(tertiles == 'Intermediate burden', "2", '3')))) %>%
  mutate(tertiles_fungal_n = as.numeric(if_else(tertiles_fung == "Low fungal burden", "1", if_else(tertiles_fung == 'Intermediate fungal burden', "2", '3')))) %>%
  mutate(tertiles_alpha_n = as.numeric(if_else(tertiles_alpha == "Low diversity", "1", if_else(tertiles_alpha == 'Intermediate diversity', "2", '3')))) %>%
  mutate(event_n = as.numeric(if_else(event == "Intubated", "0", if_else(event == "Extubated", "1", "2"))))
attach(df)
cif1 <- cuminc(ftime=time,fstatus=event_n,group = tertiles_n)
ggcompetingrisks(fit = cif1, multiple_panels = T,
                 xlab = "Days since start of mechanical ventilation",
                 ylab = "% of patients",
                 palette = c("#1b849e", "#bea93b"),
                 #orientation = "reverse",
                 legend = "top", legend.title="",  xlim = c(0, 60))
cif2 <- cuminc(ftime=time,fstatus=event_n,group = tertiles_fungal_n)
ggcompetingrisks(fit = cif2, multiple_panels = T,
                 xlab = "Days since start of mechanical ventilation",
                 ylab = "% of patients",
                 palette = c("#1b849e", "#bea93b"),
                 #orientation = "reverse",
                 legend = "top", legend.title="",  xlim = c(0, 60))               
cif3 <- cuminc(ftime=time,fstatus=event_n,group = tertiles_alpha_n)
ggcompetingrisks(fit = cif3, multiple_panels = T,
                 xlab = "Days since start of mechanical ventilation",
                 ylab = "% of patients",
                 palette = c("#1b849e", "#bea93b"),
                 #orientation = "reverse",
                 legend = "top", legend.title="",  xlim = c(0, 60))                
 rm(cif1, cif2, cif3)                
```

Univariable and multivariable competing risk regression models (Table E3)
```
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ age, data=extdat,weight= fgwt)) # univariable; change age for sex, bmi, bal_sofa_total, bal_bact, bal_aspergillosis, antibiotics_before_ICU_admission, antibiotic_exposure_score, antifungal
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ log16s+age+sex+bmi+bal_sofa_total+bal_bact+bal_aspergillosis+antibiotics_before_ICU_admission+antibiotic_exposure_score+antifungal,  
      data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ log18s+age+sex+bmi+bal_sofa_total+bal_bact+bal_aspergillosis+antibiotics_before_ICU_admission+antibiotic_exposure_score+antifungal,  
      data=extdat,weight= fgwt))
rm(extdat, alpha)
```

Subdistribution hazard ratios for mortality and cause-specific hazard ratios for extubation and mortality (Table E4)
```
summary(coxph(Surv(time, event) ~ log16s, data = df, id=sample_id)) # change log16s for log18s, age, sex, bmi, bal_sofa_total, bal_bact, bal_aspergillosis, antibiotics_before_ICU_admission, antibiotic_exposure_score, antifungal

extdat <- finegray(Surv(time, event) ~ ., data=df, etype="Deceased")
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ log16s, data=extdat,weight= fgwt)) # change log16s for log18s, age, sex, bmi, bal_sofa_total, bal_bact, bal_aspergillosis, antibiotics_before_ICU_admission, antibiotic_exposure_score, antifungal
```

Two sensitivity analyses were performed to test particular assumptions of above-described models. 
First, adjudication of clinical outcomes at 60 days following sample collection (instead of 60 days following intubation) to account for differences in time between intubation and sample collection. 

```
df <- df %>%
  dplyr::select(-vfd60, -event, -time, -event_n) %>%
  mutate(vfd60_collection = as.factor(if_else(ventfreedays60_bal == "0", "Deceased or Intubated", "Extubated"))) %>%
  mutate(event = if_else(mortality_d60_bal == "1", "Deceased",  
                         if_else(ventfreedays60_bal == "0", "Intubated", "Extubated"))) %>% 
  mutate(event = as.factor(event)) %>%
  mutate(event = fct_relevel(event, "Intubated", "Extubated", "Deceased")) %>%
  mutate(time = if_else(mechvent_since_bal > 60, 60, mechvent_since_bal)) %>%
  mutate(time = as.numeric(time))
df$event <- relevel(df$event, ref = "Intubated")
df$vfd60_collection <- relevel(df$vfd60_collection, ref = "Extubated")

comparisons <- list(c("Extubated", "Deceased or Intubated"))

df %>%
  ggplot(aes(x = vfd60_collection, y = `16s`, fill = vfd60_collection)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+
  scale_y_continuous(trans=log_epsilon_trans(epsilon=0.5),limits=c(0,750))+ 
  geom_hline(linetype = 2, yintercept = 1.32) +
  geom_hline(linetype = 2, yintercept = 0.293) +
  scale_colour_manual(values = c("#297568","#752936"))+
  scale_fill_manual(values = c("#297568","#752936"))+
  theme_bw()+
  xlab("")+
  ylab("Bacterial burden")+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+
  theme(legend.position = "none")

df %>%
  ggplot(aes(x = vfd60_collection, y = `18s_pseudocount`, fill = vfd60_collection)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+
  scale_y_continuous(trans=log_epsilon_trans(epsilon=0.1),limits=c(0,5000))+  
  geom_hline(linetype = 2, yintercept = 0.63) +
  geom_hline(linetype = 2, yintercept = 0.1625) +
  scale_colour_manual(values = c("#297568","#752936"))+
  scale_fill_manual(values = c("#297568","#752936"))+
  theme_bw()+
  xlab("")+
  ylab("Fungal burden")+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+
  theme(legend.position = "none")

df %>%
  ggplot(aes(x = vfd60_collection, y = Shannon, fill = vfd60_collection)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+ # add dots for each observation
  geom_hline(linetype = 2, yintercept = 2.167) +
  geom_hline(linetype = 2, yintercept = 1.5873) +
  scale_colour_manual(values = c("#297568","#752936"))+
  scale_fill_manual(values = c("#297568","#752936"))+
  theme_bw()+
  xlab("")+
  ylab("Shannon diversity")+
  theme(legend.position = "none")
```
```
extdat <- finegray(Surv(time, event) ~ ., data=df, etype="Extubated")

summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ log16s, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ tertiles, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ log18s, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ tertiles_fung, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ Shannon, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ tertiles_alpha, data=extdat,weight= fgwt))
rm(extdat)

df <- df %>%
  mutate(event_n = as.numeric(if_else(event == "Intubated", "0", if_else(event == "Extubated", "1", "2"))))

attach(df)
cif1 <- cuminc(ftime=time,fstatus=event_n,group = tertiles_n)
ggcompetingrisks(fit = cif1, multiple_panels = T,
                 xlab = "Days since first BAL",
                 ylab = "% of patients",
                 palette = c("#1b849e", "#bea93b"),
                 #orientation = "reverse",
                 legend = "top", legend.title="",  xlim = c(0, 60))
cif2 <- cuminc(ftime=time,fstatus=event_n,group = tertiles_fungal_n)
ggcompetingrisks(fit = cif2, multiple_panels = T,
                 xlab = "Days since first BAL",
                 ylab = "% of patients",
                 palette = c("#1b849e", "#bea93b"),
                 #orientation = "reverse",
                 legend = "top", legend.title="",  xlim = c(0, 60))               
cif3 <- cuminc(ftime=time,fstatus=event_n,group = tertiles_alpha_n)
ggcompetingrisks(fit = cif3, multiple_panels = T,
                 xlab = "Days since first BAL",
                 ylab = "% of patients",
                 palette = c("#1b849e", "#bea93b"),
                 #orientation = "reverse",
                 legend = "top", legend.title="",  xlim = c(0, 60))                
rm(cif1, cif2, cif3)  
```

Second, we extended our follow up and assessed clinical outcomes at 90 days following sample collection

```
df <- df %>%
  dplyr::select(-vfd60_collection, -event, -time, -event_n) %>%
  mutate(vfd90 = as.factor(if_else(ventfreedays90_bal == "0", "Deceased or Intubated", "Extubated"))) %>%
  mutate(event = if_else(mortality_d90_bal == "1", "Deceased",  
                         if_else(ventfreedays90_bal == "0", "Intubated", "Extubated"))) %>% 
  mutate(event = as.factor(event)) %>%
  mutate(event = fct_relevel(event, "Intubated", "Extubated", "Deceased")) %>%
  mutate(time = if_else(mechvent_length > 90, 90, mechvent_length)) %>%
  mutate(time = as.numeric(time))
df$event <- relevel(df$event, ref = "Intubated")
df$vfd90 <- relevel(df$vfd90, ref = "Extubated")

comparisons <- list(c("Extubated", "Deceased or Intubated"))

df %>%
  ggplot(aes(x = vfd90, y = `16s`, fill = vfd90)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+
  scale_y_continuous(trans=log_epsilon_trans(epsilon=0.5),limits=c(0,750))+ 
  geom_hline(linetype = 2, yintercept = 1.32) +
  geom_hline(linetype = 2, yintercept = 0.293) +
  scale_colour_manual(values = c("#297568","#752936"))+
  scale_fill_manual(values = c("#297568","#752936"))+
  theme_bw()+
  xlab("")+
  ylab("Bacterial burden")+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+
  theme(legend.position = "none")

df %>%
  ggplot(aes(x = vfd90, y = `18s_pseudocount`, fill = vfd90)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+
  scale_y_continuous(trans=log_epsilon_trans(epsilon=0.1),limits=c(0,5000))+  
  geom_hline(linetype = 2, yintercept = 0.63) +
  geom_hline(linetype = 2, yintercept = 0.1625) +
  scale_colour_manual(values = c("#297568","#752936"))+
  scale_fill_manual(values = c("#297568","#752936"))+
  theme_bw()+
  xlab("")+
  ylab("Fungal burden")+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+
  theme(legend.position = "none")

df %>%
  ggplot(aes(x = vfd90, y = Shannon, fill = vfd90)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+ # add dots for each observation
  geom_hline(linetype = 2, yintercept = 2.167) +
  geom_hline(linetype = 2, yintercept = 1.5873) +
  scale_colour_manual(values = c("#297568","#752936"))+
  scale_fill_manual(values = c("#297568","#752936"))+
  theme_bw()+
  xlab("")+
  ylab("Shannon diversity")+
  theme(legend.position = "none")
```
```
extdat <- finegray(Surv(time, event) ~ ., data=df, etype="Extubated")

summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ log16s, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ tertiles, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ log18s, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ tertiles_fung, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ Shannon, data=extdat,weight= fgwt))
summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ tertiles_alpha, data=extdat,weight= fgwt))
rm(extdat)

df <- df %>%
  mutate(event_n = as.numeric(if_else(event == "Intubated", "0", if_else(event == "Extubated", "1", "2"))))

attach(df)
cif1 <- cuminc(ftime=time,fstatus=event_n,group = tertiles_n)
ggcompetingrisks(fit = cif1, multiple_panels = T,
                 xlab = "Days since first BAL",
                 ylab = "% of patients",
                 palette = c("#1b849e", "#bea93b"),
                 #orientation = "reverse",
                 legend = "top", legend.title="",  xlim = c(0, 90))
cif2 <- cuminc(ftime=time,fstatus=event_n,group = tertiles_fungal_n)
ggcompetingrisks(fit = cif2, multiple_panels = T,
                 xlab = "Days since first BAL",
                 ylab = "% of patients",
                 palette = c("#1b849e", "#bea93b"),
                 #orientation = "reverse",
                 legend = "top", legend.title="",  xlim = c(0, 90))               
cif3 <- cuminc(ftime=time,fstatus=event_n,group = tertiles_alpha_n)
ggcompetingrisks(fit = cif3, multiple_panels = T,
                 xlab = "Days since first BAL",
                 ylab = "% of patients",
                 palette = c("#1b849e", "#bea93b"),
                 #orientation = "reverse",
                 legend = "top", legend.title="",  xlim = c(0, 90))                
rm(cif1, cif2, cif3)  
rm(df)
```


## Step 5 - Lung microbiota community composition is associated with successful extubation in COVID-19-related ARDS
We visualized lung microbial composition of the first available BAL sample per patient (n=114) by principal coordinate analysis
```
p1 <- microbiome::transform(P, "compositional")
set.seed(88)
wunifrac <- phyloseq::distance(p1, method = "wunifrac") 
adonis2((wunifrac ~ vfd60), by = 'margin', data = m, permutations =9999) 
#permutest(betadisper(wunifrac, m$vfd60), permutations = 9999)
adonis2((wunifrac ~ vfd60+age+sex+bmi+bal_sofa_total+bal_bact+bal_aspergillosis+antibiotics_before_ICU_admission+antibiotic_exposure_score+antifungal), by = 'margin', data = m, permutations =9999) 

adonis2((wunifrac ~ log10(as.numeric(`16s`))), by = 'margin', data = m, permutations =9999) 
adonis2((wunifrac ~ vfd60 + log10(as.numeric(`16s`))), by = 'margin', data = m, permutations =9999)

ord.wuni <- ordinate(p1, method = 'PCoA', distance = wunifrac) 
wuni_df <- as.data.frame(as.matrix(ord.wuni$vectors)) 
wuni_df$sample <- row.names(wuni_df) 
wuni_df <- wuni_df %>%
  left_join(m) %>%
  select(Axis.1, Axis.2, sample, vfd60)
centroids <- aggregate(cbind(Axis.1,Axis.2)~ vfd60, data= wuni_df, mean)
wuni_df <- merge(wuni_df,centroids,by="vfd60",suffixes=c("",".centroid"))

ggplot(wuni_df, aes(Axis.1, Axis.2, color = vfd60)) +
  geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, color= vfd60), alpha = 0.3)+ 
  geom_point(data=wuni_df,aes(color=vfd60),size=5,alpha=1.0) +
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, label=c("Deceased or Intubated", "Extubated")), size=10, fill = ggplot2::alpha(c("white"),0.76)) +
  scale_colour_manual(values = c("#752936","#44aa99"))+
  theme_bw() +
  xlab("Axis 1") +
  ylab("Axis 2") + 
  theme(legend.position = "none")
```

To identify individual bacterial taxa associated with successful extubation, biplot and rank abundance analyses were used
```
df.fam <- plot_ordination(p1, ord.wuni, type="species", color="Family", justDF = T) %>% 
  select(Family, Axis.1, Axis.2)%>%
  group_by(Family) %>%
  summarise_at(c("Axis.1", "Axis.2"), mean, na.rm=T)

df.fam <- df.fam %>%
  mutate(x = as.numeric(Axis.1)) %>%
  mutate(y = as.numeric(Axis.2)) %>%
  mutate(length = ((x*x)+(y*y))) %>%
  select(-x, -y) 

ggplot() +
  geom_segment(data = df.fam, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), arrow=arrow(), size=0.8) + #lines per Family
  geom_text(data = df.fam, aes(x = Axis.1, y = Axis.2, label = Family), size = 6, color="#1b849e")+ #names per Family
  theme_bw() + 
  theme(legend.position = "none")+
  xlab("Axis 1") + #Label X-axis
  ylab("Axis 2") + #Label Y-axis
  labs(color="")
  
rm(p1, df.fam, wuni_df, ord.wuni, centroids)
```

```
tax <- get.tax(P)
ps_gen <- phyloseq::tax_glom(P, "Genus", NArm = TRUE)
taxgen <- get.tax(ps_gen) %>%
  select(-Species)
for (i in 1:6){taxgen[,i] <- as.character(taxgen[,i])} 
taxgen[is.na(taxgen)] <- "unspecified" 
tax_table(ps_gen)<-set.tax(taxgen) 
topN <- 15 
most_abundant_gen <- as.data.frame(sort(taxa_sums(ps_gen), TRUE)[1:topN])
most_abundant_gen$otu <- row.names(most_abundant_gen)
composition <- get.otu.melt(ps_gen, filter.zero =F) %>%
  mutate(Genus = if_else(otu %in% most_abundant_gen$otu, Genus, "Others")) %>%
  mutate(Genus = fct_relevel(Genus, "Massilia", "Staphylococcus", "Streptococcus", "Ralstonia",
                             "Prevotella", "Pseudarthrobacter", "Porphyromonas", "Enterococcus",
                             "Klebsiella","Methylobacterium","Alloprevotella", "Atopobium",  
                             "Lactobacillus", "Peptostreptococcus","Bifidobacterium", "Others")) %>%
  mutate(Phylum = if_else(Phylum == "Firmicutes", "Firmicutes", 
                          if_else(Phylum == "Bacteroidetes", "Bacteroidetes", 
                                  if_else(Phylum == "Proteobacteria", "Proteobacteria", 
                                          if_else(Phylum == "Actinobacteria", "Actinobacteria", "Other"))))) %>%
  mutate(Phylum = if_else(Genus == "Others", "Other", Phylum)) %>%
  mutate(Phylum = fct_relevel(Phylum, "Proteobacteria", "Firmicutes", "Bacteroidetes", 
                              "Actinobacteria", "Other")) %>%
  mutate(vfd60 = fct_relevel(vfd60, "Extubated", "Deceased or Intubated"))

comp_vfd <- composition %>%
  group_by(Phylum, Genus, vfd60) %>%
  dplyr::summarize(mean = mean(pctseqs), sd = sd(pctseqs))
getPalette <- c("#849e1b",  "#9e351b", "#1b849e", "#bea93b", "#808080")

ggplot(comp_vfd, aes(x=Genus, y=mean, fill=Phylum))+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(~vfd60) +
  scale_fill_manual(values = getPalette, na.value="grey50")+
  theme_bw(base_size=15)+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Mean relative abundance")+ 
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")

p.genus <- composition %>%
  group_by(Genus) %>%
  summarize(p = wilcox.test(pctseqs ~ vfd60)$p.value)
p.genus$BH = p.adjust(p.genus$p, method = "BH")

rm(comp_vfd, composition, df.fam, p1, most_abundant_gen, p.genus, ps_gen, tax, taxgen, topN, getPalette)
```

In addition, a DESeq2 model and random forest analyses were used to identify differentially abundant taxa. 
```
ps.deseq <- tax_glom(P, "Genus")
ps.deseq <- core(ps.deseq, detection = 1, prevalence = 10/100, include.lowest = T) # prevalence of 10%

gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

dsq <- phyloseq_to_deseq2(ps.deseq,~vfd60 ) 
geoMeans <- apply(counts(dsq), 1, gm_mean)
dsq <- estimateSizeFactors(dsq, geoMeans = geoMeans) 
dsq <- DESeq(dsq,  fitType="local")    

res <- results(dsq, cooksCutoff = FALSE, pAdjustMethod = "BH" )
deseq <- res[which(res$padj < 0.05), ]  #adjusted p-value <0.05
deseq <- cbind(as(deseq, "data.frame"), as(tax_table(ps.deseq)[rownames(deseq), ], "matrix"))

deseq <- deseq %>% 
  select(Genus,log2FoldChange) %>%
  group_by(Genus) %>%
  summarise_at(c("log2FoldChange"), sum, na.rm=T)%>%
  mutate(group=ifelse(log2FoldChange<0, "Deceased or Intubated", "Extubated")) 

ggplot(deseq, aes(x=reorder(Genus,log2FoldChange), y=log2FoldChange, fill=group), 
              stat="identity", color= "black")+
  geom_bar(stat = "identity") + 
  coord_flip() +
  ylab("log 2-Fold Change")+
  xlab("") +
  ggtitle("Differentially abundant Genera for Extubated vs. Deceased/Intubated") +
  theme_bw() +
  scale_fill_manual(values = c("#752936","#44aa99"))+
  theme(legend.title = element_blank(),
        legend.position = "none") 
        
rm(deseq, dsq, geoMeans, res, ps.deseq)
```
```
ps.rf <- P
ps.rf <- tax_glom(ps.rf, "Genus") 
ps.rf <- microbiome::transform(ps.rf, "compositional")
predictors <- (otu_table(ps.rf)) 
response <- as.factor(sample_data(ps.rf)$vfd60) 
rf.data <- data.frame(response, predictors)
tax <- get.tax(ps.rf)

set.seed(88)
classify <- randomForest(response~., data = rf.data, ntree = 501, na.action=na.exclude,
                         importance = TRUE, proximities = TRUE)

imp <- importance(classify) 
imp <- data.frame(predictors = rownames(imp), imp)
imp$predictors <- factor(imp$predictors, levels = imp$predictors) 
imp <- imp %>%
  mutate(otu = predictors) %>%
  left_join(tax) %>%
  select(MeanDecreaseAccuracy, Family, Genus, otu) %>%
  group_by(Genus)%>%
  summarize(MeanDecreaseAccuracy = sum(MeanDecreaseAccuracy))

imp <- arrange(imp, desc(MeanDecreaseAccuracy)) 
imp.10 <- head(imp, n=10)
ggplot(imp.10, aes(x = reorder(Genus, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#1b849e") +
  coord_flip() +
  ggtitle("Most discriminating Genera for Extubated vs. Deceased/Intubated") +
  theme_bw()

rm(ps.rf, predictors, response, rf.data, tax, classify, imp, imp.10)
```
```
ps.rf <- P
ps.rf <- tax_glom(ps.rf, "Genus")
ps.rf <- microbiome::transform(ps.rf, "compositional")

# Prepara dataframes for analysis
predictors <- (otu_table(ps.rf)) 
response <- as.numeric(sample_data(ps.rf)$ventfreedays60_intub)
rf.data <- data.frame(response, predictors)
tax <- get.tax(ps.rf)

set.seed(88)
classify <- randomForest(response~., 
                         data = rf.data, ntree = 501,na.action=na.exclude,
                         importance = TRUE, proximities = TRUE)

imp <- classify$importance
imp <- data.frame(predictors = rownames(imp), imp) 
imp$predictors <- factor(imp$predictors, levels = imp$predictors) 
imp.10 <- imp %>%
  mutate(otu = predictors) %>%
  left_join(tax) %>%
  select(X.IncMSE, Family, Genus) %>%
  group_by(Genus)%>%
  summarize(X.IncMSE = sum(X.IncMSE))
imp.10 <- arrange(imp.10, desc(X.IncMSE)) 
imp.10 <- head(imp.10, n=10)

ggplot(imp.10, aes(x = reorder(Genus, X.IncMSE), y = X.IncMSE, fill = Genus)) +
  geom_bar(stat = "identity", fill = "#1b849e") +
  coord_flip() +
  ggtitle("Most discriminating Genera for number of ventilator-free days") +
  theme_bw()

rm(ps.rf, predictors, response, rf.data, tax, classify, imp, imp.10)
```

Outcome assessment at 60 and 90 days following sample collection yielded similar associations between clinical outcomes and lung microbiota community composition. 

```
m <- m %>%
  mutate(vfd60_bal = if_else(ventfreedays60_bal == 0, "Deceased or Intubated", "Extubated")) %>%
  mutate(vfd90_bal = if_else(ventfreedays90_bal == 0, "Deceased or Intubated", "Extubated"))
sample_data(P) <- set.samp(m)
p1 <- microbiome::transform(P, "compositional")

set.seed(88)
wunifrac <- phyloseq::distance(p1, method = "wunifrac") 
adonis2((wunifrac ~ vfd60_bal), by = 'margin', data = m, permutations =9999) 
ord.wuni <- ordinate(p1, method = 'PCoA', distance = wunifrac) 

wuni_df <- as.data.frame(as.matrix(ord.wuni$vectors)) 
wuni_df$sample <- row.names(wuni_df) 
wuni_df <- wuni_df %>%
  left_join(m) %>%
  select(Axis.1, Axis.2, sample, vfd60_bal)
centroids <- aggregate(cbind(Axis.1,Axis.2)~ vfd60_bal, data= wuni_df, mean) 
wuni_df <- merge(wuni_df,centroids,by="vfd60_bal",suffixes=c("",".centroid")) 

ggplot(wuni_df, aes(Axis.1, Axis.2, color = vfd60_bal)) +
  geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, color= vfd60_bal), alpha = 0.3)+ 
  geom_point(data=wuni_df,aes(color=vfd60_bal),size=5,alpha=1.0) + 
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, label=c("Deceased or Intubated", "Extubated")), size=10, fill = ggplot2::alpha(c("white"),0.76)) +
  scale_colour_manual(values = c("#752936","#44aa99"))+
  theme_bw() +
  xlab("Axis 1") +
  ylab("Axis 2") + 
  theme(legend.position = "none")
rm(wuni_df, ord.wuni, centroids)
```
```
adonis2((wunifrac ~ vfd90_bal), by = 'margin', data = m, permutations =9999)
ord.wuni <- ordinate(p1, method = 'PCoA', distance = wunifrac)
wuni_df <- as.data.frame(as.matrix(ord.wuni$vectors))
wuni_df$sample <- row.names(wuni_df)
wuni_df <- wuni_df %>%
  left_join(m) %>%
  select(Axis.1, Axis.2, sample, vfd90_bal)
centroids <- aggregate(cbind(Axis.1,Axis.2)~ vfd90_bal, data= wuni_df, mean) 
wuni_df <- merge(wuni_df,centroids,by="vfd90_bal",suffixes=c("",".centroid"))

ggplot(wuni_df, aes(Axis.1, Axis.2, color = vfd90_bal)) +
  geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, color= vfd90_bal), alpha = 0.3)+ 
  geom_point(data=wuni_df,aes(color=vfd90_bal),size=5,alpha=1.0) + 
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, label=c("Deceased or Intubated", "Extubated")), size=10, fill = ggplot2::alpha(c("white"),0.76)) +
  scale_colour_manual(values = c("#752936","#44aa99"))+
  theme_bw() +
  xlab("Axis 1") + 
  ylab("Axis 2") + 
  theme(legend.position = "none")
rm(wuni_df, ord.wuni, centroids, wunifrac)
```



## Step 6 - Lung microbiota alterations correspond with secondary pulmonary infections
First, we aimed to verify if secondary bacterial and fungal infections were associated with a higher bacterial and fungal burden

```
cultures <- metadata %>%
  mutate(bal_bact_pathogen_other = as.character(bal_bact_pathogen_other)) %>%
  
  # Staphylococcus
  mutate(bal_bact_staph = if_else(bal_bact_pathogen == "5", "Staphylococcus spp", 
                                  if_else(bal_bact_second_pathogen == "5", "Staphylococcus spp", "n"))) %>%
  # Klebsiella
  mutate(bal_bact_kleb = if_else(bal_bact_pathogen == "3", "Klebsiella spp", 
                                 if_else(bal_bact_second_pathogen == "3", "Klebsiella spp", "n"))) %>%
  mutate(bal_bact_kleb = if_else(bal_bact_pathogen_other %like% "Klebsiella", "Klebsiella spp", bal_bact_kleb)) %>%
  mutate(bal_bact_kleb = if_else(bal_bact_pathogen_other %like% "klebsiella", "Klebsiella spp", bal_bact_kleb)) %>%
  # Serratia
  mutate(bal_bact_serratia = if_else(bal_bact_pathogen_other == "Serratia marcescens", "Serratia marcescens", 
                                     if_else(bal_bact_second_pathogen_other == "Serratia marcescens", "Serratia marcescens", "n"),  # if something else
                                     if_else(bal_bact_second_pathogen_other == "Serratia marcescens", "Serratia marcescens", "n"))) %>%  # if NA
  mutate(bal_bact_pathogen = as.character(bal_bact_pathogen)) %>%
  mutate(bal_bact_second_pathogen = as.character(bal_bact_second_pathogen)) %>%
  mutate(bal_bact_pathogen = if_else(bal_bact_pathogen_other == "Serratia marcescens", "7", bal_bact_pathogen, bal_bact_pathogen)) %>%
  mutate(bal_bact_second_pathogen = if_else(bal_bact_second_pathogen_other == "Serratia marcescens", "7", bal_bact_second_pathogen, bal_bact_second_pathogen)) %>%
  mutate(bal_bact = as.character(bal_bact)) %>%
  mutate(bal_bact = if_else(bal_bact_pathogen_other == "gisten", "0", bal_bact, bal_bact)) %>%   # Remove 'gisten' (=yeasts) as a positive bacterial culture
  select(sample, sample_id, patient_id, bal_nr, 
         bal_bact, bal_bact_pathogen, bal_bact_pathogen_other,
         bal_bact_second, bal_bact_second_pathogen, bal_bact_second_pathogen_other, 
         bal_bact_staph, bal_bact_kleb, bal_bact_serratia,
         bal_aspergillosis, `16s`, `18s`)

comparisons <- list(c("0", "1"), c("0","2"), c("1", "2"))
colours <- c("#849e1b", "#1b849e", "#9e351b")

cultures %>%
  filter(!is.na(bal_bact)) %>%
  ggplot(aes(x=bal_bact, y=`16s`, fill=bal_bact))+ 
  geom_boxplot(outlier.shape = NA, alpha = .20) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+ 
  theme_cowplot(12)+  
  scale_fill_manual(values = colours) + 
  scale_color_manual(values = colours) +
  xlab("")+ # x-axis text
  ylab("Bacterial DNA Burden")+ 
  scale_y_continuous(trans=log_epsilon_trans(epsilon=10),limits=c(0,2200))+  
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+ 
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(),
        legend.title = element_blank())+
  labs(fill = "")

comparisons <- list(c("0", "1"))
colours <- c("#849e1b", "#9e351b")
cultures %>%
  filter(!is.na(bal_aspergillosis)) %>%
  ggplot(aes(x=bal_aspergillosis, y=`18s`, fill=bal_aspergillosis))+ 
  geom_boxplot(outlier.shape = NA, alpha = .20) + 
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+ 
  theme_cowplot(12)+  
  scale_fill_manual(values = colours) + 
  scale_color_manual(values = colours) + 
  xlab("")+ 
  ylab("Fungal DNA Burden")+ 
  scale_y_continuous(trans=log_epsilon_trans(epsilon=10),limits=c(0,3000))+  
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+ 
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(),
        legend.title = element_blank())+
  labs(fill = "")
```

Next, we asked whether bacterial ventilator-associated pneumonia was associated with elevated relative abundances of the causative pathogens in BAL fluid. 
```
cultures <- cultures %>%
  filter(bal_bact == "2")
psculture <- ps
sample_data(psculture) <- set.samp(cultures)
psculture <- microbiome::transform(psculture, "compositional")
psculture <- phyloseq::tax_glom(psculture, "Genus", NArm = TRUE) 
cult_gen <- get.otu.melt(psculture, filter.zero =F) 

# Staphylococcus
staph <- cult_gen %>%
  mutate(culture.staph = if_else(bal_bact_staph == "Staphylococcus spp", "Positive", "Negative", "Negative")) %>%
  filter(Genus == "Staphylococcus") %>%
  group_by(sample_id, culture.staph, bal_bact_second, Genus) %>%
  dplyr::summarize(sum = sum(pctseqs))

comparisons <- list(c("Negative", "Positive"))
colours <- c("#849e1b", "#9e351b")

ggplot(staph, aes(x=culture.staph, y=sum, fill=culture.staph))+ 
  geom_boxplot(outlier.shape = NA, alpha = .20) + 
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+ 
  theme_cowplot(12)+  
  scale_fill_manual(values = colours) +
  scale_color_manual(values = colours) + 
  xlab("")+ # x-axis text
  ylab("Relative abundance")+ 
  ggtitle("Staphylococcus")+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+ 
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(),
        legend.title = element_blank())+
  labs(fill = "")

# Klebsiella
kleb <- cult_gen %>%
  mutate(culture.kleb = if_else(bal_bact_kleb == "Klebsiella spp", "Positive", "Negative", "Negative")) %>%
  filter(Genus == "Klebsiella") %>%
  group_by(sample_id, culture.kleb, bal_bact_second, Genus) %>%
  dplyr::summarize(sum = sum(pctseqs))

ggplot(kleb, aes(x=culture.kleb, y=sum, fill=culture.kleb))+ 
  geom_boxplot(outlier.shape = NA, alpha = .20) + 
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+ 
  theme_cowplot(12)+  
  scale_fill_manual(values = colours) + 
  scale_color_manual(values = colours) +
  ggtitle("Klebsiella")+
  xlab("")+ # x-axis text
  ylab("Relative abundance")+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+ 
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        legend.title = element_blank())+
  labs(fill = "")

# Serratia
serr <- cult_gen %>%
  mutate(culture.serr = if_else(bal_bact_serratia == "Serratia marcescens", "Positive", "Negative", "Negative")) %>%
  filter(Genus == "Serratia") %>%
  group_by(sample_id, culture.serr, bal_bact_second, Genus) %>%
  dplyr::summarize(sum = sum(pctseqs))

ggplot(serr, aes(x=culture.serr, y=sum, fill=culture.serr))+ 
  geom_boxplot(outlier.shape = NA, alpha = .20) + 
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.15)+ 
  theme_cowplot(12)+  
  scale_fill_manual(values = colours) + 
  scale_color_manual(values = colours) + 
  ggtitle("Serratia")+
  xlab("")+ 
  ylab("Relative abundance")+ 
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", size=6)+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        legend.title = element_blank())+
  labs(fill = "")
rm(cultures, cult_gen, comparisons, colours, staph, kleb, serr, psculture)
```


## Step 7 - Lung bacterial and fungal burden are correlated with alveolar inflammatory cytokine responses
We evaluated if microbial burdens were correlated with alveolar cytokine responses, as measured in both initial and follow-up samples. 

```
cytokinesbal <- read_excel("Documents/PhD/Covid-19/Data Github/cytokines-covid-ards.xlsx") %>% # Cytokine data
  mutate(sample_id = as.character(sample_id))
cytokinesbal <- left_join(metadata, cytokinesbal)

cytokinesbal <- cytokinesbal %>%
  mutate(`18s_pseudocount` = `18s`+0.009018145) # add pseudocount (50% of lowest value) to enable log transformations
  
library(formulaic)
cytokinenames <- colnames(cytokinesbal)[42:104]
cytokinenames <- add.backtick(cytokinenames, include.backtick = "as.needed")
cytokines <- as.matrix(cytokinesbal[,42:104])
```
```
# Bacterial
res.bb.bal = data.frame(variable = "bacterialburden", 
                        cytokinename=cytokinenames, spearmanrho=0, pvalue=2, padj=2, stringsAsFactors = F)
for(i in 1:63){
  ct = cor.test(cytokinesbal$`16s`, cytokines[,i], method="spearman")
  res.bb.bal[i, "cytokinename"] = cytokinenames[i]
  res.bb.bal[i, "pvalue"] = ct$p.value
  res.bb.bal[i, "spearmanrho"] = ct$estimate
}

#cor.test((cytokinesbal$`16s`), cytokinesbal$TNFa, method="spearman") #check to see if values are correct

res.bb.bal$padj = p.adjust(res.bb.bal$pvalue, method = "BH")
res.bb.bal <- res.bb.bal%>%
  mutate(Sign = if_else(spearmanrho > 0, "Positive", "Negative")) %>%
  mutate(Sign = if_else(padj > 0.05, "Neutral", Sign))
table(res.bb.bal$Sign)

res.bb.bal %>%
  ggplot(aes(x = reorder(cytokinename, -spearmanrho), y = spearmanrho, fill = 'blue')) +
  geom_bar(stat = "identity", fill = "#1b849e") +
  geom_text(aes(y = -.25, x = reorder(cytokinename, -spearmanrho), label = round(padj,3)))+
  coord_flip() +
  ggtitle("Correlation between BAL biomarkers and bacterial burden") +
  xlab("")+
  theme_bw()
```

```
# Fungal
res.fb.bal = data.frame(variable = "fungalburden", 
                        cytokinename=cytokinenames, spearmanrho=0, pvalue=2, padj=2, stringsAsFactors = F)

for(i in 1:63){
  ct = cor.test(cytokinesbal$`18s`, cytokines[,i], method="spearman")
  res.fb.bal[i, "cytokinename"] = cytokinenames[i]
  res.fb.bal[i, "pvalue"] = ct$p.value
  res.fb.bal[i, "spearmanrho"] = ct$estimate
}

cor.test((cytokinesbal$`18s`), cytokinesbal$TNFa, method="spearman") #check to see if values are correct

res.fb.bal$padj = p.adjust(res.fb.bal$pvalue, method = "BH")
res.fb.bal <- res.fb.bal%>%
  mutate(Sign = if_else(spearmanrho > 0, "Positive", "Negative")) %>%
  mutate(Sign = if_else(padj > 0.05, "Neutral", Sign))

table(res.fb.bal$Sign)

res.fb.bal %>%
  ggplot(aes(x = reorder(cytokinename, -spearmanrho), y = spearmanrho, fill = 'blue')) +
  geom_bar(stat = "identity", fill = "#1b849e") +
  geom_text(aes(y = -.31, x = reorder(cytokinename, -spearmanrho), label = round(padj,3)))+
  coord_flip() +
  ggtitle("Correlation between BAL biomarkers and fungal burden") +
  xlab("")+
  theme_bw()
```

Next, we focused on cytokines that are commonly elevated in ARDS or COVID-19, and serve as markers of (alveolar) inflammation in these conditions
```
cytokinesbal %>%
  ggplot(aes(x = `16s`, y = TNFa)) +
  geom_point(size=2.5, colour = "#000000") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  theme_bw()+
  scale_x_continuous(trans=log_epsilon_trans(epsilon=5),limits=c(0,400))+  
  scale_y_continuous(trans=log_epsilon_trans(epsilon=10),limits=c(0,3500))+  
  theme(legend.position = "none") +
  ylab("TNF-alpha")+
  xlab("Bacterial burden")
cor.test((cytokinesbal$`16s`), cytokinesbal$TNFa, method="spearman") 

cytokinesbal %>%
  ggplot(aes(x = `16s`, y = `IL-6`)) +
  geom_point(size=2.5, colour = "#000000") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  theme_bw()+
  scale_x_continuous(trans=log_epsilon_trans(epsilon=5),limits=c(0,400))+  
  scale_y_continuous(trans=log_epsilon_trans(epsilon=100),limits=c(0,47000))+  
  theme(legend.position = "none") +
  ylab("IL-6")+
  xlab("Bacterial burden")
cor.test((cytokinesbal$`16s`), cytokinesbal$`IL-6`, method="spearman") 

cytokinesbal %>%
  ggplot(aes(x = `16s`, y = `IL-1b`)) +
  geom_point(size=2.5, colour = "#000000") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  theme_bw()+
  scale_x_continuous(trans=log_epsilon_trans(epsilon=5),limits=c(0,400))+  
  scale_y_continuous(trans=log_epsilon_trans(epsilon=100),limits=c(0,16000))+  
  theme(legend.position = "none") +
  ylab("IL-1 beta")+
  xlab("Bacterial burden")
cor.test((cytokinesbal$`16s`), cytokinesbal$`IL-1b`, method="spearman") 

cytokinesbal %>%
  ggplot(aes(x = `16s`, y = `TGF-a`)) +
  geom_point(size=2.5, colour = "#000000") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  theme_bw()+
  scale_x_continuous(trans=log_epsilon_trans(epsilon=10),limits=c(0,400))+  
  scale_y_continuous(trans=log_epsilon_trans(epsilon=5),limits=c(0,200))+  
  theme(legend.position = "none") +
  ylab("TGF alpha")+
  xlab("Bacterial burden")
cor.test((cytokinesbal$`16s`), cytokinesbal$`TGF-a`, method="spearman") 
```
```
cytokinesbal %>%
  ggplot(aes(x = `18s_pseudocount`, y = TNFa)) +
  geom_point(size=2.5, colour = "#000000") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  theme_bw()+
  scale_x_continuous(trans=log_epsilon_trans(epsilon=5),limits=c(0,1800))+  
  scale_y_continuous(trans=log_epsilon_trans(epsilon=10),limits=c(0,3500))+  
  theme(legend.position = "none") +
  ylab("TNF-alpha")+
  xlab("Fungal burden")
cor.test((cytokinesbal$`18s`), cytokinesbal$TNFa, method="spearman") 

cytokinesbal %>%
  ggplot(aes(x = `18s_pseudocount`, y = `TGF-a`)) +
  geom_point(size=2.5, colour = "#000000") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  theme_bw()+
  scale_x_continuous(trans=log_epsilon_trans(epsilon=10),limits=c(0,1800))+  
  scale_y_continuous(trans=log_epsilon_trans(epsilon=5),limits=c(0,200))+  
  theme(legend.position = "none") +
  ylab("TGF alpha")+
  xlab("Fungal burden")
cor.test((cytokinesbal$`18s`), cytokinesbal$`TGF-a`, method="spearman") 

cytokinesbal %>%
  ggplot(aes(x = `18s_pseudocount`, y = `IL-12p70`)) +
  geom_point(size=2.5, colour = "#000000") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  theme_bw()+
  scale_x_continuous(trans=log_epsilon_trans(epsilon=5),limits=c(0,1800))+  
  scale_y_continuous(trans=log_epsilon_trans(epsilon=2),limits=c(0,350))+  
  theme(legend.position = "none") +
  ylab("IL-12 p70")+
  xlab("Fungal burden")
cor.test((cytokinesbal$`18s`), cytokinesbal$`IL-12p70`, method="spearman") 

cytokinesbal %>%
  ggplot(aes(x = `18s_pseudocount`, y = `IL-17A`)) +
  geom_point(size=2.5, colour = "#000000") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  theme_bw()+
  scale_x_continuous(trans=log_epsilon_trans(epsilon=5),limits=c(0,1800))+  
  scale_y_continuous(trans=log_epsilon_trans(epsilon=20),limits=c(0,90))+  
  theme(legend.position = "none") +
  ylab("IL-17A")+
  xlab("Fungal burden")
cor.test((cytokinesbal$`18s`), cytokinesbal$`IL-17A`, method="spearman") 
```

## SessionInfo
```
sessionInfo()
```
```
## R version 4.1.1 (2021-08-10)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Monterey 12.2.1
## 
## Matrix products: default
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] nl_NL.UTF-8/nl_NL.UTF-8/nl_NL.UTF-8/C/nl_NL.UTF-8/nl_NL.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] formulaic_0.0.8             openxlsx_4.2.5              splitstackshape_1.4.8      
##  [4] data.table_1.14.2           ggpmisc_0.4.5               ggpp_0.4.3                 
## [7] ggrepel_0.9.1               scales_1.1.1                cmprsk_2.2-11              
## [10] survminer_0.4.9             survival_3.3-1              readxl_1.3.1               
## [13] randomForest_4.7-1          DESeq2_1.32.0               SummarizedExperiment_1.22.0
## [16] Biobase_2.52.0              MatrixGenerics_1.4.3        matrixStats_0.61.0         
## [19] GenomicRanges_1.44.0        GenomeInfoDb_1.28.4         IRanges_2.26.0             
## [22] S4Vectors_0.30.2            BiocGenerics_0.38.0         ggpubr_0.4.0               
## [25] cowplot_1.1.1               RColorBrewer_1.1-2          microbiome_1.14.0          
## [28] vegan_2.5-7                 lattice_0.20-45             permute_0.9-7              
## [31] yingtools2_0.0.1.60         forcats_0.5.1               stringr_1.4.0              
## [34] dplyr_1.0.8                 purrr_0.3.4                 readr_2.1.2                
## [37] tidyr_1.2.0                 tibble_3.1.6                ggplot2_3.3.5              
## [40] tidyverse_1.3.1             phyloseq_1.36.0            
```
