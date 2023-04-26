# load packages
library(qiime2R)
library(microbiome)
library(dplyr)
library(ggplot2)
library(vegan)
library(nlme)

# set global parameters
theme_set(theme_classic())
sc <- scale_color_manual(values = c("AR" = "#C34554", 
                              "AR-FMT" = "#55B35E", 
                              "PR" = "#81B1F7")) 

# functions
# calculates values for stats bars for plots
data.summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m, ymin=ymin, ymax=ymax))
}

# download files into phyloseq object

setwd("~/Kohl Lab/Whooping Cranes/QIIME processing/qza_files")

# fecal_qiime_metadata <- read.delim("fecal_qiime_metadata.txt")
# fecal_qiime_metadata$day <- factor(fecal_qiime_metadata$day)
# fecal_qiime_metadata$year <- factor(fecal_qiime_metadata$year)

# fecal_qiime_metadata <- read.delim("fecal_qiime_metadata.txt") %>% 
#   mutate(day = as.factor(day)) %>% 
#   mutate(year = as.factor(year))

WC <-qza_to_phyloseq(
  features="table-pfilt-no-mitochloro.qza",
  tree="rooted-tree.qza",
  taxonomy="taxonomy.qza",
  metadata = "fecal_qiime_metadata.txt"
)

# add factor variables
sample_data(WC)$day_factor <- as.factor(sample_data(WC)$day)
sample_data(WC)$year_factor <- as.factor(sample_data(WC)$year)

### DATA EXPLORATION AND CLEANUP

# data summaries
summarize_phyloseq(WC)
sample_sums(WC)

# rarefy
# one sample at 8113 the rest above 21285
# one blank is very high so need to remove blanks apart from rarefying (and decontam at some point)
WC.rarelow <- rarefy_even_depth(WC, sample.size=8113, rngseed=14, verbose=TRUE)
WC.rarehigh <- rarefy_even_depth(WC, sample.size=21285, rngseed=14, verbose=TRUE)

# remove blanks
WC.rarelow.noblank <- subset_samples(WC.rarelow, sample.type == "fecal")

# subset only 2019 and 2022 (where there are all three groups)
WC.1922 <- subset_samples(WC, year == "2019" | year == "2022")
WC.rare.1922 <- rarefy_even_depth(WC.1922, sample.size=8113, rngseed=14, verbose=TRUE)

# subset only day 14
WC.d14 <- subset_samples(WC, day == 14)
WC.rare.d14 <- rarefy_even_depth(WC.d14, sample.size=min(sample_sums(WC.d14)), rngseed=14, verbose=TRUE)

### ALPHA DIVERSITY ANALYSIS

# alpha diversity plots
plot_richness(WC.rarelow.noblank, x="group") # maybe a slight trend towards PR higher?
plot_richness(WC.rarelow.noblank, x="year", color = "group") # pretty similar
plot_richness(WC.rarelow.noblank, x="sex") # similar
plot_richness(WC.rarelow.noblank, x="day", color="group") # no obvious trends

# plots for only 2019 and 2022 data
sample_data(WC.rare.1922)$shannon <- microbiome::alpha(WC.rare.1922, index = "shannon")
plot_richness(WC.rare.1922, x="group") # maybe still PR a little higher
plot_richness(WC.rare.1922, x="year", color = "group") # FMT maybe tends to be on the low side within a given year
plot_richness(WC.rare.1922, x="day", color = "group")

# shannon plot over time by group
# for full dataset
alpha <- microbiome::alpha(WC.rarelow.noblank)

meta.rarelow.noblank <- meta(WC.rarelow.noblank) %>% 
  mutate(shannon = alpha$diversity_shannon) %>% 
  mutate(observed = alpha$observed)

####PLOTS FOR PRESENTATION

meta.rarelow.noblank %>% 
  filter(day != "21") %>% 
  ggplot(aes(x=day, y=shannon, group = group)) + sc +
  geom_point(aes(color = group), position = position_dodge(width = 2), shape = 1, size = 3) +
  stat_summary(fun.y=mean, geom = "line", size = 1, aes(color = group)) +
  stat_summary(fun.data = data.summary, geom = "errorbar", size = 1, width = 2, na.rm = TRUE, aes(color = group)) +
  stat_summary(fun.y=mean, geom = "point", size = 3, aes(color = group)) +
  theme_classic(base_size = 15) +
  ylab("shannon diversity")

# observed richness with all points
meta.rarelow.noblank %>% 
  filter(day != "21") %>% 
  ggplot(aes(x=day, y=observed, group = group)) + sc +
  geom_point(aes(color = group), position = position_dodge(width = 2), shape = 1, size = 3) +
  stat_summary(fun.y=mean, geom = "line", size = 1, aes(color = group)) +
  stat_summary(fun.data = data.summary, geom = "errorbar", size = 1, width = 2, na.rm = TRUE, aes(color = group)) +
  stat_summary(fun.y=mean, geom = "point", size = 3, aes(color = group)) +
  theme_classic(base_size = 15) +
  ylab("observed richness")
# observed richness but no AR-FMT
meta.rarelow.noblank %>% 
  filter(day != "21") %>% 
  ggplot(aes(x=day, y=observed, group = group, alpha = group)) + sc +
  scale_alpha_manual(values = c(1, 0, 1)) +
  geom_point(aes(color = group), position = position_dodge(width = 2), shape = 1, size = 3) +
  stat_summary(fun.y=mean, geom = "line", size = 1, aes(color = group)) +
  stat_summary(fun.data = data.summary, geom = "errorbar", size = 1, width = 2, na.rm = TRUE, aes(color = group)) +
  stat_summary(fun.y=mean, geom = "point", size = 3, aes(color = group)) +
  theme_classic(base_size = 15) +
  ylab("observed richness")
# observed richness but no AR-or AR
meta.rarelow.noblank %>% 
  filter(day != "21") %>% 
  ggplot(aes(x=day, y=observed, group = group, alpha = group)) + sc +
  scale_alpha_manual(values = c(0, 0, 1)) +
  geom_point(aes(color = group), position = position_dodge(width = 2), shape = 1, size = 3) +
  stat_summary(fun.y=mean, geom = "line", size = 1, aes(color = group)) +
  stat_summary(fun.data = data.summary, geom = "errorbar", size = 1, width = 2, na.rm = TRUE, aes(color = group)) +
  stat_summary(fun.y=mean, geom = "point", size = 3, aes(color = group)) +
  theme_classic(base_size = 15) +
  ylab("observed richness")


# lm for shannon
data_for_lm <- meta.rarelow.noblank %>% 
  filter(day!= "21") %>% 
  mutate(bird.name = C(factor(bird.name, sum))) %>% 
  mutate(group = C(factor(group, sum)))
model <- lme(shannon ~ group * day_factor, random =~1|bird.code, data = data_for_lm)
summary(model)

# repeat for 2019/2022
alpha.1922 <- microbiome::alpha(WC.rare.1922)
meta.rare.1922 <- meta(WC.rare.1922) %>% 
  mutate(shannon = alpha.1922$diversity_shannon)
meta.rare.1922 %>% 
  ggplot(aes(x=day, y=shannon, group = group)) +
  geom_point(aes(color = group)) +
  stat_summary(fun.y=mean, geom = "point", color = "black") +
  stat_summary(fun.y=mean, geom = "line", aes(color = group)) +
  stat_summary(fun.data = data.summary, geom = "pointrange", size = 0.3, na.rm = TRUE, aes(color = group))

#### BETA DIVERSITY ANALYSIS

# beta diversity without blanks or day 21
WC.rarelow.noblank.no21 <- subset_samples(WC.rarelow.noblank, day_factor != "21")
BrayDist.NoBlanks.No21 <- distance(WC.rarelow.noblank.no21, method = "bray")
BrayOrd.NoBlanks.No21 <- ordinate(WC.rarelow.noblank.no21, "NMDS", distance = BrayDist.NoBlanks.No21)
BrayPlot.NoBlanks.No21 <- plot_ordination(WC.rarelow.noblank.no21, BrayOrd.NoBlanks.No21, justDF = TRUE)

# by group
BrayPlot.NoBlanks.No21 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) + sc +
  geom_point(size=4, aes(color = group)) +
  labs(shape = "day")
# by group and day
BrayPlot.NoBlanks.No21 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) + sc +
  geom_point(size=4, aes(color = group, shape = day_factor)) +
  labs(shape = "day")
# by group but only PR and AR
BrayPlot.NoBlanks.No21 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, alpha = group)) + sc +
  scale_alpha_manual(values = c(1, 0, 1)) +
  geom_point(size=4, aes(color = group)) +
  labs(shape = "day")
# by group but only PR
BrayPlot.NoBlanks.No21 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, alpha = group)) + sc +
  scale_alpha_manual(values = c(0, 0, 1)) +
  geom_point(size=4, aes(color = group)) +
  labs(shape = "day")


meta.rarelow.noblank.No21 <- meta(WC.rarelow.noblank.no21)
adonis2(BrayDist.NoBlanks.No21 ~ group * day_factor, data = meta.rarelow.noblank.No21) # groups significantly different but replicates from each bird included

mod <- betadisper(BrayDist.NoBlanks.No21, sample_data(WC.rarelow.noblank.no21)$group, type = c("median","centroid"), bias.adjust = FALSE,
           sqrt.dist = FALSE, add = FALSE)
anova(mod)
(mod.HSD <- TukeyHSD(mod))
plot(mod, ellipse=TRUE, hull=FALSE)
boxplot(mod)

mod$distances %>% ggplot(aes(x=group, y=distances))

# by day
BrayPlot.NoBlanks %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(size=3, aes(color = day_factor))
adonis2(BrayDist.NoBlanks ~ day, data = meta.NoBlanks) # significant by day

# by bird
BrayPlot.NoBlanks %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(size=3, aes(color = bird.name))
adonis2(BrayDist.NoBlanks ~ bird.name, data = meta.NoBlanks) # significant

# beta diversity without blanks for only 2019 and 2022
BrayDist.NoBlanks.1922 <- distance(WC.rare.1922, method = "bray")
BrayOrd.NoBlanks.1922 <- ordinate(WC.rare.1922, "NMDS", distance = BrayDist.NoBlanks.1922)
BrayPlot.NoBlanks.1922 <- plot_ordination(WC.rare.1922, BrayOrd.NoBlanks.1922, justDF = TRUE)
meta.NoBlanks.1922 <- meta(WC.rare.1922)

# by group
BrayPlot.NoBlanks.1922 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(size=3, aes(color = group))
adonis2(BrayDist.NoBlanks.1922 ~ group, data = meta.NoBlanks.1922) # significant but much less than before

# by day
BrayPlot.NoBlanks.1922 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(size=3, aes(color = day))
adonis2(BrayDist.NoBlanks.1922 ~ day, data = meta.NoBlanks.1922) # very significant

# by bird
BrayPlot.NoBlanks.1922 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(size=3, aes(color = bird.name))
adonis2(BrayDist.NoBlanks.1922 ~ bird.name, data = meta.NoBlanks.1922) # not significant

# beta diversity without blanks within a single day (day 14)
BrayDist.NoBlanks.d14 <- distance(WC.rare.d14, method = "bray")
BrayOrd.NoBlanks.d14 <- ordinate(WC.rare.d14, "NMDS", distance = BrayDist.NoBlanks.d14)
BrayPlot.NoBlanks.d14 <- plot_ordination(WC.rare.d14, BrayOrd.NoBlanks.d14, justDF = TRUE)
meta.NoBlanks.d14 <- meta(WC.rare.d14)

# by group
BrayPlot.NoBlanks.d14 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(size=3, aes(color = group))
adonis2(BrayDist.NoBlanks.d14 ~ group, data = meta.NoBlanks.d14) # significant

# by year
BrayPlot.NoBlanks.d14 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(size=3, aes(color = year))
adonis2(BrayDist.NoBlanks.d14 ~ year, data = meta.NoBlanks.d14) # significant

## CLEARLY, GROUP AND YEAR AND DAY ALL HAVE AN EFFECT SO THE BEST COURSE OF ACTION 
## MIGHT BE A MODEL TRYING TO DETECT THE EFFECT OF GROUP WITH YEAR AS A RAND VAR AND DAY AS EITHER?
## SUBJECT AS A RAND


## Things to do:
# See if I can run decontam and figure out what was in the blanks
  # requires making new Sample_or_Control TF column
  # https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
# Factor the year and day variables
# Run a model?? Look into using Maaslin
