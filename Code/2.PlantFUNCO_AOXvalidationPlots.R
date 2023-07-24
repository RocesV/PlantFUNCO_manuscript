###################################
#    TAIOTACTICAL Mutant Assays   #
###################################

### Load libs

library(ggpubr)
library(car)

### AOX assay

## Root Phenotype Analyses

AOX <- read.delim("clipboard", header = T)
shapiro.test(AOX$Root_Length)
shapiro.test(AOX$Hypocotil_Length)
shapiro.test(AOX$Ratio)
AOX_C <- AOX[which(AOX$Condition == "Mock"),]
AOX_S <- AOX[which(AOX$Condition == "AA"),]

shapiro.test(AOX_C$Root_Length)
shapiro.test(AOX_C$Hypocotil_Length)
shapiro.test(AOX_C$Ratio)

shapiro.test(AOX_S$Root_Length)
shapiro.test(AOX_S$Hypocotil_Length)
shapiro.test(AOX_S$Ratio)

leveneTest(AOX_C$Root_Length ~ AOX_C$Genotype)
leveneTest(AOX_S$Root_Length ~ AOX_S$Genotype)

model<-aov(Root_Length~Genotype, data=AOX_C)
out <- HSD.test(model, "Genotype", group = TRUE, console = TRUE)

model<-aov(Root_Length~Genotype, data=AOX_S)
out <- HSD.test(model, "Genotype", group = TRUE, console = TRUE)


leveneTest(AOX_C$Hypocotil_Length ~ AOX_C$Genotype)
leveneTest(AOX_S$Hypocotil_Length ~ AOX_S$Genotype)

model<-aov(Hypocotil_Length~Genotype, data=AOX_C)
out <- HSD.test(model, "Genotype", group = TRUE, console = TRUE)

model<-aov(Hypocotil_Length~Genotype, data=AOX_S)
out <- HSD.test(model, "Genotype", group = TRUE, console = TRUE)


leveneTest(AOX_C$Ratio ~ AOX_C$Genotype)
leveneTest(AOX_S$Ratio ~ AOX_S$Genotype)

model<-aov(Ratio~Genotype, data=AOX_C)
out <- HSD.test(model, "Genotype", group = TRUE, console = TRUE)

model<-aov(Ratio~Genotype, data=AOX_S)
out <- HSD.test(model, "Genotype", group = TRUE, console = TRUE)

my_comparisons = list(c("Col-0", "AOX1A"),c("Col-0", "AOX1C"),c("Col-0", "AOX1D"))
AOX_C$Genotype <- factor(AOX_C$Genotype, levels = c("Col-0", "AOX1A", "AOX1C", "AOX1D"))
ggboxplot(AOX_C, x = "Genotype", y = "Root_Length",
          color = "Genotype", palette = "npg", add = "jitter", size = 1.5)+ 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", ref.group = "Col-0")+
  stat_compare_means(label.y = 5, method = "anova") + ggtitle("Mock Root Length")

AOX_S$Genotype <- factor(AOX_S$Genotype, levels = c("Col-0", "AOX1A", "AOX1C", "AOX1D"))
ggboxplot(AOX_S, x = "Genotype", y = "Root_Length",
          color = "Genotype", palette = "npg", add = "jitter", size = 1.5)+ 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", ref.group = "Col-0")+
  stat_compare_means(label.y = 5, method = "anova") + ggtitle("Antimycin A Root Length")

AOX$Genotype <- factor(AOX$Genotype, levels = c("Col-0", "AOX1A", "AOX1C", "AOX1D"))
ggboxplot(AOX, x = "Genotype", y = "Root_Length",
          color = "Genotype", palette = "npg", add = "jitter", size = 1.5, facet.by = "Condition")+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", ref.group = "Col-0")+
  stat_compare_means(label.y =2.2, method = "kruskal.test") + ggtitle("Root Length")


ggboxplot(AOX, x = "Genotype", y = "Hypocotil_Length",
          color = "Genotype", palette = "npg", add = "jitter", size = 1.5, facet.by = "Condition")+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", ref.group = "Col-0")+
  stat_compare_means(label.y = 0.4, method = "kruskal.test") + ggtitle("Hypocotyl Length")

ggboxplot(AOX, x = "Genotype", y = "Hypocotil_Length",
          color = "Genotype", palette = "npg", add = "jitter", size = 1.5, facet.by = "Condition")+ 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", ref.group = "Col-0")+
  stat_compare_means(label.y = 0.4, method = "anova") + ggtitle("Hypocotyl Length")

my_comparisons = list(c("Col-0", "AOX1C"),c("Col-0", "AOX1D"))
ggboxplot(AOX, x = "Genotype", y = "Ratio",
          color = "Genotype", palette = "npg", add = "jitter", size = 1.5, facet.by = "Condition")+ 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", ref.group = "Col-0")+
  stat_compare_means(label.y = 0.4, method = "anova") + ggtitle("Ratio Root/Hypocotyl")


## Cotyledon DAB (ROS) analyses

AOX.dab <- read.delim("clipboard", header = T)
AOX.dab$Genotype <- factor(AOX.dab$Genotype, levels = c("Col-0", "AOX1A", "AOX1C", "AOX1D"))

shapiro.test(AOX.dab$CTCF[which(AOX.dab$Condition == "Mock")])
leveneTest(AOX.dab$CTCF[which(AOX.dab$Condition == "PEGxHeat")] ~ AOX.dab$Genotype[which(AOX.dab$Condition == "PEGxHeat")])


my_comparisons = list(c("Col-0", "AOX1A"),c("Col-0", "AOX1C"),c("Col-0", "AOX1D"))
ggbarplot(AOX.dab, x = "Genotype", y = "CTCF",
          color = "Genotype", palette = "npg", facet.by = "Condition", add = "mean_se", position = position_dodge(0.8)) +
stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", ref.group = "Col-0")+
stat_compare_means(label.y = 8e+05, method = "kruskal.test") + ggtitle("DAB")

ggbarplot(AOX.dab, x = "Condition", y = "CTCF",
          color = "Condition", palette = "jama", facet.by = "Genotype", add = "mean_se", position = position_dodge(0.8)) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", ref.group = "Col-0")+
  stat_compare_means(label.y = 2e+05, method = "kruskal.test") + ggtitle("DAB")


