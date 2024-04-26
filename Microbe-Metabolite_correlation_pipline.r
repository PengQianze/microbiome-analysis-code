
#########################目水平微生物组与代谢物关联分析
#########################Association analysis of microbiome and metabolites at the order level
#所有代谢物变异系数 rda
#Coefficient of variation for all metabolites

library(vegan)
library(ggplot2)

taxonomy=read.csv("RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata
taxonomy[is.na(taxonomy)]<-0

idx = rownames(metadata) %in% colnames(taxonomy)
metadata = metadata[idx, , drop = F]
taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]

#RDA CCA判断

#judgment of RDA and CCA 
otu.helli=decostand(t(taxonomy),method = "hellinger")
decorana(otu.helli)
#Axis Lengths的第一轴的大小

#如果大于4.0,就应选CCA（基于单峰模型，典范对应分析）
#如果在3.0-4.0之间，选RDA和CCA均可
#如果小于3.0, RDA的结果会更合理（基于线性模型，冗余分析）

#analysisgroup3="CFA"
#analysisgroup4="LPA"

colorsF <-c(Pseudomonadales="#F4BD27",Enterobacteriales="#55C1D7",Bacteroidales="#7BCBDF",
Burkholderiales="#B3C84D",Xanthomonadales="#4DFAA1",Rhizobiales="#EEEAB4",Sphingomonadales="#FADB75",
Sphingobacteriales="#A3F4B7",Flavobacteriales="#BAA5E6",
Rhodospirillales="#BDD2C7",Other="#FFFFFF",
HCA="#007808",FLA="#00AC0B",LPA="#009A3E",CFA="#13ED6A",SPLT="#38DF8E",SPA="#1ADFA4")

corrmetabo<-c("LPA","HCA","FLA","CFA","SPLT","SPA")
sampFile = as.data.frame(metadata[, corrmetabo], row.names = row.names(metadata))

###RDA分析
###RDA Redundancy Analysis
otu_rda_all <- rda(t(taxonomy) ~ LPA+HCA+FLA+CFA+SPLT+SPA, data = sampFile)

####代谢物与微生物组相关性检验
####Correlation test between metabolites and microbiome
otu_rda_all_envfit <- envfit(otu_rda_all,sampFile,permu=999)
cor_metabo <- cbind(as.data.frame(otu_rda_all_envfit$vectors$r),as.data.frame(otu_rda_all_envfit$vectors$pvals))
colnames(cor_metabo) = c("r", "p")
cor_metabo <-rownames_to_column(cor_metabo, var = "metabo")

p14 <- ggplot(cor_metabo,aes(fill=metabo, y = r,x =reorder(metabo,-r)))+
geom_bar(color="#BAB6B6",stat = 'identity', width = 0.9,size=1)+
geom_text(aes(y = r+0.05, label =paste("p =",p, sep = ""),colour=metabo) ,size = 4, fontface = "bold") +
labs(x = '', y = '')+
xlab("Metabolite factor")+
ylab(expression(R^"2"))+
theme_classic()+
scale_color_manual(values=colorsF)+
scale_fill_manual(values=colorsF)
p14

ggsave("metabo_effection_plot.pdf", p14, width=6, height=5, units="in")
ggsave("metabo_effection_mhatonplot.jpg",p14, width=6, height=5, units="in")

####代谢组对微生物组变异的贡献度
####Contribution of metabolome to microbiome variation

variance_all = (otu_rda_all$CCA$tot.chi/otu_rda_all$tot.chi)
perm_otu_rda = anova.cca(otu_rda_all, permutations = 1000, parallel = 4)
p.val = perm_otu_rda[1, 4]

rda.scaling1 <- summary(otu_rda_all, scaling = 1)

eig1 = rda.scaling1$concont$importance

points = as.data.frame(rda.scaling1$sites)
tax<-c("Pseudomonadales","Xanthomonadales")
tax1<-t(taxonomy)
tax1<-as.data.frame(tax1[, tax], row.names = row.names(tax1))
sampFiletax =cbind(sampFile,tax1)
points = cbind(sampFiletax, points)

p10 = ggplot(points) + 
  geom_point(aes(x =  RDA1 , y = RDA2, fill=Pseudomonadales),shape=21,alpha = 0.8, size =4, stroke = 0.1,color = "#BAB6B6") + 
  labs(x = paste("RDA 1 (", format(100 * eig1[2,1], digits = 4), "%)", sep = ""), 
       y = paste("RDA 2 (", format(100 * eig1[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_all, digits = 3), " % of variance; P = ", format(p.val, digits = 2), sep = "")) +
  theme_classic() + 
  theme(text = element_text(family = "sans", size = 7))+
  scale_fill_gradient2(low="#B5CBE2",high="#F4BD27",mid="#FBFBFB",midpoint=40)
#  scale_colour_gradient2(low="#B5CBE2",high="#E8BE74",mid="#FBFBFB",midpoint=40)
p10

####添加解释变量
####Add explanatory variable

otu_rda_env = as.data.frame(rda.scaling1$biplot)

cor_metabo2<-cor_metabo[order(cor_metabo$r,decreasing = TRUE),]
cor_metabo2<-cor_metabo2[1:6,]
otu_rda_env2<-otu_rda_env[cor_metabo2$metabo,]

otu_rda_env2 <-rownames_to_column(otu_rda_env2, var = "metabo")

p11<-p10+geom_segment(data=otu_rda_env2,aes(x=0,y=0,xend=RDA1*4,yend=RDA2*4,colour=metabo),
,size=0.2,linetype=1,arrow=arrow(angle = 35,length=unit(0.3,"cm"))) +
geom_text(data=otu_rda_env2,aes(x=RDA1*4,y=RDA2*4,label=metabo,colour=metabo),size=3.5,
hjust=(1-sign(otu_rda_env2$RDA1))/2,angle=(180/pi)*atan(otu_rda_env2$RDA2/otu_rda_env2$RDA1))+
geom_vline(xintercept = 0,lty="dashed",color = 'black', size = 0.2)+
geom_hline(yintercept = 0,lty="dashed",color = 'black', size = 0.2)+
scale_color_manual(values=colorsF)+
theme_classic() 
p11

####响应代谢组的微生物变量，rda中添加最受影响的前3个响应变量
####Microbial variables in response to the metabolite, and the top three most affected response variables were added to RDA

rda<-rda(sampFile,t(taxonomy))
otu_cca_envfit <- envfit(rda,t(taxonomy),permu=999)

cor_com1 <- cbind(as.data.frame(otu_cca_envfit$vectors$r),as.data.frame(otu_cca_envfit$vectors$pvals))
colnames(cor_com1) = c("r", "p")
cor_com1 <-rownames_to_column(cor_com1, var = "microbiome")
cor_com2 <- cor_com1[cor_com1$p<0.05,]

cor_com3<-cor_com2[order(cor_com2$r,decreasing = TRUE),]
cor_com3<-cor_com3[1:10,]

cor_comTAX<-cor_com3[1:3,]
taxname<-cor_comTAX$microbiome
otu_rda_species = rda.scaling1$species
otu_rda_species = as.data.frame(otu_rda_species[taxname,])
otu_rda_species <-rownames_to_column(otu_rda_species, var = "microbiome")

p12<-p11+geom_segment(data=otu_rda_species,aes(x=0,y=0,xend=RDA1/12,yend=RDA2/12,colour=microbiome)
,size=0.2,linetype=1,arrow=arrow(angle = 35,length=unit(0.3,"cm"))) +
geom_text(data=otu_rda_species,aes(x=RDA1/12,y=RDA2/12,label=microbiome,colour=microbiome),size=3.5,
hjust=(1-sign(otu_rda_species$RDA1))/2,angle=(180/pi)*atan(otu_rda_species$RDA2/otu_rda_species$RDA1))+
theme_classic() 
#scale_color_manual(values=colorsF)
#scale_fill_manual(values=colorsO)
p12

ggsave("rda_plot.pdf", p12, width=8, height=5, units="in")
ggsave("rda_plot.jpg",p12, width=8, height=5, units="in")

####微生物组中各变量与代谢物组相关性检验
####Correlation test between variables in microbiome and metabolite group

p13 <- ggplot(cor_com3,aes(fill=microbiome, x = r,y =reorder(microbiome,r)))+
geom_bar(color= "#BAB6B6",stat = 'identity', width = 0.9,size=1)+
geom_text(aes(x = r+0.1, label =paste("p =",p, sep = ""),colour=microbiome) ,size = 5, fontface = "bold") +
labs(x = '', y = '')+
ylab("Microbiome factor")+
xlab(expression(R^"2"))+
theme_classic()+
scale_color_manual(values=colorsF)+
scale_fill_manual(values=colorsF)
p13

ggsave("microbiome_impact_plot.pdf", p13, width=7, height=5, units="in")
ggsave("microbiome_impact_plot.jpg",p13, width=7, height=5, units="in")


################################单一代谢物对微生物组的解释度
################################Variance of the microbiome explantied by a single metabolite
############################################################HCA
indiv_rda_all <- rda(t(taxonomy) ~ HCA, data = sampFile)

variance_indiv = (indiv_rda_all$CCA$tot.chi/indiv_rda_all$tot.chi)
perm_indiv_rda = anova.cca(indiv_rda_all, permutations = 1000, parallel = 4)
p.val2 = perm_indiv_rda[1, 4]

indiv_rda.scaling1 <- summary(indiv_rda_all, scaling = 1)

eig2 = indiv_rda.scaling1$cont$importance

points2 = as.data.frame(indiv_rda.scaling1$sites)
tax<-c("Pseudomonadales","Xanthomonadales")
tax1<-t(taxonomy)
tax1<-as.data.frame(tax1[, tax], row.names = row.names(tax1))
sampFiletax =cbind(sampFile,tax1)
points2 = cbind(sampFiletax, points2)

p15 = ggplot(points2) + 
  geom_point(aes(x =  RDA1 , y = PC1, fill=HCA),shape=21,alpha = 0.8, size =4, stroke = 0.1, color = "#BAB6B6") + 
  labs(x = paste("RDA 1 (", format(100 * eig2[2,1], digits = 4), "%)", sep = ""), 
       y = paste("PC 1 (", format(100 * eig2[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_indiv, digits = 3), " % of variance; P = ", format(p.val2 , digits = 2), sep = "")) +
  theme_classic() + 
  theme(text = element_text(family = "sans", size = 7))+
  scale_fill_gradient2(low="#B5CBE2",high="#007808",mid="#FBFBFB",midpoint=1.5)
p15

ggsave("HCA_effection_plot.pdf", p15, width=5, height=5, units="in")
ggsave("HCA_effection_plot.jpg",p15, width=5, height=5, units="in")

################################################################LPA
indiv_rda_all <- rda(t(taxonomy) ~ LPA, data = sampFile)

variance_indiv = (indiv_rda_all$CCA$tot.chi/indiv_rda_all$tot.chi)
perm_indiv_rda = anova.cca(indiv_rda_all, permutations = 1000, parallel = 4)
p.val2 = perm_indiv_rda[1, 4]

indiv_rda.scaling1 <- summary(indiv_rda_all, scaling = 1)

eig2 = indiv_rda.scaling1$cont$importance

points2 = as.data.frame(indiv_rda.scaling1$sites)
tax<-c("Pseudomonadales","Xanthomonadales")
tax1<-t(taxonomy)
tax1<-as.data.frame(tax1[, tax], row.names = row.names(tax1))
sampFiletax =cbind(sampFile,tax1)
points2 = cbind(sampFiletax, points2)

p16 = ggplot(points2) + 
  geom_point(aes(x =  RDA1 , y = PC1, fill=LPA),shape=21,alpha = 0.8, size =4, stroke = 0.1, color = "#BAB6B6") + 
  labs(x = paste("RDA 1 (", format(100 * eig2[2,1], digits = 4), "%)", sep = ""), 
       y = paste("PC 1 (", format(100 * eig2[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_indiv, digits = 3), " % of variance; P = ", format(p.val2 , digits = 2), sep = "")) +
  theme_classic() + 
  theme(text = element_text(family = "sans", size = 7))+
  scale_fill_gradient2(low="#B5CBE2",high="#009A3E",mid="#FBFBFB",midpoint=1.5)
#  scale_fill_gradient(low="#B5CBE2",high="#009A3E")
p16

ggsave("LPA_effection_plot.pdf", p16, width=5, height=5, units="in")
ggsave("LPA_effection_plot.jpg",p16, width=5, height=5, units="in")

###########################################################FLA
indiv_rda_all <- rda(t(taxonomy) ~ FLA, data = sampFile)

variance_indiv = (indiv_rda_all$CCA$tot.chi/indiv_rda_all$tot.chi)
perm_indiv_rda = anova.cca(indiv_rda_all, permutations = 1000, parallel = 4)
p.val2 = perm_indiv_rda[1, 4]

indiv_rda.scaling1 <- summary(indiv_rda_all, scaling = 1)

eig2 = indiv_rda.scaling1$cont$importance

points2 = as.data.frame(indiv_rda.scaling1$sites)
tax<-c("Pseudomonadales","Xanthomonadales")
tax1<-t(taxonomy)
tax1<-as.data.frame(tax1[, tax], row.names = row.names(tax1))
sampFiletax =cbind(sampFile,tax1)
points2 = cbind(sampFiletax, points2)

p17 = ggplot(points2) + 
  geom_point(aes(x =  RDA1 , y = PC1, fill=FLA),shape=21,alpha = 0.8, size =4, stroke = 0.1, color = "#BAB6B6") + 
  labs(x = paste("RDA 1 (", format(100 * eig2[2,1], digits = 4), "%)", sep = ""), 
       y = paste("PC 1 (", format(100 * eig2[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_indiv, digits = 3), " % of variance; P = ", format(p.val2 , digits = 2), sep = "")) +
  theme_classic() + 
  theme(text = element_text(family = "sans", size = 7))+
  scale_fill_gradient2(low="#B5CBE2",high="#00AC0B",mid="#FBFBFB",midpoint=0.3)
p17

ggsave("FLA_effection_plot.pdf", p17, width=5, height=5, units="in")
ggsave("FLA_effection_plot.jpg",p17, width=5, height=5, units="in")

###############################################################CFA
indiv_rda_all <- rda(t(taxonomy) ~ CFA, data = sampFile)

variance_indiv = (indiv_rda_all$CCA$tot.chi/indiv_rda_all$tot.chi)
perm_indiv_rda = anova.cca(indiv_rda_all, permutations = 1000, parallel = 4)
p.val2 = perm_indiv_rda[1, 4]

indiv_rda.scaling1 <- summary(indiv_rda_all, scaling = 1)

eig2 = indiv_rda.scaling1$cont$importance

points2 = as.data.frame(indiv_rda.scaling1$sites)
tax<-c("Pseudomonadales","Xanthomonadales")
tax1<-t(taxonomy)
tax1<-as.data.frame(tax1[, tax], row.names = row.names(tax1))
sampFiletax =cbind(sampFile,tax1)
points2 = cbind(sampFiletax, points2)

p18 = ggplot(points2) + 
  geom_point(aes(x =  RDA1 , y = PC1, fill=CFA),shape=21, size =4, stroke = 0.1,color = "#BAB6B6") + 
  labs(x = paste("RDA 1 (", format(100 * eig2[2,1], digits = 4), "%)", sep = ""), 
       y = paste("PC 1 (", format(100 * eig2[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_indiv, digits = 3), " % of variance; P = ", format(p.val2 , digits = 2), sep = "")) +
  theme_classic() + 
  theme(text = element_text(family = "sans", size = 7))+
  scale_fill_gradient2(low="#B5CBE2",high="#13ED6A",mid="#FBFBFB",midpoint=0.5)
p18

ggsave("CFA_effection_plot.pdf", p18, width=5, height=5, units="in")
ggsave("CFA_effection_plot.jpg",p18, width=5, height=5, units="in")

##################################################################SPLT
indiv_rda_all <- rda(t(taxonomy) ~ SPLT, data = sampFile)

variance_indiv = (indiv_rda_all$CCA$tot.chi/indiv_rda_all$tot.chi)
perm_indiv_rda = anova.cca(indiv_rda_all, permutations = 1000, parallel = 4)
p.val2 = perm_indiv_rda[1, 4]

indiv_rda.scaling1 <- summary(indiv_rda_all, scaling = 1)

eig2 = indiv_rda.scaling1$cont$importance

points2 = as.data.frame(indiv_rda.scaling1$sites)
tax<-c("Pseudomonadales","Xanthomonadales")
tax1<-t(taxonomy)
tax1<-as.data.frame(tax1[, tax], row.names = row.names(tax1))
sampFiletax =cbind(sampFile,tax1)
points2 = cbind(sampFiletax, points2)

p19 = ggplot(points2) + 
  geom_point(aes(x =  RDA1 , y = PC1, fill=SPLT),shape=21,alpha = 0.8,  size =4, stroke = 0.1,color = "#BAB6B6") + 
  labs(x = paste("RDA 1 (", format(100 * eig2[2,1], digits = 4), "%)", sep = ""), 
       y = paste("PC 1 (", format(100 * eig2[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_indiv, digits = 3), " % of variance; P = ", format(p.val2 , digits = 2), sep = "")) +
  theme_classic() + 
  theme(text = element_text(family = "sans", size = 7))+
  scale_fill_gradient2(low="#B5CBE2",high="#38DF8E",mid="#FBFBFB",midpoint=0.2)
p19

ggsave("SPLT_effection_plot.pdf", p19, width=5, height=5, units="in")
ggsave("SPLT_effection_plot.jpg",p19, width=5, height=5, units="in")

###############################################################SPA
indiv_rda_all <- rda(t(taxonomy) ~ SPA, data = sampFile)

variance_indiv = (indiv_rda_all$CCA$tot.chi/indiv_rda_all$tot.chi)
perm_indiv_rda = anova.cca(indiv_rda_all, permutations = 1000, parallel = 4)
p.val2 = perm_indiv_rda[1, 4]

indiv_rda.scaling1 <- summary(indiv_rda_all, scaling = 1)

eig2 = indiv_rda.scaling1$cont$importance

points2 = as.data.frame(indiv_rda.scaling1$sites)
tax<-c("Pseudomonadales","Xanthomonadales")
tax1<-t(taxonomy)
tax1<-as.data.frame(tax1[, tax], row.names = row.names(tax1))
sampFiletax =cbind(sampFile,tax1)
points2 = cbind(sampFiletax, points2)

p20 = ggplot(points2) + 
  geom_point(aes(x =  RDA1 , y = PC1, fill=SPA),shape=21,alpha = 0.8, size =4, stroke = 0.1, color = "#BAB6B6") + 
  labs(x = paste("RDA 1 (", format(100 * eig2[2,1], digits = 4), "%)", sep = ""), 
       y = paste("PC 1 (", format(100 * eig2[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_indiv, digits = 3), " % of variance; P = ", format(p.val2 , digits = 2), sep = "")) +
  theme_classic() + 
  theme(text = element_text(family = "sans", size = 7))+
  scale_fill_gradient2(low="#B5CBE2",high="#1ADFA4",mid="#FBFBFB",midpoint=0.17)
p20

ggsave("SPA_effection_plot.pdf", p20, width=5, height=5, units="in")
ggsave("SPA_effection_plot.jpg",p20, width=5, height=5, units="in")

#单一代谢物与相关物种相关性分析
#Correlation analysis between single metabolite and related species

library(ggplot2)

tax<-c("Pseudomonadales","Xanthomonadales","Burkholderiales","Sphingomonadales","Sphingobacteriales")
tax2<-t(taxonomy)
tax2<-as.data.frame(tax2[, tax], row.names = row.names(tax2))
sampFiletax2 =cbind(sampFile,tax2)

theme1<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=12),
        axis.line = element_line(color = "black",size = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        # axis.text.x = element_blank(),
        # axis.ticks.x=element_blank(),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
        # legend.position = "none",
       panel.spacing = unit(0,"lines")
       )

##########HCA_Pseudomonadales
p21 <-ggplot(sampFiletax2,aes(HCA,Pseudomonadales)) +
geom_point(aes(),pch=21,size =4, stroke = 0.1,color= "#BAB6B6",alpha=0.8,fill="#F4BD27")+
#scale_x_continuous(limits = c(1.3, 4.1), breaks = seq(1.3, 4.1, 0.5))+
geom_smooth(method = "lm", formula = NULL,linewidth=0.2,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 100,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+theme1
p21

ggsave("HCA_Pseudomonadales_plot.pdf", p21, width=6.0, height=8.2, units="in")
ggsave("HCA_Pseudomonadales_plot.jpg",p21, width=6.0, height=8.2,units="in")

##########FLA_Pseudomonadales
p22 <-ggplot(sampFiletax2,aes(FLA,Pseudomonadales)) +
geom_point(aes(),pch=21,size =4, stroke = 0.1,color= "#BAB6B6",alpha=0.8,fill="#F4BD27")+
#scale_x_continuous(limits = c(1.3, 4.1), breaks = seq(1.3, 4.1, 0.5))+
geom_smooth(method = "lm", formula = NULL,linewidth=0.2,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 100,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+theme1
p22

ggsave("FLA_Pseudomonadales_plot.pdf", p22,  width=5.9, height=8.2, units="in")
ggsave("FLA_Pseudomonadales_plot.jpg",p22,  width=5.9, height=8.2, units="in")

##########LPA_Pseudomonadales
p23 <-ggplot(sampFiletax2,aes(LPA,Pseudomonadales)) +
geom_point(aes(),pch=21,size =4, stroke = 0.1,color= "#BAB6B6",alpha=0.8,fill="#F4BD27")+
#scale_x_continuous(limits = c(1.3, 4.1), breaks = seq(1.3, 4.1, 0.5))+
geom_smooth(method = "lm", formula = NULL,linewidth=0.2,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 100,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+theme1
p23

ggsave("LPA_Pseudomonadales_plot.pdf", p23, width=5.9, height=8.2, units="in")
ggsave("LPA_Pseudomonadales_plot.jpg",p23, width=5.9, height=8.2, units="in")

##########CFA_Pseudomonadales
p24 <-ggplot(sampFiletax2,aes(CFA,Pseudomonadales)) +
geom_point(aes(),pch=21,size =4, stroke = 0.1,color= "#BAB6B6",alpha=0.8,fill="#F4BD27")+
#scale_x_continuous(limits = c(1.3, 4.1), breaks = seq(1.3, 4.1, 0.5))+
geom_smooth(method = "lm", formula = NULL,linewidth=0.2,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 100,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+theme1
p24

ggsave("CFA_Pseudomonadales_plot.pdf", p24,  width=5.9, height=8.2, units="in")
ggsave("CFA_Pseudomonadales_plot.jpg",p24,  width=5.9, height=8.2, units="in")

##########SPLT_Pseudomonadales
p25 <-ggplot(sampFiletax2,aes(SPLT,Pseudomonadales)) +
geom_point(aes(),pch=21,size =4, stroke = 0.1,color= "#BAB6B6",alpha=0.8,fill="#F4BD27")+
#scale_x_continuous(limits = c(1.3, 4.1), breaks = seq(1.3, 4.1, 0.5))+
geom_smooth(method = "lm", formula = NULL,linewidth=0.2,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 100,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+theme1
p25

ggsave("SPLT_Pseudomonadales_plot.pdf", p25,  width=5.9, height=8.2, units="in")
ggsave("SPLT_Pseudomonadales_plot.jpg",p25,  width=5.9, height=8.2, units="in")

##########SPA_Pseudomonadales

p26 <-ggplot(sampFiletax2,aes(SPA,Pseudomonadales)) +
geom_point(aes(),pch=21,size =4, stroke = 0.1,color= "#BAB6B6",alpha=0.8,fill="#F4BD27")+
#scale_x_continuous(limits = c(1.3, 4.1), breaks = seq(1.3, 4.1, 0.5))+
geom_smooth(method = "lm", formula = NULL,linewidth=0.2,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 100,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+theme1
p26

ggsave("SPA_Pseudomonadales_plot.pdf", p26,  width=5.9, height=8.2, units="in")
ggsave("SPA_Pseudomonadales_plot.jpg",p26,  width=5.9, height=8.2, units="in")

#######################################TOP5ORDER_METABO_COR
library(reshape2)

tax<-c("Pseudomonadales","Xanthomonadales","Burkholderiales","Sphingomonadales","Sphingobacteriales")
tax2<-t(taxonomy)
tax2<-as.data.frame(tax2[, tax], row.names = row.names(tax2))
sampFiletax2 =cbind(sampFile,tax2)

longsampFiletax <- melt(sampFiletax2,id.vars = -c(1:6),variable.name = 'metabo')
names(longsampFiletax)[ncol(longsampFiletax)] <- 'content'
longsampFiletax$metabo = factor(longsampFiletax$metabo, levels = c("HCA","FLA","LPA","CFA","SPLT","SPA"))

#longsampFiletax <- melt(longsampFiletax,id.vars = -c(1:5),variable.name = 'micro')
#names(longsampFiletax)[ncol(longsampFiletax)] <- 'RA'


colorsF <-c(Pseudomonadales="#F4BD27",Enterobacteriales="#55C1D7",Bacteroidales="#7BCBDF",
Burkholderiales="#B3C84D",Xanthomonadales="#4DFAA1",Rhizobiales="#EEEAB4",Sphingomonadales="#FADB75",
Sphingobacteriales="#A3F4B7",Flavobacteriales="#BAA5E6",
Rhodospirillales="#BDD2C7",Other="#FFFFFF",
HCA="#007808",FLA="#00AC0B",LPA="#009A3E",CFA="#13ED6A",SPLT="#38DF8E",SPA="#1ADFA4")

#############Xanthomonadales
p27 <-ggplot(longsampFiletax,aes(content,Xanthomonadales)) +
geom_point(aes(),pch=21,size=1.5,color= "#BAB6B6",alpha=0.8,fill="#4DFAA1")+
geom_smooth(method = "lm", formula = NULL,linewidth=0.5,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 25,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+facet_grid(~ metabo, scales = "free", space = "fixed")+theme_classic()
p27

ggsave("Metabo_Xanthomonadales_plot.pdf", p27, width=15, height=4, units="in")
ggsave("Metabo_Xanthomonadales_plot.jpg",p27, width=15, height=4, units="in")

#############Burkholderiales
p28 <-ggplot(longsampFiletax,aes(content,Burkholderiales)) +
geom_point(aes(),pch=21,size=1.5,color= "#BAB6B6",alpha=0.8,fill="#B3C84D")+
geom_smooth(method = "lm", formula = NULL,linewidth=0.5,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 35,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+facet_grid(~ metabo, scales = "free", space = "fixed")+theme_classic()
p28

ggsave("Metabo_Burkholderiales_plot.pdf", p28, width=15, height=4, units="in")
ggsave("Metabo_Burkholderiales_plot.jpg",p28, width=15, height=4, units="in")

#############Sphingomonadales
p29 <-ggplot(longsampFiletax,aes(content,Sphingomonadales)) +
geom_point(aes(),pch=21,size=1.5,color= "#BAB6B6",alpha=0.8,fill="#FADB75")+
geom_smooth(method = "lm", formula = NULL,linewidth=0.5,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 12,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+facet_grid(~ metabo, scales = "free", space = "fixed")+theme_classic()
p29

ggsave("Metabo_Sphingomonadales_plot.pdf", p29, width=15, height=4, units="in")
ggsave("Metabo_Sphingomonadales_plot.jpg",p29, width=15, height=4, units="in")

#############Sphingobacteriales
p30 <-ggplot(longsampFiletax,aes(content,Sphingobacteriales)) +
geom_point(aes(),pch=21,size=1.5,color= "#BAB6B6",alpha=0.8,fill="#B3C84D")+
geom_smooth(method = "lm", formula = NULL,linewidth=0.5,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 9,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+facet_grid(~ metabo, scales = "free", space = "fixed")+theme_classic()
p30

ggsave("Metabo_Sphingobacteriales_plot.pdf", p30, width=15, height=4, units="in")
ggsave("Metabo_Sphingobacteriales_plot.jpg",p30, width=15, height=4, units="in")



###################################################PAL_SNP_microbiome_correlation
#PAL_family_gene_structure
library(ggplot2)
library(gggenes)
library (dplyr)
library(ggpubr)

allPALgene2=read.csv("PAL_SNP.csv", header=T,  comment.char="", stringsAsFactors=F)

allPALgene<-subset(allPALgene2, lableing == "YES" )

#allPALgene<-allPALgene2

allPALgene$label = factor(allPALgene$label, levels = c("OsPAL07" ,"OsPAL04","OsPAL03", 
"OsPAL01","OsPAL06","OsPAL05", "OsPAL08","OsPAL02" ))

p49 <-ggplot(allPALgene, aes(xmin = start, xmax = end,  y = label)) +
geom_gene_arrow(fill="#F4BD27") +
geom_subgene_arrow(data = allPALgene,
                     aes(xmin = 0, xmax = end+500,
                         y =label,
                         fill = pvalue,color = pvalue,
                         xsubmin =position, xsubmax = position+20), alpha = 0.7) +
scale_colour_gradient(low="#B5CBE2",high="#007808")+
scale_fill_gradient(low="#B5CBE2",high="#007808")+
theme_genes()+labs( y = '')+
geom_subgene_label(data = subset(allPALgene, lableing == "YES" ),
    aes(xsubmin = position, xsubmax = position+2000, label = name),
    min.size = 0,grow = TRUE)
p49

ggsave(paste0("PAL_SNP.pdf"), p49, width=6, height=4, units="in")
ggsave(paste0("PAL_SNP.jpg"), p49, width=6, height=4, units="in")

#The microbiome effect size of PAL-SNP that significantly associated metabolite
library(vegan)
library(ggplot2)
library(ggprism)
library(reshape2)
library(tidyverse)

colors <-c(HAP1="#97BBE1",HAP2 ="#E7B35A",Unassigned ="#00FF7F")
level =c("HAP2","HAP1")
analysisgroup="Haplotype"
shapes<-c(HAP1=21,HAP2 =22,NS =23,leaf=24,CK=25)

metadata=read.csv("HAP_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata
taxonomy[is.na(taxonomy)]<-0

idx = rownames(metadata) %in% colnames(taxonomy)
metadata = metadata[idx, , drop = F]
taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]

data <- vegdist(t(taxonomy),method = "bray")
write.table(as.matrix(data), 'microbiobray.txt', sep = '\t', col.names = NA, quote = FALSE)#输出距离矩阵
otu_matrix<-read.table("microbiobray.txt", header=T, row.names=1, sep="\t", comment.char="")


for (col_name in colnames(metadata)) {
  if (col_name == "Group"){
    next
  }
PCOAanalysisgroup = col_name
sampFile = as.data.frame(metadata[, c(PCOAanalysisgroup)], row.names = row.names(metadata))
colnames(sampFile)[1] = "group"

#非限制性排序及组间显著性分析
#sampFile = as.data.frame(metadata[, "group"], row.names = row.names(metadata))
pcoa = cmdscale(otu_matrix, k = 3, eig = T)
points = as.data.frame(pcoa$points)
eig = pcoa$eig
points = cbind(points, sampFile)
colnames(points) = c("x", "y", "z", "group")

#显著性
adonis_table = adonis2(otu_matrix ~ group, data = sampFile, permutations = 999)
adonis_pvalue = adonis_table$`Pr(>F)`[1]
adonis_R2=adonis_table$`R2`[1]
adonis_R2pvalue = paste(col_name,"group","R2:", adonis_R2, "pvalue:", adonis_pvalue, sep = "\t")
write.table(adonis_R2pvalue, file = "microbio_pcoa_adonis2test.txt", append = TRUE, 
            sep = "\t", quote = F, row.names = F, col.names = F)

pcoa.fit <- envfit(pcoa ~ group, data=sampFile, perm=999)
envfitR2 = paste(col_name,"group R2:", pcoa.fit$factors$r, "Pvalue:", pcoa.fit$factors$pvals,sep = "\t")
write.table(envfitR2,"microbio_pcoa_envfit.txt",append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

points$group = factor(points$group, levels = level )
#可视化
p51 = ggplot(points, aes(x = x, y = z, shape = group,fill = group)) + 
  labs(x = paste("PCo 1 (", format(100 * eig[1]/sum(eig), digits = 4), "%)", sep = ""), 
       y = paste("PCo 3 (", format(100 * eig[3]/sum(eig), digits = 4), "%)", sep = ""), color = "group")+
  scale_fill_manual(values = colors)+
  scale_colour_manual(values = colors)+
  scale_shape_manual(values = shapes)
  #+ scale_fill_gradient(colors = c("Green", "Blue", "Red"),limits = c(- 10, 200))
  #+ scale_colour_distiller(palette = "Spectral")
  #+ scale_fill_distiller(palette = "Spectral")

p51 = p51+ geom_point(alpha = 0.6, size =4, stroke = 0.1, color = "#BAB6B6") + theme_classic() + 
  theme(text = element_text(family = "sans", size = 7))+
  #stat_ellipse(level = 0.60)+ 
  ggtitle(paste(col_name," p = ", pcoa.fit$factors$pvals,  " R2 = ", pcoa.fit$factors$r, sep = ""))
p51

ggsave(paste0(col_name,".pdf"), p51, width=5, height=5, units="in")
ggsave(paste0(col_name,".jpg"), p51, width=5, height=5, units="in")
}


##########correlation of p-value and microbiome effect size 

library(ggplot2)

hcagene_microbiome=read.csv("gwas_pvalue.csv", header=T, comment.char="", stringsAsFactors=F)

theme1<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=12),
        axis.line = element_line(color = "black",size = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        # axis.text.x = element_blank(),
        # axis.ticks.x=element_blank(),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
        # legend.position = "none",
       panel.spacing = unit(0,"lines")
       )

p52 <-ggplot(hcagene_microbiome,aes(microbiome_effect_size,HCA_association)) +
geom_point(aes(),pch=21,size =8, stroke = 0.1,color= "#BAB6B6",alpha=0.8,fill="#007808",position=position_jitter(0.001))+
#scale_x_continuous(limits = c(1.3, 4.1), breaks = seq(1.3, 4.1, 0.5))+
geom_smooth(method = "lm", formula = 'y ~ x',linewidth=0.2,se=T,color="blue",linetype="dashed",aes(group=1))+
stat_cor(label.y = 3.8,aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"),group=1),color="black",method = "pearson",
           label.x.npc = "left")+theme1
p52

p53<-p52+ggrepel::geom_text_repel(aes(label=name),hcagene_microbiome,
    size = 3, #注释文本的字体大小
    box.padding = 0.1, #字到点的距离
    point.padding = 0.3, #字到点的距离，点周围的空白宽度
  #  min.segment.length = 0.5, #短线段可以省略
    segment.color = "black", #segment.colour = NA, 不显示线段
    show.legend = F)
p53

ggsave("HCA_hcagene_microbiome.pdf", p53, width=6.5, height=5, units="in")
ggsave("HCA_hcagene_microbiome.jpg",p53, width=6.5, height=5,units="in")