
##########################################################loding packages and data
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(ggprism)
library(reshape2)
library(tidyverse)
library(ggpmisc)
library(plyr)
library(rstatix)

colorsw <-c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
  "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96","#FF00FF","#EE82EE","#FFA500",
  "#009ACD","#AB82FF","#00FFFF","#8B668B","#CD4F39","#8EE5EE","#FFDEAD")
colors <-c(INDICA="#B5CBE2",JAPONICA ="#E8BE74",Unassigned ="#00FF7F",
           leaf="#FFF68F",CK="#FFA500")
level =c("INDICA","JAPONICA")
topN = 20
analysisgroup="Group"
shapes<-c(INDICA=21,JAPONICA =22,NS =23,leaf=24,CK=25)

metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

###############################################################代谢物分析
###############################################################metabolite analysis 

#进化树热图数据
#Heatmap data of evolutionary tree 

analysisgroup2 = "Group1"

corrmetabo<-c("HCA","SPLT","CFA","LPA","FLA","SPA")
corrmetaboa1<-c(analysisgroup2,corrmetabo)
sampFile3 = as.data.frame(metadata[, corrmetaboa1], row.names = row.names(metadata))
mat_mean = aggregate(sampFile3[, 2:ncol(sampFile3)], by = sampFile3[1], FUN = mean)

#代谢物在组间间差异分析
#Metabolite difference analysis between groups

analysisgroup3 = "Group"
corrmetabo<-c("HCA","SPLT","CFA","LPA","FLA","SPA")
corrmetaboa1<-c(analysisgroup3,corrmetabo)
sampFile3 = as.data.frame(metadata[, corrmetaboa1], row.names = row.names(metadata))
df <- sampFile3 %>%mutate(Group=factor(Group,levels =level))
summetabo<-melt(sampFile3,id.vars=c("Group"))

summetabolite<-ddply(summetabo,c("Group","variable"),summarise,allmean=mean(value,na.rm=TRUE),sd1=sd(value,na.rm = TRUE),
                                                                n1=sum(!is.na(value)),se1=sd1/sqrt(n1))
write.csv(as.data.frame(summetabolite),"summetabolite.csv")

theme<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=12),
        axis.line = element_line(color = "black",linewidth = 0.4),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=10,angle=0),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
        legend.position = "none",
       panel.spacing = unit(0,"lines")
       )


#所有代谢物组间差异all-pvalue###校正#adjust_pvalue(method = "bonferroni")
#p-value for all metabolite groups#adjust_pvalue (method = "bonferroni")

dfn<- df[c(corrmetabo, analysisgroup3)]
dfn<-melt(dfn,id.vars=analysisgroup3)

names(dfn) <- c( "Group","metabolite","Content")
stat.test<- dfn%>%group_by(metabolite) %>% t_test(Content ~ Group)%>%adjust_pvalue(method = "bonferroni")%>%
add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
write_delim(as.data.frame(stat.test), paste0(analysisgroup3,".txt"))

#4HCA差异 5:4 PORTRAIT
#Difference analysis of 4HCA ; 5:4 PORTRAIT

col_name = "HCA"
dfn <- df[c(col_name, analysisgroup3)]
names(dfn) <- c("Content", "Group")
stat.test<- dfn %>% t_test(Content ~ Group)%>%
add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
write_delim(as.data.frame(stat.test), paste0(col_name,".txt"))

p1 <-ggplot(dfn,aes(Group,Content)) +
geom_boxplot(aes(alpha = 1,colour=Group),alpha=1, width = 0.7,outlier.shape = NA, outlier.size = 0,fill= "transparent")+
#geom_jitter(aes(fill=Group),col="black",alpha=0.9, position = position_jitter(0.2))+
geom_point(aes(fill=Group,alpha=0.7,colour =Group),alpha=0.8,size =0.8, shape = 21, stroke = 0.5,position = position_jitter(0.25))+
scale_fill_manual(values = colors)+
scale_colour_manual(values = colors)+
stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F)+
scale_size_continuous(range=c(1,3)) +  
ggtitle(paste(col_name," content"))+theme
p1

ggsave(paste0(col_name,".pdf"), p1, width=4, height=5, units="in")
ggsave(paste0(col_name,".png"), p1, width=4, height=5, units="in")

#SPLT差异 5:3 PORTRAIT
#Difference analysis of SPLT ; 5:4 PORTRAIT
col_name = "SPLT"
dfn <- df[c(col_name, analysisgroup3)]
names(dfn) <- c("Content", "Group")
stat.test<- dfn %>% t_test(Content ~ Group)%>%
add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
write_delim(as.data.frame(stat.test), paste0(col_name,".txt"))

p2 <-ggplot(dfn,aes(Group,Content)) +
geom_boxplot(aes(alpha = 1,colour=Group),alpha=1, width = 0.7,outlier.shape = NA, outlier.size = 0,fill= "transparent")+
#geom_jitter(aes(fill=Group),col="black",alpha=0.9, position = position_jitter(0.2))+
geom_point(aes(fill=Group,alpha=0.7,colour =Group),alpha=0.8,size =0.8, shape = 21, stroke = 0.5,position = position_jitter(0.25))+
scale_fill_manual(values = colors)+
scale_colour_manual(values = colors)+
stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F)+
scale_size_continuous(range=c(1,3)) +  
ggtitle(paste(col_name," content"))+theme
p2


ggsave(paste0(col_name,".pdf"), p2, width=4, height=5, units="in")
ggsave(paste0(col_name,".png"), p2, width=4, height=5, units="in")

#CFA差异 5:3 PORTRAIT
#Difference analysis of CFA ; 5:4 PORTRAIT

col_name = "CFA"
dfn <- df[c(col_name, analysisgroup3)]
names(dfn) <- c("Content", "Group")
stat.test<- dfn %>% t_test(Content ~ Group)%>%
add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
write_delim(as.data.frame(stat.test), paste0(col_name,".txt"))

p3 <-ggplot(dfn,aes(Group,Content)) +
geom_boxplot(aes(alpha = 1,colour=Group),alpha=1, width = 0.7,outlier.shape = NA, outlier.size = 0,fill= "transparent")+
#geom_jitter(aes(fill=Group),col="black",alpha=0.9, position = position_jitter(0.2))+
geom_point(aes(fill=Group,alpha=0.7,colour =Group),alpha=0.8,size =0.8, shape = 21, stroke = 0.5,position = position_jitter(0.25))+
scale_fill_manual(values = colors)+
scale_colour_manual(values = colors)+
stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F)+
scale_size_continuous(range=c(1,3)) +  
ggtitle(paste(col_name," content"))+theme
p3

ggsave(paste0(col_name,".pdf"), p3, width=4, height=5, units="in")
ggsave(paste0(col_name,".png"), p3, width=4, height=5, units="in")

#LPA差异 5:3 PORTRAIT
#Difference analysis of LPA ; 5:4 PORTRAIT

col_name = "LPA"
dfn <- df[c(col_name, analysisgroup3)]
names(dfn) <- c("Content", "Group")
stat.test<- dfn %>% t_test(Content ~ Group)%>%
add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
write_delim(as.data.frame(stat.test), paste0(col_name,".txt"))

p4 <-ggplot(dfn,aes(Group,Content)) +
geom_boxplot(aes(alpha = 1,colour=Group),alpha=1, width = 0.7,outlier.shape = NA, outlier.size = 0,fill= "transparent")+
#geom_jitter(aes(fill=Group),col="black",alpha=0.9, position = position_jitter(0.2))+
geom_point(aes(fill=Group,alpha=0.7,colour =Group),alpha=0.8,size =0.8, shape = 21, stroke = 0.5,position = position_jitter(0.25))+
scale_fill_manual(values = colors)+
scale_colour_manual(values = colors)+
stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F)+
scale_size_continuous(range=c(1,3)) +  
ggtitle(paste(col_name," content"))+theme
p4

ggsave(paste0(col_name,".pdf"), p4, width=4, height=5, units="in")
ggsave(paste0(col_name,".png"), p4, width=4, height=5, units="in")

#FLA差异 5:4 PORTRAIT
#Difference analysis of FLA ; 5:4 PORTRAIT

col_name = "FLA"
dfn <- df[c(col_name, analysisgroup3)]
names(dfn) <- c("Content", "Group")
stat.test<- dfn %>% t_test(Content ~ Group)%>%
add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
write_delim(as.data.frame(stat.test), paste0(col_name,".txt"))

p5 <-ggplot(dfn,aes(Group,Content)) +
geom_boxplot(aes(alpha = 1,colour=Group),alpha=1, width = 0.7,outlier.shape = NA, outlier.size = 0,fill= "transparent")+
#geom_jitter(aes(fill=Group),col="black",alpha=0.9, position = position_jitter(0.2))+
geom_point(aes(fill=Group,alpha=0.7,colour =Group),alpha=0.8,size =0.8, shape = 21, stroke = 0.5,position = position_jitter(0.25))+
scale_fill_manual(values = colors)+
scale_colour_manual(values = colors)+
stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F)+
scale_size_continuous(range=c(1,3)) +  
ggtitle(paste(col_name," content"))+theme
p5

ggsave(paste0(col_name,".pdf"), p5, width=4, height=5, units="in")
ggsave(paste0(col_name,".png"), p5, width=4, height=5, units="in")

#SPA差异 5:4 PORTRAIT
#Difference analysis of SPA ; 5:4 PORTRAIT

col_name = "SPA"
dfn <- df[c(col_name, analysisgroup3)]

names(dfn) <- c("Content", "Group")
stat.test<- dfn %>% t_test(Content ~ Group)%>%adjust_pvalue(method = "bonferroni")%>%
add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
write_delim(as.data.frame(stat.test), paste0(col_name,".txt"))

p6 <-ggplot(dfn,aes(Group,Content)) +
geom_boxplot(aes(alpha = 1,colour=Group),alpha=1, width = 0.7,outlier.shape = NA, outlier.size = 0,fill= "transparent")+
#geom_jitter(aes(fill=Group),col="black",alpha=0.9, position = position_jitter(0.2))+
geom_point(aes(fill=Group,alpha=0.7,colour =Group),alpha=0.8,size =0.8, shape = 21, stroke = 0.5,position = position_jitter(0.25))+
scale_fill_manual(values = colors)+
scale_colour_manual(values = colors)+
stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F)+
scale_size_continuous(range=c(1,3)) +  
ggtitle(paste(col_name," content"))+theme
p6

ggsave(paste0(col_name,".pdf"), p6, width=4, height=5, units="in")
ggsave(paste0(col_name,".png"), p6, width=4, height=5, units="in")

################代谢物PCOA
################PCoA of metabolite

analysisgroup3 = "Group"
corrmetabo<-c("LPA","HCA","FLA","CFA","SPLT","SPA")

mtbcoaFile = as.data.frame(metadata[,corrmetabo], row.names = row.names(metadata))

#####################################################metabopcoa
data <- vegdist(as.matrix(mtbcoaFile),method = "bray")
write.table(as.matrix(data), 'metabobray.txt', sep = '\t', col.names = NA, quote = FALSE)#输出距离矩阵
otu_matrix<-read.table("metabobray.txt", header=T, row.names=1, sep="\t", comment.char="")

sampFile = as.data.frame(metadata[, c(analysisgroup3)], row.names = row.names(metadata))
colnames(sampFile)[1] = "group"

pcoa = cmdscale(otu_matrix, k = 3, eig = T)
points = as.data.frame(pcoa$points)
eig = pcoa$eig
points = cbind(points, sampFile)
colnames(points) = c("x", "y", "z", "group")

#显著性
adonis_table = adonis2(otu_matrix ~ group, data = sampFile, permutations = 999)
adonis_pvalue = adonis_table$`Pr(>F)`[1]
adonis_R2=adonis_table$`R2`[1]
adonis_R2pvalue = paste("group","R2:", adonis_R2, "pvalue:", adonis_pvalue, sep = "\t")
write.table(adonis_R2pvalue, file = "adonis2test.txt", append = TRUE, 
            sep = "\t", quote = F, row.names = F, col.names = F)

pcoa.fit <- envfit(pcoa ~ group, data=sampFile, perm=999)
envfitR2 = paste("group R2:", pcoa.fit$factors$r, "Pvalue:", pcoa.fit$factors$pvals,sep = "\t")
write.table(envfitR2,"metaboenvfit.txt",append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

#可视化
p7 = ggplot(points, aes(x = y, y =x, color = group,shape = group,fill = group)) + 
  labs(y = paste("PCo 1 (", format(100 * eig[1]/sum(eig), digits = 4), "%)", sep = ""), 
       x = paste("PCo 2 (", format(100 * eig[2]/sum(eig), digits = 4), "%)", sep = ""), color = "group")+
  scale_fill_manual(values = colors)+
  scale_colour_manual(values = colors)+
  scale_shape_manual(values = shapes)

p7 = p7+ geom_point(alpha = 0.6, size = 2,stroke=0.1) + theme_classic() + 
  theme(text = element_text(family = "sans", size = 7))+stat_ellipse(level = 0.68)+ 
  ggtitle(paste( " P = ", pcoa.fit$factors$pvals,  " R2 = ", pcoa.fit$factors$r, sep = ""))
p7

ggsave(paste0("metpcoa.pdf"), p7, width=5, height=5, units="in")
ggsave(paste0("metapcoa.png"), p7, width=5, height=5, units="in")

#############代谢物相关性矩阵
#############Metabolite correlation matrix
#相关性all_smooth and correlation

library(GGally)

p9<-ggpairs(mtbcoaFile, columns = 1:ncol(mtbcoaFile),
aes(alpha = 0.3),cardinality_threshold = 1000,
lower = list(continuous = "smooth"))
p9
ggsave(paste0("all_metabocor.pdf"), p9, width=5, height=5, units="in")
ggsave(paste0("all_metabocor.png"), p9, width=5, height=5, units="in")

###############################################################2.microbiome analysis

################################alpha多样性
################################alpha diversity
alphadiversity<-read.csv("vegan.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

idx = rownames(metadata) %in% rownames(alphadiversity)
metadata1 = metadata[idx, , drop = F]
alpha = alphadiversity[rownames(metadata1),]
sampFile2 = as.data.frame(metadata1[, c(analysisgroup)], row.names = row.names(metadata1))
colnames(sampFile2)[1] = "Group"
alpha_dif = cbind(sampFile2, alpha)
sumalpha =melt(alpha_dif, id.vars = c("Group"))

sumalphadiver<-ddply(sumalpha,c("Group","variable"),summarise,allmean=mean(value,na.rm=TRUE),sd1=sd(value,na.rm = TRUE),
                                                                n1=sum(!is.na(value)),se1=sd1/sqrt(n1))
write.csv(as.data.frame(sumalphadiver),"sumalphadiver.csv")

df <- alpha_dif%>%mutate(Group=factor(Group,levels =level))

for (col_name in colnames(df)) {
  if (col_name == "Group"){
    next
  }
filename <- paste0(col_name,".pdf")
dfn <- df[c({{col_name}}, "Group")]
names(dfn) <- c("alpha", "Group")
stat.test<- dfn %>% t_test(alpha ~ Group)%>%
add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
write_delim(as.data.frame(stat.test), paste0(col_name,".txt"))

p <-ggplot(dfn,aes(Group,alpha)) +
#aes()对数据进行定义，函数内的参数为对函数定义
#geom_violin(aes(fill=Group),alpha = 0.5,width = 0.5) +
geom_boxplot(aes(alpha = 1,colour=Group),alpha=0, width = 0.7,outlier.shape = NA, outlier.size = 0,fill= "transparent")+
#geom_jitter(aes(fill=Group),col="black",alpha=0.9, position = position_jitter(0.2))+
geom_point(aes(fill=Group,alpha=0.7,colour =Group),alpha=0.8,size =0.8, shape = 21, stroke = 0.5,position = position_jitter(0.25))+
scale_fill_manual(values = colors)+
scale_colour_manual(values = colors)+
stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F)+
scale_size_continuous(range=c(1,3)) +  
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=12),
        axis.line = element_line(color = "black",size = 0.4),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=10,angle=-90),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
        legend.position = "none",
       panel.spacing = unit(0,"lines")
       )+
       labs(title = col_name)
  coord_cartesian(clip = "off")
ggsave(paste0(filename), p, width=4, height=5, units="in")
ggsave(paste0(col_name,".png"), p, width=4, height=5, units="in")
}


#######################################################PCOA of order level

taxonomy=read.csv("RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata
taxonomy[is.na(taxonomy)]<-0

idx = rownames(metadata) %in% colnames(taxonomy)
metadata = metadata[idx, , drop = F]
taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]

data <- vegdist(t(taxonomy),method = "bray")
write.table(as.matrix(data), 'microbiobray.txt', sep = '\t', col.names = NA, quote = FALSE)#输出距离矩阵
otu_matrix<-read.table("microbiobray.txt", header=T, row.names=1, sep="\t", comment.char="")

PCOAanalysisgroup="Group"
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
adonis_R2pvalue = paste("group","R2:", adonis_R2, "pvalue:", adonis_pvalue, sep = "\t")
write.table(adonis_R2pvalue, file = "microbio_pcoa_adonis2test.txt", append = TRUE, 
            sep = "\t", quote = F, row.names = F, col.names = F)

pcoa.fit <- envfit(pcoa ~ group, data=sampFile, perm=999)
envfitR2 = paste("group R2:", pcoa.fit$factors$r, "Pvalue:", pcoa.fit$factors$pvals,sep = "\t")
write.table(envfitR2,"microbio_pcoa_envfit.txt",append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

#可视化
p8 = ggplot(points, aes(x = x, y = y, color = group,shape = group,fill = group)) + 
  labs(x = paste("PCo 1 (", format(100 * eig[1]/sum(eig), digits = 4), "%)", sep = ""), 
       y = paste("PCo 2 (", format(100 * eig[2]/sum(eig), digits = 4), "%)", sep = ""), color = "group")+
  scale_fill_manual(values = colors)+
  scale_colour_manual(values = colors)+
  scale_shape_manual(values = shapes)
  #+ scale_fill_gradient(colors = c("Green", "Blue", "Red"),limits = c(- 10, 200))
  #+ scale_colour_distiller(palette = "Spectral")
  #+ scale_fill_distiller(palette = "Spectral")

p8 = p8+ geom_point(alpha = 0.6, size = 0.6) + theme_classic() + 
  theme(text = element_text(family = "sans", size = 7))+stat_ellipse(level = 0.60)+ 
  ggtitle(paste( " p = ", pcoa.fit$factors$pvals,  " R2 = ", pcoa.fit$factors$r, sep = ""))
p8

ggsave("microPCOA12.pdf", p8, width=5, height=5, units="in")
ggsave("microPCOA12.jpg", p8, width=5, height=5, units="in")


###################################################目水平堆叠图
###################################################Stacking diagram at order level

taxonomy=read.csv("RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata
taxonomy[is.na(taxonomy)]<-0

idx = rownames(metadata) %in% colnames(taxonomy)
metadata = metadata[idx, , drop = F]
taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]

colorsO <-c(Pseudomonadales="#F4BD27",Enterobacteriales="#55C1D7",Bacteroidales="#7BCBDF",
Clostridiales="#A5D9E6",Burkholderiales="#B3C84D",Xanthomonadales="#4DFAA1",
Rhizobiales="#EEEAB4",Lactobacillales="#A3F4CA",Sphingomonadales="#FADB75",
Rhodospirillales="#BDD2C7",Other="#E3E3E3")

level =c("INDICA","JAPONICA")
topN = 11
analysisgroup="Group"

mean_sort = as.data.frame(taxonomy[(order(-rowSums(taxonomy))), ])
idx = grepl("unassigned|unclassified|unknown", rownames(mean_sort), ignore.case = T)
mean_sort = rbind(mean_sort[!idx, ], mean_sort[idx, ])
other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(topN - 1), ]
mean_sort = rbind(mean_sort, other)
rownames(mean_sort)[topN] = c("Other")
amplic2 = t(mean_sort)

analysisgroup2 = "Group1"
sampFile = as.data.frame(metadata[, c(analysisgroup,analysisgroup2 )], row.names = row.names(metadata))
colnames(sampFile)[2] = "Group1"
colnames(sampFile)[1] = "Group"

mat_t2 = merge(sampFile, amplic2, by = "row.names")
mat_t2 = mat_t2[, c(-1)]
mat_mean = aggregate(mat_t2[, 3:ncol(mat_t2)], by = mat_t2[2], FUN = mean)
write.csv(mat_mean,"RA_mean.csv")

mean_sort2 = merge(mat_mean ,sampFile, by = "Group1")
mean_sort2 <- unique(mean_sort2)

data_all = as.data.frame(melt(mean_sort2, id.vars = c("Group1","Group")))
data_all$variable = factor(data_all$variable, levels = rownames(mean_sort))
data_all$value = as.numeric(data_all$value)

#x轴按丰度排序
mat_mean$Group1 <- reorder(mat_mean$Group1, -mat_mean$Pseudomonadales)
data_all <- data_all %>%mutate(Group1=factor(Group1,levels = levels(mat_mean$Group1)))
#data_all <- unique(data_all)
#x轴按自定义
#data_all <- data_all %>%mutate(variable=factor(variable,levels =level))

p7 = ggplot(data_all, aes(x = Group1, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "fill", width = 1,colour ="white",linewidth=0.001) + 
 # geom_bar(stat = "identity", position = "stack", width = 0.9,colour ="white",linewidth=0.001) + 
  scale_y_continuous(labels = scales::percent) +
  xlab("Groups") + ylab("Relative abundance (%)") + theme_classic() + 
  theme(text = element_text(family = "sans", size = 5),legend.position="top", legend.direction = "horizontal" ,
        axis.text.y = element_text(color="black",size=5),
        axis.text.x = element_text(color="black",size=5,angle=45))+
  scale_fill_manual(values = colorsO)+
  facet_grid(~ Group, scales = "free", space = "free")+guides(fill = guide_legend(nrow = 1, bycol = TRUE))
p7

ggsave("microstack.pdf", p7, width=10, height=4, units="in")
ggsave("microstack.jpg", p7, width=10, height=4, units="in")

#############################################方差分析
#############################################Analysis of variance
#######genus_t_test

taxonomy=read.csv("GenusRA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata
taxonomy[is.na(taxonomy)]<-0

idx = rownames(metadata) %in% colnames(taxonomy)
metadata = metadata[idx, , drop = F]
taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]

sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
colnames(sampFile)[1] = "Group"
level =c("INDICA","JAPONICA")

dat1 <- as.data.frame(merge(t(taxonomy),sampFile,by='row.names'))
df1 <-dat1[,-1]
dfngenus = as.data.frame(melt(df1, id.vars = c("Group")))

genusstat.test<- dfngenus %>% group_by(variable)%>% t_test(value ~ Group)%>%
add_significance(p.col = "p")
write.csv(as.data.frame(genusstat.test),"genus-ttest.csv")

sumgenus<-ddply(dfngenus,c("Group","variable"),summarise,allmean=mean(value,na.rm=TRUE),sd1=sd(value,na.rm = TRUE),
                                                                n1=sum(!is.na(value)),se1=sd1/sqrt(n1))
write.csv(as.data.frame(sumgenus),"sumgenus.csv")


#######order_t_test
taxonomy=read.csv("RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata
taxonomy[is.na(taxonomy)]<-0

idx = rownames(metadata) %in% colnames(taxonomy)
metadata = metadata[idx, , drop = F]
taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]

sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
colnames(sampFile)[1] = "Group"
level =c("INDICA","JAPONICA")

dat1 <- as.data.frame(merge(t(taxonomy),sampFile,by='row.names'))
df1 <-dat1[,-1]
dfnorder = as.data.frame(melt(df1, id.vars = c("Group")))

orderstat.test<- dfnorder %>% group_by(variable)%>% t_test(value ~ Group)%>%
add_significance(p.col = "p")
write.csv(as.data.frame(orderstat.test),"orderstat-test.csv")

sumorder<-ddply(dfnorder,c("Group","variable"),summarise,allmean=mean(value,na.rm=TRUE),sd1=sd(value,na.rm = TRUE),
                                                                n1=sum(!is.na(value)),se1=sd1/sqrt(n1))
write.csv(as.data.frame(sumorder),"sumorder.csv")


#######otu_t_test
taxonomy=read.csv("otu.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata
taxonomy[is.na(taxonomy)]<-0

idx = rownames(metadata) %in% colnames(taxonomy)
metadata = metadata[idx, , drop = F]
taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]

sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
colnames(sampFile)[1] = "Group"
level =c("INDICA","JAPONICA")

dat1 <- as.data.frame(merge(t(taxonomy),sampFile,by='row.names'))
df1 <-dat1[,-1]
dfnorder = as.data.frame(melt(df1, id.vars = c("Group")))

orderstat.test<- dfnorder %>% group_by(variable)%>% t_test(value ~ Group)%>%
add_significance(p.col = "p")
write.csv(as.data.frame(orderstat.test),"otustat-test.csv")

sumorder<-ddply(dfnorder,c("Group","variable"),summarise,allmean=mean(value,na.rm=TRUE),sd1=sd(value,na.rm = TRUE),
                                                                n1=sum(!is.na(value)),se1=sd1/sqrt(n1))
write.csv(as.data.frame(sumorder),"sumotur.csv")
##################TOP10目中属水平差异分析图
##################Analysis of genus level differences in TOP10 orders

genusdif=read.csv("Genusdif.csv", header=T,  comment.char="", stringsAsFactors=F)

df2=genusdif%>%filter(-log10(p) > 7)

shapes<-c(INDICA=21,JAPONICA =22,NS =25,leaf=24,CK=25)
colors <-c(INDICA="#B5CBE2",JAPONICA ="#E8BE74",NS ="#E3E3E3",
           leaf="#FFF68F",CK="#FFA500")

colorsO <-c(Pseudomonadales="#F4BD27",Enterobacteriales="#55C1D7",Bacteroidales="#7BCBDF",
Clostridiales="#A5D9E6",Burkholderiales="#B3C84D",Xanthomonadales="#4DFAA1",
Rhizobiales="#EEEAB4",Lactobacillales="#A3F4CA",Sphingomonadales="#FADB75",
Rhodospirillales="#BDD2C7",Other="#FFFFFF")

genusdif$Order = factor(genusdif$Order, levels = c("Pseudomonadales","Enterobacteriales","Bacteroidales",
"Clostridiales","Burkholderiales","Xanthomonadales","Rhizobiales","Lactobacillales","Sphingomonadales","Rhodospirillales"))

p9 <- ggplot(data = genusdif,aes(x = Order, y =  -log10(p))) + 
  geom_point(aes(fill = enrichment,color= enrichment,shape= enrichment, size =abs(enrichment_index) ), alpha=1,stroke = 0.5,position = position_jitter(0.45)) +
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  scale_fill_manual(values=colors)+
#  scale_x_continuous("Order") +
#  scale_y_continuous("-log10 (p-value)",limits = c(0,12),breaks = seq(0,12,4))+
  geom_hline(yintercept=c(1.3),lty=2,col="black",lwd=0.5) + theme_classic() + 
  theme(text = element_text(family = "sans", size = 5),legend.position="top", legend.direction = "horizontal" ,
        axis.text.y = element_text(color="black",size=5),
        axis.text.x = element_text(color="black",size=5,angle=45))
p9

  p10<-p9+ggrepel::geom_text_repel( aes(label=variable),df2,
    size = 1, #注释文本的字体大小
    box.padding = 0, #字到点的距离
    point.padding = 0.3, #字到点的距离，点周围的空白宽度
  #  min.segment.length = 0.5, #短线段可以省略
    segment.color = "black", #segment.colour = NA, 不显示线段
    show.legend = F)
  p10

ggsave("mhatonplot.pdf", p10, width=10, height=3, units="in")
ggsave("mhatonplot.jpg",p10, width=10, height=3, units="in")
