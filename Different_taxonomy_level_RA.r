
# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}
# 解析参数-h显示帮助信息
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
              help="Input reads count file; such as OTU table, kraken2 taxonomy counts table [default %default]"),
   make_option(c("-t", "--taxonomy"), type="character", default="result/taxonomy.txt",
              help="Input taxonomy file; such as OTU taxonomy table [default %default]"),
  make_option(c("-v", "--Visualization"), type="numeric", default=0,
              help="Visualizatio; default 0 not to visualization [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result/RA/",
              help="Output file dirction; Output Relative abundance dirction [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))


# 显示输入输出确认是否正确
print(paste("The input feature table is ", opts$input,  sep = ""))
# print(paste("Normalized filename: ", opts$normalize,  sep = ""))
# print(paste("Output alpha diversity: ", opts$output, sep = ""))

# suppressWarnings(dir.create("alpha/"))

# 1.3 安装CRAN来源常用包
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list <- c("dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

otu<-read.table(opts$input,sep = "\t",head=T,row.names=1,quote = "")
tax<-read.table(opts$taxonomy,head=T,row.names=1,quote = "")

row_names <- as.vector(rownames(otu))
otutax <- tax[row_names,]
otu_tax <- cbind(otu, tax[row_names,])

#kingdom
col_num <- ncol(otu)
tax_col_num <- ncol(tax)
for (i in 1:tax_col_num) {
group <- otu_tax[,(ncol(otu)+i)]
sums <- rowsum(otu_tax[,(1:col_num)],group)
write.csv(sums,file=paste(opts$output,names(tax)[i],".csv"))
col_sum <- colSums(sums)
sums_RA <- sums
for (j in 1:ncol(sums)) {
for (r in 1:nrow(sums)) {
sample_sum<-as.vector(col_sum[j])
RA<-((sums_RA[r,j]/sample_sum)*100)
sums_RA[r,j]<-RA
}
}
write.csv(sums_RA,file=paste(opts$output,names(tax)[i],"RA.csv"))
}

