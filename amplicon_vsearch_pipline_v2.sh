mkdir rawdata
cp raw_reads/* rawdata/
mkdir result
seqkit stat rawdata/*.fastq.gz > result/seqkit.txt
mkdir temp
rename 's/\-/_/g'  rawdata/*.fastq.gz
rename 's/\.R/\./'  rawdata/*.fastq.gz 

## 2. 序列合并和重命名 reads merge and rename
#方法1.for循环顺序处理
time for i in `tail -n+2 metadata.txt|cut -f1`;do
    vsearch --fastq_mergepairs rawdata/${i}.1.fastq.gz --reverse rawdata/${i}.2.fastq.gz \
    --fastqout temp/${i}.merged.fq --relabel ${i}.
    done &&

# 检查最后一个文件前10行中样本名
    #head temp/`tail -n+2 metadata.txt | cut -f 1 | tail -n1`.merged.fq | grep ^@

### 2.3 改名后序列整合 integrate renamed reads
time  cat temp/*.merged.fq > temp/all.fq &&
seqkit stat temp/all.fq > result/allseqkit.txt &&
    #head -n 6 temp/all.fq|cut -c1-60


## 3. 切除引物与质控 Cut primers and quality filter

    # 左端10bp标签+19bp上游引物V5共为29，右端V7为18bp下游引物
    # Cut barcode 10bp + V5 19bp in left and V7 18bp in right
    time vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 25 --fastq_stripright 15 \
      --fastq_maxee_rate 0.01 \
      --fastaout temp/filtered.fa &&
    # 查看文件了解fa文件格式
    #head temp/filtered.fa


### 4.1 序列去冗余 Dereplicate

    # 并添加miniuniqusize最小为8或1/1M，去除低丰度噪音并增加计算速度
    # -sizeout输出丰度, --relabel必须加序列前缀更规范, 1s
    vsearch --derep_fulllength temp/filtered.fa \
      --minuniquesize 10 --sizeout --relabel Uni_ \
      --output temp/uniques.fa &&

### 4.2 聚类OTU/去噪ASV Cluster or denoise
# 97%聚类OTU，适合大数据/ASV规律不明显/reviewer要求
  # 重命名relabel、按相似id=97%聚类，不屏蔽qmask
  # 记录输入sizein和输出频率sizeout
    vsearch --cluster_size temp/uniques.fa  \
     --relabel OTU_ --id 0.97 \
     --qmask none --sizein --sizeout \
     --centroids temp/otus.fa &&

    sed -i 's/;.*//' temp/otus.fa &&

### 4.3 基于参考去嵌合 Reference-based chimera detect


    # silva_16s_v123.fa
    vsearch --uchime_ref temp/otus.fa \
      -db /home/hppiserver/amplicondb/silva_16s_v123.fa \
      --nonchimeras result/raw/otus.fa &&

    sed -i 's/\r//g' result/raw/otus.fa


### 5.1 生成特征表

    # id(0.97)：97%相似度比对序列
    time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 0.97 --threads 48 \
    	--otutabout result/raw/otutab.txt &&

    # vsearch结果windows用户删除换行符^M校正为标准Linux格式
    sed -i 's/\r//' result/raw/otutab.txt &&

    # csvtk统计表行列
    # 这里一定看好列数，是不是等于你的样品数；如果不等，一般是样品命名存在问题，具体看上面解释
    csvtk -t stat result/raw/otutab.txt &&


### 5.2 物种注释，且/或去除质体和非细菌 Remove plastid and non-Bacteria

    # SILVA数据库(silva_16s_v123.fa)更好注释真核、质体序列，极慢耗时3h起
    # 置信阈值通常0.6/0.8，vserch最低0.1/usearch可选0输出最相似物种注释用于观察潜在分类
 
    time vsearch --sintax result/raw/otus.fa \
      --db /home/hppiserver/amplicondb/silva_16s_v123.fa  \
      --sintax_cutoff 0.1 \
      --tabbedout result/raw/otus.sintax &&
    #less -SN result/raw/otus.sintax
    sed -i 's/\r//' result/raw/otus.sintax
    
    #原始特征表行数筛选
   wc -l result/raw/otutab.txt
    #R脚本选择细菌古菌(真核)、去除叶绿体、线粒体并统计比例；输出筛选并排序的OTU表
    #输入为OTU表result/raw/otutab.txt和物种注释result/raw/otus.sintax
    #输出筛选并排序的特征表result/otutab.txt和
    #统计污染比例文件result/raw/otutab_nonBac.txt和过滤细节otus.sintax.discard
    Rscript /home/hppiserver/script/otutab_filter_nonBac.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonBac.stat \
      --discard result/raw/otus.sintax.discard &&

    # 筛选后特征表行数

    cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
    
    #过滤特征表对应序列
    time usearch -fastx_getseqs result/raw/otus.fa \
        -labels result/otutab.id -fastaout result/otus.fa &&
    #过滤特征表对应序列注释
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/raw/otus.sintax result/otutab.id \
        > result/otus.sintax &&


### 5.3 等量抽样标准化

    # Normlize by subsample

    #使用vegan包进行等量重抽样，输入reads count格式Feature表result/otutab.txt
    #可指定输入文件、抽样量和随机数，输出抽平表result/otutab_rare.txt和多样性alpha/vegan.txt
    mkdir -p result/alpha
    Rscript /home/hppiserver/script/otutab_rare.R --input result/otutab.txt \
      --depth 100 --seed 1 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt &&
    time usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat &&
    #cat result/otutab_rare.stat


## 8. 物种注释分类汇总

    #OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
      > result/taxonomy2.txt &&
    #less -SN result/taxonomy2.txt

    #OTU对应物种8列格式：注意注释是非整齐
    #生成物种表格OTU/ASV中空白补齐为Unassigned
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax &&
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt &&
    #less -SN result/taxonomy.txt
    #生成相对丰度表
     sed -i 's/#//g' result/otutab.txt
     mkdir -p result/RA
     Rscript /home/hppiserver/script/Different_taxonomy_level_RA.r --input result/otutab.txt \
      --taxonomy result/taxonomy.txt \
      --output result/RA/ &&

## 10. 空间清理及数据提交

    #删除中间大文件
    rm -rf temp/*.fq &&

    # 分双端统计md5值，用于数据提交
    cd rawdata
    md5sum *.1.fastq.gz  > md5sum1.txt
    md5sum *.2.fastq.gz  > md5sum2.txt
    paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | sed 's/*//g' > ../result/md5sum.txt
    rm md5sum*
    cd ..
    cat result/md5sum.txt

tar -cvf result.tar.gz result
