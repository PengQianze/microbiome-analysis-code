############################################SNP phylogenetic tree
#phliy

seqboot<seqboot.par &&mv outfile seqboot.out &&
dnadist<dnadist.par &&  mv outfile dnadist.out && 
neighbor<neighbor.par && mv  outfile nei.out && mv outtree nei.tree 

less dnadist.out | tr '\n' '|'| sed 's/| / /g' | tr '|' '\n' >infile.dist.table
less nei.tree  | tr '\n' ' '|sed 's/ //g' > outtree.nwk

#生成无根树
phylip retree W U Q 

####################################microbiome phylogenetic tree

## 1. 筛选高丰度/指定的特征

    #按相对丰度0.02%筛选高丰度OTU
    usearch -otutab_trim ../otutab.txt \
        -min_otu_freq 0.0002 \
        -output otutab.txt &&
    #统计筛选OTU表特征数量
    tail -n+2 otutab.txt | wc -l

    #修改特征ID列名
    sed -i '1 s/#OTU ID/OTUID/' otutab.txt
    #提取ID用于提取序列
    cut -f 1 otutab.txt > otutab_high.id

    # 筛选高丰度菌/指定差异菌对应OTU序列
    usearch -fastx_getseqs ../otus.fa -labels otutab_high.id \
        -fastaout otus.fa &&
    head -n 2 otus.fa

## 2. 构建进化树

    # Muscle软件进行序列对齐，3s
    muscle -in otus.fa -out otus_aligned.fas &&

    # FastTree快速建树(Linux)
    fasttree -gtr -nt otus_aligned.fas > otus.nwk 
