sample=${1}
#tar -xvzf ${sample}.tar.gz
wigToBigWig ${sample} ~/WORK/2.genome/hg19.chrom.sizes ${sample}.bw
