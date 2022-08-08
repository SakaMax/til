# VCFの取り扱いまとめ
## vcfのハードフィルタリング
参照: [gatkのガイド](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering)

カレントディレクトリにjoint call した*.vcf.gzがあるとする。
1. ディレクトリの作成 `mkdir -p tmp/snp tmp/indel tmp/snp-filtered tmp/indel-filtered`
2. SNPの抽出 `parallel -j $(nproc) --result result-selectSNP 'gatk SelectVariants -V {} -select-type SNP -O tmp/snp/{}' ::: *.vcf.gz`
3. indelの抽出 `parallel -j $(nproc) --result result-selectIndel 'gatk SelectVariants -V {} -select-type INDEL -O tmp/indel/{}' ::: *.vcf.gz`
4. SNPのフィルタ 
```{bash}
cat << "EOF" > filter-snp.sh 
IN="tmp/snp/${1}"
TMP="tmp/snp-filtered/${1%%.*}.raw.vcf.gz"
OUT="tmp/snp-filtered/${1}"

gatk VariantFiltration\
    -V ${IN} \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${TMP} && \
vcftools --gzvcf ${TMP} \
  --recode --stdout  --remove-filtered-all | \
bcftools view -Oz - -o ${OUT} && \
tabix -p vcf ${OUT} && \
rm ${TMP} ${TMP}.tbi
EOF

parallel -j $(nproc) --result result-filterSNP 'bash filter-snp.sh {}' ::: *.vcf.gz
```
5. indelのフィルタ
```{bash}
cat << "EOF" > filter-indel.sh
IN="tmp/indel/${1}"
TMP="tmp/indel-filtered/${1%%.*}.raw.vcf.gz"
OUT="tmp/indel-filtered/${1}"

gatk VariantFiltration \
    -V ${IN} \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O ${TMP} && \
vcftools --gzvcf ${TMP} \
  --recode --stdout  --remove-filtered-all | \
bcftools view -Oz - -o ${OUT} && \
tabix -p vcf ${OUT} && \
rm ${TMP} ${TMP}.tbi
EOF

parallel -j $(nproc) --result result-filterIndel 'bash filter-indel.sh {}' ::: *.vcf.gz
```
