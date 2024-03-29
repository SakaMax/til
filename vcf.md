# VCFの取り扱いまとめ

## fastq -> gvcfの作成

### NCBI(SRX)からのfastqの取得

1. ncbiのAPI keyをとっておく
2. NCBI eutilsでDRRの名前を取得
3. sra toolsでsraファイルの回収
4. sra toolsでfastqの抽出

コード例

- 遺伝研のsingularityで動かす場合
- 遺伝研のbiotoolsにはsra toolsがなかった（このスクリプトを作った当時）のでNCBIのdockerイメージからsifを作成
  - [NCBIのdockerイメージ](https://hub.docker.com/r/ncbi/sra-tools)
  - [sifの作り方](https://sc.ddbj.nig.ac.jp/software/Apptainer/)
- `SRX2764015 8 32`のような引数をつけて走らせると良い
  - [これ](https://www.ncbi.nlm.nih.gov/sra/SRX2764015)のfastqが手に入る

<details>

<summary>SRX_to_fastq.sh</summary>

```bash
#$ -S /usr/bin/bash
#$ -cwd
#$ -l s_vmem=4G
#$ -l mem_req=4G
#$ -pe mpi 8
#$ -l short
# arg 1: DRX query
# arg 2: num. of thread
# arg 3: amount of memory in Gigabase

api_key="WRITE YOUR API KEY HERE"
esearch_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
efetch_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# search UID of DRX
uid=$(wget -O - "${esearch_url}?db=sra&term=${1}&usehistory=y&api_key=${api_key}" | \
    sed -n -r 's!^.*<Id>(.+)</Id>.*$!\1!p')
echo "UID in the SRA database: ${uid}"

# get DRR name
drr=$(wget -O - "${efetch_url}?db=sra&id=${uid}&usehistory=y&api_key=${api_key}"  | \
     sed -n -r 's!^.*<PRIMARY_ID>([DS]{1}RR[0-9]+)</PRIMARY_ID>.*$!\1!p')
echo "Run file: ${drr}"

echo $(pwd)

# obtain sra file
echo "prefetch"
singularity exec --bind $(pwd) ~/containers/sra-tools.sif prefetch \
 -v -p --max-size 50000000 -O . ${drr}

# move sra file to here and remove directory produced by prefetch
mv "${drr}"/"${drr}.sra" ./
rm -r "${drr}"

# extract fastq
echo "fasterq-dump"
singularity exec --bind $(pwd) ~/containers/sra-tools.sif fasterq-dump \
 -v -e ${2:-8} -m "${3:-30}G" "${drr}.sra" \
 && rm "${drr}.sra"

mv "${drr}_1.fastq" "${1}_1.fastq"
mv "${drr}_2.fastq" "${1}_2.fastq"
```

</details>

### フィルタリング

fastp

### マッピング fastq -> sam -> bam -> bam(ソート済み)

- bwa-mem
  - GATKで処理するために`-R "@RG\tID:${RGID}\tSM:${RGSM}\tPL:${RGPL}\tLB:${RGLB}"`をつけておく必要がある
  - RGIDとRGSMはDRXの番号, RGPLは`"Illumina"`, RGLBは`DRXの番号_1`などとする[^1]
- samtools view
- samtools sort

[^1]: `-R`オプションに関する[参考リンク](https://hashiyuki.hatenablog.com/entry/2016/05/07/164740)

### バリアントコール sorted bam -> gvcf

GATK関連のツールを遺伝研スパコンで動かすときはスクリプトに`#$ -v MALLOC_ARENA_MAX=2`が必要

- FixMateInformation
- MarkDuplicatesSpark
- HaplotypeCaller
  - `-ERC GVCF`を引数につける
  - 十分なスレッド数が使えるときは`--native-pair-hmm-threads 10`をつけるのを推奨（根拠となる調査がどっかにあった）

## gvcfからvcfを作る

g.vcf.gzのあるディレクトリで作業する

- GenomicsDBImport

 1. mapファイルの作成 `ls | egrep 'g.vcf.gz$' | sed -r 's/([A-Za-z0-9]+)(\.g\.vcf\.gz$)/\1\t\1\2\t\1\2.tbi/g' > cohort.map`
 2. chromosome.txt (例) `pv sample1.g.vcf.gz | gunzip -c | cut -f 1 | egrep -v -e '^#' -e 'scaffold' | sort | uniq | tee chrom.txt`
 3. GenomicsDBImport.shの作成
 4. インデックスの作成 `parallel -j $(nproc) 'tabix {}' ::: $(ls | grep g.vcf.gz)`
 5. `parallel -j $(nproc) -v --result results2 './GenomicsDBImport.sh $(sed -n {}p chromosome.txt) cohort.map'  ::: $(seq 1 20) &`

- GenotypeGVCFs
  - `parallel -j 4 -v --result result_gmax 'L=$(sed -n {}p chromosome.txt) && echo "Chromosome ${L}" && gatk --java-options "-Xmx8g" GenotypeGVCFs  -R Wm82.a4.v1/Gmax_508_v4.0.softmasked.fa -V "gendb://${L}" -O "${L}_gmax.vcf.gz"' ::: $(seq 1 20) &`
  - メモリを食うのでジョブ数とXmxはお好みで調整

mapファイル

```text
系統名 系統.g.vcf.gz
```

chromosome.txt

```text
(gvcfの#CHROM列に登場する染色体名)
```

GenomicsDBImport.sh

```bash
gatk --java-options "-Xmx30g" GenomicsDBImport \
        --genomicsdb-workspace-path ./"$1" \
        -L "$1" \
        --sample-name-map "$2"
```

既存のgenomicsDBに追加したいときは`--genomicsdb-workspace-path`を`--genomicsdb-update-workspace-path`に置き換え, `-L "$1" \`を消す

## vcfのハードフィルタリング

参照: [gatkのガイド](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering)

カレントディレクトリにjoint call した*.vcf.gzがあるとする。

1. ディレクトリの作成 `mkdir -p tmp/snp tmp/indel tmp/snp-filtered tmp/indel-filtered`
2. SNPの抽出 `parallel -j $(nproc) --result result-selectSNP 'gatk SelectVariants -V {} -select-type SNP -O tmp/snp/{}' ::: *.vcf.gz`
3. indelの抽出 `parallel -j $(nproc) --result result-selectIndel 'gatk SelectVariants -V {} -select-type INDEL -O tmp/indel/{}' ::: *.vcf.gz`
4. SNPのフィルタ `parallel -j $(nproc) --result result-filterSNP --compress --tmp-dir parallel-tmp 'bash filter-snp.sh {} > {}-snp.log 2>&1' ::: *.vcf.gz`
5. indelのフィルタ `parallel -j $(nproc) --result result-filterIndel --compress --tmp-dir parallel-tmp 'bash filter-indel.sh {} > {}-indel.log 2>&1' ::: *.vcf.gz`

<details>
<summary>filter-snp.sh</summary>

 ```bash
#!/bin/bash

# Path to files
BEAGLE="/usr/local/bin/beagle.22Jul22.46e.jar"
IN="tmp/snp/${1}"
TMP="tmp/snp/${1}_tmp.vcf"
RAW="tmp/snp-filtered/${1%%.*}.raw"
OUT="tmp/snp-filtered/${1%%.*}.out"

# Remove Spanning or overlapping deletions
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035531912-Spanning-or-overlapping-deletions-allele-
zcat "${IN}" | \
sed '/*/d' | \

# Annotate SNPs
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' > "${TMP}" && \

# Filter SNPs (hard-filter)
# see https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# and https://gatk.broadinstitute.org/hc/en-us/articles/360035531012--How-to-Filter-on-genotype-using-VariantFiltration
# Mark heterozygous genotype call as failure because our samples are pure lines.
gatk VariantFiltration\
    -V "${TMP}" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --genotype-filter-expression "isHet == 1" --genotype-filter-name "isHetFilter" \
    --genotype-filter-expression "DP < 10" --genotype-filter-name "DP10" \
    -O "${RAW}.vcf.gz" && \

# Remove SNPs and genotypes without PASS
gatk SelectVariants \
    -V "${RAW}.vcf.gz" \
    --exclude-filtered \
    --set-filtered-gt-to-nocall \
    -O "${RAW}.pass.vcf.gz" && \

# Remove SNPs with call rate < 10% and maf < 5%
plink2 --vcf "${RAW}.pass.vcf.gz" \
    --geno 0.1 \
    --maf 0.05 \
    --out "${OUT}" \
    --export vcf \
    --allow-extra-chr \
    --threads 1 && \

# Compress and index
bcftools view -Oz -o "${OUT}.vcf.gz" "${OUT}.vcf" && \
tabix -p vcf "${OUT}.vcf.gz" && \

# Impute
java -Xmx8g -jar "${BEAGLE}" \
    gt="${OUT}.vcf.gz" \
    out="${OUT}.imputed" \
    impute=true \
    gp=true \
    nthreads=1 && \

# Index imputed vcf
tabix "${OUT}.imputed.vcf.gz" && \

# Remove temp file
rm "${TMP}"

 ```


</details>

<details>
<summary>filter-indel.sh</summary>

```bash
#!/bin/bash

# Path to files
BEAGLE="/usr/local/bin/beagle.22Jul22.46e.jar"
IN="tmp/indel/${1}"
TMP="tmp/indel-filtered/${1%%.*}_tmp.vcf"
RAW="tmp/indel-filtered/${1%%.*}.raw"
OUT="tmp/indel-filtered/${1%%.*}.out"

# Remove Spanning or overlapping deletions
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035531912-Spanning-or-overlapping-deletions-allele-
zcat "${IN}" | \
sed '/*/d' | \

# Annotate SNPs
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' > "${TMP}" && \

# Filter SNPs (hard-filter)
# see https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# and https://gatk.broadinstitute.org/hc/en-us/articles/360035531012--How-to-Filter-on-genotype-using-VariantFiltration
# Mark heterozygous genotype call as failure because our samples are pure lines.
gatk VariantFiltration \
    -V "${TMP}" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --genotype-filter-expression "isHet == 1" --genotype-filter-name "isHetFilter" \
    --genotype-filter-expression "DP < 10" --genotype-filter-name "DP10" \
    -O "${RAW}.vcf.gz" && \

# Remove SNPs and genotypes without PASS
gatk SelectVariants \
    -V "${RAW}.vcf.gz" \
    --exclude-filtered \
    --set-filtered-gt-to-nocall \
    -O "${RAW}.pass.vcf.gz" && \

# Remove SNPs with call rate < 10% and maf < 5%
plink2 --vcf "${RAW}.pass.vcf.gz" \
    --geno 0.1 \
    --maf 0.05 \
    --out "${OUT}" \
    --export vcf \
    --allow-extra-chr \
    --threads 1 && \

# Compress and index
bcftools view -Oz -o "${OUT}.vcf.gz" "${OUT}.vcf" && \
tabix -p vcf "${OUT}.vcf.gz" && \

# Impute
java -Xmx8g -jar "${BEAGLE}" \
    gt="${OUT}.vcf.gz" \
    out="${OUT}.imputed" \
    impute=true \
    gp=true \
    nthreads=1 && \

# Index imputed vcf
tabix "${OUT}.imputed.vcf.gz" && \

# Remove temp file
rm "${TMP}"
```

</details>
