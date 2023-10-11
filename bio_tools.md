# bioinfo系ツールのメモ

## blast

- blastの検索結果はsamで出力できる
  - outfmt "17 SQ" (SQがないと塩基配列がsamに入らない）
  - sed 's/Query_1/(リファレンス名)/g' blast.sam > blast.sed.sam
  - あとはigvtoolでsort, index
- アメリカのオフィスアワーにでかいblastをNCBIでやると確かに重い
- -outfmt 7で出力したファイルからbitscoreが最大の部分のみを取り出す
  - `grep -v "#" | awk 'BEGIN {asv="begin"; bitscore=0;} ($1 != asv || $12 >= bitscore) {print $0; asv=$1; bitscore=$12+0;}'`

## vcf関連

- vcfのサンプル名を修正したいときは`bcftools reheader`を使う
  - `bcftools reheader -s "サンプル名を行区切りで書いたテキストファイル" -o out.vcf in.vcf`

## NCBI E-utilities

- pythonでアクセスするならBiopythonの`Entrez`を使うと楽。
  - Genbank形式のデータ: `handle = Entrez.efetch(db="nuccore", id="ACCESSION", retmode="xml")` して `Entrez.read(handle) ->list[dict]`
  - fasta形式のデータ: `handle = Entrez.efetch(db=db, id=ids, retmode="text", rettype="fasta")` して `handle.read() ->str`

## 変換

- sam/bam -> fasta `samtools fasta input.bam > input.fa`

## GAPITの引数一覧
```
function (Y = NULL, G = NULL, GD = NULL, GM = NULL, KI = NULL, 
    Z = NULL, CV = NULL, Aver.Dis = 1000, buspred = FALSE, bin.from = 10000, 
    bin.to = 10000, bin.by = 10000, cutOff = 0.05, effectunit = 1, 
    file.output = TRUE, FDRcut = FALSE, group.from = 1e+06, group.to = 1e+06, 
    group.by = 50, Geno.View.output = TRUE, h2 = NULL, inclosure.from = 10, 
    inclosure.to = 10, inclosure.by = 10, Inter.Plot = FALSE, 
    Inter.type = c("m", "q"), kinship.cluster = "average", kinship.group = "Mean", 
    kinship.algorithm = "Zhang", lmpred = FALSE, model = "MLM", 
    maxOut = 100, memo = NULL, Model.selection = FALSE, Multi_iter = FALSE, 
    Major.allele.zero = FALSE, Multiple_analysis = TRUE, num_regwas = 10, 
    N4 = FALSE, NQTN = NULL, N.sig = NULL, NJtree.group = NULL, 
    NJtree.type = c("fan", "unrooted"), output.numerical = FALSE, 
    output.hapmap = FALSE, QC.Y = FALSE, QTN.position = NULL, 
    QTNDist = "normal", r = 0.25, Random.model = TRUE, sangwich.top = NULL, 
    sangwich.bottom = NULL, seq.cutoff = NULL, SNP.P3D = TRUE, 
    SNP.effect = "Add", SNP.impute = "Middle", SNP.fraction = 1, 
    SNP.test = TRUE, SNP.MAF = 0, SNP.FDR = 1, testY = NULL, 
    plot.bin = 10^5, PCA.total = 0, PCA.col = NULL, PCA.3d = FALSE, 
    PCA.legend = NULL, PCA.View.output = TRUE, Phenotype.View = TRUE, 
    Predict.type = "GEBV", WS = c(1, 1000, 10000, 1e+05, 1e+06, 
        1e+07), WS0 = 10000) 
```

GAPITの[ソース](https://github.com/jiabowang/GAPIT/blob/f0291d2d42a3d0b5f08146033215aed56edde519/R/GAPIT.R)にある説明：

```
#' @param Y  data.frame of phenotype data where each row is a sample and each column is a trait, the first column is the sample names
#' @param G  data.frame of genotypic data in HAPMAP format
#' @param GD data.frame of genetic data in numerical format, where each row is a sample and each column is a variant.
#' @param GM a data.frame of genomic coordinates for the genetic map
#' @param KI an $NxN$ matrix of kinship coefficients
#' @param Z  an $NxN$ (for MLM) or an $NxN`$ (CMLM) matrix of index, which is made with 0 and 1 value to indicate indivdual belong to each group.
#' @param CV Covariance matrix
#' @param testY data.frame of phenotype data in testing population, where each row is a sample and each column is a trait, the first column is the sample names.
#' @param group.from integer, minimum number of group(s) to consider in CMLM
#' @param group.to integer, maximum number of group(s) to consider in CMLM
#' @param group.by integer, increment for evaluating group size in CMLM
#' @param kinship.cluster algorithm for calculating kinship centroid (options: "average", "complete", "ward", "single", "mcquitty", "median", and "centroid") 
#' @param kinship.group method for calculating group membership (options: "Mean", "Max", "Min", and "Median")
#' @param kinship.algorithm algorithm to calculate the kinship matrix (options: "VanRaden", "EMMA", "Loiselle", and "Zhang")
#' @param buspred logical, option for prediction after GWAS。
#' @param lmpred logical (vector), option for seletion of linear model prediction or (and) ABLUP.
#' @param FDRcut logical, filter pseudo QTN based on FDR cut-off in BLINK
#' @param bin.from integer, minimum number of bin(s) to consider in SUPER
#' @param bin.to integer, maximum number of bin(s) to consider in SUPER
#' @param bin.by integer, increment for evaluating bin size in SUPER
#' @param inclosure.from integer, minimum number of pesudo QTNs to consider in SUPER
#' @param inclosure.to integer, maximum number of pesudo QTNs to consider in SUPER
#' @param inclosure.by integer, increment for evaluating number of pesudo QTNs in SUPER
#' @param SNP.P3D logical, to use P3D or Not for Testing SNPs
#' @param SNP.effect genetic model for coding the SNP effect (options: "Add" (additive), "Dom", "Left", and "Right")
#' @param SNP.impute SNP imputation method (options: "Middle", "Major", and "Minor")
#' @param PCA.total integer, number of principal components to include in Q matrix (can be zero)
#' @param SNP.fraction numerical input between 0 and 1, fraction of SNPs Sampled to Estimate Kinship and PCs
#' @param SNP.MAF numerical input between 0 and 1, minor allele frequency to filter SNPs in GWAS reports
#' @param SNP.FDR numerical input between 0 and 1, false discovery rate for filtering SNPs
#' @param PCA.col list for points color in PCA plot. The total length of PCA.col should be equal to the number of individuals in the GD or G file.
#' @param PCA.3d logical, whether output 3D PCA plot.
#' @param NJtree.group numeric, set the number of clustering groups in the NJtree plot.
#' @param NJtree.type type of neighbor joining tree (options: "fan" and "unrooted")
#' @param sangwich.top Model type to run in the first iteration of SUPER, (options: "MLM", "GLM", "CMLM","Fast-LMM")
#' @param sangwich.bottom Model type to run in the last iteration of SUPER, (options: "MLM", "GLM", "CMLM","Fast-LMM")
#' @param file.output logical, whether output all result files.
#' @param cutOff numeric value, the threshold for filtering significant markers from all. It would be transfor as Bornferrni cutoff in Manhattan plots.
#' @param Model.selection logical, whether evaluate optimum number of CV file. If TRUE, all likelyhood values should be evaluated for each CV combination.
#' @param output.numerical logical,whether output numerical genotype file from HapMap file.
#' @param output.hapmap logical,whether output numerical HapMap file from numerical genotype file.
#' @param Multi_iter logical, whether add more iterations for FarmCPU, BLINK.
#' @param num_regwas numeric, the maximum number of selective significant markers into re-GWAS model.
#' @param Major.allele.zero logical, whether set major allele as 0, and minor allele as 2, if FALSE, they will be set as reverse.
#' @param Random.model logical, whether ran random model to estimate PVE values for significant markers after GWAS.
#' @param memo text, from users to remark for output files.
#' @param Inter.Plot logical, whether to output the interactive Manhattan and QQ plots.
#' @param Inter.type Interactive plot type for Manhattan and QQ plots."m" indicate manhattan plot and "q" indicate QQ plot.
#' @param WS numeric or numeric vector, the distance between detected markers and real QTN should be recognized as a real power.
#' @param WS0 numeric, the cutoff threshold for distance between markers to display in GAPIT.Genotype.Distance_R_Chro.pdf file.
#' @param Aver.Dis=1000 numeric, average display windowsize in LD decay plot,
#' @param maxOut numeric, set the number of markers in the power calculation, the top maxOut number of P values markers should be selected.
#' @param QTN.position numeric vector, set where are the QTNs' position. Its maximun values should be equal to total marker number, and its length should be equal to the NQTN.
#' @param PCA.View.output logical, whether to output the PCA view
#' @param Geno.View.output logical whether to output the Genotype analysis including MAF, heterzygosity, LD decay, and other genotype distribution output.
#' @param h2 numeric value, to set simulation phenotype heritability. It ranged from 0 to 1 means 0% to 100%.
#' @param NQTN numeric value, to set simulation number of QTN. It ranged from 1 to the total markers number.
#' @param QTNDist option for distribution of simulated QTN genetic effect in the simulation,(options: "normal" and "geometry")
#' @param effectunit numeric value, the effect unit of the first choosed marker in the simulation pheotype. default as 1
#' @param Multiple_analysis logical, whether to output the mulitple mahattan and QQ plots. default as TRUE 
#' @param model model type to run, (options: "MLM", "GLM", "CMLM", "MMLM", "SUPER", "FarmCPU", "gBLUP",  "cBLUP", and "sBLUP"
#' @param Predict.type option to display which type predicted factor again real phenotype in the GAPIT.Association.Prediction pdf file.(options: "GEBV","BLUP" and "BLUE")
#' @param SNP.test logical, whether to do GWAS or GS.
#' @param seq.cutoff numeric value, the threshold for filtering significant markers from all. It would be transfor as Bornferrni cutoff in GGS.
#' @details 
```
