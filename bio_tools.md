# blast
- blastの検索結果はsamで出力できる
  - outfmt "17 SQ" (SQがないと塩基配列がsamに入らない）
  - sed 's/Query_1/(リファレンス名)/g' blast.sam > blast.sed.sam
  - あとはigvtoolでsort, index
- アメリカのオフィスアワーにでかいblastをNCBIでやると確かに重い
- -outfmt 7で出力したファイルからbitscoreが最大の部分のみを取り出す 
  - `grep -v "#" | awk 'BEGIN {asv="begin"; bitscore=0;} ($1 != asv || $12 >= bitscore) {print $0; asv=$1; bitscore=$12+0;}'`

# vcf関連
- vcfのサンプル名を修正したいときは`bcftools reheader`を使う
  - `bcftools reheader -s "サンプル名を行区切りで書いたテキストファイル" -o out.vcf in.vcf`

# NCBI E-utilities
- pythonでアクセスするならBiopythonの`Entrez`を使うと楽。 
  - Genbank形式のデータ: `handle = Entrez.efetch(db="nuccore", id="ACCESSION", retmode="xml")` して `Entrez.read(handle) ->list[dict]`
  - fasta形式のデータ: `handle = Entrez.efetch(db=db, id=ids, retmode="text", rettype="fasta")` して `handle.read() ->str`

# 変換
- sam/bam -> fasta `samtools fasta input.bam > input.fa`
