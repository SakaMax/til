# blast
- blastの検索結果はsamで出力できる
  - outfmt "17 SQ" (SQがないと塩基配列がsamに入らない）
  - sed 's/Query_1/(リファレンス名)/g' blast.sam > blast.sed.sam
  - あとはigvtoolでsort, index
- アメリカのオフィスアワーにでかいblastをNCBIでやると確かに重い

# vcf関連
- vcfのサンプル名を修正したいときは`bcftools reheader`を使う
  - `bcftools reheader -s "サンプル名を行区切りで書いたテキストファイル" -o out.vcf in.vcf`

# NCBI E-utilities
- pythonでアクセスするならBiopythonの`Entrez`を使うと楽。 
  - efetchするなら`handle = Entrez.efetch(db="nuccore", id="ACCESSION", retmode="xml")｀ して `Entrez.read(handle)`でデータが取り出せる。

# シェル周り
  ## tr
  - tqdmのログとかに含まれる`^M`を改行に変換する: `cat Hoge.seq.e | tr '\r' '\n'`
