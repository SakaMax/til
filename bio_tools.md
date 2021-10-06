# blast
- blastの検索結果はsamで出力できる
  - outfmt "17 SQ" (SQがないと塩基配列がsamに入らない）
  - sed 's/Query_1/(リファレンス名)/g' blast.sam > blast.sed.sam
  - あとはigvtoolでsort, index
- アメリカのオフィスアワーにでかいblastをNCBIでやると確かに重い
