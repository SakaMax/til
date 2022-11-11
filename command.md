# シェルコマンド

## parallel

### ファイルの行単位で並列処理

`parallel -j $(nproc) 'echo $(sed -n {}p arguments.txt)' ::: $(seq 1 10)`

- `arguments.txt`が一行ずつ表示される
- ただし順不同に並列実行されるため出力順はファイルの上から順とは異なる

## sed

## マッチ部分を参照する

`echo "hoge.g.vcf.gz" | sed -r 's/(^[A-Za-z0-9]+)(\.g\.vcf\.gz$)/\1 \1\2/g'` → `hoge hoge.g.vcf.gz`

- 参照したい部分を`()`でくくる
- 参照するときは`\番号`で参照する

## tr

### tqdmのログとかに含まれる`^M`を改行に変換する

`cat Hoge.seq.e | tr '\r' '\n'`
