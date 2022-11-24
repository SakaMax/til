# シェルコマンド

## parallel

### ファイルの行単位で並列処理

`parallel -j $(nproc) 'echo $(sed -n {}p arguments.txt)' ::: $(seq 1 10)`

- `arguments.txt`が一行ずつ表示される
- ただし順不同に並列実行されるため出力順はファイルの上から順とは異なる

### 重い処理でエラーになって落ちるとき

`--compress --tmpdir ./parallel-tmp`オプションをつける

## pv

`pipe view`の略らしい

### 進捗を見ながら圧縮

[参考](https://qiita.com/xkumiyu/items/9df519d767f921e29d20)

<details>
<summary>tarc-pv.sh</summary>

```bash
#!/bin/bash
#$1 target
tar -cf - ${1} | pv -c --name ${1} -s $(du -sb ${1} | awk '{print $1}')
```

</details>

上記のスクリプトを使って
`./tarc-pv.sh hoge | pbzip2 -9 -c -p6 -m500 > hoge.tar.bz2`

圧縮先のデータの流れも見たいなら
`./tarc-pv.sh hoge | xz -9 -c -T 6 -M 12G | pv -c --name xz > hoge.tar.xz`

## sed

## マッチ部分を参照する

`echo "hoge.g.vcf.gz" | sed -r 's/(^[A-Za-z0-9]+)(\.g\.vcf\.gz$)/\1 \1\2/g'` → `hoge hoge.g.vcf.gz`

- 参照したい部分を`()`でくくる
- 参照するときは`\番号`で参照する

## tr

### tqdmのログとかに含まれる`^M`を改行に変換する

`cat Hoge.seq.e | tr '\r' '\n'`
