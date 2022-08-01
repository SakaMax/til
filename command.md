# ファイルの行単位で並列処理 (parallel)
`parallel -j $(nproc) 'echo $(sed -n {}p arguments.txt)' ::: $(seq 1 10)`

- `arguments.txt`が一行ずつ表示される
- ただし順不同に並列実行されるため出力順はファイルの上から順とは異なる
