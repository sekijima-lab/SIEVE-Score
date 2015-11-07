このスクリプトは、glideの出力をパースし、PCAとk-meansクラスタリングを行い、その結果を可視化・スコアリングするスクリプトです。

実行環境はMac/Linuxを想定しています。また、maestroがインストールされており、$SCHRODINGERが正しく設定されている必要があります。個別のライセンスは必要ありません。

また、glideによるドッキングの際にoutputタブの''Write per-residue interaction scores''オプションをONにする必要があります。

実行方法：

$SCHRODINGER/run python ./glide_kmeans.py -i input.maegz -o output.csv -hits [hits.txt | None] options

options ()内はデフォルト：
  -cl: K-meansのクラスタ数(5)
  -n : hits/non-hits周囲のn化合物に加点/減点(100)
  -title: グラフのタイトル(-iの値)
  -p : 加点の値(int,1)
  -m : 減点の値(int,1)
  -propose: 書き出す化合物数(1000)

  -p,-mの値は合計点が0になるように, -m -1 -p (#non-hits)/(#hits)と設定するのを推奨

Output:
1. output.csv
    スコア上位の化合物
2. output_loading.csv
    PCAの重み
3. input_cluster.maegz
    入力にクラスタ番号とスコアを追記したもの
4. input_cluster.csv
    クラスタリングとPCAの結果
4. output.png
    PCA/クラスタリングの可視化結果
5. output.log
    オプションの設定(再現用)

hits.txtの書き方
1行に1化合物を記入する
hit化合物はname,1、non-hit化合物はname,-1のように記入する
正の値のものがhit、負の値のものがnon-hitと認識され、#で始まる行はコメントとみなされる
スコアに重みを付ける場合は、name,2やname,-3など、倍率を入力する

hitの指定数にもよりますが、10万リガンドで30分程度かかります。
現バージョンでは50万リガンド以上での使用は非推奨です。