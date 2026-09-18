# Bianchi IX 複数時刻履歴・デコヒーレンス監査

## 目的

第三壁までの監査では、各時刻の A/B+/B- セクター確率質量を計算した。しかし一時刻の確率から
「どの壁を順に選んだか」の確率を作ることはできない。

今回は同じ縮約 Bianchi IX Hamiltonian と同じ2次元平方根量子化を変えず、三つの時計時刻

- s1 = 4.6（第一壁）
- s2 = 10.385（第二壁）
- s3 = 32.46（第三壁）

で壁セクター射影 P_A, P_B+, P_B- を挿入する。

第一壁では既存計算で A セクターが約98%なので、主解析は A に条件付けた9履歴

C_(A,j,k) = P_k(s3) U(s3,s2) P_j(s2) U(s2,s1) P_A(s1) U(s1,s0)

を扱う。第一壁の非A成分は捨てたことを隠さず、別の欠落質量として報告する。完全な27履歴は今回の主計算ではない。

## デコヒーレンス汎関数

各クラス演算子の branch state を

|psi_alpha> = C_alpha |psi0>

として

D(alpha,beta) = <psi_beta|psi_alpha>

を直接計算する。

D の対角要素は「履歴重み」だが、非対角要素が無視できないとき通常の加法的確率として解釈しない。
特に最終壁 k だけを指定した粗視化事象について、

|| sum_j |psi_(A,j,k)> ||^2

と

sum_j ||psi_(A,j,k)||^2

の差を additivity defect として保存する。差があれば、第二壁を細かく区別した履歴の干渉が残る。

## 数値上の重要点

長時間版は s=6,12,20 で有限箱を拡大する。従来の regrid_state は全波束を再正規化するが、
これを各 branch に適用すると branch の相対振幅を破壊する。

この監査では同じ cubic interpolation を使うが、再格子化は線形・非正規化のまま適用する。
全9 branch を最後に足し戻した状態と、中間射影を入れずに A branch を進めた状態の fidelity を検査する。

主設定は従来の長時間版と同じ grid / potential cap / time step / Krylov 次数を用いる。
第三区間について Krylov 次数60の対照も取り、デコヒーレンス行列の数値依存を記録する。

## 先行研究との区別

decoherent / consistent histories では、複数時刻の class operator と decoherence functional を使い、
履歴に確率を割り当てられる条件を調べる。量子宇宙論では Wheeler-DeWitt 型模型に対する class operator
の構成も研究されている（例: Halliwell, arXiv:1108.5991; arXiv:0909.2597）。

今回の有限時計 s と逐次射影によるクラス演算子は、選んだ時間発展を持つ縮約模型に対する直接的な有限時間監査である。
timeless WDW 制約と可換な不変 class operator を構成したものではない。

## 判定

CIで実計算した summary.json を成果物として保存する。次を分離して報告する。

1. branch の対角重み
2. decoherence functional の最大非対角要素
3. 正規化 pair coherence
4. coarse/fine additivity defect
5. branch 再結合 closure
6. Krylov 次数対照
7. 第一壁で条件付けにより除外した非A質量

数値的に非対角成分が十分小さくない場合、P(A→B-→B+) のような量を通常の履歴確率と呼ばない。

## 限界

- 一様真空 Bianchi IX minisuperspace の限定模型。
- 一つの平方根量子化・L2表現・内部時計 s=-alpha。
- A/B+/B- は beta 平面の選んだ漸近壁分割。
- 有限箱、potential cap、三回の再格子化を含む。
- 第一壁Aに条件付けた9履歴であり、完全な27履歴ではない。
- 現実のブラックホール観測、特異点解消、量子重力の検証ではない。
