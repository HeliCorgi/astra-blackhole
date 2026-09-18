# Bianchi IX 複数時刻履歴・デコヒーレンス監査

## 問い

一時刻の A/B+/B- 質量ではなく、第一・第二・第三壁をまたぐ有限時計の class operator

C_(A,j,k)=P_k(s3) U(s3,s2) P_j(s2) U(s2,s1) P_A(s1) U(s1,s0)

を作り、branch state |psi_alpha>=C_alpha|psi0> と

D(alpha,beta)=<psi_beta|psi_alpha>

を直接計算する。時刻は s1=4.6, s2=10.385, s3=32.46。

既存第一壁計算で A の一時刻質量が約98%なので、主解析は A に条件付けた9履歴を扱う。
第一壁の残り約2%は除外質量として別に報告し、完全27履歴を計算したとはしない。

## 確率と呼ぶ条件

D の対角成分は履歴重みであり、非対角成分が無視できない限り通常の加法的確率とは呼ばない。
最終壁 k の粗視化確率

|| sum_j |psi_(A,j,k)> ||^2

と細分履歴の対角和

sum_j ||psi_(A,j,k)||^2

の差も additivity defect として保存する。

## 数値上の修正：独立Lanczosを棄却

最初の試作では第二壁の3 branchを別々の単一ベクトルLanczosで第三壁まで進めた。
各Lanczos空間が入力branchに依存するため、数値時間発展が厳密には同じ線形演算子にならない。

試作では全branchを足し戻した状態と未射影状態の相対ノルム差が約2.39%になり、
調べたい干渉量と同程度だった。この試作結果は科学的判定から棄却する。

本版では、3 branch全体から共通の block-Krylov 空間を各時間刻みで構成する。
投影された Hermitian 行列 H=Q^* A Q の同じ matrix function

Q exp(+i ds sqrt(H)/hbar) Q^*

を全branchへ作用させる。これにより、表現されたbranch span内で同一の線形近似時間発展を使う。

s=6,12,20 の再格子化も branch ごとに再正規化せず、同じ cubic interpolation を線形に適用する。
main block dimension と大きい control dimension を比較し、D行列の差を保存する。

## 先行研究との区別

decoherent / consistent histories では class operator と decoherence functional によって履歴確率が成立する条件を調べる。
量子宇宙論への応用として Halliwell, arXiv:1108.5991 および arXiv:0909.2597 などがある。

今回の構成は選んだ内部時計 s と時間依存平方根Hamiltonianを持つ有限時間模型での逐次射影。
timeless Wheeler-DeWitt制約と可換な不変class operatorを構成したものではない。

## 判定項目

- 9履歴の対角重み
- max |D(alpha,beta)|, alpha!=beta
- max |Re D(alpha,beta)|
- sqrt(Daa Dbb) で規格化したpair coherence
- coarse/fine additivity defect
- 最終射影の再結合closure
- branch Gram matrixの変化（主に非ユニタリな数値regridの影響）
- block-Krylov dimension control
- 第一壁で条件付けから除外した非A質量

## 限界

一様真空 Bianchi IX minisuperspace、一つの平方根量子化、L2表現、内部時計 s=-alpha に限定。
A/B+/B- は選んだbeta平面分割。有限箱、potential cap、cubic regridを含む。
現実のブラックホール観測、特異点解消、量子重力の検証ではない。


## 実行結果（2026-09-18）

最初の独立 single-vector Lanczos 版は棄却した。branch ごとに異なる Krylov 空間を作ると「同じ U を各 branch に作用させる」という class operator の線形性を数値的に壊し、branch 再結合と未射影状態の相対ノルム差が約 0.0239 になったためである。

公開採用版は3 branch を同じ block-Krylov 空間で同時に伝播する。主値 block dimension=102、対照=114。102次元主値では：

- 第一壁 A 重み: 0.9797417699（非A 0.0202582301）
- max |D(alpha,beta)|, alpha!=beta: 0.0194181091
- max |Re D(alpha,beta)|: 0.0043431962
- max normalized pair coherence: 0.179989172
- max normalized real pair: 0.089656130
- 対角和: 0.995267558
- 最終射影の branch 再結合相対ノルム誤差: 0
- branch Gram の s2→s3 相対変化: 0.006379632（主に三回の数値 regrid の系統を含む）
- 102→114 の最大 D 行列差: 0.00501324
- 102→114 の最大非対角差: 0.00190454

最大複素非対角は数値対照差より約10倍大きい。従って、**この選んだ9履歴について厳密な decoherence を採用しない**。一方、弱い consistency に直接効く実部は regrid 系統と同程度なので、弱い条件の成立／不成立をこの計算だけで確定しない。

第三時刻の coarse/fine 加法性欠損は A=+0.0000505, B+=-0.007367, B-=+0.012049。完全にゼロではないが、その一部は弱い実部判定と同じ数値系統の影響を受ける。

A条件付きの対角重み上位は：

| 履歴 | 量子 branch 対角重み | 古典 Wigner 条件付き頻度 |
|---|---:|---:|
| A→B-→B+ | 0.324336 | 0.353748 |
| A→B-→A | 0.265352 | 0.171792 |
| A→B-→B- | 0.150691 | 0.213469 |
| A→B+→B- | 0.114997 | 0.090978 |

量子列の値は decoherence 不成立のため通常の確率ではない。古典列は4096本の決定論的 Hamilton 軌道の固定三時刻ラベル頻度で、ds=.01→.005 でラベル変更は0だった。両列を「確率差」として採点しない。

### 再現

```sh
python -m unittest discover -s tests -p 'test_bianchi_ix_histories.py' -v
python research/run_bianchi_ix_histories.py --out artifacts/bianchi-ix-histories --main-maxdim 102 --control-maxdim 114

python -m unittest discover -s tests -p 'test_bianchi_ix_classical_histories.py' -v
python research/run_bianchi_ix_classical_histories.py --out artifacts/bianchi-ix-classical-histories
```

GitHub Actions の採用実行は量子 run 35328641671、古典 run 35328749788。CI成功は数値再現性の記録であり、この量子化・時計・壁分割の物理的正しさを保証しない。
