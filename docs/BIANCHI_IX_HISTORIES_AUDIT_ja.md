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
