# Bianchi IX: A vs {B+,B-} の "ever entered region" CAP 監査

## 目的

固定三時刻の逐次射影ではなく、「第一壁で A に条件付けた後、第三壁時刻までに B 型領域へ一度でも入ったか」を問う。

現行の内部時計 s=-alpha と平方根 Hamiltonianを保った有限時計版として、B={B+,B-} に複素吸収ポテンシャルを置く。

G_eff = G - i V0 F_B

no-B branch はこの非Hermitian発展、ever-B branch は unrestricted U との差

C_everB |psi_A> = U|psi_A> - C_noB|psi_A>

で定義する。

Halliwellの量子宇宙論では、領域へ入らないclass operatorを複素ポテンシャルのS行列として作り、timeless Hamiltonian constraintと可換なoperatorへする（arXiv:0909.2597, 1108.5991）。今回の実装はその**有限内部時計版の数値監査**であり、constraint-compatible timeless class operatorを実装したとはしない。

## 領域

A/B±の三つの正の漸近壁項には共通因子 exp(-4s) があるため、A と B={B+,B-} の境界は beta 平面で時間非依存になる。

B iff

3 beta_+ + sqrt(3)|beta_-| > 0.

吸収窓 F_B は境界からの符号付き距離にcompact cubic smoothstepを使い、深いAで0、深いBで1とする。widthは物理定数ではなく数値/粗視化regulator。

## 数値法

第一壁 s=4.6 でhard A-sectorへ射影し、その条件付き状態をunit normにする。既存長時間版と同じ4有限箱を使い、s=6,12,20でcubic regridする。ただしclass branchの相対振幅を保つためregrid後の正規化はしない。

各時間刻みで unrestricted state と no-B state を同じ block-Krylov subspaceへ入れ、Hamiltonian unitary substepに**同一の近似作用素**を使う。no-B列だけ前後に exp[-V0 F_B ds/(2hbar)] を作用させるStrang CAP stepとする。

最終2 branchについて

D(alpha,beta)=<psi_beta|psi_alpha>/||U psi_A||^2

を計算する。diagonal weightを通常の確率と呼ぶのはoffdiagonalが十分小さく、V0/width/時間刻み/Krylov依存が制御された場合だけ。

## 最初のpilot

CIで以下を別jobとして走らせる。

- weak: V0=0.025, width=0.30, maxdim=66
- main: V0=0.050, width=0.30, maxdim=66
- strong: V0=0.100, width=0.30, maxdim=66
- krylov: V0=0.050, width=0.30, maxdim=78

V0は典型初期generator scale ~1より小さい範囲から開始する。plateauがなければclass operatorの数値を採用せず、吸収reflection/有限区間依存を追加監査する。

古典対照は同じ4096 Sobol Wigner初期集団をs=4.6でA条件付けし、連続積分中にhard B領域へ一度でも入った割合とfirst-entry時刻分位点をds=.01/.005で比較する。

## 限界

- 一様真空Bianchi IX minisuperspace。
- 一つの平方根量子化と内部時計。
- finite-clock CAPでありtimeless WDW class operatorではない。
- CAP strength/widthはregulator。
- 有限箱、potential cap、cubic regridを継承。
- 「壁へ衝突した」厳密な幾何学的交差数ではなく、選んだB領域へのentryを問う。
- 実在ブラックホールの観測、特異点解消、量子重力検証ではない。
