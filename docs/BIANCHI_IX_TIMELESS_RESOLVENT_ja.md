# Bianchi IX timeless stationary resolvent / T-operator pilot

## 目的

直前の finite-window S-matrix pilot では、有限 Dirichlet 箱のまま
unphysical parameter window を広げると、`[C,S_T]` と successive window 差が
むしろ増大した。これは Halliwell の無限 scattering limit を有限反射箱で近似する
設計が不適切である可能性を示す負の対照として保存する。

次段階では時間窓を伸ばさず、stationary resolvent へ切り替える。

Halliwell の complex-potential class operatorでは、領域へ入らない class operator
はその領域に局在した complex potential の S 行列として構成され、constraint と
可換であることが重要である（arXiv:0909.2597, 1108.5991）。

本 pilot は **S 行列そのものをまだ構成しない**。まず有限箱上の outgoing
resolvent / T-operator の数値核が安定化する見込みを検査する。

## 定義

前段と同じ新しい二階 constraint

[
C=P_s^2-A(s)
]

と、同じ B={B+,B-} 領域の非負関数 `V=V0 F_B` を使う。

[
H_{eff}=C-iV,qquad
G_0(z)=(z-C)^{-1},qquad
G_V(z)=(z-C+iV)^{-1},
]

[
z=E+ieta,qquad E=0,quad eta>0.
]

`eta` は outgoing resolvent を有限格子で正則化する数値 regulator であり、
物理パラメータではない。

摂動 `U=-iV` に対し有限行列で

[
G_V=G_0+G_0UG_V,
]

[
T(z)=U+UG_VU
]

を使う。最初の式は実装の代数的 negative control として直接検査する。

## reference state

前段と同じ発想で、near-zero constraint mode の部分空間から有限箱 edge mass を
小さくする packet を作る。これは物理状態の選択則ではなく、境界反射を診断から
できるだけ分離するための数値 reference。

## scan

主 grid は前段の `WDWGridSpec` を継承する。

- `E=0`
- `eta = .10, .05, .025`
- `V0 = .025, .05`

各点で

- Dyson identity residual
- free resolvent norm
- `||T psi||`
- Born項 `||U psi||` に対する比
- near-zero shellへの T-action 投影率
- interacting resolvent response の norm
- response の outer one-cell edge mass

を保存する。

さらに `eta` を半減したときの

[
rac{|T_{eta/2}psi-T_etapsi|}{|T_{eta/2}psi|}
]

を記録する。

## 判定

この pilot だけで class operator や履歴確率へ進まない。

次段階へ進む最低条件は：

1. Dyson residual が数値丸めレベルである。
2. `eta` を下げたとき T-action の successive change が低下する傾向を持つ。
3. response の edge mass が支配的でない。
4. その後、独立な box / potential-cap scan でも同じ傾向が残る。

条件2–4を満たさなければ、有限 Dirichlet resolvent も負の対照として固定し、
outgoing/PML 型境界または別の continuum scattering discretization へ進む。

## 限界

- 一様真空 Bianchi IX minisuperspace の新しい二階 WDW constraint 量子化。
- 既存 finite-clock square-root 量子化とは同一視しない。
- finite Dirichlet box と potential cap を継承。
- `eta`, shell projection, Euclidean grid norm は数値診断。
- induced/Rieffel physical inner product は未実装。
- decoherence functional と physical history probability は未計算。
- 現実のブラックホール観測、特異点解消、量子重力検証ではない。
