# Kantowski–Sachs量子監査：clock / 同一状態KG対応 / ordering / mass-curvature domain

更新: 2026-09-19

## 結論

Schwarzschild interior → Kantowski–Sachs の古典bridgeは既にPASSした。
今回はその先を4つの独立gateへ分けた。

| gate | 結果 |
|---|---|
| KS clock | PARTIAL |
| 同一状態 L2 positive-frequency ↔ KG | PARTIAL |
| factor ordering | FAIL — ORDERING-SENSITIVE |
| quantum mass / Kretschmann ordering-domain | PARTIAL |
| inverse-radius regulator | FAIL |

したがって、ブラックホール特異点の量子的解消／持続はまだ判定しない。

## 1. KS clock audit

repoのclassical constraintは

[
C=p_T^2-p_x^2-e^{-2x}=0.
]

正周波数Schrödinger方程式

[
ihpartial_Tpsi=Hpsi,qquad
H=sqrt{p_x^2+e^{-2x}}
]

はcanonical momentumでは

[
p_T=-H
]

の枝に対応する。

Poisson bracketは

[
{T,C}=2p_T.
]

有限 (x) では

[
H>|p_x|
]

なので (p_T
eq0)。従ってこのbranch上でTはclassical global clockとして単調。

### xはglobal clockではない

[
{x,C}=-2p_x.
]

既存classical solutionには (p_x=0) の反射点があり、そこでx-clockは停止する。
これは「xの反射」を「面積半径のbounce」と読んではいけない理由とも一致する。

### areal-radius clock

[
Y=x+T=-log(4r)
]

とすると

[
{Y,C}=2(p_T-p_x)=-2(H+p_x)<0
]

で、有限radiusのSchwarzschild内部branchでは単調。

正準変換

[
Y=T+x,quad X=x,quad
p_Y=p_T,quad p_X=p_x-p_T
]

ではconstraintは

[
-p_X^2-2p_Xp_Y-e^{-2X}=0,
]

従って

[
p_Y=-rac{p_X^2+e^{-2X}}{2p_X}.
]

ここで量子化には (p_X^{-1}) が必要になる。
classical branchでは有限radiusで (p_X>0) だがhorizon limitで0へ近づくため、
domain / extensionを決めずにT-clockとY-clockの量子同値性は主張しない。

従ってclock gateは **PARTIAL**。

## 2. 同じ物理状態をL2とKGで対応させる

二階WDWは

[
(partial_T^2+H^2/h^2)Psi=0.
]

positive-frequency sectorのKG inner productは

[
(Psi,Phi)_{KG}
=
rac{2}{h}langlePsi,HPhiangle_{L2}.
]

L2 state (chi) と同じ物理状態を表すKG field amplitudeは

[
Psi=
sqrt{rac h2},H^{-1/2}chi.
]

逆写像は

[
chi=
sqrt{rac2h},H^{1/2}Psi.
]

observableも同時に

[
O_{KG}=S O_{L2}S^{-1},
qquad
S=sqrt{rac h2}H^{-1/2}
]

と写す必要がある。

有限箱で h=.4,.2,.1 を検査すると、

- KG norm
- L2へのmap-back
- mapped x expectation
- mapped r expectation

は全て概ね (10^{-14}) 以下で一致した。

一方、**同じfield amplitudeをそのまま (|Psi|^2dx) で読むnegative control**では
(langle xangle) が最大約0.297ずれた。

従って「同じ波動関数を同じ確率密度で読む」のは同一状態比較ではない。

### continuum domain

有限箱では (H>0) なので (H^{-1/2}) は有界。
右端を16→24→32へ広げると最小energyは

[
0.04326	o0.02791	o0.02059
]

へ下がり、(H^{-1/2}) の最大scaleは増加する。

continuumのzero-energy thresholdを含む完成空間では
(H^{-1/2}) のdomainを別途指定する必要がある。
したがってinner-product gateも **PARTIAL**。

Chiba–Matsui–MurataのKG再現で得たscalar pullbackには有限boost timeで非零定数裾があり、
そのfieldをそのままL2(dx) stateとはできない。
今回の (H^{pm1/2}) mapは、この以前のnegative controlと矛盾しない。

## 3. factor ordering

明示的stress-test familyとして

[
A_q=
-h^2e^{-qx}partial_x(e^{qx}partial_x)+e^{-2x}
]

を採用した。

これは

[
L^2(e^{qx}dx)
]

でformal symmetric。
unitary map (g=e^{qx/2}f) によりflat L2では

[
A_q^{flat}
=
A_0+rac{h^2q^2}{4}.
]

従って全qでclassical principal symbolは同じ

[
p_x^2+e^{-2x},
]

差は (O(h^2))。

同じflat-representation初期stateを使い、q=-2,-1,0,1,2で比較した。
既存のbox mean-x数値許容差は (3	imes10^{-4})。

最大差は

[
max|Deltalangle xangle|
simeq0.526
]

(h=.4, |q|=2)。

h=.2でも

- |q|=1: 約0.0344
- |q|=2: 約0.1328

h=.1でも

- |q|=1: 約0.00789
- |q|=2: 約0.0314

で、数値誤差より十分大きい。

従って**この宣言したfamily内ではORDERING-SENSITIVE**。

ただしこれは全ての物理的に許されるorderingの分類ではない。
最近のflat-minisuperspace研究にはpath-integral measureとorderingを同時に扱うと
制限されたprescription群が物理的に同値になり得るという議論もある
(Franken et al., arXiv:2512.23656)。
従って次段階は「任意qを全部物理的候補」とするのではなく、
physical inner product / path-integral / symmetry principleで許容familyを狭めること。

## 4. classical massからquantum massへ

古典bridgeで得た

[
mu=
rac14e^{x-T}p_T(p_T-p_x)
]

に対し、局所ordering family

[
widehat M_d
=
rac14e^{x-T}
[hat p_T^2-hat p_That p_x
+2dh(hat p_T-hat p_x)]
]

を調べた。

SymPyで厳密に

[
[widehat M_d,widehat C]
=
-rac{ih}{2}e^{x-T}
(hat p_T+2dh)widehat C
]

を確認した。

従ってこのfamilyはconstraint kernelを保つ。

しかしflat kinematical L2でformal adjointを要求すると

- (d_x) coefficientから (operatorname{Im}d=1/4)
- (d_T) coefficientから (operatorname{Im}d=3/4)
- identity termから (operatorname{Re}d=0,\ operatorname{Im}d=1/2)

が同時に要求され、解がない。

これは重要だが、**physical KG inner productでのno-goではない**。
flat kinematical L2をphysical inner productと同一視してはいけない。

### constraint kernel と正周波数sectorは別

さらに同じ物理状態のKG表現へ (\widehat M_d) を作用させ、

[
(\hat p_T+H)\widehat M_d\Psi
]

を直接測った。これは (widehat M_dPsi) が現在採用する
(p_T=-H) のpositive-frequency sectorに残るかを見る検査。

(d=0,,i/4,,i/2,,3i/4) を有限箱dx=.04で試すと、
relative branch residualの最良値でも約

[
0.9913
]

だった。

従って、**full WDW constraint kernelを保つことは、選択したpositive-frequency
Hilbert sectorを保つことを意味しない**。現在の局所 (M_d) familyをそのまま
T-clock量子理論のmass observableとして採用しない。

これは非局所なDirac observableやphysical KG inner product上の別orderingを否定しない。

### finite-box candidates

left / symmetric / six-permutation Weyl candidatesも検査した。

dx=.02で

- symmetric candidate mass expectation drift: 約4.79%
- Weyl6: 約1.69%

だった。

driftはdx精密化で低下しているため、これだけをcontinuum非保存の証明にはしない。

## 5. domain: massとKretschmann

古典Schwarzschildでは

[
K=48rac{mu^2}{r^6}.
]

しかし既存tail auditではq=2,4,6のpositive exponential momentが
finite regulator removalに失敗している。

特に

[
|e^xpsi|^2
=
int e^{2x}|psi|^2dx
]

はmass orderingの一部で必要になり得るためq=2 domain問題に直結する。

またbare

[
langle r^{-6}angle
]

はq=6 exponential momentと同じ問題を持つ。

ただしmass operatorを先に作用させるnonlocal/symmetric orderingではtail cancellationやdomainが変わり得る。
従って「全ての (48widehat{mu^2/r^6}) が発散する」とは結論しない。

現在の正しい判定は:

- selected (widehatmu): **なし**
- selected quantum Kretschmann: **なし**
- domain gate: **PARTIAL**
- inverse-radius regulator: **FAIL**

## 6. 現在のBH量子claim

古典Schwarzschild/KS bridgeはPASS。

しかし量子側は

- clock: PARTIAL
- physical inner product/domain: PARTIAL
- factor ordering: FAIL / sensitive
- curvature regulator: FAIL
- semiclassical common-state limit: PENDING

なので、特異点解消／持続の物理解釈はまだ昇格しない。

## 7. 次の狭い課題

1. physical KG Hilbert spaceを保つmass Dirac observableを、局所family以外も含めて構成する。
2. path-integral/symmetry principleからKS ordering familyを狭める。
3. Y-clockの (p_X^{-1}) domainを定義し、同一Dirac stateをT/Y clockで比較する。
4. 選ばれた (widehatmu) からpositive quadratic-form型Kretschmann候補を作る。
5. tail domainとregulator removalをその**選ばれたoperator**で再監査する。

この順番を飛ばして有限箱の曲率値を「量子ブラックホール特異点の答え」とはしない。


## 背景文献

- Ashtekar, Tate, Uggla, arXiv:gr-qc/9302026, 9302027 — minisuperspaceのDirac quantization、observable、deparametrizationの古典的背景。
- Mostafazadeh, arXiv:gr-qc/0205049, 0306003 — Klein–Gordon型方程式のpositive-definite Hilbert-space structureとobservable表現の背景。
- Franken et al., arXiv:2512.23656 — flat minisuperspaceでpath-integral measureとorderingを同時に扱い、許容orderingの物理的等価性を議論する最近の結果。今回のKS ordering stress testへその結論を自動適用したものではない。
