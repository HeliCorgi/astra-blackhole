# Schwarzschild interior → Kantowski–Sachs classical CAS bridge

## 目的

Bianchi IX 監査からブラックホール特異点の本筋へ戻るため、既存の
`research/wdw_model.py` が使う真空 Kantowski–Sachs WDW constraintを、
Schwarzschild内部から古典的に接続する。

確認対象は:

1. Schwarzschild内部metricがKS ansatzへ入ること
2. ADM+GHY縮約の符号・係数
3. Misner変数のHamiltonian constraint
4. repoの `x,T` 変数への正準変換
5. 面積半径 `r` の同一性
6. Schwarzschild mass parameterを与えるclassical Dirac observable

である。

このbridgeは**古典的な正準構造の照合**であり、量子時計、L2内積、正周波数枝、
factor ordering、量子Kretschmann演算子、特異点解消を証明しない。

## 1. Schwarzschild内部

幾何学的質量を

[
mu=G M_{m phys}
]

とする。(0<ho<2mu) では

[
F(ho)=rac{2mu}{ho}-1>0
]

として、Schwarzschild metricは

[
ds^2=
-F^{-1}dho^2
+F,dchi_S^2
+ho^2 dOmega_2^2.
]

内部では (ho) が時間的、元のSchwarzschild時刻方向が空間的になる。

非コンパクト空間方向の座標規格化を明示するため

[
chi_S=lambdachi,qquad lambda>0
]

と置くと、KS形式

[
ds^2=-N^2dt^2+A^2dchi^2+R^2dOmega_2^2
]

へ

[
t=ho,qquad
N=F^{-1/2},qquad
A=lambdasqrt F,qquad
R=ho
]

で入る。

(lambda) は非コンパクト空間座標の規格化であり、新しい物理定数ではない。

## 2. KS ADM縮約

空間metric

[
h_{ij}dx^idx^j=A^2dchi^2+R^2dOmega_2^2
]

をMaxima/ctensorとCadabraで独立に計算し、

[
{}^{(3)}R=rac{2}{R^2}
]

を得た。

Misner変数は既存WDW監査と同じ

[
A=e^{sqrt3eta},
qquad
R=e^{-sqrt3(Omega+eta)}.
]

zero shiftで

[
K^chi{}_chi=rac{dot A}{NA},
qquad
K^	heta{}_	heta=K^phi{}_phi=rac{dot R}{NR}.
]

従って両CASで

[
K_{ij}K^{ij}-K^2
=
rac{6}{N^2}
(-dotOmega^2+doteta^2)
]

を確認した。

Einstein–Hilbert作用に時間境界のGHY項を含め、空間の非コンパクト方向を
fiducial interval (L_0) に制限すると共通因子は (L_0/(4G))。
この共通因子で作用全体を割ったnormalized reduced Lagrangianは

[
L=
rac{6e^{-2sqrt3Omega-sqrt3eta}}{N}
(-dotOmega^2+doteta^2)
+2Ne^{sqrt3eta}.
]

fiducial factorを落とすことは物理的な質量規格化の消去ではない。
以下のcanonical momentumはこのnormalized actionに対するもの。

## 3. canonical momenta とconstraint

[
P_Omega=
-rac{12e^{-2sqrt3Omega-sqrt3eta}}{N}dotOmega,
]

[
P_eta=
rac{12e^{-2sqrt3Omega-sqrt3eta}}{N}doteta.
]

Legendre transformは

[
H=
-rac{N e^{2sqrt3Omega+sqrt3eta}}{24}
left[
P_Omega^2-P_eta^2+48e^{-2sqrt3Omega}
ight].
]

従ってclassical constraint surfaceは

[
C_{m M}
=
P_Omega^2-P_eta^2+48e^{-2sqrt3Omega}=0.
]

これは既存WDW文書が出発点として使っていたconstraintと一致する。
今回は文献式を入力して合わせたのではなく、ADM縮約から再取得した。

## 4. repoの (x,T) 変数への正準変換

既存コードの定義

[
x=sqrt3Omega-log4,
qquad
T=sqrt3eta
]

に対し、canonical one-formを保存するmomentumは

[
p_x=rac{P_Omega}{sqrt3},
qquad
p_T=rac{P_eta}{sqrt3}.
]

両CASで

[
{x,p_x}=1,qquad
{T,p_T}=1
]

を確認した。

代入すると

[
C_{m M}
=
3left(p_x^2-p_T^2+e^{-2x}ight),
]

よってrepoの

[
C_{m repo}
=
p_T^2-p_x^2-e^{-2x}
]

とは

[
C_{m repo}=-rac13 C_{m M}
]

の関係にある。

非零定数倍なのでclassical constraint surfaceは同じ。
量子側でもkernelだけを見る限り全体定数倍は同じ方程式だが、
これはfactor orderingやinner productの一意性を意味しない。

## 5. Schwarzschild解上のmomentum

Schwarzschild interior embeddingを上のnormalized canonical variablesへ代入すると

[
p_T=-4mulambda,
]

[
p_x=4lambda(mu-ho).
]

これらを (C_{m repo}) に代入するとCAS上で厳密に0。

特に

[
Eequiv -p_T=4mulambda
]

で、repoのreduced classical Hamiltonian

[
E=sqrt{p_x^2+e^{-2x}}
]

と一致する。

標準Schwarzschild空間座標規格化は (lambda=1)。
一般の (lambda) は空間座標のrescalingを保持したもの。

## 6. 面積半径 — shared singularity variable

repoで使ってきた

[
r=rac14 e^{-x-T}
]

へSchwarzschild embeddingを代入すると

[
r=ho
]

が厳密に成り立つ。

したがって既存WDWコードが追跡していた `r` は、
このclassical bridge上では**Schwarzschildの面積半径そのもの**である。

classical singularityは

[
r=ho	o0.
]

これは重要な改善だが、
(langle rangle>0) や有限時計区間の正規化だけを
「量子特異点解消」とする理由にはならない。

## 7. mass Dirac observable

repo variablesだけでSchwarzschildの幾何学的質量を再構成するclassical量として

[
mu_D=
rac14 e^{x-T}p_T(p_T-p_x)
]

を得る。

Schwarzschild embeddingでは

[
mu_D=mu.
]

さらにcanonical Poisson bracketを直接計算すると

[
{mu_D,C_{m repo}}
=
-rac{p_T}{2}e^{x-T}C_{m repo}.
]

従ってconstraint surface上で

[
{mu_D,C_{m repo}}approx0.
]

つまり (mu_D) はこの縮約classical systemで弱い意味のDirac observable。

既存 `classical_path` に対する数値比較でもmass drift最大値は
約 (3.28	imes10^{-16})。

## 8. Schwarzschild radius relation

repoのclassical solution

[
x(T)=logcosh(T+u_0)-log E
]

に対して

[
mu=rac14 E e^{u_0}
]

となり、

[
rac{r}{2mu}
=
rac{1}{1+e^{2(T+u_0)}}.
]

既存 `tests/test_wdw_geometry.py` のSchwarzschild radius relationは
このbridgeの特殊化である。

独立checkerでは複数初期値に対して最大差
(4.44	imes10^{-16})。

## 9. 4次元曲率との接続

古典Schwarzschildでは

[
K_{m Schw}
=
R_{abcd}R^{abcd}
=
rac{48mu^2}{r^6}.
]

今回確定した (r) と (mu_D) により、classical levelでは

[
K_{m Schw}
=
rac{48mu_D^2}{r^6}
]

と書ける。

ただしこれは**量子演算子をまだ定義していない**。
今後必要なのは

- (widehat{mu}) と (widehat r) のordering
- (r^{-6}) のquadratic-form/domain
- physical inner product
- clock dependence
- self-adjointness / observable interpretation

である。

従って今回のPASSを「量子Kretschmannが有限／発散」とは読まない。

## 10. CASとrepo照合

初回bridge run:

- workflow: `Schwarzschild to KS classical CAS bridge`
- run: `35430090228`
- job: `105862991584`
- artifact: `10580970845`
- artifact SHA-256:
  `90d15a1fe9feafa0cc7434ff67f6a8bba1c926edce77eae99d49a45fcc851271`

backend:

- Maxima/ctensor 5.46.0: PASS, 約0.37 s
- Cadabra 2.5.14: PASS, 約2.50 s
- `wdw_model.py` comparison: PASS

repo numerical residual:

- classical constraint: (6.66	imes10^{-16})
- mass drift: (3.28	imes10^{-16})
- radius identity: 0
- Schwarzschild radius relation: (4.44	imes10^{-16})

## 11. gate判定

`known_limits.schwarzschild_ks_bridge` はこのclassical scopeでは **PASS** とする。

これにより「Bianchi IXからBHへ戻る古典bridge」は閉じる。

ただし次のblocking itemsは残る:

1. KS clock robustness
2. physical / induced inner product
3. operator domain
4. factor ordering
5. (widehat{mu}), (widehat K) のquantum observable定義
6. regulator / boundary stability
7. semiclassical correspondence

従って `PHYSICAL INTERPRETATION NOT IDENTIFIED` は維持する。

## 文献との照合

Kantowski–Sachs Misner metricと通常のWDW constraint

[
P_Omega^2-P_eta^2+48e^{-2sqrt3Omega}=0
]

は既存 `WDW_AUDIT_ja.md` が参照する López-Domínguez, Obregón,
Zacarías, arXiv:0910.4408 の通常GR sectorと整合する。

Schwarzschild内部で時空がKantowski–Sachs形式になること、
および古典SchwarzschildのKretschmann
(48mu^2/r^6) は標準GRの既知関係として照合した。

今回のCAS derivationはそれらをrepoのcanonical normalizationと変数へ
接続する役割を持つ。
