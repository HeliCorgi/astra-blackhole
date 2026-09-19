# KS mass affine-covariance / ordering-selection obstruction 監査

更新: 2026-09-20

## 結論

PR #20 で positive-frequency KS Hamiltonian

\[
H=\sqrt{p_x^2+e^{-2x}},\qquad H\ge0
\]

の continuum KG completion と \(I+V\ge0\) は確認したが、

\[
q_M[\psi]=\frac14
\langle H e^{X/2}\psi,(I+V)H e^{X/2}\psi\rangle
\]

の closability / associated self-adjoint operator / ordering uniqueness は未解決のまま残った。

今回、「古典Schwarzschild massのaffine則を量子でも厳密に課せばorderingを選べるか」を切り分けた。

| 対象 | 判定 |
|---|---|
| 古典 \((u,E)\) 正準変数とmass affine bracket | PASS |
| null座標 \((U,W)\) でのconstraint / mass再表現 | PASS |
| positive-frequency sectorでの exact affine covariance による自己共役mass選択 | **不可能** |
| 現行 \(M_0\) の exact affine commutator | 負の対照: O(1) residual |
| \(q_M\) closability / closed form / associated operator | **PENDING（変更なし）** |

重要なのは、これは「mass operatorが存在しない」という一般no-goではない。
**半有界なpositive-frequency \(H\) と、非零のpositive self-adjoint \(M\) に、古典のaffine scalingを全実時刻で厳密に要求することが両立しない**、という選択原理のno-goである。

従って今後、current \(M_0\) のclosabilityを調べる際に「exact affine covarianceを満たさないから棄却する」という基準は使わない。

## 1. 模型・範囲

- 真空 Schwarzschild interior / Kantowski–Sachs WDW minisuperspace
- dimensionless repo conventions
- retained: \(x,p_x\)、selected positive-frequency \(T\) branch
- omitted: matter, inhomogeneous modes, environment
- current mass candidate:
  \[
  M_0=\frac14e^{X/2}H(I+V)He^{X/2},
  \qquad V=\frac{i}{h}[H,X]
  \]
- 実在ブラックホール、特異点解消／持続、完全量子重力の判定ではない

## 2. 古典affine構造

\[
E=\sqrt{p_x^2+e^{-2x}},
\qquad
u=\operatorname{arsinh}(p_x e^x)
\]

と置く。

SymPyでcanonical Poisson bracketを直接計算し、

\[
\boxed{\{u,E\}=1}
\]

を確認した。

基準時計 \(T=0\) のclassical massは

\[
\mu_0=\frac14e^xE(E+p_x).
\]

\[
e^u=p_xe^x+\sqrt{1+p_x^2e^{2x}}
=e^x(E+p_x)
\]

なので、

\[
\boxed{\mu_0=\frac{E}{4}e^u}.
\]

従って

\[
\boxed{\{E,\mu_0\}=-\mu_0},
\qquad
\{\mu_0,E\}=+\mu_0.
\]

これは classical reduced dynamics で \(\mu_D(T)=e^{-T}\mu_0\) がDirac massとして一定になることと整合する。

## 3. null座標への再表現

以前のChiba照合でも使ったnull minisuperspace coordinatesを、速度演算子 \(V\) と区別するためここでは \(U,W\) と書く:

\[
U=-\frac14e^{-x-T},\qquad
W=\frac14e^{-x+T}.
\]

canonical one-formから

\[
p_T=-Up_U+Wp_W,
\qquad
p_x=-Up_U-Wp_W.
\]

従って

\[
p_T^2-p_x^2=-4UWp_Up_W.
\]

また \(e^{-2x}=-16UW\) なのでrepo constraint

\[
C=p_T^2-p_x^2-e^{-2x}
\]

は

\[
\boxed{
C=-4UW(p_Up_W-4)
}.
\]

constraint surfaceでは \(p_Up_W=4\)。

Schwarzschild massは

\[
\boxed{
\mu=\frac18p_Tp_W
}
\]

となり、on shellで

\[
\mu=\frac18(Wp_W^2-4U).
\]

## 4. exact quantum affine covariance を課すと何が起こるか

古典則から、量子でも

\[
U(t)MU(t)^\dagger=e^{-t}M,
\qquad
U(t)=e^{-iHt/h}
\]

を全 \(t\in\mathbb R\) で厳密に要求することを考える。

仮定:

1. \(H\) はself-adjointかつsemibounded。現行positive-frequency sectorでは \(\sigma(H)=[0,\infty)\)。
2. \(M\) は非零のpositive self-adjoint operator。
3. 上のexact covarianceを全実 \(t\) で満たす。

\(M\) にkernelがあればそのsupportへ制限する。covarianceによりsupportは \(U(t)\)-invariantであり、その上で \(M\) はinjective。
spectral calculusで \(Q=\log M\) を定義すると

\[
U(t)QU(t)^\dagger=Q-tI.
\]

さらに \(W(s)=e^{isQ}\) とすれば

\[
U(t)W(s)=e^{-ist}W(s)U(t).
\]

従ってStoneの一意性から

\[
\boxed{
W(s)^\dagger H W(s)=H+hsI
}.
\]

しかしunitary conjugationはspectrumを保存するので

\[
\sigma(H)=\sigma(H)+hs
\]

が全 \(s\in\mathbb R\) で必要になる。

非空のsemibounded spectrumは全実平行移動に不変ではない。したがって矛盾。

\[
\boxed{
H\ge0\ \text{のselected sector上では、非零positive self-adjoint }M
\text{ に古典affine scalingを厳密には課せない。}
}
\]

これはfinite box数値に依存しない。

## 5. 現行 \(M_0\) はexact affine operatorではない

現行candidateについて、finite-box matrixで

\[
[H,M_0]\stackrel{?}{=}-ihM_0
\]

を負の対照として直接測った。

設定:

- \(h=.2\)
- left=-4
- dx=.08
- right=8,12,16
- repoと同じpositive fourth-order Dirichlet discretization

結果:

| right | relative residual |
|---:|---:|
| 8 | 1.6208417247 |
| 12 | 1.6208419879 |
| 16 | 1.6208419864 |

commutatorのanti-Hermiticity residualは \(4.4\times10^{-13}\)〜\(9.6\times10^{-13}\)。

したがって、このfinite-box candidateはexact affine commutatorを満たさない。
これは前節のno-goと矛盾しない。

**このO(1) residualをmass candidateの失格理由にはしない。**
current candidateはpositive-frequency Hilbert space内のpositive formとrelational transport

\[
M_D(T)=U(T)M_0U(T)^\dagger
\]

を選んだもので、exact affine representationを宣言していない。

## 6. Cavaglià–de Alfaro–Filippov 1995 との比較

Cavaglià, de Alfaro, Filippov, arXiv:gr-qc/9508062 は別のSchwarzschild minisuperspace canonical quantizationで、gauge-invariant quantities \(I,J,N\) がaffine algebraを作ることを用いてphysical measureを決めている。

そのgauge-fixed positive-definite Hilbert spaceでは、classical massに対応するHermitian \(\hat J\) がhalf-line supportのためself-adjoint extensionを持たず、一方 \(\hat J^2\) はboundary conditionを選んでself-adjointにできる、と報告している。

これは今回のrepo \(M_0\) と同一operatorではない。canonical variables、gauge fixing、physical Hilbert spaceが異なるため、同論文のdeficiency indicesを \(M_0\) へ移植しない。

ただし「Schwarzschild massのquantizationではsupportとaffine structureがself-adjointnessを支配し得る」という独立の警告として重要で、今回のspectral no-goと方向は整合する。

## 7. mass closabilityについて今回わかったこと／わからないこと

今回わかったこと:

- classical massには明確なaffine bracketがある
- exact affine covarianceはpositive-frequency semibounded \(H\) 上のpositive self-adjoint mass選択原理として使えない
- current \(M_0\) がexact affine operatorでないことは数値的にも確認できる

まだ未確定:

\[
q_M[\psi]
=\frac14
\langle H e^{X/2}\psi,(I+V)H e^{X/2}\psi\rangle
\]

の

- closability
- closed form completion
- representation theoremで得るassociated self-adjoint operator
- ordering uniqueness

従って mass_continuum_form は **PENDINGを維持**する。

次の解析では、exact affine symmetryによる一意化ではなく、

\[
T_M=(I+V)^{1/2}H e^{X/2}
\]

そのもののclosability、またはそのadjoint domainを直接調べる。

## 8. 再現

~~~sh
python -m unittest discover -s tests -p 'test_ks_affine_mass_obstruction.py' -v

python research/ks_affine_mass_obstruction_audit.py \
  --out artifacts/ks-affine-mass/summary.json

python tools/check_ks_affine_mass_obstruction_results.py \
  --artifact artifacts/ks-affine-mass/summary.json \
  --baseline research/ks_affine_mass_obstruction_results/summary.json
~~~

ローカル独立実行では4 unit tests、audit本体、regression checkerがPASSした。
GitHub Actionsでclean repository上の同じ計算を再実行する。

Lean sourceは変更しない。

## 9. claim boundary

今回の正しい主張は:

- selected positive-frequency quantizationで **exact affine covarianceを自己共役massのordering-selection principleにはできない**
- current \(M_0\) のcontinuum closabilityは依然PENDING
- selected explicit-radius Kretschmann 3 formのFAILも変更しない

主張しないもの:

- 「量子massは存在しない」
- 「current \(M_0\) は非closable」
- 「Cavagliàらのmass operatorとrepo \(M_0\) は同じ」
- 「ブラックホール特異点が解消／持続した」
