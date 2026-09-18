# Bianchi IX GR reduction CAS audit

## 目的

`research/bianchi_ix_model.py` と `research/quantum_bianchi_ix_model.py` が使う
Bianchi IX の縮約 Hamiltonian を、実装コードから独立な3次元幾何から再導出し、
符号・係数・constraint・`s=-alpha` 後の generator を照合する。

この監査は**古典縮約規約の整合性**を対象とする。量子化、時計、物理内積、
operator domain、factor ordering、regulator removal、ブラックホール特異点解消を
証明しない。

## 実装側の参照式

repo が採用する規約は

\[
C=
\frac12 e^{-3\alpha}
\left(-p_\alpha^2+p_+^2+p_-^2\right)
+e^\alpha V(\beta_+,\beta_-)=0,
\]

\[
V=\frac16\left[
e^{-8\beta_+}
+2e^{4\beta_+}\left(\cosh(4\sqrt3\beta_-)-1\right)
-4e^{-2\beta_+}\cosh(2\sqrt3\beta_-)
\right].
\]

内部変数は \(s=-\alpha\) とし、

\[
W(s,\beta_+,\beta_-)=2e^{-4s}V,
\qquad
G_s=-\sqrt{p_+^2+p_-^2+W}.
\]

## 共通の独立3-geometry

SU(2) の標準Euler角 one-form \(\sigma_i\) を使い、

\[
\omega_i=\frac12\sigma_i
\]

と規格化する。この規格化では等方極限
\(a=b=c=e^\alpha\) の空間scalar curvatureが

\[
{}^{(3)}R=6e^{-2\alpha}
\]

になる。

対角Bianchi IX metricを

\[
h=a^2\omega_1^2+b^2\omega_2^2+c^2\omega_3^2
\]

とし、Misner変数を

\[
a=e^{\alpha+\beta_++\sqrt3\beta_-},\quad
b=e^{\alpha+\beta_+-\sqrt3\beta_-},\quad
c=e^{\alpha-2\beta_+}
\]

と置く。

両backendはrepositoryの `potential()` から曲率式を読み戻さず、このEuler角metricから
独立にscalar curvatureを計算する。

## Backend 1: Cadabra 2.5.14

Cadabra 2.5.14 の `cadabra2-cli` 内で、SU(2)左不変frameから空間曲率を再導出する。

初期実装は一般 (a,b,c) のEuler角座標metricを直接逆行列化し、全Christoffel/Ricci成分を三角関数込みで簡約していた。複数runでは成功したが、main run 35388086412 で `Run independent Cadabra Bianchi IX reduction` が約20分間完了せずworkflow timeoutでcancelされた。数式のFAILではなく、generic symbolic simplificationの実行時間不安定性だった。

そこで同じ3-geometryを座標展開せず、直交左不変frameで計算する。規格化 (omega_i=sigma_i/2) のdual frameを (e_i) とすると、向きの全体符号を除いて

[
[e_2,e_3]=2e_1,qquad
[e_3,e_1]=2e_2,qquad
[e_1,e_2]=2e_3.
]

metric

[
h=a^2omega_1^2+b^2omega_2^2+c^2omega_3^2
]

のorthonormal frame (E_1=e_1/a, E_2=e_2/b, E_3=e_3/c) では

[
[E_2,E_3]=rac{2a}{bc}E_1,
]

およびcyclic permutationとなる。

Cadabraプロセス内のscalar bridgeを使い、structure constantsからKoszul公式

[
2Gamma_{ij}{}^k
=
c_{ij}{}^k-c_{jk}{}^i+c_{ki}{}^j
]

でLevi-Civita接続を構成し、

[
R(E_i,E_j)E_k
=

abla_i
abla_jE_k
-
abla_j
abla_iE_k
-
abla_{[E_i,E_j]}E_k
]

を代数的に縮約する。

閉形式を入力として使わずに得たscalar curvatureを、最後に

[
{}^{(3)}R=
rac{2[-a^4-b^4-c^4+2(a^2b^2+b^2c^2+c^2a^2)]}
{a^2b^2c^2}
]

と比較する。追加negative controlとして等方極限 (a=b=c=r) が

[
{}^{(3)}R=6/r^2
]

になることも厳密に確認する。

Misner変数へ代入後は従来どおり

[
{}^{(3)}R+12e^{-2alpha}V=0
]

を検査し、ADM kinetic termとLegendre transformへ進む。

このframe pathは、座標Christoffelの三角関数簡約を避けるため、CIの実行時間を物理式とは無関係なCAS simplification heuristicに依存させにくい。workflowではCadabra実行自体に180秒の外部timeoutも置き、予期しないhangを20分後まで待たない。

## Backend 2: Maxima / ctensor 5.46.0

第二backendはSymPy系を避け、Ubuntu 24.04の Maxima 5.46.0 と `ctensor` を使う。

`cas/maxima/bianchi_ix_reduction.mac` は同じEuler角metricをMaxima matrixとして作り、

- `cmetric()`
- `christof(false)`
- `ricci(false)`
- `scurvature()`

で3次元曲率を計算する。

その後のMisner変数代入、ADM kinetic algebra、canonical momenta、Legendre transformも
Maxima自身で行う。

採用branch run 35387048247 では以下がすべて厳密に0だった。

- `ASTRA_MAXIMA_R3_CLOSED_DELTA`
- `ASTRA_MAXIMA_R3_MISNER_DELTA`
- `ASTRA_MAXIMA_KINETIC_DELTA`
- `ASTRA_MAXIMA_PA_DELTA`
- `ASTRA_MAXIMA_PP_DELTA`
- `ASTRA_MAXIMA_PM_DELTA`
- `ASTRA_MAXIMA_CONSTRAINT_DELTA`

従ってMaximaも独立に

\[
{}^{(3)}R=-12e^{-2\alpha}V
\]

とrepo規約のconstraintを再取得した。

## ADM normalization

shiftを0とし、

\[
K_{ij}=\frac{1}{2N}\dot h_{ij}
\]

を採用する。全体符号を逆にしても
\(K_{ij}K^{ij}-K^2\) は変わらない。

共通の空間体積/\((2\kappa)\)因子を落とし、縮約ADM作用を12で割る
repoのcanonical normalizationでは

\[
L=
\frac{e^{3\alpha}}{2N}
(-\dot\alpha^2+\dot\beta_+^2+\dot\beta_-^2)
-Ne^\alpha V.
\]

canonical momentaは

\[
p_\alpha=-\frac{e^{3\alpha}}{N}\dot\alpha,\qquad
p_+=\frac{e^{3\alpha}}{N}\dot\beta_+,\qquad
p_-=\frac{e^{3\alpha}}{N}\dot\beta_-.
\]

Legendre transformから両backendで

\[
C=
\frac12e^{-3\alpha}
(-p_\alpha^2+p_+^2+p_-^2)+e^\alpha V
\]

を再取得した。

## repo実装との比較

CadabraとMaximaの結果はそれぞれ独立checkerで

- `research/bianchi_ix_model.py::potential`
- `wall_term_and_gradient`
- `rhs`
- `research/quantum_bianchi_ix_model.py::_wall_vec`

へ照合する。

採用runの結果:

| backend | symbolic reduction | max wall relation diff | wall/gradient diff | reduced Hamilton rhs diff |
|---|---|---:|---:|---:|
| Cadabra 2.5.14 | PASS | \(2.22\times10^{-21}\) | 0 | 0 |
| Maxima/ctensor 5.46.0 | PASS | \(2.33\times10^{-21}\) | 0 | 0 |

従って、**repoが明示するBianchi IX古典縮約規約は2つの独立CAS backendと一致した**。

physics audit の `GR_REDUCTION` gate は **PASS** とする。

## xAct

`cas/xact/bianchi_ix_reduction.wl` も同じEuler角metric、Christoffel/Ricci、
ADM kinetic、Legendre transformを再計算するsourceとして保持する。

ただしrepository CIにはlicensed Wolfram Engine/xAct runtimeを設定していないため、
xActは未実行。現在はCadabra+Maximaで二独立backend条件を満たしているので、
xActは任意の第三cross-checkであり、未実行であること自体は `GR_REDUCTION` PASSを
取り消さない。

## 再現

Cadabra:

~~~sh
cadabra2-cli cas/cadabra/bianchi_ix_reduction.cdb

python tools/check_bianchi_ix_gr_reduction.py \
  --reference cas/bianchi_ix_reference.json \
  --cadabra-result artifacts/gr-reduction/cadabra_result.json \
  --out artifacts/gr-reduction/comparison.json
~~~

Maxima:

~~~sh
maxima --very-quiet -b cas/maxima/bianchi_ix_reduction.mac \
  | tee artifacts/gr-reduction-maxima/maxima.log

python tools/check_bianchi_ix_gr_reduction_maxima.py \
  --maxima-log artifacts/gr-reduction-maxima/maxima.log \
  --out artifacts/gr-reduction-maxima/comparison.json
~~~

第二backend初回成功run:

- workflow: `Bianchi IX GR reduction CAS`
- run: 35387048247
- Maxima job: 105736275544
- artifact: 10564665571
- artifact SHA-256:
  `0dbe0f4ff34b614b027d9a99ebcd16741b0c62e031a7d2ab3015d811bb9736ef`

同runではCadabra jobも再実行して成功した。

## ブラックホール本筋との関係

Bianchi IXは特定のSchwarzschild内部そのものではない。
このGR reduction PASSを「BH特異点の物理が確定した」と読み替えない。

次のblocking known-limit obligationとして
`known_limits.schwarzschild_ks_bridge` を追加した。

詳細は [BLACK_HOLE_RETURN_GATE_ja.md](BLACK_HOLE_RETURN_GATE_ja.md)。

次段階ではSchwarzschild interior → Kantowski–Sachsの古典bridgeを独立CASで導出し、
既存 `wdw_model.py` のconstraint・時計・半径・曲率observableへ戻す。

## 限界

1. CadabraはCadabra 2.5.14プロセス内のscalar bridgeを使い、SU(2)左不変frameのstructure constantsからKoszul/Riemannを代数的に構成する。
2. Maximaは別CASの`ctensor`を使うため第二backendとして独立扱いした。
3. xActはsourceのみで未実行。
4. one-form規格化、ADM sign convention、overall canonical normalizationは明示的に固定した規約であり、規約の一意性を証明していない。
5. 古典縮約の一致から、量子factor ordering、clock、physical inner product、operator domainは決まらない。
6. timeless scattering regulator gateはFAILのままで、履歴確率は未採用。
7. BH singularity interpretationはSchwarzschild/KS return gate以降で別に判定する。
