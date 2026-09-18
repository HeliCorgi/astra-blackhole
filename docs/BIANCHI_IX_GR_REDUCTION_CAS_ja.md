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

Cadabra 2.5.14 の `cadabra2-cli` 内でEuler角座標metricを直接構成した。
Cadabraプロセス内のSymPy scalar backendを用いた明示的3×3 component calculationで
Christoffel、Ricci tensor、Ricci scalarを再計算した。

得られた閉形式は

\[
{}^{(3)}R=
\frac{2[-a^4-b^4-c^4+2(a^2b^2+b^2c^2+c^2a^2)]}
{a^2b^2c^2}.
\]

Misner変数へ代入すると厳密に

\[
{}^{(3)}R+12e^{-2\alpha}V=0
\]

となった。

さらに

\[
K_{ij}K^{ij}-K^2
=
\frac{6}{N^2}
\left(-\dot\alpha^2+\dot\beta_+^2+\dot\beta_-^2\right)
\]

とLegendre transformを確認した。

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

1. CadabraはCadabra 2.5.14プロセス内のSymPy scalar backendを使うcomponent計算。
2. Maximaは別CASの`ctensor`を使うため第二backendとして独立扱いした。
3. xActはsourceのみで未実行。
4. one-form規格化、ADM sign convention、overall canonical normalizationは明示的に固定した規約であり、規約の一意性を証明していない。
5. 古典縮約の一致から、量子factor ordering、clock、physical inner product、operator domainは決まらない。
6. timeless scattering regulator gateはFAILのままで、履歴確率は未採用。
7. BH singularity interpretationはSchwarzschild/KS return gate以降で別に判定する。
