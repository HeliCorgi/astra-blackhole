# Bianchi IX GR reduction CAS audit

## 目的

`research/bianchi_ix_model.py` と `research/quantum_bianchi_ix_model.py` が使う
Bianchi IX の縮約 Hamiltonian を、実装コードから独立な3次元幾何から再導出し、
符号・係数・constraint・`s=-alpha` 後の generator を照合する。

この監査は**古典縮約規約の整合性**を対象とする。量子化、時計、物理内積、
operator domain、factor ordering、regulator removal の正しさを証明しない。

## 実装側の参照式

repo が採用する規約は

[
C=
rac12 e^{-3alpha}
left(-p_alpha^2+p_+^2+p_-^2ight)
+e^alpha V(eta_+,eta_-)=0,
]

[
V=rac16left[
e^{-8eta_+}
+2e^{4eta_+}left(cosh(4sqrt3eta_-)-1ight)
-4e^{-2eta_+}cosh(2sqrt3eta_-)
ight].
]

内部変数は (s=-alpha) とし、

[
W(s,eta_+,eta_-)=2e^{-4s}V,
qquad
G_s=-sqrt{p_+^2+p_-^2+W}.
]

## 独立な3-geometry

SU(2) の標準Euler角 one-form (sigma_i) を使い、

[
omega_i=rac12sigma_i
]

と規格化した。この規格化では等方極限
(a=b=c=e^alpha) の空間scalar curvatureが

[
{}^{(3)}R=6e^{-2alpha}
]

になる。

対角Bianchi IX metricを

[
h=a^2omega_1^2+b^2omega_2^2+c^2omega_3^2
]

とし、Misner変数を

[
a=e^{alpha+eta_++sqrt3eta_-},quad
b=e^{alpha+eta_+-sqrt3eta_-},quad
c=e^{alpha-2eta_+}
]

と置く。

Cadabra 2.5.14 の `cadabra2-cli` 内でEuler角座標metricを直接構成し、
CadabraのSymPy scalar backendを用いた明示的3×3 component calculationで
Christoffel、Ricci tensor、Ricci scalarを再計算した。

得られた閉形式は

[
{}^{(3)}R=
rac{2[-a^4-b^4-c^4+2(a^2b^2+b^2c^2+c^2a^2)]}
{a^2b^2c^2}.
]

Misner変数へ代入するとCAS上で厳密に

[
{}^{(3)}R+12e^{-2alpha}V=0
]

となった。

## ADM kinetic term

shiftを0とし、

[
K_{ij}=rac{1}{2N}dot h_{ij}
]

を採用する。全体符号を逆にしても
(K_{ij}K^{ij}-K^2) は変わらない。

対角固有値を使うとCASで

[
K_{ij}K^{ij}-K^2
=
rac{6}{N^2}
left(-dotalpha^2+doteta_+^2+doteta_-^2ight)
]

を確認した。

共通の空間体積/(2kappa)因子を落とし、縮約ADM作用を12で割る
repoのcanonical normalizationを使うと

[
L=
rac{e^{3alpha}}{2N}
(-dotalpha^2+doteta_+^2+doteta_-^2)
-Ne^alpha V.
]

canonical momentaは

[
p_alpha=-rac{e^{3alpha}}{N}dotalpha,qquad
p_+=rac{e^{3alpha}}{N}doteta_+,qquad
p_-=rac{e^{3alpha}}{N}doteta_-.
]

Legendre transformから

[
C=
rac12e^{-3alpha}
(-p_alpha^2+p_+^2+p_-^2)+e^alpha V
]

を再取得した。

## repo実装との比較

`tools/check_bianchi_ix_gr_reduction.py` はCAS outputを読み、

- `research/bianchi_ix_model.py::potential`
- `wall_term_and_gradient`
- `rhs`
- `research/quantum_bianchi_ix_model.py::_wall_vec`

と独立に比較する。

採用runでは:

- Cadabra symbolic zero checks: 全て0
- `max_wall_relation_abs = 2.2234614865425384e-21`
- classical/quantum wall gradient差: 0
- reduced Hamilton rhs差: 0

だった。

従って、**repoが明示するBianchi IX古典縮約規約は、実行済みCadabra backendと一致した**。

## xAct

`cas/xact/bianchi_ix_reduction.wl` に同じEuler角metric、Christoffel/Ricci、
ADM kinetic、Legendre transformを再計算するWolfram/xAct sourceを保存した。

ただしrepository CIにはlicensed Wolfram Engine/xAct runtimeを設定していないため、
xAct sourceは**未実行**である。設定ファイルが存在することを成功実行とは扱わない。

このため physics audit の `GR_REDUCTION` gate は:

- Cadabra backend: PASS
- xAct backend: PENDING
- aggregate: **PARTIAL**

とする。

## 再現

```sh
# CIはCadabra 2.5.14 noble x86_64 .debを公開SHA-256で固定する。
cadabra2-cli cas/cadabra/bianchi_ix_reduction.cdb

python tools/check_bianchi_ix_gr_reduction.py \
  --reference cas/bianchi_ix_reference.json \
  --cadabra-result artifacts/gr-reduction/cadabra_result.json \
  --out artifacts/gr-reduction/comparison.json
```

採用branch run:

- workflow: `Bianchi IX GR reduction CAS`
- run: 35384482998
- job: 105728047567
- artifact: 10563541351
- artifact SHA-256:
  `521a576583e5cc831d5758b36409c3fc35a015fe9555c7bc92416caf72f37183`

## 限界

1. Cadabra実行はCadabra 2.5.14プロセス内のSymPy scalar backendを使ったcomponent計算。
   Cadabraの抽象 `evaluate` tensor pathそのものではない。
2. xActはsourceを保存しただけで、まだ実行していない。
3. 空間one-form規格化、ADM sign convention、overall canonical normalizationを明示的に固定している。
   規約の一意性を証明したものではない。
4. 古典縮約の一致から、量子factor orderingやphysical inner productは決まらない。
5. 現在のtimeless scattering regulator gateはFAILのままで、履歴確率は未採用。
