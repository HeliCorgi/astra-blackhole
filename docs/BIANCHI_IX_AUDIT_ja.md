# Bianchi IX / Mixmaster へ対称性を一段外す監査

従来のKantowski–Sachs (KS) は1自由度の可積分模型で、Mixmaster型カオスを表せない。そこで真空Bianchi IXのMisner異方性 `(beta+, beta-)` と共役運動量を残し、体積変数 `alpha` を内部時刻として古典方程式を積分した。

今回の連続Bianchi IX軌道は少なくとも2回の明瞭な壁反射を示した。ポテンシャルが十分小さいKasner区間で抽出した `u` は

`1.7112429 -> 1.4056553 -> 2.4651063`

となり、標準BKL写像 `u>=2ならu-1, 1<u<2なら1/(u-1)` の予測に対する相対差は約 `2.38e-4`, `1.67e-5` だった。前回のGauss/BKL写像だけの対照から一段進み、**連続Bianchi IXの壁反射とBKLカオス写像を数値的に接続した**。

有限の2反射区間そのものでは、1e-8の近接初期条件が指数的に離れる証拠は出ない。抽出した最初のuをBKL写像へ延長すると、1e-12の差は35遷移目で0.1を超える。同一初期値では完全に同じ列を返すので、これは確率ノイズではなく決定論的感度の対照である。

## 方程式

採用規約は

`C = 1/2 exp(-3 alpha)(-p_alpha^2+p_+^2+p_-^2)+exp(alpha)V=0`

`V=(1/6)[exp(-8 beta+)+2 exp(4 beta+)(cosh(4 sqrt(3) beta-)-1)-4 exp(-2 beta+)cosh(2 sqrt(3) beta-)]`。

`alpha` を時刻にすると `H=-p_alpha=sqrt(p_+^2+p_-^2+2 exp(4 alpha)V)`。本コードでは特異側へ進む `s=-alpha` を使う。指数壁の浮動小数点オーバーフローを避けるため `2 exp(-4s)V` を指数項へ展開して評価するが、壁を切り詰めたり反射条件を手で入れたりしない。

規約の出典: Phys. Rev. D 109, 044038 (2024) accepted manuscript, Eq. (5)-(7), https://link.aps.org/accepted/10.1103/PhysRevD.109.044038 。Bianchi IX連続流とBKL写像の関係: Imponente & Montani, arXiv:gr-qc/0401086。一般的なcosmological billiards: Damour, Henneaux, Nicolai, arXiv:hep-th/0212256。

## 軌道・カオス・量子を分ける

- KS: 古典軌道は一本で可積分。
- Bianchi IX: 各初期値には古典軌道があるが、壁反射列はBKL写像へ近づき、長い列では初期値感度が強い。**軌道がないのではなく、決定論的に予測が難しくなる**。
- 確率過程: 法則にノイズを加える。今回未実装。
- 量子Bianchi IX: 一般には波動関数/状態を基本にし、一本の古典軌道を基本変数にしない。**今回は同じBianchi IX Hamiltonianの量子時間発展をまだ実装していない。** リポジトリの既存WDWは概念対照にはなるが、別の縮約模型なので同一条件比較には数えない。

量子Mixmasterには先行研究がある。Lecian, Montani, Moriconi, Phys. Rev. D 88, 103511 (2013), arXiv:1311.6004 はpolymer表現で波束の広がりと古典カオスの変化を扱う。Bergeron et al., Phys. Rev. D 92, 061302 (2015), arXiv:1501.02174 は別の量子化でBianchi IXを扱う。これらを今回再現したとはしない。

## 限界と次

これは空間一様なBianchi IX宇宙模型であり、現実のSchwarzschild内部そのものではない。BKLが一般的な時空特異点近傍の局所モデルになるという仮説と、特定ブラックホール深部の物理は分ける。

次は**同じBianchi IX reduced Hamiltonianを量子化した2次元異方性波束**を実装し、波束ピーク・分散/干渉・Kasner遷移列を平均へ潰したときの情報損失を同じ初期データで比較する。固定三角形の玩具ビリヤードを「量子Bianchi IX」と呼ぶことは避ける。
