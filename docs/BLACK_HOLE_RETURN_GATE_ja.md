# BLACK HOLE RETURN GATE — Bianchi IX 監査からBH特異点へ戻る条件

更新: 2026-09-19

## 結論

戻る経路はある。ただし Bianchi IX そのものを「ブラックホール内部」と読み替えない。

このリポジトリのブラックホール本筋は、すでに実装済みの**真空 Kantowski–Sachs /
Schwarzschild 内部 WDW 系列**にある。Bianchi IX は、その一段外側で

- 異方性自由度を増やしたときの古典壁反射
- 量子波束の分岐
- histories / class operator
- timeless constraint
- regulator / boundary 問題
- GR reduction の符号・係数監査

を壊れやすい条件で検査するために使った。

したがって Bianchi IX の結果だけから「BH特異点が解消した／しない」と結論しない。
BH singularity claim へ戻る前に、以下の bridge を blocking gate とする。

## 現在のBH側の基盤

### 1. Kantowski–Sachs WDW

[WDW_AUDIT_ja.md](WDW_AUDIT_ja.md) では

[
a=e^{sqrt3eta},qquad
r=e^{-sqrt3(Omega+eta)}
]

とし、変数変換後

[
widehat CPsi=
[h^2partial_T^2-h^2partial_x^2+e^{-2x}]Psi=0,
qquad
r=rac14e^{-x-T}
]

を使っている。

一つの追加選択として

[
ihpartial_Tpsi
=
sqrt{-h^2partial_x^2+e^{-2x}}psi
]

と L2(dx) 内積を採用した。

この系列では有限Tで半径反転を確認しておらず、有限Tの正規化を
特異点回避とは扱っていない。

### 2. 非有界曲率関連量

[WDW_TAIL_AUDIT_ja.md](WDW_TAIL_AUDIT_ja.md) では逆半径モーメントの
有限値への収束が不合格になった。

ただし計算した J6 は完全な量子Kretschmann演算子ではない。
同文書自身が、次に必要なのは

- 質量を含む4次元曲率の関係的演算子
- 二次形式 / operator domain
- 時計・内積選択への頑健性

だと明記している。

### 3. 別時計・KG形式との比較

[CHIBA_REPRODUCTION_ja.md](CHIBA_REPRODUCTION_ja.md) では、
球対称真空の一セクターでKG形式を再計算し、同じ解・同じ内部切片・同じ流束を
保てば座標変換後も期待値が一致することを確認した。

一方、そのスカラー波動関数をそのまま以前の L2(dx) 状態へ読み替えられないことも
確認した。従って「時計を変えただけ」として状態・内積を混ぜない。

## Blocking bridge

machine-readable obligation:

`known_limits.schwarzschild_ks_bridge`

を追加し、以下を満たすまで BH singularity interpretation を昇格しない。

### A. Schwarzschild interior → Kantowski–Sachs の古典bridge

Schwarzschild内部 metric から Kantowski–Sachs ansatz への写像を明示し、

- lapse / signature
- 面積半径 r
- mass parameter M
- canonical variables
- Hamiltonian constraint
- boundary term
- Poisson bracket
- 特異点 r→0 の位置

を固定する。

既存KSコードへ「似ている式」を合わせるのではなく、
独立CASから係数・符号を比較する。

### B. 既存KS WDW constraintとの照合

Schwarzschild/KS bridge の古典constraintから、
現在使っている

[
C=p_T^2-p_x^2-e^{-2x}=0
]

へ至る正準変換・定数規格化・枝選択を記録する。

ここが一致しなければ、既存WDW系列を「Schwarzschild内部量子化」として
上位claimへ使わない。

### C. 同じ物理状態 / 内積の比較

少なくとも

- L2 positive-frequency branch
- KG / induced-inner-product sector

について、同じ古典境界データから対応する状態を定義する。

別のGaussianを置き直して「時計依存」と比較しない。

### D. 4次元 singularity observable

「特異点回避」の判定対象を先に固定する。

候補には

- 関係的な面積半径 r
- 4次元 Kretschmann scalar
- mass observable を含む曲率
- geodesic / relational-time endpoint

があるが、有限箱の (langle rangle) や x の反射だけでは代用しない。

曲率演算子では operator / quadratic-form domain を明示する。

### E. Regulator / clock / ordering gate

BH側でも Bianchi IX と同じ監査を通す。

- clock
- factor ordering
- physical inner product
- finite box / grid / cutoff
- boundary condition
- semiclassical / Schwarzschild limit

のどれかで結論が変われば、その依存ラベルを残す。

## Bianchi IX の役割

Bianchi IX はここで捨てるのではない。

BH singularity近傍を一般化したときに

- KSの一自由度結論が異方性を増やしても残るか
- historiesの定義がconstraint-compatibleか
- regulator removalが可能か

を見る**stress test**として使う。

しかし最終的なBH claimは、Schwarzschild/KS bridgeを通したobservableへ戻して
再評価する。

## 次の実装順

1. Bianchi IX GR reduction を2独立CAS backendで閉じる。
2. **Schwarzschild interior → Kantowski–Sachs classical CAS bridge** を作る。
3. 既存 `wdw_model.py` のconstraint変数へ正準変換を照合する。
4. KS側 clock / inner-product / factor-ordering obligationsを具体化する。
5. 質量を含む4次元曲率observableとdomainを定義する。
6. その後にだけ、BH singularity avoidance / persistence を評価する。

この順序により、Bianchi IX監査が自己目的化してブラックホール本筋から離れることを防ぐ。
