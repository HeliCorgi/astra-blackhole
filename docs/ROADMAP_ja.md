# 研究段階と次の切り分け

2026-09-16。実装済みと未実施を分ける。

## 1. 角方向分布の省略：限定モデルで実施済み

[方向分布の監査](ANGULAR_AUDIT_ja.md)を追加した。同じ初期エネルギー運動量テンソル、粒子数、全径方向スペクトル、指定した低次モーメントを保ち、P4/P6/P8の角方向成分だけを変えた。正値性・重力制約・初期一致、曲率の時間微分の差、質量ゼロと対照を検査した。

これは一様なKantowski–Sachs系、有限区間の結果。一般的な空間非一様崩壊や現実のブラックホールに対する完了宣言ではない。

## 2. 衝突と反応時間：質量ゼロ・定数RTAで実施済み

[衝突監査](COLLISION_AUDIT_ja.md)で、粒子数・エネルギーに整合する局所平衡、運動量非依存の定数緩和率を採用した。質量ゼロのエネルギー積分した角分布と時空を連立し、7緩和率×3方向成分×2符号を同じ固有時間で比較した。初期状態の区別の縮小と完全流体参照との差を分離した。粗い角度展開の負分布を棄却し、数値的に未分解の差も明示した。

二体散乱積分、運動量・密度依存の衝突率、質量ありの粒子混合は未実施。モデル依存の検査として残る。RTAが自然の微視的相互作用を説明したとはしない。

### 2a. 最小状態変数の三者比較：今回実施

[閉じ方の監査](CLOSURE_AUDIT_ja.md)。密度のみA、密度と圧力差B、さらにc4を残すCを、同じ90条件・時間幅[0,.30]で比較した。別幾何学60条件で最小の総合合格候補はA=12、B=30、C=11、該当なし=7。保持変数の絶対的な最小数ではなく、三つのゼロ打ち切り候補についての結果である。

全270低次軌道中38件では、有限モーメントが正の測度と両立しても、文字通りの打ち切り多項式は負になる。情報不足と再構成の不適切さを区別する必要がある。

### 2b. 同じ変数数で正値性を保つ閉じ方：限定モデルで実施済み

[正の角分布による比較](POSITIVE_CLOSURE_AUDIT_ja.md)。前回90条件・同じ時間区間で、角エントロピーを最大化する指数分布から上位係数を代数的に決めた。ρ・ΠのBは60→87/90、ρ・Π・c4のCは82→90/90で総合合格。旧Cの8不合格を変数追加なしで解消した。ただしK誤差はBの24条件、Cの31条件で増えた。新しい未使用テストではなく、前回条件を再使用した比較。全fの熱力学エントロピーを最大化したとの主張はしない。

次は未使用の初期幾何学・混合分布を先に固定して再評価し、その後に事前警告を進める。特にBの圧力誤差が残る3条件や、境界近くのモーメント、低い緩和率での適用限界を保留せず調べる。モデル全体を変更した検証、二粒子相関や量子効果とは区別する。

### 2c. 適用範囲の事前警告：未実施

閉じ方の比較後に、簡略モデル自身の時間尺度比・圧力差・履歴だけで、基準との誤差を予告できるかを検査する。現在の表は採点後の適用範囲であり、未来を見ずに判定する警告器ではない。学習と評価を分け、判断不能を残す。

## 3. 二粒子相関：未実施

一粒子分布を一致させ、二粒子の相関だけが異なる状態を比較する。平均値の違いと分散・相関の違いを分ける。現在のVlasov/RTAはいずれも一粒子分布の記述で、二粒子相関は未実装。

## 4. 量子応力の揺らぎ：未実施

量子状態、再正規化、粗視化、量子応力の反作用を明示した別のモデルが必要。量子補正という名前の任意関数を置いただけで自然界の量子効果を計算したとはしない。

## 5. 真空幾何学の量子化：限定したWDW枝で実施

[WDW監査](WDW_AUDIT_ja.md)。T=log(a)を時計に選び、平方根HamiltonianとL2内積で幾何学の波束を進めた。5状態、格子精密化・領域拡大・古典Wigner集団との比較を実施。有限区間で半径の幾何平均の反転は出ず、xの反射との取り違えを否定する対照を保存した。これは従来のRTA模型の量子化ではなく、真空の別モジュール。

次は時計／演算子順序の選択への依存と、裾の収束を含む量子曲率の定義を検査する。量子曲率、特異点終端、全自由度への誤差評価は未実施。古典閉包の未使用条件での評価と事前警告の課題は取り消さない。

## 検証基盤

従来33件に正の閉じ方15件を追加して計48件。従来の `reproduce.py` は全単体テストと三つの従来監査を実行する。第四の三者比較は `research/run_closure_audit.py` と専用CIで実行する。保存済みの集計結果と全条件の再計算を比較し、詳細な各行は成果物へ残す。個々の行を保存値と照合する回帰検査とは区別する。

第五の比較は `research/run_positive_closure_audit.py` と独立したCI設定で実行する。CIでは構造・ソース指紋と保存集計値を回帰比較し、全450行はZIP成果物に保存する。個別行の比較は元の配布ZIPまたはCI成果物を用いて別に行う。GitHub上の実行済み状況は検証記録で確認する。

Lean 4.19.0の既存8命題、公理の検査方針、成功記録は維持する。今回の閉じ方・実積分・曲率時間発展のLean形式化は未実施。CI通過と自然界での妥当性は別。CIはpush・PR・手動起動のみで、定期研究や自動コミットは行わない。

## 共通の判定

何を一致させ、何を変え、どの予測が区別不能になるかを明記する。差が出ない対照も保存する。数値誤差・測定誤差・モデルの省略を分ける。追加仮説を使う場合は、実観測の検証と混同しない。

WDW用15件のテストと専用CIを追加。量子計算の再現性と、採用した時計・内積の物理的妥当性は別に扱う。Leanの既存8命題は変更していない。

## 量子幾何学の補足：裾の監査を実施

[WDWの裾と非有界モーメント](WDW_TAIL_AUDIT_ja.md)を追加。元の時間発展を変えず、q=2,4,6の逆半径モーメント、初期丸め誤差、領域と観測量の切断、低エネルギー漸近形を比較した。有限値への収束は不合格。連続モデルの発散に関する結論は、明示した閾値漸近展開の条件に依存する。4次元曲率演算子そのものや測地線完全性の結論ではない。

次は質量を含む関係的な4次元曲率演算子の順序・二次形式の定義域を定め、時計・内積選択への依存を検査する。都合よく裾や低エネルギー成分を除去して特異点回避とはしない。古典閉包の未使用条件テストや事前警告の課題は未完のまま保持する。既存Lean命題は変更なし。新しい数値・スペクトル解析を形式証明済みとはしない。

## 文献照合の次段階：同じ解の再現を実施

[Chibaらの限定セクターの再現](CHIBA_REPRODUCTION_ja.md)。境界データ・Cartesian時計・KG形式を合わせた四条件の差分／積分比較と、同じ内部切片の座標変換を実施。同じ規則を保持すれば流束と期待値が一致することを検査した。スカラー引戻しの定数裾から、以前のL2 Gaussian状態と単純には同一視できないことも確認した。著者の生データの再現と、全領域・全枝の量子同値性は未実施。

次は、正周波数の選択と物理内積をそろえた状態／観測量の表現変換を明示する。別のGaussianを置き直して同じ状態の比較とはしない。既存の古典閉包・曲率演算子・相関の未実施課題は維持する。


## 6. Bianchi IX / Mixmaster：軌道・量子波束・履歴を同じ縮約模型で比較

古典側では、Kantowski--Sachsの1自由度可積分系から一段対称性を外した一様真空 Bianchi IX を実装した。連続方程式の最初の二壁反射から抽出した Kasner u は BKL 写像と一致し、長期初期値鋭敏性は BKL 対照で確認した。有限二反射区間だけをカオスの証明とはしない。

量子側では同じ古典主記号に対する一つの2次元平方根量子化を有限箱で実装し、第一壁から第三壁まで段階的に箱を拡大した。一本の古典中心軌道からのずれの大部分は同じ初期 Wigner 集団の古典分布でも生じる一方、第三壁では A/B+/B- 一時刻セクターがほぼ三方向へ広がる。これは一時刻分布であり履歴確率ではない。

### 6a. 三時刻 class operator / decoherence functional：限定模型で実施

[複数時刻履歴監査](BIANCHI_IX_HISTORIES_AUDIT_ja.md)で s=4.6,10.385,32.46 の逐次射影を使い、第一壁 A に条件付けた9 branch state と D(alpha,beta) を直接計算した。branch別の独立 Lanczos が同一線形作用素を与えない失敗を保存して棄却し、共通 block-Krylov へ変更した。

102次元主計算では max |D_off|=0.019418、最大規格化 pair coherence=0.179989。114次元対照との差は非対角で最大0.001905。この細分履歴集合に厳密なデコヒーレンスは採用しない。そのため A→B-→B+ の対角重み0.324336を通常の量子履歴確率とは呼ばない。弱い実部整合性については、cubic regrid による branch Gram 相対変化約0.00638と同程度の系統が残るため保留する。

古典 Wigner 4096軌道では第一壁 A 条件付きの A→B-→B+ 固定三時刻頻度は0.353748、ds=.01→.005でラベル変更0。量子対角重みと古典頻度は異なる種類の量であり、量子側のデコヒーレンスなしに確率差とはしない。

次は、(1) 第二壁情報を二分割／より粗い領域へまとめたとき何段階までデコヒーレンスが改善するか、(2) cubic regrid の代替と領域拡大で弱い整合性の系統を下げられるか、(3) 固定時刻セクターではなく「領域へ入る」class operator を用いる timeless / constraint-compatible 構成との差を調べる。環境や非一様自由度を追加したデコヒーレンスは別模型として扱う。

### 6b. 第二壁粗視化：全5集合分割を監査

検証済み9×9 D 行列を再利用し、第二壁 A/B+/B- の全5集合分割をblock和で評価した。A / {B+,B-} の二分割では max |D_off| が0.004799まで下がるが、102→114差を二項和へ伝播した対照尺度0.003809の1.26倍にとどまるため厳密デコヒーレンスは保留。B+単独・B-単独の二分割は対照尺度の5倍以上の複素干渉が残る。第二壁を完全に捨てた最終時刻だけの三分岐は直交射影なので自明に対角であり、動的デコヒーレンスとは数えない。

次は cubic regrid の系統を下げて A/{B±} 二分割を再判定し、その後に固定時刻射影から「領域へ入る」constraint-compatible class operatorへ進む。


### 6c. ever-entered region：finite-clock CAP pilot

固定三時刻の射影から離れ、第一壁A条件付け後にB={B+,B-}領域へ一度でも入ったかを、G_eff=G-iV0 F_B の複素吸収ポテンシャルで監査した。unrestricted/no-entryは各時間刻みで同じblock-Krylov作用素を共有し、regridは線形・非正規化とした。

V0=.025→.075でstrict off-diagonalは.1105→.0073へ低下するが、.085で.0136、.10で.0232へ再増大する。V0=.075,maxdim90でも|Doff|=.00578が残る。width=.20/.30/.45でstrict量は近いがRe Doffは符号を変える。広いregulator plateauがないためstrict decoherenceもweak consistencyも採用せず、branch対角重みをever-entry確率とは呼ばない。

古典Wigner 4096軌道ではA条件付きever-B頻度0.998729、ds=.01/.005一致。量子対角重みと確率比較しない。

次はfinite-clock CAPを負の対照として固定し、timeless WDW constraintと可換なcomplex-potential S-matrix class operatorを実装する。必要ならCAP reflectionをfluxで独立診断する。

### 6d. timeless前段：平方根branchと二階WDW constraintの非同値性

現行branch `i hbar d_s psi=G(s)psi`, `G=-sqrt(A(s))` について、naive constraint `C=(i hbar d_s)^2-A(s)` を作用させると `C psi=i hbar(d_s G)psi` が残る。既存packetで相対mismatchは s=3.0で0.02105、s=4.6で0.07566、s=5.5で0.09042。delta=.02/.01/.005とKrylov 60/72/84で安定した。

よってHalliwell型の二階timeless WDW constraintへ進む場合、それは現行finite-clock平方根量子化の単なる表現変更ではなく新しい量子化選択として扱う。次は小さい3D minisuperspace格子で二階constraintとcomplex-potential S-matrixを構成し、まずclass operatorのconstraint commutatorとfinite-window収束を監査する。physical induced inner product上の確率はその後。


### 6e. timeless S-matrix small-grid pilot：有限Dirichlet箱は不合格

新しい二階constraint C=P_s^2-A(s) を有限3D minisuperspace格子へ置き、B領域のcomplex potentialについて finite-window interaction-picture S_T を計算した。V0=0ではS_T=Iと[C,S_T]=0の負の対照が1e-12級で通る。

near-zero固有modeの境界依存を下げるため、20個のnear-zero modeからedge mass最小packetを作成し、edge massを0.0765まで下げた。それでもV0=.025でcommutator relativeはT=.25/.5/1.0に対し0.280/0.356/0.434と増加し、successive S_T差も0.0467→0.0906へ増加。V0=.05/.10ではさらに悪化。

よってこの有限Dirichlet scattering設計をHalliwellのT→∞ class operator近似として採用しない。次は外側境界をoutgoing/absorbing化するか、Green/resolventによるscattering operatorへ切り替える。induced physical inner productとdecoherence functionalはconstraint-compatibleなscattering収束を確認した後に実装する。


### 6f. timeless stationary resolvent / T-operator pilot：有限Dirichlet resolventも不合格

finite-window S-matrix の反射箱問題を避けるため、同じ二階constraint \(C=P_s^2-A(s)\) に対して
\(G_0(z)=(z-C)^{-1}\), \(G_V(z)=(z-C+iV)^{-1}\), \(T=U+UG_VU\), \(U=-iV\) の
stationary pilotへ切り替えた。20 near-zero modeからedge mass最小packetを作り、edge massは0.07652。

Dyson identity residualは \(10^{-14}\) 級で通る一方、\(\eta=.10/.05/.025\) でfree resolvent normは
9.93/19.47/36.42へ増大。T-action successive changeもV0=.025で0.137→0.211、V0=.05で
0.234→0.331と増加した。response edge massは約0.078–0.081で急増しないが、
\(\eta\to0^+\) の安定域は確認できない。

よってこの有限Dirichlet resolvent/T-operatorをHalliwell型on-shell S-matrix/class operatorの
近似として採用しない。finite-window S-matrixとstationary resolventの二つの負の対照を保存する。
次はoutgoing/PML outer boundaryまたはcontinuum spectral densityを扱えるscattering discretizationを実装し、
box/cap/boundary依存とconstraint commutationを先に監査する。induced physical inner productと
decoherence functionalはその後。


### 6g. outer absorbing-layer control：境界依存は改善するが未収束

finite-window S-matrixとstationary Dirichlet resolventの負の対照を受け、
B領域のclass-operator用complex potentialとは別に、有限箱外側の反射を抑える
数値absorbing layer \(W_{out}\) をbackgroundへ共通に入れた。

3-cell layerではouter strengthを0→.20とすると、V0=.025の \(\eta\) successive
T-action changeが0.137/0.211→0.076/0.081へ低下し、\(\eta=.025\) のbackground
one-cell edge massも0.0841→0.0244へ低下した。したがってDirichlet外壁が前段の
不安定性へ寄与していた可能性は高い。

一方、3-cell profileはreference packetの23.3%と重なる。1/2/3-cell対照では
1-cellでも改善は残るが \(\eta\to0\) 方向でT-action差が再増大する。
\(\eta=.025,\gamma=.20\) のlayer 1→2 / 2→3差はV0=.05で0.117/0.098、
outer strength依存も約0.10以上残る。

よってこのouter CAPを都合のよいstrength/layerへ固定してHalliwell型on-shell
class operatorとはしない。次はtrue PML / exterior complex scaling / continuum
spectral densityを扱うscattering discretizationへ進み、incoming/outgoing fluxなどで
reflectionを独立測定する。constraint-compatibleなscattering limitを確認するまで
induced physical inner productとdecoherence functionalは実装しない。


### 6h. physics audit gate / claim compiler：基盤を実装

量子化選択の「正解」をLeanへ決めさせるのではなく、GR reduction、Lean semantics、clock、inner product/domain、factor ordering、regulator、known limit、claim compilerを別gateへ分解する。

`PHYSICS_AUDIT.md` と `audit/physics_obligations.json` をsource of truthとし、blocking gateにPENDING/PARTIAL/FAILが残れば物理解釈の昇格を禁止する。現在はregulator gateがFAIL、known-limitがPARTIAL、CAS/clock/domain/orderingがPENDING。

Leanには3つの抽象意味論定理を追加するが、continuous operatorのdomain/self-adjointnessや量子化の物理的正しさは証明対象外。xAct/Cadabraによる独立GR reductionはCAS gateの次段階としてPENDINGのまま残す。


### 6i. Bianchi IX GR reduction CAS：Cadabra実行済み、第二backend待ち

数値実装が採用するBianchi IXの古典縮約を、SU(2) Euler角の3-metricから独立再導出するcheckerを追加した。Cadabra 2.5.14を公開SHA-256で固定し、Cadabraプロセス内のSymPy scalar backendで明示的3×3 Christoffel/Ricci計算を行った。

確認済み:
- ({}^{(3)}R=-12e^{-2\alpha}V)
- (K_{ij}K^{ij}-K^2=6N^{-2}(-\dot\alpha^2+\dot\beta_+^2+\dot\beta_-^2))
- repo規約のnormalized Lagrangian / canonical momenta
- (C=\frac12e^{-3\alpha}(-p_\alpha^2+p_+^2+p_-^2)+e^\alpha V)
- (W=2e^{-4s}V), (G_s=-\sqrt{p_+^2+p_-^2+W})
- classical/quantum wall実装とgradient、reduced Hamilton equations

Cadabra backendはPASS。一方、xAct/xCoba sourceは用意したがlicensed Wolfram runtimeがCIにないため未実行。したがって physics audit のGR_REDUCTIONはPENDINGからPARTIALへ進めるが、完全PASSにはしない。

次はxActを実際に走らせる実行環境を接続するか、同等に独立した第二CAS/解析certificateを追加する。その後もclock/domain/ordering/regulatorの各gateは別に残る。


### 6j. 第二独立CAS backend：Maxima/ctensorでGR_REDUCTIONをPASS

Cadabra実行が内部でSymPy scalar backendを使うため、第二backendに同じSymPy系を重ねず、
Ubuntu 24.04のMaxima 5.46.0 + ctensorを採用した。

同じSU(2) Euler角metricからMaxima ctensorで `cmetric` / `christof` / `ricci` /
`scurvature` を実行し、Misner変数代入、ADM kinetic、canonical momenta、
Legendre transformもMaxima自身で計算した。

symbolic deltaは全項0。repo実装比較でもwall relation差は最大約2.33e-21、
classical/quantum wall gradient差0、reduced Hamilton rhs差0。
Cadabra 2.5.14とMaxima/ctensor 5.46.0の二独立実行が同じconstraintへ到達したため、
Bianchi IXの `GR_REDUCTION` gateをPASSへ昇格する。xActは任意の第三cross-checkとして
sourceのみ保持し、未実行であることを明示する。

### 6k. ブラックホール本筋へ戻すblocking gate

Bianchi IXはSchwarzschild内部そのものではなく、異方性を増やしたsingularity stress test。
この系列だけを続けてBH特異点の主張へすり替えないため、
`known_limits.schwarzschild_ks_bridge` をblocking obligationとして追加した。

次は [BLACK_HOLE_RETURN_GATE_ja.md](BLACK_HOLE_RETURN_GATE_ja.md) に従い、

1. Schwarzschild interior → Kantowski–Sachs の古典写像を独立CASで導出
2. 既存KS WDW `C=p_T^2-p_x^2-e^{-2x}=0` への正準変換・規格化を照合
3. L2 positive-frequency branch とKG/induced-inner-product sectorで同じ物理状態を対応付け
4. 質量を含む4次元曲率observableとdomainを固定
5. clock / ordering / regulator gateをKS側でも通す

の順で本来のBH singularity questionへ戻る。

Bianchi IXのtimeless scattering regulator gateがFAILである事実は保持し、
BH側で同じ数値設計を無条件に移植しない。


### 6l. Schwarzschild interior → Kantowski–Sachs classical bridge：PASS

Bianchi IXからBH特異点の本筋へ戻すblocking gateを実物化した。

Schwarzschild内部
[
F=2mu/ho-1,quad
N=F^{-1/2},quad
A=lambdasqrt F,quad
R=ho
]
をKS ansatzへ入れ、Maxima/ctensorとCadabraで独立にADM縮約した。

確認済み:

- (^{(3)}R=2/R^2)
- (K_{ij}K^{ij}-K^2=6N^{-2}(-dotOmega^2+doteta^2))
- normalized ADM Lagrangian / canonical momenta
- (P_Omega^2-P_eta^2+48e^{-2sqrt3Omega}=0)
- (x=sqrt3Omega-log4, T=sqrt3eta) がcanonical
- repo constraint (p_T^2-p_x^2-e^{-2x}=0) はMisner constraintの (-1/3) 倍
- repoの (r=rac14e^{-x-T}) はSchwarzschild areal radius (ho)
- (mu_D=rac14e^{x-T}p_T(p_T-p_x)) はSchwarzschild massへ一致
- ({mu_D,C}=-(p_T/2)e^{x-T}C) なのでconstraint surface上でDirac observable

初回CAS run 35430090228 はMaxima約0.37s、Cadabra約2.50sでPASS。
既存 `wdw_model.classical_path` との比較もconstraint/mass/radius/Schwarzschild relationが (10^{-15}) 級以下で一致。

`known_limits.schwarzschild_ks_bridge` はPASSへ昇格。

次のBH主線:

1. KS clock audit — (T) の単調性、branch sector、代替clock
2. L2 positive-frequency branch とKG/induced-inner-product sectorで同じ状態を対応付け
3. factor-ordering family
4. mass observableと (K=48mu^2/r^6) の量子ordering/domain
5. regulator/boundary/semiclassical gate
6. その後にのみsingularity avoidance/persistenceを評価


### 6m. KS black-hole quantum gates：clock PARTIAL / ordering sensitive / mass-curvature未選定

Schwarzschild/KS古典bridgeの後段を4 gateへ分けて実装した。

1. **clock**
   - (\{T,C\}=2p_T)、repo branch (p_T=-H) では有限xで0にならずclassical T-clockはPASS。
   - (Y=x+T=-\log(4r)) もclassicalには単調。
   - xは (p_x=0) の反射点でglobal clockにならない。
   - Y-clock deparametrizationは (p_X^{-1}) を要求するため量子domain未解決。aggregate PARTIAL。

2. **same-state L2 ↔ KG**
   - positive-frequency KG sectorで
     [
     \Psi_{KG}=\sqrt{h/2}\,H^{-1/2}\chi_{L2}
     ]
     を用いるとfinite boxでnorm・mapped x/r期待値が約1e-14で一致。
   - 同じscalar fieldをそのまま (|\Psi|^2dx) と読むnegative controlは最大約0.297のx期待値差。
   - continuum zero-energy thresholdで (H^{-1/2}) domain未解決。PARTIAL。

3. **factor ordering**
   - (A_q=-h^2e^{-qx}\partial_x(e^{qx}\partial_x)+e^{-2x}) を
     (L^2(e^{qx}dx)) 上のformal-symmetric familyとして監査。
   - flat representationでは (A_q=A_0+h^2q^2/4)、同じclassical principal symbol、差はO(h²)。
   - q=-2..2で最大 (Delta\langle x\rangle\approx0.526)、既存box error threshold 3e-4を大幅に超える。
   - **ORDERING-SENSITIVE** と判定。

4. **mass / Kretschmann**
   - classical (mu=\frac14e^{x-T}p_T(p_T-p_x)) からconstraint-kernel preserving local ordering familyを構成。
   - そのfamily内ではflat kinematical L2 formal-adjoint条件が、d_x / d_T / identity coefficientからそれぞれ異なる Im(d) を要求して衝突。
   - full WDW constraint kernelを保ってもpositive-frequency sectorは保たず、d=0,i/4,i/2,3i/4の最良branch residualは約0.991。
   - finite-box symmetric/Weyl candidateのmass expectation driftはdx=.02で約4.79% / 1.69%だが精密化で低下するためcontinuum no-goとはしない。
   - q=2,4,6 exponential-tail regulatorは既存監査で不合格。量子 (widehat\mu), (widehat K) は未選定。PARTIAL。

次は、positive-frequency physical KG Hilbert spaceを保つnonlocal/Dirac mass observableを構成し、path-integral/symmetry原理でordering familyを狭めた後、Y-clockのinverse-momentum domainとcommon-state semiclassical limitを監査する。


### 6n. positive-frequency mass Dirac candidate → Y-clock → Kretschmann

局所 `M_d` familyがpositive-frequency sectorを保たなかったため、selected sector内部からmass observableを作り直した。

[
V=(i/h)[H,X],qquad
M_0=rac14e^{X/2}H(1+V)He^{X/2},
qquad
M_D(T)=U(T)M_0U(T)^dagger.
]

finite boxではpositive/Hermitian、physical-KG adjointness、mass expectation保存がPASS。
hとgridを同時精密化し、同じGaussian Wigner ensembleのclassical mass平均と比較するとrelative errorは
0.203→0.0523→0.0170→0.00826 と低下した。packet-center比較は別diagnosticとして保持する。

ただしcontinuum quadratic-form closureとordering uniquenessは未証明なのでmass gateはPARTIAL。

同じDirac stateを (Y=x+T=-log(4r)) のnull sliceへpull backすると、
KG norm差最大2.93e-4、same mass flux差最大4.29e-6で一致。
一方deparametrized Y Hamiltonianは (p_X^{-1}) を含み、right=12→32で
min p_X eigenvalueは1.92e-5→1.42e-7、inverse normはoperator levelで増大する。
従ってsame-state slice比較はPASS、standalone Y-time domainはPENDING。

最後に選んだmass candidateを

- (48|r^{-3}M_Dpsi|^2)
- (48|M_Dr^{-3}psi|^2)
- (48|r^{-3/2}M_Dr^{-3/2}psi|^2)

へ入れた。T=2、right=12→32で3 formともlog10値が約62–66 decades増え、
mass-applied threshold overlapもright=24, dx=.04/.02で約0.3724/0.3730と非零にresolveされた。
従って**この3つのexplicit-radius positive formではKretschmann regulator/domain gateはFAIL**。

これは全てのnonlocal relational curvature observableへのno-goではない。
次はmass formのcontinuum closability / associated self-adjoint operator、
Y-clock zero-mode domain、そして物理的に動機づけられたnonlocal curvature候補を狭く検討する。


### 6o. KS continuum threshold / physical-KG completion：KG gateをPASSへ更新

positive-frequency KS Hamiltonian `H^2=-h^2 d_x^2+e^(-2x)` を continuum Liouville spectrum で監査した。完全規格化一般化固有関数は `phi_k=pi^(-1)sqrt(2k sinh(pi k)) K_(ik)(e^(-x)/h)`, `E=hk`。
zero thresholdでは `phi_k/k -> sqrt(2/pi) K0(e^(-x)/h)`。現行 h=.2 Gaussian のthreshold coefficientは `|B0|=0.0143884187466872` で3つのdomain/grid controlに相対約2.1e-15で安定し、`c(k)=B0 k+O(k^3)`。従ってreference packetは ordinary-L2 の `D(H^-1/2)` と `D(H^-1)` に入り、generic nonzero-B0 stateは `H^-3/2` でlog threshold obstructionへ達する。

positive-frequency KG spaceは `D(H^1/2)` を `||Psi||_KG^2=(2/h)||H^1/2 Psi||^2` で完成して定義でき、`S=sqrt(h/2)H^-1/2` はL2からこの完成空間へのunitary mapへ延長できる。したがって `inner_product.ks_continuum_completion` は **PASS** へ更新。ordinary-L2上で `H^-1/2` がunboundedである事実は保持し、全KG amplitudeをordinary-L2関数とは扱わない。

mass側では `V=(i/h)[H,X]` が `HV+VH=2P` を満たし、`|P|<=H` から `-I<=V<=I`、従って `I+V>=0` をcontinuumで得る。ただしfinite-boxの min eig(I+V)は right=8/12/16 で 2.79e-4 / 4.69e-5 / 1.37e-5 と0へ向かう。よってstrict positive gapを使ったmass form closureは採用しない。

次の主課題は (1) `q_M=(1/4)<H e^(X/2)psi,(I+V)H e^(X/2)psi>` のclosability / closed form / associated self-adjoint operator、(2) same Dirac stateについてのY-clock `p_X^-1` continuum domain、(3) その後のnonlocal relational curvature候補。Kretschmann explicit-radius 3 formのFAILは変更しない。


### 6p. KS mass affine-covariance selection obstruction

PR #20後のmass form closability課題に対し、まずclassical massのaffine構造を固定した。\(E=\sqrt{p_x^2+e^{-2x}}\), \(u=\operatorname{arsinh}(p_xe^x)\) で \(\{u,E\}=1\)、\(\mu_0=(E/4)e^u\)、\(\{E,\mu_0\}=-\mu_0\)。null variables \(U=-e^{-x-T}/4, W=e^{-x+T}/4\) ではconstraintは \(C=-4UW(p_Up_W-4)\)、massは \(\mu=p_Tp_W/8\)。

positive-frequency Hilbert spaceでsemibounded \(H\ge0\) と非零positive self-adjoint \(M\) に exact covariance \(e^{-iHt/h}Me^{iHt/h}=e^{-t}M\) を全実tで要求すると、\(Q=\log M\) からWeyl relationが生じ、unitary conjugationで \(H\mapsto H+hs\) を全実sについて要求する。これは \(\sigma(H)=[0,\infty)\) と両立しない。従ってexact affine covarianceをordering一意化の条件には採用しない。

現行finite-box \(M_0\) の [H,M0]=-ihM0 residualはright=8/12/16で約1.620842とO(1)のまま。これは負の対照で、mass closability FAILではない。Cavaglià–de Alfaro–Filippov (gr-qc/9508062) の別Schwarzschild量子化でもhalf-line supportがmass自己共役性を制約するが、そのdeficiency indicesをrepo operatorへ移植しない。

**次は** \(T_M=(I+V)^{1/2}H e^{X/2}\) のclosability / adjoint-domainを直接解析する。mass_continuum_form=PENDING, Y-clock inverse-pX=PENDING, selected explicit-radius Kretschmann=FAILは変更しない。
