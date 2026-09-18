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
