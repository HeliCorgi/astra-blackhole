# Bianchi IX timeless stationary resolvent / T-operator pilot

## 目的

直前の finite-window S-matrix pilot では、有限 Dirichlet 箱のまま
unphysical parameter window を広げると、\([C,S_T]\) と successive window 差が
むしろ増大した。この失敗を負の対照として保存し、時間窓を伸ばさない stationary
resolvent へ切り替える。

Halliwell の complex-potential class operator では、領域へ入らない class operator
をその領域に局在した complex potential の S 行列として構成し、constraint と
可換な量を得ることが要件になる（arXiv:0909.2597, 1108.5991）。

本 pilot は **S 行列そのものをまだ構成しない**。有限箱上の outgoing-resolvent
regulator を使って、T-operator の数値核に \(\eta\to0^+\) の安定化傾向があるかを
先に監査する。

## 模型と数値設定

前段と同じ、新しい二階 WDW constraint 量子化

\[
C=P_s^2-A(s)
\]

を \(q=(s,\beta_+,\beta_-)\) の有限3次元 minisuperspace に置く。変数は
\(s,\beta_+,\beta_-\) の3つだけを保持し、物質、環境、空間非一様自由度は含めない。
数値は既存模型の無次元化規約を継承し、実天体の SI 単位へ同定しない。

主格子は：

- \(s\in[2,6]\), \(\beta_+\in[-3.5,3.5]\), \(\beta_-\in[-4,4]\)
- interior points: \(12\times14\times14=2352\)
- \(\hbar=0.2\)
- potential absolute cap = 8
- B 境界 smooth width = 0.30

B={B+,B-} 領域には従来と同じ \(V=V_0F_B\ge0\) を使う。

\[
H_{\mathrm{eff}}=C-iV,\qquad
G_0(z)=(z-C)^{-1},\qquad
G_V(z)=(z-C+iV)^{-1},
\]

\[
z=E+i\eta,\qquad E=0,\quad \eta>0.
\]

\(\eta\) は outgoing resolvent を有限格子で正則化する数値 regulator であり、
物理パラメータではない。

摂動 \(U=-iV\) に対して有限行列で

\[
G_V=G_0+G_0UG_V,\qquad
T(z)=U+UG_VU
\]

を使う。最初の Dyson identity は実装の代数的 negative control として直接検査する。

## reference packet

20個の near-zero constraint mode の部分空間で、有限箱の outer one-cell projector
を最小化する superposition を数値 reference とした。これは物理状態の選択則ではない。

採用 packet:

- edge mass (1 cell): 0.07652061
- capped-region mass: \(1.78\times10^{-5}\)
- constraint energy mean: -0.01048549
- constraint energy spread: 0.00581048
- \(\|C\psi\|/\|\psi\|\): 0.01198780

## scan

- \(E=0\)
- \(\eta=.10,.05,.025\)
- \(V_0=.025,.05\)

各点で Dyson residual、free resolvent norm、\(\|T\psi\|\)、Born 項との比、
near-zero shell 投影率、interacting response norm、response edge mass を保存した。

## 実行結果

Dyson identity residual は全点で \(3.6\times10^{-15}\) から
\(1.53\times10^{-14}\) で、有限行列実装の代数的対照は通った。\(V_0=0\) では
T-action norm は0。

一方、\(\eta\) を小さくすると free resolvent norm は
9.93 → 19.47 → 36.42 と増大し、T-action の successive change も減らない。

| V0 | eta .10→.05 | eta .05→.025 |
|---:|---:|---:|
| .025 | 0.13749 | 0.21063 |
| .050 | 0.23375 | 0.33086 |

response の outer one-cell edge mass は約0.078–0.081で支配的には増えていないが、
このことだけでは有限箱の離散 near-zero poles の影響を除外できない。T-action 自体も
\(\eta\) 低下で安定値へ近づく傾向を示していない。

**結論：この finite-Dirichlet stationary resolvent/T-operator pilot でも、
\(\eta\to0^+\) の採用可能な安定域を確認できなかった。Halliwell型の on-shell
S-matrix / class operator として採用しない。**

これは Halliwell formalism の否定ではない。現在の finite box、離散 spectrum、
potential cap、および boundary 条件を含む数値設計の負の結果である。

## 次

次は Dirichlet 箱を固定したまま regulator を増やすのではなく、outer boundary を
outgoing/PML 型に変更するか、continuum spectral density を明示的に扱える scattering
discretizationへ進む。そこで box / cap / boundary scan と on-shell constraint
commutation を先に確認する。

induced/Rieffel physical inner product、decoherence functional、履歴確率は、
constraint-compatible な scattering limit を確認した後に実装する。

## 再現

```sh
python -m unittest discover -s tests -p 'test_bianchi_ix_timeless_resolvent.py' -v
python research/run_bianchi_ix_timeless_resolvent.py \
  --out artifacts/timeless-resolvent \
  --check-against research/bianchi_ix_timeless_resolvent_results/summary.json
```

採用 CI run 35377267600 では4テスト、pilot生成が成功した。artifact
`bianchi-ix-timeless-resolvent` (ID 10560925377, SHA-256
`7dcce55c35bd89bc44a1f44208d78e84a3918d32de96f1157b42c4895f1622b8`) を保存した。
この成功は数値再現性の記録であり、物理的妥当性の証明ではない。

## 限界

- 一様真空 Bianchi IX minisuperspace の新しい二階 WDW constraint 量子化。
- 既存 finite-clock square-root 量子化とは同一視しない。
- finite Dirichlet box と potential cap を継承。
- \(\eta\), shell projection, Euclidean grid norm は数値診断。
- reference packet は edge mass を小さくする数値選択で、物理状態処方ではない。
- induced/Rieffel physical inner product は未実装。
- decoherence functional と physical history probability は未計算。
- 現実のブラックホール観測、特異点解消、量子重力検証ではない。
