# Bianchi IX timeless resolvent: outer absorbing-layer pilot

## 目的

finite-window S-matrix と stationary resolvent の両方で、有限 Dirichlet 箱から
Halliwell 型 scattering limit へ収束する挙動を確認できなかった。

次の切り分けでは、B領域に置く class-operator 用 complex potential とは別に、
外側境界の反射だけを減らす数値 absorbing layer を導入する。

この outer absorber は **物理ポテンシャルではない**。また、本実装を PML
(perfectly matched layer) とは呼ばない。有限 Dirichlet 壁は吸収層のさらに外に残る。

## 定義

前段と同じ二階 constraint \(C=P_s^2-A(s)\) と B 領域 potential
\(V_B=V_0F_B\) を用いる。

outer one-cell から内側へ3セルの smooth profile \(F_{out}\) を作り、

\[
W_{out}=\gamma F_{out}\ge0
\]

とする。background と B-absorbing 系の両方に同じ outer absorber を入れる：

\[
H_{bg}=C-iW_{out},
\]

\[
G_{bg}(z)=(z-C+iW_{out})^{-1},
\]

\[
G_{full}(z)=(z-C+iW_{out}+iV_B)^{-1}.
\]

\(U_B=-iV_B\) とすれば

\[
G_{full}=G_{bg}+G_{bg}U_BG_{full},
\qquad
T_B=U_B+U_BG_{full}U_B.
\]

したがって outer absorber の効果を background に共通化し、B領域の T-action だけを
比較する。

## scan

主格子は前段と同じ2352次元。

- \(E=0\)
- \(\eta=.10,.05,.025\)
- \(V_0=.025,.05\)
- outer strength \(\gamma=0,.05,.10,.20\)
- outer layer = 3 cells
- reference = 20 near-zero modes の edge-minimized packet

判定は「ある1点の \(\gamma\) で改善したか」ではなく、

1. \(\eta\) 半減時の T-action change が複数の \(\gamma\) で低下するか
2. T-action が \(\gamma\) 自体へ強く依存しないか
3. background/response の edge・outer-layer support が制御されるか

を見る。

条件を満たさない場合、outer CAPを都合のよい強さへ調整して採用しない。

## 未実施

- 独立な incoming/outgoing flux による reflection 測定
- 真の PML / exterior complex scaling
- continuum spectral density
- induced physical inner product
- decoherence functional / history probability

## 再現

~~~sh
python -m unittest discover -s tests -p 'test_bianchi_ix_timeless_outer_cap.py' -v
python research/run_bianchi_ix_timeless_outer_cap.py --out artifacts/timeless-outer-cap
~~~

CI成功は数値再現性だけを示し、constraint-compatible class operator の物理的妥当性を
保証しない。


## 実行結果（2026-09-19）

模型は一様真空 Bianchi IX minisuperspace の新しい二階 constraint 量子化で、
保持する変数は \(s,\beta_+,\beta_-\) の3つ。物質、環境、空間非一様自由度は
含めない。既存コードの無次元化規約を継承し、実天体の SI 単位へ同定しない。

主格子は \(12\times14\times14=2352\) interior points、\(\hbar=0.2\)、
\(s\in[2,6]\)、\(\beta_+\in[-3.5,3.5]\)、\(\beta_-\in[-4,4]\)。
potential absolute cap=8、B境界width=.30。reference stateは20 near-zero modes
からouter one-cell massを最小化したpacketで、edge mass=0.0765206、
\(\|C\psi\|/\|\psi\|=0.0119878\)。

### 3-cell outer layer の主scan

Dyson identity residualは全点で丸め誤差級、\(V_0=0\) の T-action は0。
外側吸収を強めると、有限Dirichlet resolventで悪化していた \(\eta\) 感度と
edge supportは系統的に低下した。

| V0 | outer strength \(\gamma\) | T change .10→.05 | T change .05→.025 |
|---:|---:|---:|---:|
| .025 | 0 | .13749 | .21063 |
| .025 | .05 | .10822 | .13983 |
| .025 | .10 | .09214 | .10840 |
| .025 | .20 | .07608 | .08132 |
| .050 | 0 | .23375 | .33086 |
| .050 | .05 | .18485 | .22244 |
| .050 | .10 | .15721 | .17237 |
| .050 | .20 | .12876 | .12809 |

\(\eta=.025\) では \(\gamma=0\to.20\) によりbackground one-cell edge massが
0.08410→0.02441へ低下した。V0=.025 の interacting response edge massも
0.08104→0.02769へ低下する。したがって有限箱の外側反射が前段の不安定性へ
寄与していた可能性は高い。

ただし、\(\eta=.025\) の T-action は outer strength 自体にも依存する。
V0=.025で \(\gamma=0\to.05/.05\to.10/.10\to.20\) の相対変化は
0.1177/0.0639/0.0689、V0=.05では0.1762/0.1007/0.1113。
広い \(\gamma\)-independent plateau とは言えない。

### layer幅対照

3-cell layerはreference packetの0.2330と重なる。改善が境界反射だけでなく
reference stateそのものの減衰を含む可能性を切り分けるため、1/2/3-cellを比較した。

- 1-cell profile overlap: 0.07652
- 2-cell profile overlap: 0.14559
- 3-cell profile overlap: 0.23300

狭い1-cellでも \(\gamma=.20\) により前段より改善するが、V0=.025の
T successive changeは0.1098→0.1460、V0=.05では0.1859→0.2274と、
\(\eta\) を下げる側で再び増える。

さらに \(\eta=.025,\gamma=.20\) ではlayer 1→2 / 2→3 のT-action差が
V0=.025で0.0741/0.0631、V0=.05で0.1172/0.0983残る。
したがって3-cell・\(\gamma=.20\) の比較的平坦な値だけを選んで
「収束した」と判定することはしない。

**結論：数値outer CAPはedge supportと \(\eta\) 感度を明確に改善するが、
strength・layer幅に対するregulator-independentなscattering limitを確認できない。
このpilotからHalliwell型on-shell class operatorを採用しない。**

これはHalliwell formalismの否定ではなく、現在のfinite grid・Dirichlet wall・
outer CAP・potential capを含む数値設計についての結果である。decoherence
functional、induced/Rieffel physical inner product、履歴確率はまだ計算しない。

## 次

outer CAPを強くするscanはここで止める。次は真のPML、exterior complex scaling、
またはcontinuum spectral densityを明示的に扱うscattering discretizationへ進む。
その前後でincoming/outgoing fluxまたは別の散乱診断によりreflectionを独立測定し、
box / cap / boundary依存とconstraint commutationを先に確認する。

## 再現と検証記録

~~~sh
python -m unittest discover -s tests -p 'test_bianchi_ix_timeless_outer_cap.py' -v
python research/run_bianchi_ix_timeless_outer_cap.py --out artifacts/timeless-outer-cap
python research/run_bianchi_ix_timeless_outer_cap_layer_scan.py --out artifacts/timeless-outer-cap
python research/check_bianchi_ix_timeless_outer_cap_results.py \
  --artifact-dir artifacts/timeless-outer-cap \
  --baseline research/bianchi_ix_timeless_outer_cap_results/summary.json
~~~

採用前pilotのCI run 35378231923、artifact 10561071426、
SHA-256 \`32034f50fdba0184c9e82d41a133c040aa43a411b48f25c4f933ae4f0d8617ef\`
をbaseline生成元として保存した。公開版CIではこの保存summaryへ数値回帰する。
CI成功は数値再現性の記録であり、物理的妥当性の証明ではない。
