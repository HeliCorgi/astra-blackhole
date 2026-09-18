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
