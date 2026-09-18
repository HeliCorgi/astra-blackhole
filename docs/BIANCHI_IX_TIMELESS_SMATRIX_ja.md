# Bianchi IX timeless complex-potential S-matrix: small-grid pilot

## 目的

finite-clock CAP監査では、内部時計sを選んだまま「B領域へ一度でも入った」を構成したが、CAP強度に対する広いplateauを得られなかった。

Halliwellのdecoherent-histories量子宇宙論では、WDW型constraint H Psi=0 に対し、領域Deltaへ入らないclass operatorをDeltaに局在したcomplex potentialのS-matrixとして構成する。重要な要件はclass operatorがconstraintと可換になること。arXiv:0909.2597 と 1108.5991 を参照。

前段のfactorization監査で、既存の時間依存平方根branchはnaive二階WDW constraintの単なる枝ではないことを確認した。したがって本pilotは**新しい量子化**として扱う。

## constraint

小さい有限3D minisuperspace box q=(s,beta+,beta-) にDirichlet差分を置く。

C = P_s^2 - A(s)
  = hbar^2 L_s - hbar^2(L_+ + L_-) - W(s,beta).

Wは数値pilotのため絶対値capを置く。これは既存finite-clock模型のpotential regulatorと同じではなく、新しい数値regulator。

B領域は従来と同じ

3 beta_+ + sqrt(3)|beta_-| > 0

で、complex potential V=V0 F_B を置く。

## finite-window S-matrix

無限のunphysical parameter時間を直接扱えないので、t1=-T,t2=Tで

S_T = exp(+i C T/hbar)
      exp[-i(C-iV)2T/hbar]
      exp(+i C T/hbar)

を計算する。

Halliwell型の正しいtimeless class operatorはT→∞のscattering limitでconstraintと可換になる。したがってpilotの最初の合格条件は確率ではなく：

1. V0=0でS_T=Iとなること
2. CがHermitianでnear-zero modeを作れること
3. T増加で [C,S_T] の作用が減る傾向
4. S_T自体のT依存が減ること
5. box/cap依存を後で下げられる見込み

を確認すること。

## 何をまだしないか

- induced/Rieffel physical inner product
- zero-constraint continuum delta normalization
- D(noB,everB)の物理確率解釈
- infinite-window極限の主張

Euclidean grid normは数値診断だけ。

このpilotが収束しなければ、有限Dirichlet boxでscattering S-matrixを近似する設計自体を負の対照として保存し、absorbing/outgoing boundaryまたは直接resolvent法へ切り替える。


## Pilot実行結果

最初のnear-zero固有modeは境界1セルに約27.9%の質量を持ったため、そのままの判定を棄却した。20個のnear-zero modeの部分空間で境界projectorを対角化し、edge massを最小化したlow-energy packetへ変更した。採用packetは：

- edge mass (1 cell): 0.07652
- capped-region mass: 1.78e-5
- constraint energy mean: -0.01049
- constraint energy spread: 0.00581

V0=0では S_T=I の誤差4.4e-14以下、commutator relative 1.2e-12以下で負の対照は通った。

しかしV0>0ではTを増やしてもconstraint可換性へ近づかない。例としてV0=.025：

| T | norm² | constraint spread | commutator relative |
|---:|---:|---:|---:|
| .25 | .9259 | .00935 | .2799 |
| .50 | .8597 | .01147 | .3562 |
| 1.00 | .7458 | .01434 | .4336 |

さらにS_Tのsuccessive window差も .25→.5 で0.0467、.5→1.0で0.0906へ増える。V0=.05/.10では悪化はさらに大きい。

**結論：この有限Dirichlet箱によるscattering S-matrix近似はHalliwellのT→∞ limitへ収束している証拠がなく、constraint-compatible class operatorとして不採用。**

これはHalliwell formalism自体の否定ではない。有限箱が反射境界を持つため、無限scattering問題の数値近似として不適切である可能性が高い。次はoutgoing/absorbing outer boundary、またはresolvent/scattering formulationへ変更し、このpilotを負の対照として保持する。
