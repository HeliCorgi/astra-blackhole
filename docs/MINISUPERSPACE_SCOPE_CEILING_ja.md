# Minisuperspace hard scope ceiling

更新: 2026-09-20

## 結論

このrepositoryのKS / Bianchi IX / WDW量子計算には、ordering、regulator、clock、inner product、operator domainよりさらに外側にある制約がある。

**最初に一般相対論の無限自由度の場をhomogeneous minisuperspaceへ縮約し、その後で量子化している。**

したがって、仮に現在のminisuperspace WDW理論についてclock、physical inner product、factor ordering、self-adjoint operator domain、regulator removal、decoherent histories、semiclassical limitがすべて数学的に完成しても、それだけから full quantum general relativity が完成したとは結論しない。

これは未解決gateではなく、**hard scope ceiling** である。内部gateがすべてPASSしても消えない。

## 1. 何を先に捨てているか

full GRのcanonical configurationは局所場 g_ij(x) であり、一般には無限個の空間自由度を持つ。

このrepositoryでは先に対称性縮約する。

- Kantowski–Sachs: homogeneous spherical-interior sector
- Bianchi IX: homogeneous anisotropic sector

従って、残るconfiguration variablesは有限個である。

この縮約により、少なくとも次は量子理論のdynamical degrees of freedomとして存在しない。

- local gravitational-wave modes
- generic inhomogeneous metric modes
- generic inhomogeneous matter modes
- spatially local mode coupling
- local backreaction between homogeneous and discarded inhomogeneous sectors

Bianchi IXがKSより自由度を増やしても、homogeneous minisuperspaceである点は変わらない。

## 2. full local constraint algebraも存在しない

full canonical GRでは空間各点に H_perp(x), H_i(x) というHamiltonian / spatial-diffeomorphism constraintsがあり、そのPoisson algebraはmetric-dependent structure functionsを含む。

minisuperspaceではhomogeneityを課した時点で、この無限個のlocal constraintsは有限個のreduced constraintsへ縮退する。

したがってrepoで reduced WDW constraintを正しく量子化した、reduced constraint kernelをoperatorが保った、Leanで抽象的constraint-kernel preservationを証明した、reduced class operatorがconstraintと可換だったとしても、full quantum constraint algebraのclosure/anomaly問題を解いたことにはならない。

その問題自体がこのreduced modelには存在しない。

## 3. ordering / regulatorとの違い

ordering、regulator、clock、domainは、現在選んだminisuperspace量子模型の**内部**での不確定性である。

minisuperspace truncationはそれより一段外側で、full GR → symmetry-reduced finite-dimensional model → quantization の最初の矢印に入る仮定である。

従って、ordering dependenceを消した、regulator-independent limitを得た、physical Hilbert spaceを完成した、unique self-adjoint mass/curvature observableを選んだとしても、discarded local degrees of freedomの影響が小さいことは自動では示されない。

## 4. claim compilerで永久に禁止する昇格

machine-readable auditでは scope_ceiling.kind = MINISUPERSPACE, hard = true とする。

claim compilerは内部gateと独立に MINISUPERSPACE-SCOPED を付ける。

このlabelは、全blocking gateが将来PASSしても消えない。

少なくとも以下への昇格は禁止する。

- FULL QUANTUM GR
- ANOMALY-FREE QUANTUM CONSTRAINT ALGEBRA
- FULL LOCAL BLACK-HOLE QUANTUM DYNAMICS

一方、selected KS minisuperspace quantizationでのoperator-domain result、selected Bianchi IX minisuperspace quantizationでのhistory/decoherence result、selected symmetry-reduced modelでのclock dependence、selected reduced theoryでのregulator resultのようなscoped claimは、内部gateが十分に通れば可能である。

## 5. 次にfull theoryへ近づけるなら

このscope ceilingを下げるには、現在のminisuperspace内部監査をさらに精密化するだけでは足りない。別のmodel extensionが必要になる。

1. KS/Bianchi background上に少数のinhomogeneous perturbative modesを追加する。
2. homogeneous sectorとinhomogeneous modesのbackreactionを明示する。
3. midi-superspaceへ進み、少なくとも1つの空間座標依存を残す。
4. local Hamiltonian/diffeomorphism constraintsとそのquantum algebraを扱う。
5. truncationを拡張したとき、現在のmass/history/curvature conclusionsが安定かを比較する。

これは現在の研究を無効化するものではなく、**現在の結果がどのmodel classまで届くかを固定する境界条件**である。

## 6. 現在の正しい表現

強すぎる表現:

- WDW量子重力を完成した
- full quantum GRで特異点を解決した
- 量子重力のconstraint algebraを検証した

許される表現:

- 選んだhomogeneous minisuperspace WDW量子化の内部で、指定したclock / inner product / ordering / domainについて結果を得た
- reduced constraintに対する数学的・数値的整合性を監査した
- full local degrees of freedomとfull constraint algebraはmodel外であり、結論をそこへ昇格しない
