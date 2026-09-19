# KS physical mass / Y-clock / Kretschmann follow-up audit

更新: 2026-09-20

## 結論

Schwarzschild/Kantowski–Sachs量子監査の次段として、

1. positive-frequency physical-KG Hilbert spaceを保つmass Dirac candidate
2. 同じDirac stateを使った \(T\)-slice と \(Y=x+T=-\log(4r)\)-slice の比較
3. そのmass candidateを実際に入れた \(48\mu^2/r^6\) のpositive quadratic-form ordering

を順番に監査した。

結果は:

| 対象 | 判定 |
|---|---|
| finite-box positive-frequency mass Dirac candidate | PASS |
| physical-KG adjointness / conserved mass expectation | PASS |
| mass continuum quadratic-form closure | PENDING |
| mass ordering uniqueness | PENDING |
| same-Dirac-state \(T\leftrightarrow Y\) KG flux | PASS |
| standalone \(Y\)-time \(p_X^{-1}\) domain | PENDING |
| selected-mass Kretschmann positive forms | FAIL |
| Kretschmann regulator removal | FAIL |

従って「量子mass operatorが完全に確定した」でも、
「量子Kretschmannが一般に発散すると証明した」でもない。

一方、以前の「mass operatorをまだ一つも具体化できていない」という状態からは進み、
**positive-frequency sector内で使える明示的finite-box Dirac candidate**を固定した。

## 1. 出発点

古典Schwarzschild/KS bridgeで

\[
C=p_T^2-p_x^2-e^{-2x}=0
\]

および

\[
\mu
=
\frac14e^{x-T}p_T(p_T-p_x)
\]

を得ている。

repoのpositive-frequency conventionは

\[
ih\partial_T\psi=H\psi,
\qquad
H=\sqrt{-h^2\partial_x^2+e^{-2x}},
\]

したがってcanonical momentumでは

\[
p_T=-H
\]

のbranch。

以前試した局所 \(\widehat M_d\) familyはfull WDW constraint kernelを保ったが、
選択したpositive-frequency sectorを保たなかった。
そこで今回、full WDW上の局所orderingからではなく、
**positive-frequency physical Hilbert spaceそのものからDirac observableを構成する**
方針へ切り替えた。

## 2. positive-frequency mass candidate

\(X\) をmultiplication by \(x\) とし、

\[
V=\frac{i}{h}[H,X]
\]

をHeisenberg velocityとする。

classical principal symbolでは

\[
V_{\rm cl}=\frac{p}{H}.
\]

そこでreference clock \(T=0\) で

\[
\boxed{
M_0=
\frac14
e^{X/2}H(1+V)He^{X/2}
}
\]

を定義した。

classical symbolは

\[
\frac14e^x H^2
\left(1+\frac pH\right)
=
\frac14e^xH(H+p),
\]

すなわちSchwarzschild bridgeで得たclassical massと一致する。

### positive form

finite boxのreference case \(h=.2, dx=.04\) では

\[
\lambda_{\min}(1+V)
\simeq1.67\times10^{-6}>0,
\]

また構成した \(M_0\) matrixの最小固有値も

\[
\lambda_{\min}(M_0)
\simeq3.14\times10^{-4}>0.
\]

matrix Hermiticity residualは約

\[
3.18\times10^{-13}.
\]

この有限次元regularizationではpositive self-adjoint matrixとして機能する。

これはcontinuum \(1+V\ge0\) やclosed quadratic formの証明ではない。

## 3. Dirac / relational transport

reference-clock operatorを

\[
\boxed{
M_D(T)=U(T)M_0U(T)^\dagger
},
\qquad
U(T)=e^{-iHT/h}
\]

でtransportする。

この定義ではselected positive-frequency Hilbert spaceの内部で完結する。

同じSchrödinger stateを

\[
\chi(T)=U(T)\chi(0)
\]

とすると

\[
\langle\chi(T),M_D(T)\chi(T)\rangle
=
\langle\chi(0),M_0\chi(0)\rangle
\]

は有限次元では構成上保存される。

実測した最大relative driftは \(h=.2\) caseで約

\[
1.34\times10^{-13}.
\]

これは、以前の局所 \(\widehat M_d\) がpositive-frequency sectorを保たなかった問題を
回避する。

## 4. physical KG representation

前段監査で

\[
S=\sqrt{\frac h2}H^{-1/2}
\]

によりL2 positive-frequency stateとKG field amplitudeを対応させた。

今回massも

\[
\boxed{
M_{KG}(T)=S M_D(T)S^{-1}
}
\]

と同時に写す。

finite boxでは

- L2とKGのmass expectation差: 約 \(10^{-14}\) 以下
- KG adjoint bilinear residual: reference \(h=.2\) で約 \(1.25\times10^{-13}\)
- mass expectation conservation: 約 \(10^{-13}\)

となった。

従って**finite-box same-state physical-KG representationではPASS**。

ただしcontinuumで \(H^{-1/2}\) 自体のdomain問題が既に残っているため、
これだけでcontinuum physical observableとは呼ばない。

## 5. semiclassical check

mass candidateのsemiclassical controlでは、packet中心だけではなく、量子packetと同じGaussian Wigner ensembleにclassical mass Dirac observableを作用させた平均を比較対象にする。packet中心massは有限\(h\)でspreadを無視するため、別diagnosticとして残す。\(h\) を下げると同時にgridも精密化して比較した。

| h | dx | relative error vs same classical Wigner ensemble |
|---:|---:|---:|
| .4 | .08 | 0.2028 |
| .2 | .04 | 0.05227 |
| .1 | .02 | 0.01697 |
| .05 | .01 | 0.00826 |

同一Wigner ensembleとの差は単調に低下した。packet-centerだけとの比較も 0.5076→0.1711→0.0720→0.0350 と低下するが、primary semiclassical comparatorには採用しない。

従ってこのcandidateについては、**same-state semiclassical correspondenceの一部に肯定的evidence**がある。

ただしcurvature observableまで含むcommon-state WKB/HJ gateではないため、
project全体の `SEMICLASSICAL LIMIT PASSED` は出さない。

## 6. finite-box regulator diagnostic

\(h=.2,dx=.04\) のreference packetでmass expectationを右box端だけ変えると

| right | mass |
|---:|---:|
| 12 | 0.01990091416 |
| 16 | 0.01990091637 |
| 24 | 0.01990091712 |

で、moderate box rangeでは安定している。

したがってmass expectationそのものには、
今回の範囲でKretschmannほどのbox explosionは見えない。

## 7. continuum mass domainはまだPENDING

candidate formは概念的には

\[
q_M[\psi]
=
\frac14
\langle
H e^{X/2}\psi,\,
(1+V)H e^{X/2}\psi
\rangle.
\]

candidate coreとして

- smooth compact support
- finite energy
- \(e^{X/2}\psi\) が上のcomposed form domainに入るstate

を想定できる。

transported domainは

\[
D(M_D(T))=U(T)D(M_0)
\]

と置ける。

しかし未証明なのは

- continuum formのclosability
- closed formへのcompletion
- associated self-adjoint operator
- orderingの一意性

である。

したがってmass gateは **PARTIAL**。

KuchařのSchwarzschild canonical reductionではmassがembedding-independentなcanonical
variableとして残ることが知られている（Phys. Rev. D 50, 3961; arXiv:gr-qc/9403003）。
またAshtekar–Tate–UgglaはminisuperspaceでDirac observablesとdeparametrizationを
明示してから量子解釈を行う枠組みを論じている（arXiv:gr-qc/9302027）。
今回の \(M_D\) はそれらの一般原理と整合するよう設計したが、
同一operatorを文献から導出したという意味ではない。

## 8. 同じDirac stateをY-clock sliceへ移す

areal-radius clock

\[
Y=T+x=-\log(4r),\qquad X=x
\]

を使う。

canonical transformは

\[
p_Y=p_T,
\qquad
p_X=p_x-p_T,
\]

constraintは

\[
-p_X^2-2p_Xp_Y-e^{-2X}=0.
\]

\(Y=const\) sliceはminisuperspace上でnull slice。

ここで**別のGaussianをY-clock theoryへ置き直していない**。
同じpositive-frequency WDW solutionを

\[
T=Y-X
\]

へpull backしてKG fluxを比較した。

### same-state flux

\(h=.2,right=16,dx=.04\) で

- max KG norm difference:
  \[
  2.93\times10^{-4}
  \]
- max same-mass-observable flux difference:
  \[
  4.29\times10^{-6}
  \]
- max imaginary flux artifact:
  \[
  5.14\times10^{-7}
  \]

だった。

従ってfinite boxの同一Dirac state / 同一observable比較では
**T-sliceとY-sliceは整合**。

これは「独立に量子化したT-clock theoryとY-clock theoryがunitarily equivalent」
という証明ではない。

## 9. standalone Y-time theoryのzero-mode問題

constraintから形式的に

\[
p_Y
=
-\frac{p_X^2+e^{-2X}}{2p_X}
\]

なので、Y-time Hamiltonianは \(p_X^{-1}\) を必要とする。

positive finite-box \(p_X\) formを調べると:

| right | min pX eigenvalue | max inverse eigenvalue |
|---:|---:|---:|
| 12 | 1.92e-5 | 5.21e4 |
| 16 | 4.24e-6 | 2.36e5 |
| 24 | 5.58e-7 | 1.79e6 |
| 32 | 1.42e-7 | 7.05e6 |

box拡大でzero thresholdへ向かう。

一方、今回のpacket固有の

\[
\|p_X^{-1}\psi\|
\]

は約174→175で、この有限box系列では発散していない。

従って現在の判定は:

- same-Dirac-state slice comparison: PASS
- standalone Y-time operator: **PENDING**
- zero-mode / continuum domain: **PENDING**

である。

## 10. selected mass candidateをKretschmannへ入れる

古典Schwarzschildでは

\[
K=48\frac{\mu^2}{r^6}.
\]

今回、\(M_D(T)\) を使って3つのpositive quadratic formを宣言した。

### A. mass then radius

\[
K_{MR}[\psi]
=
48\|r^{-3}M_D\psi\|^2.
\]

### B. radius then mass

\[
K_{RM}[\psi]
=
48\|M_Dr^{-3}\psi\|^2.
\]

### C. balanced

\[
K_S[\psi]
=
48\|r^{-3/2}M_Dr^{-3/2}\psi\|^2.
\]

全てclassical principal symbolは \(48\mu^2/r^6\) に対応するが、
operator domainは異なる。

## 11. box removal結果

\(h=.2,dx=.04,T=2\) でrightを12→16→24→32へ広げた。

### mass then radius

\[
\log_{10}K:
30.46\to41.59\to64.98\to92.58.
\]

約62.1 decades増加。

### radius then mass

\[
38.38\to51.14\to77.40\to104.13.
\]

約65.8 decades増加。

### balanced

\[
37.86\to50.62\to76.88\to103.62.
\]

約65.8 decades増加。

3つとも単調に増え、plateauを示さない。

従って今回のdeclared positive formsでは
**regulator removal = FAIL**。

## 12. massを先に作用させてもthreshold tailは消えない

既存tail auditでは、低energy thresholdとのnonzero overlapが
固定Tでのpolynomial tailを生成することが示唆されていた。

今回 \(M_0\psi\) 自身のthreshold overlapを測ると:

| right | dx | |B_mass| |
|---:|---:|---:|
| 16 | .04 | 0.0208 |
| 24 | .04 | 0.3724 |
| 24 | .02 | 0.3730 |

right=24ではgrid refinementで約0.373に安定し、
zeroではないことが数値的にresolveされた。

従って少なくとも

\[
r^{-3}M_D\psi
\]

型では「mass operatorを先に作用させればthreshold componentが消えて
exponential radius weightが救われる」という挙動は見えていない。

## 13. domain判定

今回の3つについて:

### \(48\|r^{-3}M_D\psi\|^2\)

mass-applied stateにnonzero threshold overlapが残るため、
既存のpolynomial-tail解析を適用する限り \(r^{-3}M_D\psi\) はL2 domainから外れる。

### \(48\|M_Dr^{-3}\psi\|^2\)

そもそも \(r^{-3}\psi\) がmass operatorへ入る前に
既存のpositive exponential moment domain問題へ入る。

### balanced form

\(r^{-3/2}\psi\) がintermediate form domainへ入らない同種の問題を持つ。

従って**この3 orderingはFAIL**。

ただし、explicit \(r^{-3}\) multiplicationとは異なる
genuinely nonlocal relational curvature observableがtail cancellationを持つ可能性までは
除外していない。

「全ての量子Kretschmannが必ず発散する」というno-go theoremではない。

## 14. machine-readable gate

詳細結果:

- `research/ks_mass_dirac_results/summary.json`
- `verification/ks_mass_dirac_verified_run.json`
- `audit/ks_black_hole_obligations.json`

今回のgate分解:

- classical/same-state clock slicing: PASS
- independently deparametrized Y-clock \(p_X^{-1}\) domain: PENDING
- finite-box L2↔KG same-state map: PASS
- continuum KG completion: PENDING
- factor ordering family: FAIL / ORDERING-SENSITIVE
- mass Dirac candidate: PARTIAL
- selected-mass positive Kretschmann forms: FAIL
- curvature regulator: FAIL
- common-state semiclassical: PARTIAL

claim labelsは変わらない:

- `MODEL-INTERNAL NUMERICAL RESULT`
- `ORDERING-SENSITIVE`
- `REGULATOR-UNSTABLE`
- `PHYSICAL INTERPRETATION NOT IDENTIFIED`

## 15. 次に残ったもの

順序として次は:

1. \(M_0\) positive quadratic formのcontinuum closability / associated self-adjoint operatorを解析する。
2. Y-clock \(p_X^{-1}\) のzero-mode domainを、今回のsame-Dirac-state sectorに対して定義する。
3. explicit-radius orderingではないrelational curvature observableの候補が物理的に動機づけられるか検討する。
4. それでもdomain-safe \(K\) が得られなければ、このminisuperspace/quantizationでのcurvature singularity persistenceを「scoped negative result」として評価する。

有限boxで巨大なKを得たことだけを、自然界のblack-hole singularityの証明とは扱わない。


## 16. 2026-09-20 continuum threshold / KG completion follow-up

[KS continuum threshold/domain audit](KS_CONTINUUM_DOMAIN_AUDIT_ja.md) で、positive-frequency Hamiltonianのcontinuum Liouville spectrumを明示した。

- normalized continuum mode: `phi_k=pi^(-1)sqrt(2k sinh(pi k))K_(ik)(e^(-x)/h)`, `E=hk`
- threshold: `phi_k/k -> sqrt(2/pi)K0(e^(-x)/h)`
- h=.2 reference Gaussian: `|B0|=0.0143884187466872`, domain/grid control relative spread約2.05e-15
- therefore `c(k)=B0 k+O(k^3)` and the packet lies in ordinary-L2 `D(H^-1/2)` and `D(H^-1)`; generic nonzero-B0 states reach the logarithmic threshold obstruction at `H^-3/2`
- positive-frequency KG spaceをenergy norm completionとして定義すると、`S=sqrt(h/2)H^-1/2` はL2からその完成空間へのunitary mapへ延長できる

従って **continuum positive-frequency KG completionはPASS** へ更新する。これはordinary-L2で `H^-1/2` がboundedになったという意味ではない。

mass candidateについては、`HV+VH=2P` と `|P|<=H` から `-I<=V<=I`, `I+V>=0` をcontinuumで得た。しかしfinite-boxの min eig(I+V)はbox拡大で0へ向かうため、coercive gapはない。したがって `q_M` のclosability / associated self-adjoint operatorは引き続き **PENDING**。

「次に残ったもの」の優先順位は変えず、まずmass form closabilityを解析し、その次にY-clock `p_X^-1` のstate-specific continuum domainへ進む。
