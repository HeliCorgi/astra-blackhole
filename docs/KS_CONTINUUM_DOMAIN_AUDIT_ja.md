# KS continuum threshold / physical-KG completion / mass-form domain 監査

更新: 2026-09-20

## 結論

Schwarzschild/Kantowski–Sachs の positive-frequency sector

[
H^2=-h^2partial_x^2+e^{-2x},qquad Hge0
]

について、finite box だけに依存していた zero-energy threshold と物理内積の判定を continuum で切り分けた。

| 対象 | 判定 |
|---|---|
| exact Liouville continuum spectrum / threshold normalization | PASS |
| positive-frequency physical-KG Hilbert completion | PASS |
| reference Gaussian の ordinary-L2 (H^{-1/2}) domain | PASS |
| continuum velocity (-Ile Vle I), (I+Vge0) | PASS |
| mass quadratic form closability / associated self-adjoint operator | PENDING |
| standalone Y-time (p_X^{-1}) continuum domain | PENDING |

重要なのは、**(H^{-1/2}) が ordinary (L^2) 上の有界演算子でないことと、positive-frequency KG Hilbert space が定義できないことは同義ではない**、という点である。KG側は energy norm の完成として定義でき、repo の same-state map はその完成へ unitary に延長できる。

一方、mass candidate

[
q_M[psi]=rac14langle H e^{X/2}psi,(I+V)H e^{X/2}psiangle
]

では、中央因子 (I+V) のcontinuum positivityまでは示せるが、full weighted/nonlocal formの closability はまだ示していない。有限箱での正の最小固有値を continuum coercivity と読み替えない。

## 1. 模型・範囲

- 真空 Kantowski–Sachs minisuperspace
- dimensionless repo conventions
- retained: (x) と positive-frequency (T) branch
- omitted: matter, inhomogeneous modes, environment
- 実在ブラックホール観測、特異点解消／持続、完全量子重力の判定ではない

## 2. exact Liouville spectrum

(q=-x-log h) と置くと

[
rac{H^2}{h^2}=-partial_q^2+e^{2q}.
]

Liouville quantum mechanics の完全規格化連続固有関数を repo 変数へ戻すと

[
phi_k(x)=rac1pisqrt{2ksinh(pi k)},
K_{ik}(e^{-x}/h),qquad k>0,
]

[
Hphi_k=hk,phi_k.
]

Anderson–Nilsson–Pope–Stelle, arXiv:hep-th/9401007 の Eq. (2.2.9) の規格化を使う。spectrum は ([0,infty)) の連続スペクトルで、normalizable zero mode はない。

(k	o0) では

[
sqrt{2ksinh(pi k)}/pi
=sqrt{2/pi},k+O(k^3),
]

かつ (K_{ik}(z)=K_0(z)+O(k^2)) なので

[
oxed{rac{phi_k(x)}k	osqrt{rac2pi}K_0(e^{-x}/h).}
]

## 3. reference packet の threshold coefficient

現行 mass/Y-clock audit と同じ

[
h=.2,quad x_0=2,quad p_0=-1,quad sigma=sqrt h
]

の Gaussian に対し

[
B_0=sqrt{rac2pi}int K_0(e^{-x}/h)psi(x),dx
]

を直接積分した。

| left | right | dx | (|B_0|) |
|---:|---:|---:|---:|
| -4 | 8 | .01 | 0.0143884187466872 |
| -6 | 10 | .005 | 0.0143884187466872 |
| -8 | 12 | .0025 | 0.0143884187466873 |

相対spreadは約 (2.05	imes10^{-15})。従ってこの状態では

[
c(k)=B_0 k+O(k^3),qquad B_0
e0.
]

よって (H^{-s}psi) のzero-threshold寄与は

[
int_0^epsilon k^{2-2s}dk
]

で決まり、genericな (B_0
e0) 状態では

[
oxed{s<3/2}
]

が局所可積分条件になる。reference packet は ordinary (L^2) の (D(H^{-1/2})) と (D(H^{-1})) に入るが、(H^{-3/2}) で logarithmic threshold obstruction に達する。

## 4. positive-frequency KG Hilbert space

positive-frequency KG inner productは

[
(Psi,Phi)_{KG}=rac2hlanglePsi,HPhiangle_{L^2}.
]

ordinary (L^2) の (D(H^{1/2})) を

[
|Psi|_{KG}^2=rac2h|H^{1/2}Psi|_2^2
]

で完成した空間を (mathcal H^+_{KG}) とする。

repoで使っている

[
S=sqrt{h/2},H^{-1/2}
]

は ordinary (L^2) 上では unbounded である。しかし spectral representation では

[
|Schi|_{KG}=|chi|_2
]

なので、定義域上のisometryから完成によって

[
oxed{S:L^2	omathcal H^+_{KG}}
]

へ unitary に延長できる。

したがって、以前の「continuumで (H^{-1/2}) がunboundedなのでphysical KG completionはPENDING」という理由は強すぎた。**continuum positive-frequency KG completion自体はPASSへ更新する。**

注意: 完成後のすべてのKG amplitudeが ordinary (L^2(dx)) 関数であるとは限らない。ordinary-L2 observableをそのまま読むnegative controlは引き続き不採用。

## 5. continuum velocity の positivity

[
P=-ihpartial_x,qquad H^2=P^2+e^{-2X},
]

[
V=rac{i}{h}[H,X].
]

共通core上で

[
HV+VH=rac{i}{h}[H^2,X]=2P.
]

従って Sylvester/semigroup 表現

[
V=2int_0^infty e^{-tH}P e^{-tH},dt
]

を使える。(H^2ge P^2) と square-root のoperator monotonicityから (|P|le H)。また (ker H={0})。よって

[
|langlepsi,Vpsiangle|
le2int_0^infty
langle e^{-tH}psi,H e^{-tH}psiangle dt
=|psi|^2.
]

したがって bounded self-adjoint extension について

[
oxed{-Ile Vle I,qquad I+Vge0.}
]

これはfinite boxの固有値が正だったことより強いcontinuum statementである。

## 6. finite-box negative control: gapは閉じる

repoと同じpositive fourth-order Dirichlet discretizationで、別に

[
P_{comm}=rac{i}{2h}[H^2,X]
]

を作り、(HV+VH=2P_{comm}) を検査した。

| right | Sylvester relative residual | min eig((I+V)) |
|---:|---:|---:|
| 8 | 5.67e-14 | 2.786e-4 |
| 12 | 8.99e-14 | 4.685e-5 |
| 16 | 1.29e-13 | 1.365e-5 |

(V) は数値的にも ([-1,1]) に入り、同時に min eig((I+V)) はbox拡大で0へ向かう。これはasymptotic left-moving/null channelと整合する。

従って finite box の strict positivity gap を使ってmass form closureを証明する方針は採用しない。

## 7. mass quadratic form: まだPENDING

現在のcandidate:

[
q_M[psi]=rac14
langle H e^{X/2}psi,(I+V)H e^{X/2}psiangle.
]

今回確定したのは:

- (I+V) はcontinuumで bounded positive
- (C_c^infty(mathbb R)) はcandidate coreとして有限値を与える
- finite-box positivityはcontinuum gapの証拠ではない

未確定:

- full formのclosability
- closed form completion
- representation theoremで得られるassociated self-adjoint operator
- ordering uniqueness

特に (I+V) がnull channelで退化するため、単純なcoercive lower boundから (H e^{X/2}) のgraph normへ還元できない。ここを次の主課題にする。

## 8. Y-clockへの含意

standalone Y-time Hamiltonian

[
K_Y=rac{p_X^2+e^{-2X}}{2p_X}
]

が (p_X^{-1}) を含む問題も、同じasymptotic null channelに関係する。operatorとして inverse はunboundedになり得る一方、特定stateがdomainに入る可能性は別問題。

現行finite-box packetでは (|p_X^{-1}psi|) がright=12→32で約174→175に留まるため、次は**same Dirac stateのcontinuum spectral densityに対して直接domain条件を出す**。operator normの発散だけでstateを棄却しない。

## 9. 文献との位置づけ

- Anderson et al., arXiv:hep-th/9401007: Liouville Hamiltonianのfully-normalized (K_{ik}) eigenfunction。今回のKS変数へ座標変換してthreshold解析に使用。
- Mostafazadeh, arXiv:math-ph/0209014 / gr-qc/0205049: Klein–Gordon型方程式のsolution spaceにpositive-definite Hilbert structureを構成する一般的背景。今回の具体的energy completionはrepoのpositive-frequency sectorから直接定義したもの。
- Kuchař, arXiv:gr-qc/9403003: Schwarzschild massがembedding-independent canonical variableとして残る古典/Dirac量子化の背景。今回の (M_0) orderingをその論文から導出したわけではない。

## 10. 再現

```sh
python -m unittest discover -s tests -p 'test_ks_continuum_domain.py' -v
python research/ks_continuum_domain_audit.py \
  --out artifacts/ks-continuum/summary.json
python tools/check_ks_continuum_domain_results.py \
  --artifact artifacts/ks-continuum/summary.json \
  --baseline research/ks_continuum_domain_results/summary.json
```

この監査はLean sourceを変更しない。CI successは解析の入力・数値対照の再現記録であり、mass form closabilityや特異点の物理解釈を証明しない。
