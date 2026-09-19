# PHYSICS_AUDIT — 物理主張の機械監査

更新: 2026-09-19

この文書は「どの量子化が自然界で正しいか」を自動決定するものではない。目的は、専門家が通常行う監査を **obligation / gate / kill 条件**へ分解し、数理証明・CAS・数値計算・物理解釈を混同しないようにすること。

対象はまず Bianchi IX timeless class-operator 系列。将来は Kantowski–Sachs、WDW曲率、古典閉包にも同じ schema を適用する。

## 原則

Lean が証明するのは、採用した定義・仮定から結論が論理的に従うこと。Lean build が通っても、量子化、時計、内積、operator domain、factor ordering、regulator、自然界への同定が正しいことは意味しない。

有限行列が Hermitian であることと、連続作用素の self-adjointness は分ける。CI成功は再現性の証拠であり、物理的妥当性の証明ではない。

## Gate / kill ledger

| gate | 現在 | promotion rule | kill / downgrade rule |
|---|---|---|---|
| GR reduction / CAS | PASS | 2独立CASで縮約Hamiltonian・constraint・符号・係数を再導出 | どちらかのbackend不一致ならFAILへ戻す |
| Lean semantics | PASS | factorization residual / branch recombination / constraint-kernel preservationをLean build | build/axiom audit不合格なら ALGEBRAICALLY VERIFIED を出さない |
| clock | PENDING | {T,C}、単調性、branch一意性、複数clock比較 | 結論がclockで変われば CLOCK-DEPENDENT |
| inner product / domain | PENDING | 保存則、対称性、domain、extensionを明示 | 有限格子Hermiticityだけでは通さない |
| factor ordering | PENDING | principal symbol、symmetry、formal symmetry、semiclassical limit、ordering scan | 結論がorderingで変われば ORDERING-SENSITIVE |
| regulator | FAIL | box/cap/grid/eta/CAP/PML/Krylovを物理量から分離し安定域を確認 | 現在のtimeless scatteringは regulator-independent limit 未確認 |
| known limits | PENDING（一部PARTIAL） | WKB/HJ、Schwarzschild/KS、current、free case、soluble toy model | 一部だけでは SEMICLASSICAL LIMIT PASSED を出さない |
| claim compiler | IMPLEMENTED | obligations JSONのみからclaim labelを生成 | 人手で上位claim labelを上書きしない |

現在の regulator gate は、finite-window Dirichlet S-matrix、stationary resolvent、outer absorbing-layer control の三段階を経ても strength / layer / eta dependence が残るため FAIL。これは Halliwell formalism の否定ではなく、現在の数値散乱設計の不採用を意味する。

## Claim labels

claim compiler が扱うラベル:

- `ALGEBRAICALLY VERIFIED`
- `MODEL-INTERNAL NUMERICAL RESULT`
- `CLOCK-DEPENDENT`
- `ORDERING-SENSITIVE`
- `REGULATOR-UNSTABLE`
- `SEMICLASSICAL LIMIT PASSED`
- `PHYSICAL INTERPRETATION NOT IDENTIFIED`

`ALGEBRAICALLY VERIFIED` は Lean semantic gate の範囲だけを指す。プロジェクト全体や量子重力の物理的正しさを指さない。

`PHYSICAL INTERPRETATION NOT IDENTIFIED` は、blocking obligation に PENDING / PARTIAL / FAIL が一つでも残れば自動付与する。

## Lean semantic checker

Lean 4.19.0 のCI run 35381665961で `lake build` と `Audit.lean` が完了し、新規3定理はいずれも公理依存なしと確認した。これは以下の抽象定理だけに対する `ALGEBRAICALLY VERIFIED` である。


`lean/AstraBlackhole/PhysicsAudit.lean` は次の3点だけを証明する。

1. time-dependent factorization の residual が、branch equation・product rule・square relationから従う。
2. branch recombination が成立している状態へ **同じ** evolution operatorを作用させれば、recombination equality が保存される。
3. constraint と可換で zero state を保つ operator は constraint kernel を保つ。

これらは抽象的な等式定理で、微分作用素のdomain、self-adjointness、WDW量子化の選択を証明しない。

## CAS reduction checker

`cas/README.md` に外部CAS result manifestの契約を定義する。xAct/xTensor と Cadabra は候補backendだが、このrepoのclean-clone CI依存にはまだ入れない。

CAS gateをPASSにするには、少なくとも backend/version、入力action/metric/ADM convention、保持/削除自由度、reduced Lagrangian/Hamiltonian/constraint、boundary term、Poisson bracket、係数表、source/output hash、独立照合結果を保存する。

Cadabra 2.5.14 と Maxima/ctensor 5.46.0 の2独立backendを実行し、Bianchi IXの3-curvature、ADM kinetic term、Hamiltonian constraint、reduced generatorがrepo実装と一致した。したがって Bianchi IX の `GR_REDUCTION` gate はPASS。xAct sourceは保存済みだがWolfram/xAct runtime未設定のため、任意の第三cross-checkとして未実行のまま残す。詳細は `docs/BIANCHI_IX_GR_REDUCTION_CAS_ja.md`。

## Black-hole return gate

Bianchi IXはブラックホール内部そのものではない。BH singularity interpretationへ戻る経路を曖昧にしないため、`known_limits.schwarzschild_ks_bridge` をblocking obligationとして追加する。

Schwarzschild内部→Kantowski–Sachsの古典写像と既存KS WDW constraintへの正準変換はMaxima/ctensorとCadabraで再導出し、`known_limits.schwarzschild_ks_bridge` はPASSした。共有変数として areal radius `r=(1/4)exp(-x-T)` がSchwarzschildの面積半径そのものに一致し、classical mass Dirac observableも固定した。詳細は `docs/SCHWARZSCHILD_KS_BRIDGE_ja.md`。量子状態/内積・4次元曲率演算子・clock/order/domainは別gateとして残る。

## Numerical regulator checker

既存の保存結果とverification JSONをevidenceとして使う。現在は outer absorbing-layer control が数値再現済みだが、eta、outer strength、layer width への依存が残るため regulator gate は FAIL。

望ましい値へ regulator を合わせてPASSにはしない。

## Machine-readable obligations

`audit/physics_obligations.json` がsource of truth。`tools/check_physics_audit.py` が schema、evidence、gate状態を検査し、claim labelsを生成する。

```sh
python tools/check_physics_audit.py \
  --self-test \
  --obligations audit/physics_obligations.json \
  --claims-out artifacts/physics-audit/claims.json
```

Lean変更時はAGENTS.mdに従い、必ず実際に:

```sh
cd lean
lake build
lake env lean -DwarningAsError=true Audit.lean
```

を実行し、axiom reportを保存する。

## 現在未達

- clock Poisson bracket / monotonicity の機械監査
- timeless induced physical inner product
- continuum operator domain / self-adjoint extension監査
- factor-ordering family scan
- WKB/Hamilton–Jacobiまで含むknown-limit gate
- true PML / exterior complex scaling / continuum scattering limit

従って、現時点で timeless Bianchi IX class operator に物理履歴確率を割り当てない。


## KS black-hole quantum target

Schwarzschild/Kantowski–Sachs側はBianchi IX timeless targetと別のmachine-readable targetとして
`audit/ks_black_hole_obligations.json` に分離した。

現在:

- classical Schwarzschild/KS bridge: PASS
- classical T / areal-radius clock admissibility: PASSだが、quantum multi-clock domainを含むclock gateはPARTIAL
- finite-box positive-frequency L2 ↔ KG same-state map: PASSだがcontinuum (H^{-1/2}) domain未解決のためPARTIAL
- explicit Sturm–Liouville factor-ordering family: FAIL / sensitive
- quantum mass/Kretschmann ordering-domain: PARTIAL
- inverse-radius regulator removal: FAIL
- common-state semiclassical limit: PENDING

したがってKS targetのclaim compilerは `ORDERING-SENSITIVE`,
`REGULATOR-UNSTABLE`, `PHYSICAL INTERPRETATION NOT IDENTIFIED` を保持する。
詳細は `docs/KS_QUANTUM_AUDIT_ja.md`。
