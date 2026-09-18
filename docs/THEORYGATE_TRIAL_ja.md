# TheoryGate v0.4 実戦trial：Bianchi IX ever-entered history probability

## 目的

TheoryGateを単なるサンプルではなく、`astra-blackhole` の既存evidenceへ実際に接続する。

作者推奨の

```sh
theorygate template init HISTORY_PROBABILITY \
  --model-id astra-bianchi-ix \
  --output theorygate.yaml
```

相当のtemplate生成から開始し、Lean / SymPy / robustness evidenceを実際に生成して、
`HISTORY_PROBABILITY` が昇格するかを検査した。

TheoryGateは commit
`cb33044551f2f98a03ab813e20c4caf8fe9f9743`
(v0.4.0) に固定した。

## モデル固有のrewrite

built-in `HISTORY_PROBABILITY` templateはstarter policyなので、Bianchi IX finite-clock
ever-entered問題へ合わせて義務を具体化した。

テンプレート由来の中心義務:

- finite-clock CAP class operator definition
- no-B / ever-B history partition
- physical inner product
- decoherence
- regulator stability
- boundary robustness

に加えて、canonical gravityで無視できない

- clock robustness
- factor-ordering robustness

をblocking obligationへ追加した。

この追加はTheoryGateのtemplateを「間違い」としたものではない。
TheoryGate自身のdocumentationが、templateはmodel-specificに書き換えるscaffoldingだと
明記しているため、その規則に従った。

## 実行したevidence adapters

### Lean

TheoryGate Lean adapterを実際に使い、

- `AstraBlackhole.time_dependent_factorization_residual`
- `AstraBlackhole.common_operator_preserves_recombination`

を `--no-axioms --require-pass` で監査した。

結果:

`LEAN_HISTORY_ALGEBRA = PASS`

これは抽象意味論のPASSであって、時計・量子化・物理内積の正当化ではない。

### SymPy

TheoryGate SymPy adapterで2つの恒等式を検査した。

1. Bianchi IX potentialから (W=2e^{-4s}V) の6指数wall項への展開
2. (alpha=-s, p_s=-p_alpha) 後のconstraint rescaling identity

結果:

`REDUCED_DYNAMICS_ALGEBRA = PASS`

これはCAS algebraのPASSであり、GR reductionの独立再導出そのものではない。
GR reductionは別途Cadabra+Maxima gateで監査済み。

### Regulator robustness

既存のfinite-clock CAP scan

[
V_0=.025,.040,.050,.060,.075,.085,.100
]

の (|D_{m off}|) をTheoryGate robustness adapterへ投入した。

数値閾値は後付けしていない。
既存監査が採用した「狭い最小値を選ばず、limit方向でdriftが再増大するなら不採用」
というdecisionを
`require_nonincreasing_successive_drift=true`
で表した。

結果:

`REGULATOR_STABILITY = FAIL`

reason:

`successive relative drift increases toward the requested limit`

既存の「V0=.075だけ選んで確率化しない」という判断と一致した。

### Boundary robustness

V0=.075, maxdim=78で smoothing width=.20/.30/.45 のstrict offdiagonalを比較した。

TheoryGateは最大symmetric relative difference 約0.05079を記録したが、
事前登録した許容閾値がないため:

`BOUNDARY_ROBUSTNESS = PARTIAL`

となった。

これは適切な挙動。5%程度だからPASS、とTheoryGateが勝手に決めなかった。

## 最終claim

TheoryGateの出力:

- `FINITE_CLOCK_CAP_CLASS_OPERATOR_DEFINED = PASS`
- `HISTORY_PARTITION_DEFINED = PASS`
- `LEAN_HISTORY_ALGEBRA = PASS`
- `REDUCED_DYNAMICS_ALGEBRA = PASS`
- `PHYSICAL_INNER_PRODUCT = OPEN`
- `DECOHERENCE = FAIL`
- `REGULATOR_STABILITY = FAIL`
- `BOUNDARY_ROBUSTNESS = PARTIAL`
- `CLOCK_ROBUSTNESS = OPEN`
- `ORDERING_ROBUSTNESS = OPEN`

Claims:

- **MODEL_INTERNAL_EVER_ENTERED_WEIGHT: SUPPORTED**
- **HISTORY_PROBABILITY: BLOCKED**

strongest supported claim:

`MODEL_INTERNAL_EVER_ENTERED_WEIGHT`

これは現在のrepoの物理判断と一致する。

## TheoryGate側で実戦上見えた点

### 良かった点

1. FAIL/PARTIALが「テスト失敗で消える」のではなくevidenceとして残る。
2. `--require HISTORY_PROBABILITY` は期待どおりnonzeroとなり、上位claimをCIで止められる。
3. 一方で `--require MODEL_INTERNAL_EVER_ENTERED_WEIGHT` はPASSする。
4. Lean adapterのclean provenance / axiom policyを既存Lean監査へ自然に接続できた。
5. robustness adapterが「閾値なしならPARTIAL」を守った。
6. regulatorの非単調性を、数値閾値を捏造せずFAILへできた。

### 不足

現行v0.4のCAS adapterは

- SymPy
- xAct
- Cadabra

に対応するが、`astra-blackhole` の第二独立GR backendである
**Maxima/ctensor** を直接evidence化できない。

手動evidence ingestionは可能だが、

- exact command
- Maxima version
- script hash
- exit code
- PASS marker

をTheoryGate標準形式で自動記録できない。

実戦上は dedicated Maxima adapter、または任意external-CAS用のgeneric script adapterがあるとよい。

## 検証記録

workflow:

`TheoryGate Bianchi IX history-probability trial`

run:

`35389840066`

job:

`105745359311`

artifact:

`10566075329`

artifact SHA-256:

`73f66296342c4ec21a2821c565243d23b48bb05a7c16e984241cb9e22732e97a`

workflow conclusion: success.

この成功はTheoryGateが既存evidenceを正しくscope/gate化したことを示す。
HISTORY_PROBABILITY自体はBLOCKEDであり、新しい物理的確率を確立した結果ではない。
