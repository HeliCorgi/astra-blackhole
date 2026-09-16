# Lean環境と証明の範囲

2026-09-16。ツールチェーンは `leanprover/lean4:v4.19.0` に固定。これは最新版という意味ではない。Mathlib等の外部Lakeパッケージを追加せず、まず標準ライブラリだけで小さな検証対象を作る。

## 実行

Lean公式の手順でelanを導入した後、リポジトリ内で実行する。

```sh
python3 tools/check_lean_certificates.py
cd lean
lake build
lake env lean -DwarningAsError=true Audit.lean
```

`lake build` の標準ターゲットは `AstraBlackhole`。`lean-toolchain` からバージョンが選ばれる。`Audit.lean` は8定理の依存公理を出力する。CIはその出力が全件、公理依存なしであることを確認する。証明穴、独自の公理、`native_decide` は使用しない。

## 追加した証明ソース

### 有限有理数計算：6定理

`AngularCertificates.lean` は `Std.Internal.Rat` による有理数計算を使う。偶数冪の係数列 c に対して

```
M(c,s) = sum_k c[k] / (2*(s+k)+1)
```

という**有限演算**を定義し、`P4/P6/P8` の係数列について次を `by decide` で検証する。

|係数列|ゼロになるs|次の値|
|---|---|---|
|P4|0,1|M(P4,2)=8/315|
|P6|0,1,2|M(P6,3)=16/3003|
|P8|0,1,2,3|M(P8,4)=128/109395|

この有限演算は通常の計算では `(1/2) integral[-1,1] x^(2s) P_l(x) dx` に対応する。ただし、その**実積分との対応や、一般のLegendre直交定理はまだLeanで証明していない**。曲率の時間微分、重力方程式、粒子分布の正値性もこのファイルでは証明しない。

### 観測の区別：2定理

`Observability.lean` は、観測が同じなら決定論的予測器の出力も同じになること、真の未来が異なるという仮定のもとでは両方へ完全一致できないことを証明する。このような状態が自然界に存在すると証明するものではない。

## CI

`.github/workflows/lean.yml` を古典数値計算のワークフローと分離した。関連パスのpush・PR・手動起動で動く。権限は `contents: read`。自動コミット・定期研究は行わない。直接指定するActionsはコミットSHAに固定し、Mathlib取得とLakeキャッシュは無効化した。

PythonのFractionによる独立検算、Leanビルド、依存公理検査を別々に扱う。GitHub上の実行結果が確認できるまでは、設定の追加を証明成功と言い換えない。ローカルでLeanを実行できない場合も明記する。

## 参照

- Lean公式導入手順: https://lean-lang.org/install/manual/
- Lake公式文書: https://lean-lang.org/doc/reference/latest/Build-Tools-and-Distribution/Lake/
- 固定ツールチェーン: https://github.com/leanprover/lean4/releases/tag/v4.19.0
- 固定したLean Action: https://github.com/leanprover/lean-action/tree/50fcf42d2e460296f1a34b402e990d1b24f8b596

数値計算の再現性、有限数学の証明、自然界での妥当性は別々の主張である。
