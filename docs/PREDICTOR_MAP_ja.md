# 曲率予測器の対応表

ここでいう予測器は、既知の方程式から作った過去の曲率を入力し、未来の曲率を予測する計算です。自然界の特異点の最終状態を決定する理論ではありません。

| CLI | 今回の入口 | 元の実装・位置づけ |
|---|---|---|
| `--method power` | `predictors.power_law` | Kが残り時間の逆4乗に比例すると仮定し、Kと増加率から外挿 |
| `--method polynomial` | `predictors.polynomial` | `legacy/collapse_extensions/forecast.py::forecast`。元の関数を変更せず保存 |
| `--method adaptive` | `adaptive_forecast.py::predict` | 元のLTB予測器から `fit_polynomials`、`ratio`、`predict` の本体を変更せず抽出 |
| `--method neural` | `predictors.neural` | 元の1次元潜在変数モデルのノイズ対策済み重み。元ZIPから復元後に使用 |
| `--method neural-clean` | `predictors.neural` | ノイズ対策前の重み。元ZIPから復元後に使用 |

多項式版は q=(K/K現在)^(-1/4) を、現在を基準とした一次〜三次の多項式で延長します。未来を学習に使わず、過去の窓ごとに係数を当てはめます。適応版は次数間の食い違いが大きいと先読み幅を半減しますが、食い違いの閾値は真の誤差の上限ではありません。

元の適応版の完全なシミュレーターは、保存ZIP内の `legacy/inhomogeneous_singularity_forecast/forecast.py` にあります。今回のGit版では予測関数だけを独立させ、過去の物理生成器を動かさず利用できます。

ニューラル版のコード・正規化係数・重みは `tools/restore_legacy.py` で復元します。重みの読み込みは `weights_only=True`。このオプションやハッシュ確認は、物理的妥当性を保証するものではありません。

最新の研究が使用するのは固定先読み幅の旧多項式版です。24個の重なった予測窓で採点しています。旧ニューラル版を新しいEinstein–Vlasov系へ無条件に転用したわけではありません。
