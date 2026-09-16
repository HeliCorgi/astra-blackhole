# astra-blackhole

曲率の短期予測と、平均化によって捨てられる物理状態の区別を調べる、再現可能な計算記録です。

**従来の物質研究は古典的モデル内の合成データによる試験です。追加したWDW研究も、時計と内積を限定した量子幾何学模型です。実在するブラックホール内部の観測、特異点の解消、量子重力の完成、新しい自然法則の発見を主張しません。**

## 最新：波束の裾と曲率関連量の期待値

[導出と4次元曲率との区別](docs/WDW_TAIL_AUDIT_ja.md) / [数値要約](research/wdw_tail_results/summary.json)

既存WDWの5波束を37格子・状態組で再計算しました。半径の幾何平均が安定しても、逆半径2/4/6乗の期待値には有限値への収束が見られません。低エネルギー展開は約x^-6の確率密度の裾を予測し、数値解と一致しました。初期丸め誤差、計算領域、観測量の上限を別々に検査しています。

**球面固有曲率と、古典的な曲率／質量二乗を先に規格化した量の監査です。4次元Kretschmann演算子や自然界の特異点の証明ではありません。** 漸近導出の条件と未証明の残差評価を明示し、有限箱の巨大値を収束した曲率として扱いません。

```sh
python research/run_wdw_tail_audit.py --out artifacts/wdw-tails --check-against research/wdw_tail_results/summary.json
python research/plot_wdw_tails.py --results artifacts/wdw-tails
```

## 追加：真空幾何学のWheeler–DeWitt量子化

[模型・時計・結果・限界](docs/WDW_AUDIT_ja.md) / [数値結果](research/wdw_results/summary.json)

既存の古典物質コードと別に、Kantowski–Sachs真空幾何学の一つの正周波数枝を実装しました。波束・古典点軌道・同じ初期Wigner分布の古典集団を比較しています。5状態の有限区間では半径の幾何平均は減少を続け、反転を示しませんでした。波束のx方向の反射は面積半径の反転ではありません。T=log(a)は固有時間でなく、量子曲率・特異点終端・時計や順序に依存しない結論は未検証です。

```sh
python research/run_wdw_audit.py --out artifacts/wdw --check-against research/wdw_results/summary.json
python research/plot_wdw.py --results artifacts/wdw
```

新規の独立テストは15件。古典模型とLeanの既存ソースは変更していません。CIの成功状態は実行記録で確認してください。

## 前回：同じ変数で正の角分布から補う比較

[導出・失敗・限界](docs/POSITIVE_CLOSURE_AUDIT_ja.md) / [全450行はCI成果物](../../actions) / [要約](research/positive_closure_results/summary.json)

前回と同じ90条件・同じ時間区間[0,.30]で、ρ・ΠのB、ρ・Π・c4のCを保ち、上位係数をゼロにする方法と、角エントロピーに基づく正の指数分布から補う方法を比較しました。総合合格はBで60→87/90、Cで82→90/90。変数を追加せずに旧Cの不合格8条件を解消しました。一方、K誤差はBの24条件、Cの31条件で増えています。正値性と精度は別の判定です。

**前回の条件を再使用した比較であり、新しい未使用データでの検証ではありません。** これは角度別エネルギーのエントロピーを使う近似で、全粒子分布の熱力学エントロピーを最大化したものでも、特異点の解消でもありません。独立に進める変数は増やしませんが、各時点で代数方程式を解くコストは増えます。

```sh
python research/run_positive_closure_audit.py --out artifacts/positive
python research/plot_positive_closure.py --results artifacts/positive
```

前回のローカル実行は [配布時の記録](verification/positive_closure_publication_checks.json)、今回の公開再試行は [再試行記録](verification/positive_closure_retry.json) を参照してください。CI設定の存在とGitHub上の実行成功を混同しないでください。

## 前回：最小状態変数と近似の適用範囲

[三者比較の導出・結果・限界](docs/CLOSURE_AUDIT_ja.md) / [結果要約](research/closure_results/summary.json) / [全270件の詳細はCI成果物](../../actions)

完全流体A、密度と圧力差を進めるB、高次方向成分も進めるCを、90条件・同じ固有時間幅[0,.30]で比較しました。別幾何学60条件での最小合格候補はA=12、B=30、C=11、三候補とも不合格=7。合格には曲率・圧力差の精度だけでなく、保存則と分布の妥当性も要求しています。自然界での最小変数数を決めたものではありません。

**残したモーメントが正の分布と両立することと、上位成分をゼロにした再構成が正であることは別でした。** 全270低次軌道中38件は前者を満たしても後者に失敗しました。その後、同じ変数で正値性を保つ閉じ方を比較しました（上記）。

```sh
python research/run_closure_audit.py --out artifacts/closure --check-against research/closure_results
python research/plot_closure.py --results artifacts/closure
```

この第四の監査は専用の `Minimal-state closure audit` CIで実行します。既存の `reproduce.py` は三つの従来監査と全単体テストを実行し、新しい全90条件の監査は上の別コマンドで実行します。新しいLean命題は追加していません。

## 前回：衝突と反応時間

[導出・対照・失敗・限界](docs/COLLISION_AUDIT_ja.md) / [全設定の結果要約](research/collision_results/results.json)

質量ゼロ・運動量に依存しない緩和時間近似を追加しました。粒子数とエネルギーを保ち、方向分布と時空を連立します。強い衝突で初期方向構造の影響は縮みますが、有限の反応時間による完全流体近似との差は別に残ります。二体散乱積分や量子相関の計算ではありません。

P4・t=.45の例：ν=0で状態対の曲率差2.173775%、ν=64で0.00785146%。同じν=64でも完全流体参照からの最大相対差は3.024821%（分母K_fluid）。前者は対称相対差で、両者とも学習予測器の誤差ではありません。粗い角度展開で負の分布が出た失敗と、数値的に未分解の小さな差も保存しています。

## 入口

- [今回：方向分布の省略、導出・結果・限界](docs/ANGULAR_AUDIT_ja.md)
- [今回の全ケース・数値検査・厳密式](research/angular_results/results.json)
- [前回：径方向分布の研究の導出・結果・限界](research/REPORT_ja.md)
- [曲率予測器の対応表](docs/PREDICTOR_MAP_ja.md)
- [前回の主結果・全係数](research/results/results.json)
- [前回の追加24対](research/results/cohort.json)
- [旧予測器の24窓での採点](research/results/legacy_predictor_on_new_model.json)
- [初回公開の検証範囲](verification/publication_checks.json) / [今回](verification/angular_publication_checks.json)
- [研究段階と次の課題](docs/ROADMAP_ja.md)
- [GitHub Actions](../../actions)
- [Lean環境・証明の範囲](docs/LEAN_ja.md)

## 最短の実行

Python 3.10以降を使います。仮想環境の利用を推奨します。

```sh
python -m pip install -r requirements.txt
python predict_curvature.py example_input.json --method polynomial --horizon 0.1
python predict_curvature.py example_input.json --method adaptive
python reproduce.py --out artifacts
```

`reproduce.py` は全単体テスト（従来48件とWDW用15件）、計量からの独立した幾何学検算、径方向・角方向の二つのEinstein–Vlasov監査と、今回の衝突RTA監査を実行します。過去の全研究の再学習・再実行ではありません。`--out` を指定すると今回の主結果を保存済みベースラインとも比較します。幾何学検算だけは従来どおり `research/results/symbolic_geometry.json` に書き込みます。

グラフと軌道の配列は再生成します。

```sh
python research/run_research.py
python research/plot_results.py
python research/run_angular_audit.py
python research/plot_angular.py
```

ヘッドレス環境では `MPLBACKEND=Agg` を指定してください。再実行するとJSONの実行時間・環境欄なども更新されます。独立した出力先には各監査の `--out /path/to/output` を使えます。方向監査の曲線配列・PNGは生成物で、コミットにはソースと数値結果JSONを収録しています。

## 今回：方向分布は平均圧力に見えないことがある

粒子数・全径方向スペクトル・初期密度・圧力・幾何学を一致させ、方向成分 `P4/P6/P8` だけを変えました。**初期圧力が等方的でも、粒子分布が等方的とは限りません。**

|方向分布の差|最初に異なる曲率の時間微分|質量1：t=.45の曲率差|質量0：t=.45の曲率差|
|---|---:|---:|---:|
|P4|1階|0.741751%|2.173775%|
|P6|2階|0.142426%|0.461370%|
|P8|3階|0.042004%|0.146659%|

**差は対称相対差であり、予測誤差ではありません。各行は別の状態対です。** 質量ゼロでも方向情報の影響は残りました。前回の「径方向スペクトルだけなら差が消える」対照は、方向分布固定という条件のもとで維持されています。詳しくは[今回の導出](docs/ANGULAR_AUDIT_ja.md)を参照してください。

## 前回：運動量の大きさの分布

同じ初期密度・圧力・幾何学を持つが、高次の運動量分布が異なる正の粒子分布を構成しました。圧力の時間変化に最初の省略が入り、曲率の時間微分へ伝わります。

| 初期に一致させた情報 | 最初に異なる曲率の時間微分 | t=0.45の曲率の対称相対差 |
|---|---|---:|
| 粒子数・密度・圧力・幾何学 | 1階 | 0.325256% |
| 上記と4次モーメント | 2階 | 0.095363% |
| 上記と6次モーメント | 3階 | 0.030313% |

**各行は別の状態対です。最後の列は予測誤差ではありません。** 指定した低次の初期情報では当て分けられない反例であり、精密な過去履歴を含む全観測が同じという主張ではありません。

前の無衝突の両監査は一様なKantowski–Sachs Einstein–Vlasov系です。量子相関・衝突・蒸発・回転・外部接続は含みません。GitHub Actionsは再現性を検査するもので、自然界での正しさの保証ではありません。Leanの小規模な証明ソースと独立CIを追加しました。対象は有限有理数恒等式と観測の区別に関する補題であり、重力モデル全体の形式証明ではありません。実際のビルド状況はActionsで確認してください。

## 過去の保存版とニューラル予測器

Git初回公開は、最新研究を再実行できるソース・数値結果と、解析的／多項式／適応幅の予測器を収録しました。**会話で配布したZIP全体のミラーではありません。過去の学習済み重み・全履歴・生成画像は含めていません。**

元の `curvature_predictors_and_state_audit.zip` を会話から保存したうえで、次の操作により歴史的な `legacy/` 一式をローカルへ復元できます。外部から自動ダウンロードはしません。

```sh
python tools/restore_legacy.py /path/to/curvature_predictors_and_state_audit.zip
python -m pip install -r requirements-neural.txt
python predict_curvature.py example_input.json --method neural
python verify_predictors.py
```

復元ツールは元ZIPのSHA-256を照合し、変更済みファイルを上書きしません。ニューラル版は元のLTBモデル単位で過去[-0.30,0]の41点を要求し、未来[0.10,0.20,0.30,0.40]を出力します。現在の対数曲率の微分も必要です。別のモデルや実測への有効性は未検証です。

詳しい来歴と初回公開時の変更点は [PROVENANCE.json](PROVENANCE.json) に記録しています。初回公開時のMIT表記は歴史的な記録です。現行ライセンスは次のとおりです。

## ライセンス

**Apache License 2.0 — Copyright 2026 HeliCorgi**

[LICENSE](LICENSE) / [NOTICE](NOTICE) / [変更記録](docs/LICENSE_HISTORY_ja.md)

## Leanによる小規模な証明検査

```sh
python3 tools/check_lean_certificates.py
cd lean
lake build
lake env lean -DwarningAsError=true Audit.lean
```

Lean 4.19.0に固定し、Mathlibはまだ使いません。Pythonの有理数検算だけではLeanのビルド成功を意味しません。証明対象と未形式化部分は[Leanノート](docs/LEAN_ja.md)を参照してください。

### Git公開時の保存形式（2026-09-16）

計算ソース・プロトコル・検証記録と `research/positive_closure_results/summary.json` をGitに保存する。全450行、90基準解の診断、事後診断と全時系列は再生成してGitHub ActionsのZIP成果物に保存する。元の配布ZIPの内容は変更せず保管しており、Gitへの公開はその全ファイルのミラーではない。

```sh
python tools/verify_positive_closure_artifact.py artifacts/positive
python tools/check_positive_summary.py artifacts/positive/summary.json research/positive_closure_results/summary.json --self-test
```

後者は集計値の回帰検査で、450行を個別に照合したことを意味しない。公開時の詳細照合とCIの成功は別の検証記録で報告する。新しい物理計算や未使用条件での検証を追加した公開ではない。
