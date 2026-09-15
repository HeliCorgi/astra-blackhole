# astra-blackhole

曲率の短期予測と、平均化によって捨てられる物理状態の区別を調べる、再現可能な計算記録です。

**現段階の結果は、既知の古典的モデル内の合成データによる試験です。実在するブラックホール内部の観測、特異点の解消、量子重力の完成、新しい自然法則の発見を主張しません。**

## 入口

- [研究の導出・結果・限界](research/REPORT_ja.md)
- [曲率予測器の対応表](docs/PREDICTOR_MAP_ja.md)
- [保存した主結果・全係数](research/results/results.json)
- [追加24対の結果](research/results/cohort.json)
- [旧予測器の24窓での採点](research/results/legacy_predictor_on_new_model.json)
- [公開時の検証範囲](verification/publication_checks.json)
- [次の未実施課題](docs/ROADMAP_ja.md)

## 最短の実行

Python 3.10以降を使います。仮想環境の利用を推奨します。

```sh
python -m pip install -r requirements.txt
python predict_curvature.py example_input.json --method polynomial --horizon 0.1
python predict_curvature.py example_input.json --method adaptive
python reproduce.py
```

`reproduce.py` は3件の単体テスト、計量からの独立した幾何学検算、最新のEinstein–Vlasov研究を実行します。過去の全研究の再学習・再実行ではありません。曲率履歴の入力例にも未来の答えは入っていません。

グラフと主6軌道の配列は再生成します。

```sh
python research/run_research.py
python research/plot_results.py
```

ヘッドレス環境では `MPLBACKEND=Agg` を指定してください。再実行するとJSONの実行時間・環境欄なども更新されます。独立した出力先には `python research/run_research.py --out /path/to/output` を使えます。

## 現在の結果

同じ初期密度・圧力・幾何学を持つが、高次の運動量分布が異なる正の粒子分布を構成しました。圧力の時間変化に最初の省略が入り、曲率の時間微分へ伝わります。

| 初期に一致させた情報 | 最初に異なる曲率の時間微分 | t=0.45の曲率の対称相対差 |
|---|---|---:|
| 粒子数・密度・圧力・幾何学 | 1階 | 0.325256% |
| 上記と4次モーメント | 2階 | 0.095363% |
| 上記と6次モーメント | 3階 | 0.030313% |

**各行は別の状態対です。最後の列は予測誤差ではありません。** 指定した低次の初期情報では当て分けられない反例であり、精密な過去履歴を含む全観測が同じという主張ではありません。

モデルは一様なKantowski–Sachs Einstein–Vlasov系です。質量ゼロの対照では、今回の初期等方分布に限りスペクトルの差が幾何学へ影響しない結果も保存しています。量子相関・衝突・蒸発・回転・外部接続はこの研究には含みません。

## 過去の保存版とニューラル予測器

このGit初回公開は、最新研究を再実行できるソース・数値結果と、解析的／多項式／適応幅の予測器を収録しています。**会話で配布したZIP全体のミラーではありません。過去の学習済み重み・全履歴・生成画像は、このコミットには含めていません。**

元の `curvature_predictors_and_state_audit.zip` を会話から保存したうえで、次の操作により歴史的な `legacy/` 一式をローカルへ復元できます。外部から自動ダウンロードはしません。

```sh
python tools/restore_legacy.py /path/to/curvature_predictors_and_state_audit.zip
python -m pip install -r requirements-neural.txt
python predict_curvature.py example_input.json --method neural
python verify_predictors.py
```

復元ツールは元ZIPのSHA-256を照合し、変更済みファイルを上書きしません。ニューラル版は元のLTBモデル単位で過去[-0.30,0]の41点を要求し、未来[0.10,0.20,0.30,0.40]を出力します。現在の対数曲率の微分も必要です。別のモデルや実測への有効性は未検証です。

詳しい来歴と公開時の変更点は [PROVENANCE.json](PROVENANCE.json) に記録しています。MITライセンスは初期リポジトリのものを保持しています。
