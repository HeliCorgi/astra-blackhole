# Bianchi IX timeless class operator 前段：平方根枝とWDW constraintのfactorization監査

## なぜ先にこれをするか

現行Bianchi IX量子計算は内部時計 s=-alpha を選び、

i hbar d_s psi = G(s) psi,
G(s) = -sqrt(A(s)),
A(s) = -hbar^2 Delta_beta + W(s,beta)

を直接時間発展させている。

Halliwellのtimeless decoherent historiesは、外部時間を選ばず

H Psi = 0

というHamiltonian constraintの解に対し、constraintと可換なclass operatorを作る。complex-potential法では「領域に入らない」class operatorが複素ポテンシャルのscattering/S-matrixとして構成される。Halliwellはこのclass operatorがconstraintと可換であることを重要条件としている（arXiv:0909.2597; 1108.5991）。

しかし現行平方根枝をそのまま二階WDW constraint

C_WDW = (i hbar d_s)^2 - A(s)

の一方の周波数枝と同一視してはいけない。

## 正確な差

P_s=i hbar d_s と書く。

P_s psi = G psi

なら

P_s^2 psi = P_s(G psi)
            = i hbar (d_s G) psi + G P_s psi
            = i hbar (d_s G) psi + G^2 psi.

G^2=Aなので

[(i hbar d_s)^2-A] psi
= i hbar (d_s G) psi.

したがってAがs依存するBianchi IXでは一般にゼロではない。

operator factorizationでも

(P_s-G)(P_s+G)
= P_s^2-G^2+[P_s,G]
= P_s^2-A+i hbar d_s G,

となる。逆順なら最後の符号が反転する。

## この監査

既存packetを既存平方根Hamiltonianでs=2から進め、s=3.0,4.6,5.5で

R = i hbar (d_s G) psi

を数値評価し、

||R|| / ||A psi||

を「naive second-order WDW constraintへ移したときの局所factorization mismatch」の尺度として記録する。

d_s Gは中心差分、sqrt(A)はLanczos。delta=.02/.01/.005とKrylov 60/72/84を別々に振る。

## 解釈

非ゼロなら「timeless WDWが間違い」ではない。
意味は限定的で、

- 現行の内部時計平方根量子化
- naiveな二階WDW constraint量子化

が、単なる座標変換・表現変換ではないということ。

次のHalliwell型timeless実装では、どのconstraint operatorを採るかを新しい模型選択として明示する必要がある。

候補は少なくとも二つある。

1. naive second-order WDW型 C=P_s^2-A(s)
2. 現行branchを正確に埋め込むfirst-order parametrized constraint C=P_s-G(s)

後者は現行ダイナミクスを保つが、標準的なKlein-Gordon型WDW constraintとは別物。
前者はHalliwellのminisuperspace形式に近いが、現行branchとは量子化が変わる。

## 限界

この段階ではtimeless S-matrix class operator、induced physical inner product、zero-constraint spectral stateをまだ構成しない。
目的は量子化の境界を先に固定すること。
