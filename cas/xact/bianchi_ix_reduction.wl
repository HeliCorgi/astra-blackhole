(* Bianchi IX GR reduction cross-check for xAct/Wolfram Language.
   This file is intentionally committed before execution on a licensed
   Wolfram+xAct runner. Until an executed manifest is stored, the xAct
   obligation remains PENDING.

   The coordinate metric is the same SU(2) Euler-angle metric used by the
   Cadabra backend, with omega_i=sigma_i/2. The calculation below uses
   explicit component Christoffel/Ricci formulas in Wolfram Language after
   loading xAct/xTensor and xCoba; future revisions may replace the component
   loop with CTensor/MetricCompute without changing the convention manifest.
*)

Needs["xAct`xTensor`"];
Needs["xAct`xCoba`"];

ClearAll[th, ph, ps, a, b, c, alpha, bp, bm, N, ad, bpd, bmd,
  pa, pp, pm, s];

coords={th,ph,ps};
h=1/4 {
 {a^2 Cos[ps]^2+b^2 Sin[ps]^2,
  (a^2-b^2) Cos[ps] Sin[ps] Sin[th], 0},
 {(a^2-b^2) Cos[ps] Sin[ps] Sin[th],
  (a^2 Sin[ps]^2+b^2 Cos[ps]^2) Sin[th]^2+c^2 Cos[th]^2,
  c^2 Cos[th]},
 {0,c^2 Cos[th],c^2}
};
hi=FullSimplify[Inverse[h]];

Gamma=Table[
 FullSimplify[1/2 Sum[
   hi[[i,l]] (D[h[[j,l]],coords[[k]]] + D[h[[k,l]],coords[[j]]] -
              D[h[[j,k]],coords[[l]]]),{l,3}]],
 {i,3},{j,3},{k,3}];

Riemann=Table[
 FullSimplify[
  D[Gamma[[i,j,l]],coords[[k]]]-D[Gamma[[i,j,k]],coords[[l]]] +
  Sum[Gamma[[i,k,m]] Gamma[[m,j,l]]-Gamma[[i,l,m]] Gamma[[m,j,k]],{m,3}]
 ],
 {i,3},{j,3},{k,3},{l,3}];

Ricci=Table[FullSimplify[Sum[Riemann[[i,j,i,l]],{i,3}]],{j,3},{l,3}];
R3=FullSimplify[Sum[hi[[j,l]] Ricci[[j,l]],{j,3},{l,3}]];
closed=FullSimplify[
 2(-a^4-b^4-c^4+2(a^2 b^2+b^2 c^2+c^2 a^2))/(a^2 b^2 c^2)
];
If[FullSimplify[R3-closed]=!=0,Print["R3 mismatch"];Exit[2]];

aa=Exp[alpha+bp+Sqrt[3] bm];
bb=Exp[alpha+bp-Sqrt[3] bm];
cc=Exp[alpha-2 bp];
V=(Exp[-8 bp]+2 Exp[4 bp](Cosh[4 Sqrt[3] bm]-1)-
    4 Exp[-2 bp] Cosh[2 Sqrt[3] bm])/6;
If[FullSimplify[closed/.{a->aa,b->bb,c->cc}+12 Exp[-2 alpha] V]=!=0,
   Print["Misner curvature mismatch"];Exit[3]];

db1=bpd+Sqrt[3] bmd; db2=bpd-Sqrt[3] bmd; db3=-2 bpd;
ks={(ad+db1)/N,(ad+db2)/N,(ad+db3)/N};
kin=FullSimplify[Total[ks^2]-Total[ks]^2];
If[FullSimplify[kin-6(-ad^2+bpd^2+bmd^2)/N^2]=!=0,
   Print["Kinetic mismatch"];Exit[4]];

L=Exp[3 alpha]/(2N)(-ad^2+bpd^2+bmd^2)-N Exp[alpha] V;
mom={D[L,ad],D[L,bpd],D[L,bmd]};
vel={ad->-N Exp[-3 alpha] pa,bpd->N Exp[-3 alpha] pp,bmd->N Exp[-3 alpha] pm};
C=FullSimplify[(pa ad+pp bpd+pm bmd-L)/N/.vel];
Cexpected=Exp[-3 alpha]/2(-pa^2+pp^2+pm^2)+Exp[alpha] V;
If[FullSimplify[C-Cexpected]=!=0,Print["Constraint mismatch"];Exit[5]];

Print["ASTRA_XACT_BIANCHI_IX_OK"];
Print[InputForm[R3]];
Print[InputForm[C]];
