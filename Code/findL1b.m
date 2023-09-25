
(* ***********************************

GL(3) L-function code

These programs allow one to find, refine, and evaluate L-functions.
The code is general, but development concentrates on degree-3
L-functions.

The approach is based on a smoothed approximate functional equation,
where we exploit the ability to use several test functions.

4b: program was using ever increasing memory, presumably because we
keep increasing the precision (and we were calculating the detectors
to more precision than the basic equations).  So changed to use the
same precision for the detectors (in solveeqns6f.m) and added the function 
ClearGam[] to clear the values of the Gam[] (Gamma) function in memory.

4c: changed findSumlim so that it first skips by 10, then backs up and
skips by 2.  So the number of terms needed will always be even.  This should
not change the number of unknowns needed, and will make it twice as fast.

4d: changed Lambda so that del >= 1 is the sumlim  (was del > 1)

4e: changed gg to begin at 1/2

L1b: changing the FE format.  Added functions to convert between formats.
****************************** *)

(* the new format for functional equations is
        { {{del_1,lam_1},...}, {{kappa_1,mu_1},...}, N, eps}
   and the old format was
        { {del_1/2,...,kappa_1,...}, {lam_1/2,...,mu_1,...}, {1/2,...,1,...},
           Sqrt[N]/(Sqrt[Pi]...(2Pi)...), eps, 2^(-I (mu_1 + ...)) }

   Most of the code uses the old format.
   (As of version L1b, L, Z, testandsave, and findstartingvalues use the new format.)

   Also, the oldold format omitted the final 2^(-I xxx) phase factor, mistakenly
   assuming it was always 1.
*)


FEnewtoold[fe_]:= Block[{j, reshifts, imshifts, sfactors, Qfactor, sign, phasefactor},
(*	Print["fe in new notation", fe];  *)
        reshifts = Flatten[{Table[fe[[1,j,1]]/2,{j, 1, Length[fe[[1]]]}],
                        Table[fe[[2,j,1]],{j, 1, Length[fe[[2]]]}]}];
        imshifts = Flatten[{Table[fe[[1,j,2]]/2,{j, 1, Length[fe[[1]]]}],
                        Table[fe[[2,j,2]],{j, 1, Length[fe[[2]]]}]}];
        sfactors = Flatten[{Table[1/2,{j, 1, Length[fe[[1]]]}],
                        Table[1,{j, 1, Length[fe[[2]]]}]}];
        Qfactor = Sqrt[fe[[3]]]/( Pi^(Length[fe[[1]]]/2)(2 Pi)^(Length[fe[[2]]]) );
        sign = fe[[4]];
        phasefactor = 2^(-I Sum[fe[[2,j,2]],{j, 1, Length[fe[[2]]]}]);
(*	Print["fe in old notation", {reshifts, imshifts, sfactors, Qfactor, sign, phasefactor}];  *)
        {reshifts, imshifts, sfactors, Qfactor, sign, phasefactor}
];

FEoldtonew[fe_] := Block[{j,thegammaRs,thegammaCs,theConductor,theSign,thereshifts,theimshifts,thekappas,thisdel,thislam,thiskappa,thisGammaR,thisGammaC},
  thegammaRs = {};
  thegammaCs = {};
  theConductor = fe[[4]]^2;
  theSign = fe[[5]];
  thereshifts = fe[[1]];
  theimshifts = fe[[2]];
  thekappas = fe[[3]];
  For[j = 1, j <= Length[thekappas], ++j,
   If[thekappas[[j]] == 1/2,
    thisdel = 2 thereshifts[[j]];
    thislam = 2 theimshifts[[j]];
    thisGammaR = {thisdel, thislam};
    AppendTo[thegammaRs, thisGammaR];
    theConductor *= Pi (* because we squared initially*),
    (* else, this kappa is 1*)
    thiskappa = thereshifts[[j]];
    thislam = theimshifts[[j]];
    thisGammaC = {thiskappa, thislam};
    AppendTo[thegammaCs, thisGammaC];
    theConductor *= (2 Pi)^2]
   ];
  {thegammaRs, thegammaCs, theConductor, theSign}
  ];

(* Rubinstein's g(s): test function in the approximate functional equation *)
(* gg[b_,w_]:= (1+w)^(b[[1]]) E^(-1*b[[2]] I w + b[[3]] w^2) *)
gg[b_,w_]:= (w - 1/2)^(b[[1]]) E^(-1*b[[2]] I w + b[[3]] w^2) 
       (* NB: this function is indeterminate when encountering 0^0 *)


gg[u_,b_,v_]:= (1+v-u)^(b[[1]]) E^(-1*b[[2]] I (v-u) + b[[3]] (v-u)^2)

gg[gflag_,u_,b_,v_]:=If[gflag==0,gg[b,v],gg[u,b,v]];  (* gflag ==1 is Stefan's g*)

(* Gamma factors in the functional equation: H(s)L(s)=H(1-s)L(1-s) *)
HOLD[FE_,s_,PRECIS_]:=Block[{lambdaR,lambdaI,kappa,Qa,FACTORQ},
        lambdaR=FE[[1]]; lambdaI=FE[[2]]; kappa=FE[[3]]; QQ=FE[[4]];
        FACTOR=If[Length[FE] < 6, 1, FE[[6]]];
        FACTOR (QQ^s  Product[
        Gam[kappa[[j1]] (s) + lambdaR[[j1]] + I lambdaI[[j1]],PRECIS],{j1,1,Length[kappa]}])];

Lambdaold[FE_,b_,s_,Ev_,gflag_,PRECIS_]:=Block[{QQ,omega,nu,istep,ev,del,LIMS,gs,extrafactor},
        QQ=FE[[4]]; omega=FE[[5]];
        extrafactor = If[Length[FE] < 6, 1, FE[[6]]];
        nu=Ev[[1]];
        del=Ev[[2]];
	istep=stepsizeRM[nu, PRECIS]; 
(*Print["istep = ",istep]; *)
	ev={nu,istep};
	If[del>=1,LIMS=del,
	  If[del>0,LIMS=findSumlim[FE, b, s, ev, del,gflag,PRECIS,1][[1]],
	    LIMS=findSumlim[FE, b, s, ev, -del,gflag,PRECIS,0][[1]]]];
        gg[gflag,s,b,s]^(-1) ( QQ^s *
            extrafactor Sum[(bb1[n] + I bb2[n])/n^s f1[FE,b,s,n,ev,gflag,PRECIS],{n,1,LIMS}] +
            Conjugate[extrafactor] omega QQ^(1-s) *
            Sum[(bb1[n] - I bb2[n])/n^(1-s) f2[FE,b,1-s,n,ev,gflag,PRECIS],{n,1,LIMS}])];

H[FE_,s_,PRECIS_] := HOLD[FEnewtoold[FE],s,PRECIS];
Lambda[FE_,b_,s_,Ev_,gflag_,PRECIS_]:= Expand[Lambdaold[FEnewtoold[FE],b,s,Ev,gflag,PRECIS]];
L[FE_,b_,s_,Ev_,gflag_,PRECIS_]:= Expand[Lold[FEnewtoold[FE],b,s,Ev,gflag,PRECIS]];
Z[FE_,b_,s_,Ev_,gflag_,PRECIS_]:= Expand[Zold[FEnewtoold[FE],b,s,Ev,gflag,PRECIS]];

fourpoints[thisL_, fcn_] := Block[{},
   testcase = thisL;
   testfe = testcase[[1, 2]] /. {XX[1] -> testcase[[1, 1, 1]], XX[2] -> testcase[[1, 1, 2]]};
   lhalf = evaluateLfromAp[testfe, fcn, 1/2, {3, 10^-10}, 1, 50, testcase[[1, 3]], testcase[[1, 4]]];
   lone = evaluateLfromAp[testfe, fcn, 1, {3, 10^-10}, 1, 50, testcase[[1, 3]], testcase[[1, 4]]];
   lambdahalf = evaluateLambdafromAp[testfe, fcn, 1/2, {3, 10^-10}, 1, 50, testcase[[1, 3]], testcase[[1, 4]]];
   lambdaone = evaluateLambdafromAp[testfe, fcn, 1, {3, 10^-10}, 1, 50, testcase[[1, 3]], testcase[[1, 4]]];
   lamsize = Floor[Log[10, Abs[lambdahalf]]] - 9;
   N[Chop[{{"L(1/2)", lhalf}, {"L(1)", lone}, {"Lambda(1/2)", lambdahalf}, {"Lambda(1)", lambdaone}}, 10^lamsize], 16]
   ];
fourpoints[thisL_] := fourpoints[thisL, {0, 0, 0}];

evaluateZfromAp[FE_, b_, s_, Ev_, gflag_, PRECIS_, ep_, ap_] := 
  Block[{numterms},
   If[ep[[1, 1]] > 0, numterms = NextPrime[Length[ap]] - 1,
    numterms = NextPrime[Length[ap]/2] - 1];
   thisEv = {Ev[[1]], numterms};
   rawZ = Z[FE, b, s, thisEv, gflag, PRECIS];
   theseAnSubs = subsAn[anFromAp[ep, ap, numterms]];
   rawZ /. theseAnSubs];

evaluateLfromAp[FE_, b_, s_, Ev_, gflag_, PRECIS_, ep_, ap_] :=
  Block[{numterms},
   If[ep[[1, 1]] > 0, numterms = NextPrime[Length[ap]] - 1,
    numterms = NextPrime[Length[ap]/2] - 1];
   thisEv = {Ev[[1]], numterms};
   rawZ = L[FE, b, s, thisEv, gflag, PRECIS];
   theseAnSubs = subsAn[anFromAp[ep, ap, numterms]];
   rawZ /. theseAnSubs];

evaluateLambdafromAp[FE_, b_, s_, Ev_, gflag_, PRECIS_, ep_, ap_] :=
  Block[{numterms},
   If[ep[[1, 1]] > 0, numterms = NextPrime[Length[ap]] - 1,
    numterms = NextPrime[Length[ap]/2] - 1];
   thisEv = {Ev[[1]], numterms};
   rawZ = Lambda[FE, b, s, thisEv, gflag, PRECIS];
   theseAnSubs = subsAn[anFromAp[ep, ap, numterms]];
   rawZ /. theseAnSubs];

Lold[FE_,b_,s_,Ev_,gflag_,PRECIS_]:=Lambdaold[FE,b,s,Ev,gflag,PRECIS]/HOLD[FE,s,PRECIS];

Zold[FE_,b_,s_,Ev_,gflag_,PRECIS_]:=Lambdaold[FE,b,s,Ev,gflag,PRECIS]/Abs[HOLD[FE,s,PRECIS]];

Lold[FE_,b_,s_List,Ev_,gflag_,PRECIS_]:=Table[{s[[js]], L[FE,b,s[[js]],Ev,gflag,PRECIS]},{js,1,Length[s]}];
Zold[FE_,b_,s_List,Ev_,gflag_,PRECIS_]:=Table[{s[[js]], Z[FE,b,s[[js]],Ev,gflag,PRECIS]},{js,1,Length[s]}];

SetAttributes[bb1, NHoldAll];
SetAttributes[bb2, NHoldAll];

f1[FE_,b_,s_,n_,ev_,gflag_,PRECIS_]:=Block[{lambdaR,lambdaI,kappa,QQ,omega,nu,istep,LIMI1m,LIMI1p},
        {lambdaR,lambdaI,kappa,QQ,omega}=FE;
        nu=ev[[1]];istep=ev[[2]];
{LIMI1m, LIMI1p, junk1,junk2}=findIlim[FE,b,s, n, ev,gflag,PRECIS];
        (1/(2 Pi)) istep Sum[ff1[lambdaR,lambdaI,kappa,QQ,b,s,n,y,nu,gflag,PRECIS], {y,LIMI1m,LIMI1p,istep}]]

f2[FE_,b_,s_,n_,ev_,gflag_,PRECIS_]:=Block[{lambdaR,lambdaI,kappa,QQ,omega,nu,istep,LIMI2m,LIMI2p},
        {lambdaR,lambdaI,kappa,QQ,omega}=FE;
        nu=ev[[1]];istep=ev[[2]];
{junk1,junk2,LIMI2m, LIMI2p}=findIlim[FE,b,s, n, ev,gflag,PRECIS];
(1/(2 Pi)) istep Sum[ff2[lambdaR,lambdaI,kappa,QQ,b,s,n,y,nu,gflag,PRECIS], {y,LIMI2m,LIMI2p,istep}]]

ff1[lambdaR_,lambdaI_,kappa_,QQ_,b_,s_,n_,yy_,nu_,gflag_,PRECIS_]:= Product[
        Gam[kappa[[j1]] (nu + I yy +s) + lambdaR[[j1]] + I lambdaI[[j1]],PRECIS],{j1,1,Length[kappa]}]*
        (1/(nu + I yy)) gg[gflag,s,b,s + nu + I yy] (QQ/n)^(nu + I yy);

ff2[lambdaR_,lambdaI_,kappa_,QQ_,b_,s_,n_,yy_,nu_,gflag_,PRECIS_]:= Product[
        Gam[kappa[[j1]] (nu + I yy +s) + lambdaR[[j1]] - I lambdaI[[j1]],PRECIS],{j1,1,Length[kappa]}]*
        (1/(nu + I yy)) gg[gflag,1-s,b,1 - s - nu - I yy] (QQ/n)^(nu + I yy);

(* --------  *)

f1[FE_,b_,s_,n_,ev_,gflag_,PRECIS_]:=Block[{thisFE,lambdaR,lambdaI,kappa,QQ,omega,nu,istep,LIMI1m,LIMI1p},
        thisFE=FE;
        If[Length[thisFE] == 5, AppendTo[thisFE, 1]];
        {lambdaR,lambdaI,kappa,QQ,omega,extrafactor}=thisFE;
        nu=ev[[1]];istep=ev[[2]];
{LIMI1m, LIMI1p, junk1,junk2}=findIlim[thisFE,b,s, n, ev,gflag,PRECIS];
        (1/(2 Pi)) istep Sum[ff1[lambdaR,lambdaI,kappa,QQ,b,s,n,y,nu,gflag,PRECIS,extrafactor], {y,LIMI1m,LIMI1p,istep}]]

f2[FE_,b_,s_,n_,ev_,gflag_,PRECIS_]:=Block[{lambdaR,lambdaI,kappa,QQ,omega,nu,istep,LIMI2m,LIMI2p},
        thisFE=FE;
        If[Length[thisFE] == 5, AppendTo[thisFE, 1]];
        {lambdaR,lambdaI,kappa,QQ,omega,extrafactor}=thisFE;
        nu=ev[[1]];istep=ev[[2]];
{junk1,junk2,LIMI2m, LIMI2p}=findIlim[thisFE,b,s, n, ev,gflag,PRECIS];
(1/(2 Pi)) istep Sum[ff2[lambdaR,lambdaI,kappa,QQ,b,s,n,y,nu,gflag,PRECIS,extrafactor], {y,LIMI2m,LIMI2p,istep}]]

ff1[lambdaR_,lambdaI_,kappa_,QQ_,b_,s_,n_,yy_,nu_,gflag_,PRECIS_,extrafactor_]:= Product[
        Gam[kappa[[j1]] (nu + I yy +s) + lambdaR[[j1]] + I lambdaI[[j1]],PRECIS],{j1,1,Length[kappa]}]*
        (1/(nu + I yy)) gg[gflag,s,b,s + nu + I yy] (QQ/n)^(nu + I yy);

ff2[lambdaR_,lambdaI_,kappa_,QQ_,b_,s_,n_,yy_,nu_,gflag_,PRECIS_,extrafactor_]:= Product[
        Gam[kappa[[j1]] (nu + I yy +s) + lambdaR[[j1]] - I lambdaI[[j1]],PRECIS],{j1,1,Length[kappa]}]*
        (1/(nu + I yy)) gg[gflag,1-s,b,1 - s - nu - I yy] (QQ/n)^(nu + I yy);

(* --------  *)


Gam[xx_,15]:=Gam[xx,15]=N[Gamma[xx]];

Gam[xx_,PREC_]:=Gam[xx,PREC]=N[Gamma[xx],PREC];

ClearGam[]:=(Clear[Gam];
Gam[xx_,PREC_]:=Gam[xx,PREC]=N[Gamma[xx],PREC]);

stepsizeRM[nu_, PREC_] := Floor[100 (nu 2 Pi/Log[10]/PREC)]/100;

findSumlim[FE_, b_, s_, Ev_, del_,gflag_,PRECIS_,absflag_] := (* determine how many terms
                are needed to compute L(s) to accuracy del.
                Returns a pair, where the first element is
                the number of terms needed *)
  Block[{QQ,HH,thescale}, QQ = FE[[4]];
   HH = Abs[HOLD[FE, s,PRECIS]];
   GG = Abs[gg[gflag,s, b, s]];
(*    nu=Ev[[1]];
    istep=stepsizeRM[nu, PRECIS];
    ev={nu,istep};
*)
    a1coef =
    Abs[QQ^s f1[FE, b, s, 1, Ev,gflag,PRECIS]] +
     Abs[QQ^(1 - s) f2[FE, b, 1 - s, 1, Ev,gflag,PRECIS]];
   a1coef = a1coef/(GG HH);
   If[absflag==1,thescale=1,thescale=a1coef];
   (*need the nth term to be less than del, or del*a1, depending on wheher
	absflag is 1 or not *)
   n = 10;
   While[(xx = (Abs[QQ^s f1[FE, b, s, n, Ev,gflag,PRECIS]] +
          Abs[QQ^(1 - s) f2[FE, b, 1 - s, n, Ev,gflag,PRECIS]] +
          Abs[QQ^s f1[FE, b, s, n + 1, Ev,gflag,PRECIS]] +
          Abs[QQ^(1 - s) f2[FE, b, 1 - s, n + 1, Ev,gflag,PRECIS]])/(GG HH)) >
     2*del*thescale,(* Print[n,"  ",xx,"  ",2*del*thescale]; *)
	n += 10; 
	If[n>200,Print["Error: n>200", "   ",xx, " ",FE,"  ",b,"  ",s]; Return[{5,0}]]];
   n = n - 8;
   While[(yy = (Abs[QQ^s f1[FE, b, s, n, Ev,gflag,PRECIS]] +
          Abs[QQ^(1 - s) f2[FE, b, 1 - s, n, Ev,gflag,PRECIS]] +
          Abs[QQ^s f1[FE, b, s, n + 1, Ev,gflag,PRECIS]] +
          Abs[QQ^(1 - s) f2[FE, b, 1 - s, n + 1, Ev,gflag,PRECIS]])/(GG HH)) >
     2*del*thescale, n += 2];
   {n, Log[Abs[a1coef]]}];

findIlim[FE_,b_,s_, n_, ev_,gflag_,prec_] := Block[
        {LIMi, L1m, L1p, L2m, L2p,   testf1,testf2,
          lambdaR,lambdaI,kappa,QQ,nu},
        lambdaR=FE[[1]]; lambdaI=FE[[2]]; kappa=FE[[3]]; QQ=FE[[4]];
        nu=ev[[1]];
        testf1[xx_]:=Log[10,Abs[ff1[lambdaR,lambdaI,kappa,QQ,b,s,n,xx,nu,gflag,15]]];
        testf2[xx_]:=Log[10,Abs[ff2[lambdaR,lambdaI,kappa,QQ,b,s,n,xx,nu,gflag,15]]];
  LIMi =
   10 + Floor[Sum[Abs[lambdaI[[j]]], {j, 1, Length[lambdaI]}] + Abs[Im[s]]];
  initvals1 = Table[testf1[xx], {xx, -LIMi, LIMi}];
  initvals2 = Table[testf2[xx], {xx, -LIMi, LIMi}];
  max1 = Max[initvals1];
  max2 = Max[initvals2];
  testlim = LIMi + 1;
  While[testf1[testlim] > max1 - prec, ++testlim];
  L1p = testlim;
  testlim = -LIMi - 1;
  While[testf1[testlim] > max1 - prec, --testlim];
  L1m = testlim;
  testlim = LIMi + 1;
  While[testf2[testlim] > max2 - prec, ++testlim];
  L2p = testlim;
  testlim = -LIMi - 1;
  While[testf2[testlim] > max2 - prec, --testlim];
  L2m = testlim;
  {L1m, L1p, L2m, L2p}]

(*  ****************************************  *)

convert[u_,v_]:={(u-v)/6,-(2u+v)/6}

roundfraction[frac_, prec_] := Block[{size, num},
   If[frac >= 1/1000, Return[frac],
    size = Floor[Log[frac]/Log[10] + 0.000123];  (* because Floor[Log[power of 10]/Log[10]] error *)
    num = Floor[frac/10^(size - prec)];
    num /(10^(prec - size) + 1)
    ]
   ];
        
roundfraction[frac_]:=roundfraction[frac,1];  (* 1 means 2 digits in numerator *)

(*
convert[11.761250, 18.9024]
{-1.19019, -7.07082}
*)

makeequations[FE_,eis_, glis_, svals_,Ev_,gflag_,PRECIS_] :=
Block[{v,w,j,k,sol,bvals,eqns,numeqns},
  FEtmp=FE;
  FEtmp[[2]]=eis;
  numeqns = Length[glis];
  For[j = 1, j <= Length[glis], ++j,
   For[k=1,k<=2,++k,
   v[j,k] = Expand[
       L[FEtmp,
        glis[[j,k]],
        svals[[j]], Ev,gflag,PRECIS]];
   w[j,k] = v[j,k];
    ]];
  Table[Expand[w[je,1] - w[je,2]], {je, 1, numeqns}]
(*  eqns = (eqns/.{bb1[1]->1, bb2[1]->0}); *)
  ]

makeequationsNEW[FE_,eis_, glis_, svals_,Ev_,gflag_,PRECIS_] :=
Block[{v,w,j,k,sol,bvals,eqns,numeqns},
  FEtmp=FE;
  FEtmp = (FE/.Table[XX[j]-> eis[[j]],{j,1,Length[eis]}]);
(*  FEtmp[[2]]=eis;  *)
  numeqns = Length[glis];
  For[j = 1, j <= Length[glis], ++j,
   For[k=1,k<=2,++k,
   v[j,k] = Expand[
       L[FEtmp,
        glis[[j,k]],
        svals[[j]], Ev,gflag,PRECIS]];
   w[j,k] = v[j,k];
    ]];
  Table[Expand[w[je,1] - w[je,2]], {je, 1, numeqns}]
(*  eqns = (eqns/.{bb1[1]->1, bb2[1]->0}); *)
  ]


