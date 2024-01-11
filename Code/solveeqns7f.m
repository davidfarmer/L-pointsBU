
printflag=1;
(**********************

6g:  a_{32} was previously treated as an unknown when we substitute
	the relations from the local Euler factors, but it was not
	on the list of unknowns when we solve the equations.  So we
	are adding it as an unknown when solving.  This will give an
	error check, because we can express a_{32} in terms of a_2.

     when zooming in on an eigenvalue, we have to increase the PRECISion.
	How to do this is unclear.  This version looks at the precision
	of the coefficients at the new "eigenvalue" and does the following:
	A. program failed to produce a new "eigenvalue:
		re-run with PRECIS increased by 10
	B. a new eigenvalue and coefficients were produced
		(p = precision of the coefficients)
	p < 2	re-run with PRECIS increased by 10  (that probably is too conservative)
	2 < p < 7  use new eignevalue and coefficients, and increase PRECIS by 10
	7 < p < 12  use new eignevalue and coefficients, and increase PRECIS by 5
	12 < p      use new eignevalue and coefficients, and keep the same PRECIS

6h:	Added a one-variable zoom[] function.

    TO DO:
	HIT needs to be modified:  report the distance from the center of the
	square (on the scale of the square), not just whether the medin was in
	the square.

	bb1[32], bb2[32] need to appear in numerical order (not last).  There is
	a problem when zooming and adding more coefficients.

	Need to properly implement the multiplicative relationships.

	David still thinks that we are throwing away too many digits when we solve
	and then find the interpolated new eigenvalue and coefficients.  Go through
	the steps carefully, focusing on what happens after we take the median to
	find the new eigenvalue.  That median should be set to have high accuracy,
	so that we don't lost a lot when we interpolate the coefficients.
	David tried to address this with the command   newpt=SetPrecision[newptIN,100]; 
	in  interppairs[]. It didn't work.

6i:	multsubs[sumlim]  returns the multiplicative substitutions for bb1[n],bb2[n]
		for n up to (and including) sumlim

6j:	eulersubs[3,p,cpr,cpi]  substitutions for the prime power recursions, for
		j=2,3,4, for bb1[p^j] and bb2[p^j] in terms of bb1[p],bb2[p],cpr, and cpi.
		Here cpr and cpi are the real and imaginary parts of \chi(p) occuring
		in the general good degree 3 Euler factor:
		1/(1 - a_p x + \chi (p)\overline {a_p} x^2 - \chi (p) x^3)
     *** that local factor was written wrong in the documentation of previous versions,  but the implementation is right
		And as usual, bb1[p] and bb2[p] are the real and imaginary parts of 
		the unknown a_p.  (Derivation of this function is in eulerfactors1b.nb )

	eulersubsbad[2,p_]: substitutions for the prime power recursions, for
                j=3,4,5,6 for bb1[p^j] and bb2[p^j] in terms of bb1[p],bb2[p],bb1[p^2], and bb2[p^2],
		for the most general degree 2 Euler factor:
			1/(1 - a_p x - (apsq - ap^2) x ^2) .
		(Derivation of this function is in eulerfactors1b.nb )

        eulersubsbad[1,p_]: substitutions for the prime power recursions, for
                j=2,...,10 for bb1[p^j] and bb2[p^j] in terms of bb1[p] and bb2[p],
                for the general degree 1 Euler factor:
                        1/(1 - a_p x) .
                (Derivation of this function is in eulerfactors1b.nb )

        eulersubsbad[0,p_]: substitutions for the prime power recursions, for
		j=1,...,10 for bb1[p^j] and bb2[p^j] for a missing Euler factor.
		In other words,  bb1[p^j] and bb2[p^j] are both zero.

	Describing the Euler product.  We introduce a new argument to the functions
		used to find L-functions:  EP.  This describes the Euler product.
		It has the form:   EP={good part, bad part}, where "good part"
			describes the factors at the good primes, and "bad part"
			lists the bad primes and then describes each bad local factor.
		good part:  {deg, chi}, where deg is a positive integer, and
			chi is a character.  For example:
			{1} is the trivial character
			{1,-1,-1,1,0}  is the real character mod 5
			{1,I,-I,-1,0}  is one of the non-real characters mod 5
			If deg is negative, that means real coefficients
		bad part:  {bad primes, description of bad factors}
			"bad primes" is just a (possibly empty) list of bad primes.
			"description of bad factors" currently only supports an arbitrary
				bad factor of a given degree, which is specified by giving
				the degree of the bad factor.  For example,
			{{3,5},{2,1}} says that 3 and 5 are bad, with a general degree 2
				factor at 3, and a linear factor at 5.  So a_3 and a_9 are
				unknowns, as is a_5. (the number of unknows for that prime
				is the same as the degree of the bad factor.)
		12/14/23 add another list to the bad part.  For example,
                        { {4, {1}}, { {7,11}, {3,2}, {{0,0,2},{0,1}} } }
                meaning that the reciprocal roots at 7 have size 1, 1, 1/7, and the
                reciprocal roots at 11 have size 1, 1/Sqrt[11] .

	theunknowns[ep,lim]:  list of the unknowns bb1[j],bb2[j] up to lim, for an
		Euler product ep.

	converteqnsALL: changed first argument from "lev" to "EP"
		modified function to convert the relations from the Euler
		product in 3 steps:  multiplicative, good primes, bad primes.
		Those call multsubs, goodprimepowersubs, and badprimepowersubs.

	searchgrid2: added a final argument "EP", which defaults to degree 3.
		trivial character, level 1 when omitted.
		(Same for zoom1)

6k:	converteqnsALL added another argument (absflag), which required the same modification to
	normalizeequationsA1, and temporarily made a version of
	searchgrid2 which does not normalize the detector equations.

6m:	searchzoom:  refine the eigenvalue, in the general (any degree) case.
	pairtotriple:	now handles lists of any length
        eulersubs[4,p,cpr,cpi]  substitutions for the prime power recursions,
		for j=3,4, for bb1[p^j] and bb2[p^j] in terms of
		bb1[p],bb2[p],cpr, and cpi.

6n,o:	eulersubs[5,p,cpr,cpi]

6q:	eulersubsbad[\pm 3,1,2,{p,eps}]  This is bad degree 3, it really is the
		typical good degree 3 factor that we need for prime level Paramodular forms.
	Warning:  above is not the right syntax

6s:	findsolmult: changed how precision is determined, so that we use the high precision
		coefficients, not the small coefficients

7a:     findsolall: tries to find all the useful solutions to a given set of nonlinear equations.
        (deleted in 7f)
        findsolnear: tries to find a solution close to a given vector [to-do]
        (deleted in 7f)

ToDo:	bad degree 3 case
	better names for the bad factors.
	(also do the good degree 2 case and think about degree 5)

7b:     findsolone: like findsolmult, but repeatedly searches with slightly
              different starting values, stopping after finding one solution.
              (extended scalelist in findsolone without incrementing version)

7c:     realize500 handles the case of sign being a variable (by not making the equations real), so
        we force the detectors to be real.

7d:     realize500 was not making the equations real, so the output was garbage.

7e:     normalizeequationsA1 mistakenly assumed the coefficient of bb1[1] was a number.
        Now we check for that and do nothing if it is not a number.
        We have bad prime substitutions up to p^10.  This was not reflected in theunknowns[ep_,lim_],
        which caused problems with conductor 4 once we needed 32 coefficients.

        anFromAp: make list of all Dirichlet coefficients from the list of psprime values

        Add 3rd entry to EP, encoding what we know about the bad factors (from Ralf's tables).

7f:     Make scalelist much longer and more expansive at the end.
        Have findsolone, but hard code the number of times it tries.

        Add functions for evaluating the Z-function:
        fefromdata, evaluateZfromAp, evaluateZfromApEXTRA, quickZ, quickZEXTRA, boundtailsimple
        (TODO:  fix quickZ and quickZEXTRA so that precision is from input (currently 40))

        deleted zoom1, searchgrid2

        added substitutionsFromAp which combines subsAn and anFromAp, because those were not
              designed to handle unknown EPSILONp factors.

In next version:  adjust secant method to use a smaller perturbation when the coefficiient
                  is already known with some accuracy (line 1171)
                  consider removing a_32 as an unknown

********************)

(* functions for making and solving the equations we use to search
   for GL(3) L-functions  *)

pairtotriple[pt_]:={pt[[1]],pt[[2]],-pt[[1]]-pt[[2]]};  (* eigenvalues add to zero *)
pairtotriple[pt_]:=Append[pt,-1*Sum[pt[[j]],{j,1,Length[pt]}]];

substopairs[{}]={};

substopairs[lis_]:=Block[{j}, If[Head[lis[[1]] ]==List,
	Table[substopairs[lis[[j]]],{j,1,Length[lis]}],
	Table[{lis[[j,1]],lis[[j,2]]},{j,1,Length[lis]}],
	Table[{lis[[j,1]],lis[[j,2]]},{j,1,Length[lis]}]]];
(*above and following, you need N[] inside the Depth[] because {2} has depth 2,
but {Sqrt[2]} has depth 3 *)

pairstosubs[lis_]:=Block[{j}, If[Head[lis[[1,1]] ]==List,
	Table[pairstosubs[lis[[j]]],{j,1,Length[lis]}],
        Table[lis[[j,1]]->lis[[j,2]],{j,1,Length[lis]}],
        Table[lis[[j,1]]->lis[[j,2]],{j,1,Length[lis]}]]];

substovals[lis_]:=Block[{j}, If[Head[lis[[1]] ]==List,
	Table[substovals[lis[[j]]],{j,1,Length[lis]}],
	Table[lis[[j,2]],{j,1,Length[lis]}],
	Table[lis[[j,2]],{j,1,Length[lis]}]]];

interppairs[pts_,varvals_,newptIN_]:=Block[{numvars,numpts,valstointerpolate,newval,theans,j,k},
(*
	pts:  a list of point (as ordered pairs) where we have known values
	varvals:  A list of sublists, where each sublist is a list of
		{variable, value}  pairs for the corresponding point
		  (also works if the elements are variable -> value, but the
		  function still returns an answer in the form ({variable, value} )
	newptIN: the new point where we want a guess for the {var,value} pairs
*)
If[Length[pts] != Length[varvals],Print["Length mismatch in interppairs"]; Return[]];
newpt=SetPrecision[newptIN,100];  (* this is sort-of a fudge, to prevent a low precision in the
	interpolated point from contaminating the data we are interpolating  *)
theans={};
numvars=Length[varvals[[1]]];
numpts=Length[pts];
For[j=1,j<=numvars,++j,
	valstointerpolate={};
	For[k=1,k<=numpts,++k,
	  AppendTo[valstointerpolate,Flatten[{pts[[k]],varvals[[k,j,2]]}]]];
	newval=((Fit[valstointerpolate,{1,x,y},{x,y}])/.{x->newpt[[1]],y->newpt[[2]]});
	AppendTo[theans,{varvals[[1,j,1]],newval}]
];
theans
];

interppairs[pts_,varvals_,newptIN_]:=Block[{numvars,numpts,valstointerpolate,newval,theans,j,k},
(*
        pts:  a list of points (as ordered tuples) where we have known values
        varvals:  A list of sublists, where each sublist is a list of
                {variable, value}  pairs for the corresponding point
                  (also works if the elements are variable -> value, but the
                  function still returns an answer in the form ({variable, value} )
        newptIN: the new point where we want a guess for the {var,value} pairs
*)
If[Length[pts] != Length[varvals],Print["Length mismatch in interppairs"]; Return[]];
newpt=SetPrecision[newptIN,100];  (* this is sort-of a fudge, to prevent a low precision in the
        interpolated point from contaminating the data we are interpolating  *)
theans={};
numvars=Length[varvals[[1]]];
numpts=Length[pts];
dim=Length[pts[[1]]];
For[j=1,j<=numvars,++j,
        valstointerpolate={};
        For[k=1,k<=numpts,++k,
          AppendTo[valstointerpolate,Flatten[{pts[[k]],varvals[[k,j,2]]}]]];
	If[dim==2,newval=((Fit[valstointerpolate,{1,x,y},{x,y}])/.{x->newpt[[1]],y->newpt[[2]]}),
	newval=((Fit[valstointerpolate,{1,x,y,z},{x,y,z}])/.{x->newpt[[1]],y->newpt[[2]], z->newpt[[3]]})];
   (*     newval=Apply[Interpolation[valstointerpolate],newptIN];  *)
        AppendTo[theans,{varvals[[1,j,1]],newval}]
];
theans
];

interppairs[pts_,varvals_,newptIN_]:=Block[{numvars,numpts,valstointerpolate,newval,theans,j,k,x},
(*
        pts:  a list of points (as ordered tuples) where we have known values
        varvals:  A list of sublists, where each sublist is a list of
                {variable, value}  pairs for the corresponding point
                  (also works if the elements are variable -> value, but the
                  function still returns an answer in the form ({variable, value} )
        newptIN: the new point where we want a guess for the {var,value} pairs
*)
If[Length[pts] != Length[varvals],Print["Length mismatch in interppairs"]; Return[]];
newpt=SetPrecision[newptIN,100];  (* this is sort-of a fudge, to prevent a low precision in the
        interpolated point from contaminating the data we are interpolating  *)
theans={};
numvars=Length[varvals[[1]]];
numpts=Length[pts];
dim=Length[pts[[1]]];
vars=Table[x[j],{j,1,dim}];
fitform=Flatten[{1,vars}];
For[j=1,j<=numvars,++j,
        valstointerpolate={};
        For[k=1,k<=numpts,++k,
          AppendTo[valstointerpolate,Flatten[{pts[[k]],varvals[[k,j,2]]}]]];
        If[dim==2,newval=((Fit[valstointerpolate,{1,x,y},{x,y}])/.{x->newpt[[1]],y->newpt[[2]]}),
(*
        newval=((Fit[valstointerpolate,{1,x,y,z},{x,y,z}])/.{x->newpt[[1]],y->newpt[[2]], z->newpt[[3]]})];
*)
        newval=((Fit[valstointerpolate,fitform,vars])/.Table[vars[[j]]->newpt[[j]],{j,1,Length[vars]}])];
   (*     newval=Apply[Interpolation[valstointerpolate],newptIN];  *)
        AppendTo[theans,{varvals[[1,j,1]],newval}]
];
theans
];


multipleinterppairs[pts_,varvalslist_,newpt_,numcompare_]:=Block[
{tmp,j},
(* This function is similar to interppairs, except that  varvalslist is a
list where each sublist is a list of possible {variable, value}  pairs.
So for every sublist we have to match the entries with those in the other
sublists that are closest.  Then we use those to interpolate.  Multiple
answers are returned*)
tmp=closestmatchmultiple[pts,varvalslist,numcompare];
Table[interppairs[tmp[[j,1]], tmp[[j,2]],newpt],{j,1,Length[tmp]}]
];

closestmatchmultiple[pts_,varvalslist_,numcompare_]:=Block[
{numpts,theans,jr,j,k,ptstointerpolate,valstointerpolate},
(* This function is similar to interppairs, except that  varvalslist is a
list where each sublist is a list of possible {variable, value}  pairs.
So for every sublist we have to match the entries with those in the other
sublists that are closest.  Then we use those to interpolate.  Multiple
answers are returned*)
numpts=Length[pts];
theans={};
(* first find the point with the smallest number of extries*)
smallestlen=Infinity;
smallestindex=0;
For[j=1,j<=numpts,++j,
        If[Length[varvalslist[[j]] ]< smallestlen,
                smallestindex=j;
                smallestlen=Length[varvalslist[[j]] ]]
   ];
numpairstomake=smallestlen;
For[j=1,j<=numpairstomake,++j, (*Print[j]; *)
        ptstointerpolate={pts[[smallestindex]]};
        valstointerpolate={varvalslist[[smallestindex,j]]};
(*Print["first valstointerpolate : ",valstointerpolate]; *)
        For[k=1,k<=numpts,++k, (*Print["k=",k];*) 
           If[k != smallestindex,
                AppendTo[ptstointerpolate,pts[[k]]];
                AppendTo[valstointerpolate,
        closestmatchpairs[valstointerpolate[[1]],varvalslist[[k]],numcompare]];
             ];
        ];
	(*Now we need to put the points back into their original order*)
	tmpptr=ptstointerpolate[[1]];
	tmpvvr=valstointerpolate[[1]];
	For[jr=2,jr<=smallestindex,++jr,
		ptstointerpolate[[jr-1]]=ptstointerpolate[[jr]];
		valstointerpolate[[jr-1]]=valstointerpolate[[jr]];
	    ];
	ptstointerpolate[[smallestindex]]=tmpptr;
	valstointerpolate[[smallestindex]]=tmpvvr;
        (*Print[ptstointerpolate]; *)
        (*Print[valstointerpolate]; *)
(*Return[];*)
        AppendTo[theans,{ptstointerpolate, valstointerpolate}]
];
theans
];

check4pts[pts_,answ_,detectors_]:=Block[{j(*,detectnums,detect,detectlist*)},
detectnums=Length[detectors[[1]]];

ttmmpp=closestmatchmultiple[pts,answ,4];
  detectlist = {};
alldetectors={};
dim=Length[pts[[1]]];
Print["Doing a ",dim, " dimensional search."];
For[kk=1,kk<=Length[ttmmpp],++kk,
  subdetectlist={};
  For[j=1,j<=Length[pts],++j,
     Print["j =  ",j]; 
     detect[j]=(detectors[[j]]/.ttmmpp[[kk,2,j]]);
     detect[j]=Re[detect[j]];   (* because when the sign is unknown, the equaitons are not necessarily real *)
     Print["the detector: ",detect[j]];
     AppendTo[alldetectors,{pts[[j]],detect[j]}];
  ];
  dimsubsets=Subsets[Table[j,{j,1,detectnums}],{dim}];
  For[j = 1, j <= Length[dimsubsets], ++j,
    testpts = dimsubsets[[j]];   (*Print[{j,k}]; *)
    thispt = nDsecant[pts,Table[detect[mm][[testpts]],{mm,1,Length[pts]}]];
    If[NumberQ[thispt[[1]]], 
	AppendTo[subdetectlist, thispt]
    ,
        Return[""]
    ]
(*
		nDsecant[pts,Table[detect[mm][[testpts]],{mm,1,Length[pts]}]]];
*)
  ];

  AppendTo[detectlist,subdetectlist];

]; (*for kk *)

(* since other programs appear to only use the [[1]] element of the returned list,
we are also sending back the detectors so that they can be used for something else.
The next line use to jsut be:  detectlist *)
Append[detectlist,alldetectors]
];

closestmatch[elt_,lis_,numterms_]:=Block[{},
	nelt=Take[elt,numterms];
	nlis=Table[Take[lis[[nn]],numterms],{nn,1,Length[lis]}];
	SortBy[nlis,N[Norm[#-nelt]]&][[1]]
];

closestmatchpairs[elt_, lis_, numterms_] := (* also should work for substitutions *)
  Block[{nlis},
If[numterms>Length[elt],Print["not enough terms in closestmatchpairs",elt,"ee",lis,"ll",numterms]; Return[{}]];
(*Print["in closestmatchpairs",elt,"ee",lis,"ll",numterms]; *)
nelt = elt;
   nlis = lis;
   SortBy[nlis, 
     N[Norm[Take[substovals[#1] - substovals[nelt], numterms]]] &][[1]]];

NumPtsToTest[stvls_]:=If[Length[stvls]==1,100 (*used to be 1000 *),50];
NumPtsToTest[stvls_]:=If[Length[stvls]==1,10 (*used to be 1000 *),50];


insquare[{fourpts_,ptlis_}]:=Block[{med,minx,miny,maxx,maxy},
(*
Print["insquare : ",fourpts,"  ",ptlis];
*)
If[ptlis=={},Return[False]];
minx=Min[Transpose[fourpts][[1]]];
maxx=Max[Transpose[fourpts][[1]]];
miny=Min[Transpose[fourpts][[2]]];
maxy=Max[Transpose[fourpts][[2]]];
med=Median[ptlis];
med[[1]]>minx && med[[1]]<maxx && med[[2]]>miny && med[[2]]<maxy
];

(* is this one ever used, or only the following one? *)
closepts[{eps_,ptlis_}]:=Block[{med,minx,miny,maxx,maxy,thiscp},
(*
Determine if a set of points is close together.
*)
If[ptlis=={},Return[Infinity]];
thiscp = If[eps[[1]]>0 && eps[[2]]>0,
    numpts=Length[ptlis];
    firstthird=Ceiling[numpts/5];  (* actually, first fifth!*)
    lastthird=Floor[4 numpts/5];  (* actually, last fifth!*)
    xpts=Sort[Transpose[ptlis][[1]]];
    ypts=Sort[Transpose[ptlis][[2]]];
    Max[(xpts[[lastthird]]-xpts[[firstthird]])/eps[[1]],
	    (ypts[[lastthird]]-ypts[[firstthird]])/eps[[2]] ],
    Infinity,Infinity];
If[Precision[thiscp] < 1, 1/2, thiscp]
];

closepts[epsIN_,ptlis_]:=Block[{j,numpts,thiscp},
(*
Determine if a set of points is close together.
*)
If[ptlis=={},Return[Infinity]];
If[NumberQ[epsIN],eps=Table[epsIN,{j,1,Length[ptlis[[1]]]}],
	eps=epsIN];
dim=Length[eps];
thiscp = If[Min[eps]>0,
    numpts=Length[ptlis];
    firstfifth=Ceiling[numpts/5];
    lastfifth=Floor[4 numpts/5];  
    For[j=1,j<=dim,++j,
    pt[j]=Sort[Transpose[ptlis][[j]]]
    ];
    Max[Table[(pt[j][[lastfifth]]-pt[j][[firstfifth]])/eps[[j]],
	{j,1,dim}]],
    Infinity,Infinity];
If[Precision[thiscp] < 1, 1/2, thiscp]
];


makegridpts[ll_, ur_, numx_, numy_] := Block[{},
  outlis = {};
  xstep = (ur[[1]] - ll[[1]])/(numx - 1);
  ystep = (ur[[2]] - ll[[2]])/(numy - 1);
  For[j = 1, j <= (numy), ++j,
   For[i = 1, i <= (numx), ++i,
     pt = {ll[[1]] + (i - 1)*xstep, ll[[2]] + (j - 1)*ystep};
     AppendTo[outlis, pt];
     ];
   ];
  outlis
  ];

converteqnsALL[EP_, eqns_, numterms_] := converteqnsALL[EP, eqns, numterms,0];
(* the default is to not normalize the equations *)

converteqnsALL[EP_, eqns_, numterms_,absflag_] := Block[{(*tmpeqns, tmpeqns1, tmpeqns2, tmpeqns3,tmpeqns4*)},
  tmpeqns = eqns;
  tmpeqnsA = (tmpeqns /. killunknownsGE[numterms + 1]);
  tmpeqnsB = normalizeequationsA1[tmpeqnsA,absflag]; 
(*  tmpeqns1 = realize500[tmpeqnsB]; *)
  tmpeqns1R = realize500R[tmpeqnsB];
(*  Print[N[tmpeqns1-tmpeqns1R]]; *)
  tmpeqns2 = (tmpeqns1R/.multsubs[numterms]);  (* replace a_{nm} by a_n a_m for (m,n)=1 *)
  tmpeqns3 = goodprimepowersubs[EP, tmpeqns2, numterms];
  tmpeqns4 = badprimepowersubs[EP, tmpeqns3, numterms];
  tmpeqns4noEPSILON = localsignsubs[EP, tmpeqns4];
  (tmpeqns4noEPSILON/.{bb1[1]->1,bb2[1]->0})
];
(*  tmpeqns2 = makeaALLequations[lev, tmpeqns1]]; *)
(* need to replace makeaALLequations by several steps *)

knownApsubstitutions[ep_, knownAp_, lim_] := Block[{theseunknowns},
  theseunknowns = theunknowns[ep, lim];
(*  Table[theseunknowns[[j]] -> knownAp[[j]], {j, 1, Length[knownAp]}]  *)
  Table[theseunknowns[[j]] -> knownAp[[j]], {j, 1, Length[theseunknowns]}]
  ];

anFromAp[ep_, knownAp_, lim_] := 
 Block[{theseunknowns, allAn}, 
  allAn = Flatten[Table[{bb1[j], bb2[j]}, {j, 1, lim}]];
  allAn = converteqnsALL[ep, allAn, lim];
  theseunknowns = theunknowns[ep, lim];
  allAn /. knownApsubstitutions[ep, knownAp, lim]];

subsAn[anlis_]:=Flatten[Table[{bb1[j]->anlis[[2j-1]],bb2[j]->anlis[[2j]]},{j,1,Length[anlis]/2}]];

substitutionsFromAp[ep_, ap_, numterms_] := Block[{},
  theseAnSubs = subsAn[anFromAp[ep, ap, numterms]];
(* below is bad redundancy from localsignsubs *)
  If[Length[ep[[2]]] > 2,
      If[Length[ep[[2,1]]] ==1 && Length[ep[[2,3]]] == 1, (* simplest case *)
        thebadprime = ep[[2,1,1]];
        twicetheexponent = ep[[2,3,1]];
        (* under what consitions is eps_p = Conj[a_p]/Abs[a_p] when the local factor has degree 1? *)
        theEPSsubs = {EPSILONpR -> thebadprime^(twicetheexponent/2) bb1[thebadprime], EPSILONpI -> -1*thebadprime^(twicetheexponent/2) bb2[thebadprime]};
        theEPSsubs = theEPSsubs/.theseAnSubs;  (* because EPr -> bb1[p] will leave bb1[p] as a parameter *)
        theseAnSubs = Flatten[{theEPSsubs, theseAnSubs}]
      ,
      Print["Error: unimplemented case", ep]
      ]
    ];
  theseAnSubs
]

normalizeequationsA1[eqs_,absflag_]:=Block[{aa},
(*
   Print[Table[NumberQ[Coefficient[eqs[[aa]], bb1[1]]] && Coefficient[eqs[[aa]], bb1[1]] != 0,{aa,1,Length[eqs]}]];
*)
   If[Apply[And, Table[NumberQ[Coefficient[eqs[[aa]], bb1[1]]] && Coefficient[eqs[[aa]], bb1[1]] != 0, {aa,1,Length[eqs]}]],
      Expand[eqs/(
	E^(I Arg[Coefficient[eqs,bb1[1]]]) If[absflag>0,1,
          Sum[Abs[Coefficient[eqs,bb1[j]]]+Abs[Coefficient[eqs,bb2[j]]],{j,1,61}]]) ],
      eqs]]; 

realize500R[xxlis_] :=(*Warning:assumes fewer than 500 coefficients*)
  Block[{j},
  Table[xx = xxlis[[k1]]; tt = Coefficient[xx, bb1[1]];
    realscalefactor0 = (Coefficient[xx, bb1[1]]/.Flatten[Table[{bb1[j] -> 0, bb2[j] -> 0}, {j, 1, 500}]]);
    realscalefactor = realscalefactor0 /. {EPSILONpR -> 0, EPSILONpI -> 0};
    If[realscalefactor == 0,(*Print["zero for realscalefactor"];*)
      realscalefactor = 1];
    yy = Expand[xx/realscalefactor];
    Simplify[ComplexExpand[Re[yy]], 
       Flatten[{{EPSILONpR \[Element] Reals, EPSILONpI \[Element] Reals}, Table[{bb1[p] \[Element] Reals, bb2[p] \[Element] Reals}, {p, 1, 500}]}]],
  {k1, 1, Length[xxlis]}]
];

killunknownsGE[lim_]:=Flatten[Table[{bb1[j]->0,bb2[j]->0},{j,lim,2 lim+50}]];

multsubs[sumlim_]:=Block[{j,n1,n2,n3,n4,ans={},fi},
   If[sumlim>209,Print["Error in sumlim:  products of 4 primes not implemented"]];
	For[j=6,j<=sumlim,++j,
	  fi=FactorInteger[j];
	  If[Length[fi]==1,Next[],
	    If[Length[fi]==2, 
		n1=fi[[1,1]]^fi[[1,2]];
		n2=fi[[2,1]]^fi[[2,2]];
		AppendTo[ans, bb1[n1 n2]->bb1[n1]bb1[n2]-bb2[n1]bb2[n2]];
		AppendTo[ans, bb2[n1 n2]->bb1[n1]bb2[n2]+bb2[n1]bb1[n2]],  (* if length==2 *)
		If[Length[fi]==3,
		n1=fi[[1,1]]^fi[[1,2]];
		n2=fi[[2,1]]^fi[[2,2]];
		n3=fi[[3,1]]^fi[[3,2]];
		AppendTo[ans, bb1[n1 n2 n3]->bb1[n1]bb1[n2]bb1[n3] - 
			bb1[n1]bb2[n2]bb2[n3] - bb2[n1]bb2[n2]bb1[n3] - bb2[n1]bb1[n2]bb2[n3]];
		AppendTo[ans, bb2[n1 n2 n3]-> bb2[n1]bb1[n2]bb1[n3] +
			bb1[n1]bb2[n2]bb1[n3] + bb1[n1]bb1[n2]bb2[n3] -
					bb2[n1]bb2[n2]bb2[n3] ] ] (* if length==3 *)
	]
	] (* if length not eq 1 *)
	];
   ans
];

goodprimepowersubs[EP_, eqns_, numterms_]:=Block[{eqnsTMP,p,j,deg,char,cp,badprimes},
    eqnsTMP=eqns;
    {degree,char}=EP[[1]];
    deg=Abs[degree];
    badprimes=EP[[2,1]];
    j=1;  (* j indexes the primes *)
    While[p=Prime[j]; p<=Sqrt[numterms],  (* check on the Sqrt[] *)
	If[!MemberQ[badprimes,p],
	  cp=char[[Mod[p,Length[char],1] ]];  (* Mod[a,n,1] is in [1,n] *)
          eqnsTMP = (eqnsTMP/.eulersubs[deg,p,Re[cp],Im[cp]])
    	]; (* if not a bad prime *)
	++j
    ];  (* while *)
    If[degree<0,eqnsTMP/.realsubs[numterms],eqnsTMP]
];

realsubs[numterms_]:=Table[bb2[j]->0,{j,1,numterms}];

(* substitutions for the prime powers at the bad primes. Only works for the general bad factor,
meaning that we know an upper bound on its degree and nothing else *)
badprimepowersubs[EP_, eqns_, numterms_]:=Block[{eqnsTMP,p,j,deg,char,cp,badprimes,baddegree},
    eqnsTMP=eqns;
    badprimes=EP[[2,1]];
    absbadprimes=Abs[badprimes];
    baddegs=EP[[2,2]];
    badtypes=Table["",{j,1,Length[badprimes]}];
    (* obsolete, from previous attempts at specifying details about bad primes
    If[Length[EP[[2]] ] ==2, badtypes=Table["",{j,1,Length[badprimes]}],
	badtypes=EP[[2,3]] ];
    *)
    j=1;  (* j indexes the primes *)
    While[p=Prime[j]; p<=numterms, 
        If[MemberQ[badprimes,p],
	    theposition=Position[absbadprimes,p][[1,1]] ;
     	    baddegree= (baddegs[[ theposition ]]);
     	    badtype= (badtypes[[ theposition ]]);
            eqnsTMP = (eqnsTMP/.eulersubsbad[baddegree,p,badtype])
        ]; (* if a bad prime *)
	++j
    ];  (* while *)
    eqnsTMP
];

localsignsubs[EP_, eqns_]:= Block[{eqnsTMP},
    eqnsTMP = eqns;
    If[Length[EP[[2]]] > 2,
      If[Length[EP[[2,1]]] ==1 && Length[EP[[2,3]]] == 1, (* simplest case *)
        If[Length[EP[[2,3,1]]]==1, thelocaldata = {EP[[2,3,1,1]],1}, thelocaldata = EP[[2,3,1]]];
        thebadprime = EP[[2,1,1]];
(*
        twicetheexponent = EP[[2,3,1]];
*)
        twicetheexponent = thelocaldata[[1]];
        extrafactor = thelocaldata[[2]];
        extrafactorR = Re[extrafactor];
        extrafactorI = Im[extrafactor];
        (* eps = extrafactor * conjugate(a_p = bb1[p] + I bb2[p]) , normalized to have abs val 1 *)
        thenewEPSre = Expand[thebadprime^(twicetheexponent/2) (extrafactorR bb1[thebadprime] + extrafactorI bb2[thebadprime])];
        thenewEPSim = Expand[thebadprime^(twicetheexponent/2) (-1*extrafactorR bb2[thebadprime] + extrafactorI bb1[thebadprime])];
        thesubs = {EPSILONpR -> thenewEPSre, EPSILONpI -> thenewEPSim};
        Print["epsilonsubs",thesubs];
        eqnsTMP = eqnsTMP/.thesubs
      , Print["Error: unimplemented case", EP]
      ]
    ];
    eqnsTMP
];

(* general good degree 5 Euler factor *)
eulersubs[5,p_,cpr_,cpi_]:={
bb1[p^3] -> cpr*bb1[p]^2 - bb1[p]^3 - cpr*bb1[p^2] - cpr*bb2[p]^2 + bb1[p]*(2*bb1[p^2] + 
     bb2[p]*(2*cpi + 3*bb2[p])) - cpi*bb2[p^2] - 2*bb2[p]*bb2[p^2], 
 bb2[p^3] -> bb1[p]^2*(cpi - 3*bb2[p]) - bb1[p^2]*(cpi - 2*bb2[p]) - cpi*bb2[p]^2 + bb2[p]^3 + 
   cpr*bb2[p^2] + bb1[p]* (-2*cpr*bb2[p] + 2*bb2[p^2]), 
 bb1[p^4] -> 2*cpr*bb1[p]^3 - bb1[p]^4 + bb1[p^2]^2 - cpi*bb2[p] + bb1[p^2]*(2*cpi - bb2[p])*bb2[p] + 
   2*cpi*bb2[p]^3 - bb2[p]^4 + bb1[p]^2*(bb1[p^2] + 2*bb2[p]* (cpi + 3*bb2[p])) - 
   2*cpr*bb2[p]*bb2[p^2] - bb2[p^2]^2 - bb1[p]*(cpr + 2*cpr*bb1[p^2] - 
     2*cpr*bb2[p]^2 + 2*cpi*bb2[p^2] + 2*bb2[p]*bb2[p^2]), 
 bb2[p^4] -> 2*bb1[p]^3* (cpi - 2*bb2[p]) - 2*cpr*bb2[p]^3 + 2*bb1[p^2]*bb2[p^2] - 
   bb2[p]^2*bb2[p^2] + bb1[p]^2* (-2*cpr*bb2[p] + bb2[p^2]) + bb2[p]*(cpr - 2*cpr*bb1[p^2] - 
     2*cpi*bb2[p^2]) + bb1[p]*(-cpi - 2*bb1[p^2]* (cpi - bb2[p]) + 2*cpi*bb2[p]^2 + 
     4*bb2[p]^3 + 2*cpr*bb2[p^2])
}

(* general good degree 4 Euler factor *)
eulersubs[4,p_,cpr_,cpi_]:={
 bb1[p^3] -> -bb1[p]^3 + bb1[p]*(cpr + 2*bb1[p^2] + 3*bb2[p]^2) + bb2[p]* (cpi - 2*bb2[p^2]), 
 bb2[p^3] -> -3*bb1[p]^2*bb2[p] + bb2[p]*(-cpr + 2*bb1[p^2] + bb2[p]^2) + bb1[p]* (cpi + 2*bb2[p^2]), 
 bb1[p^4] -> -cpr - bb1[p]^4 + bb1[p^2]^2 + 2*cpr*bb2[p]^2 - 
   bb1[p^2]*bb2[p]^2 - bb2[p]^4 + bb1[p]^2*(2*cpr + bb1[p^2] + 
     6*bb2[p]^2) - 2*bb1[p]*bb2[p]* bb2[p^2] - bb2[p^2]^2, 
 bb2[p^4] -> -cpi - 4*bb1[p]^3*bb2[p] + 2*bb1[p]*bb2[p]*(bb1[p^2] + 
     2*bb2[p]^2) + bb2[p]^2* (2*cpi - bb2[p^2]) + 2*bb1[p^2]*
    bb2[p^2] + bb1[p]^2*(2*cpi + bb2[p^2])
};



(* general good degree 3 Euler factor *)
eulersubs[3,p_,cpr_,cpi_]:={bb1[p^2] -> -(cpr*bb1[p]) + bb1[p]^2 - bb2[p]*(cpi + bb2[p]), 
 bb2[p^2] -> -(bb1[p]*(cpi - 2*bb2[p])) + cpr*bb2[p],
 bb1[p^3] -> cpr - 2*cpr*bb1[p]^2 + bb1[p]^3 - 2*cpr*bb2[p]^2 - 3*bb1[p]* bb2[p]^2,
 bb2[p^3] -> cpi - 2*cpi*bb2[p]^2 - bb2[p]^3 + bb1[p]^2*(-2*cpi + 3*bb2[p]),
 bb1[p^4] -> -3*cpr*bb1[p]^3 + bb1[p]^4 + bb1[p]^2*(-cpi^2 + cpr^2 + 3*cpi*bb2[p] - 
     6*bb2[p]^2) + cpr*bb1[p]* (2 + 4*cpi*bb2[p] - 3*bb2[p]^2) + bb2[p]* (-2*cpi + (cpi^2 - cpr^2)*
      bb2[p] + 3*cpi*bb2[p]^2 + bb2[p]^3),
 bb2[p^4] -> cpr*bb1[p]^2*(2*cpi - 3*bb2[p]) + bb1[p]^3*(-3*cpi + 4*bb2[p]) + 
   cpr*bb2[p]*(2 - 2*cpi*bb2[p] - 3*bb2[p]^2) +
	bb1[p]* (2*cpi + 2*(cpi^2 - cpr^2)* bb2[p] - 3*cpi*bb2[p]^2 - 4*bb2[p]^3)}

(* general good degree 2 Euler factor *)
eulersubs[2,p_,cpr_,cpi_]:={bb1[p^2] -> -cpr + bb1[p]^2 - bb2[p]^2,
 bb2[p^2] -> -cpi + 2*bb1[p]*bb2[p], 
 bb1[p^3] -> bb1[p]^3 + 2*cpi*bb2[p] - bb1[p]* (2*cpr + 3*bb2[p]^2), 
 bb2[p^3] -> -2*cpi*bb1[p] - 2*cpr*bb2[p] + 3*bb1[p]^2* bb2[p] - bb2[p]^3, 
 bb1[p^4] -> -cpi^2 + cpr^2 + bb1[p]^4 + 6*cpi*bb1[p]*bb2[p] + 3*cpr*bb2[p]^2 + bb2[p]^4 - 
   	3*bb1[p]^2*(cpr + 2*bb2[p]^2), 
 bb2[p^4] -> 2*cpi*cpr - 3*cpi*bb1[p]^2 - 6*cpr*bb1[p]* bb2[p] + 4*bb1[p]^3*bb2[p] +
	3*cpi*bb2[p]^2 - 4*bb1[p]* bb2[p]^3}

(* special bad degree 3 for Paramodular forms of squarefree level *)
eulersubsbad[-1,p_,{"paramodularp",eps_}]:= {bb1[p^2] -> -1 + eps^2/p + (eps*bb1[p])/Sqrt[p] + bb1[p]^2, 
 bb1[p^3] -> -(eps/Sqrt[p]) + (-2 + (2*eps^2)/p)*bb1[p] + (2*eps*bb1[p]^2)/Sqrt[p] + 
   bb1[p]^3,
 bb1[p^4] -> ((eps^2 - p)^2 + 2*eps*(eps^2 - 2*p)*Sqrt[p]* bb1[p] + (4*eps^2 - 3*p)*p*
     bb1[p]^2 + 3*eps*p^(3/2)* bb1[p]^3 + p^2*bb1[p]^4)/p^2, 
 bb2[p] -> 0,
 bb2[p^2] -> 0,
 bb2[p^3] -> 0, 
 bb2[p^4] -> 0, 
 bb2[p^5] -> 0, 
 bb2[p^6] -> 0, 
 bb2[p^7] -> 0, 
 bb2[p^8] -> 0, 
 bb2[p^9] -> 0, 
 bb2[p^10] -> 0}

(* general bad degree 2 euler factor *)
eulersubsbad[2,p_,""]:={bb1[p^3] -> -bb1[p]^3 + bb1[p]*(2*bb1[p^2] + 3*bb2[p]^2) - 2*bb2[p]* bb2[p^2],
 bb2[p^3] -> -3*bb1[p]^2*bb2[p] + 2*bb1[p^2]*bb2[p] + bb2[p]^3 + 2*bb1[p]*bb2[p^2],
 bb1[p^4] -> -bb1[p]^4 + bb1[p^2]^2 - bb1[p^2]*bb2[p]^2 - bb2[p]^4 + bb1[p]^2*(bb1[p^2] + 
     6*bb2[p]^2) - 2*bb1[p]*bb2[p]* bb2[p^2] - bb2[p^2]^2, 
 bb2[p^4] -> -4*bb1[p]^3*bb2[p] + 2*bb1[p]*bb2[p]*(bb1[p^2] + 2*bb2[p]^2) + bb1[p]^2* bb2[p^2] +
	(2*bb1[p^2] - bb2[p]^2)*bb2[p^2], 
 bb1[p^5] -> -2*bb1[p]^3*bb1[p^2] + 6*bb1[p]^2*bb2[p]*bb2[p^2] -
	2*bb2[p]*(3*bb1[p^2] + bb2[p]^2)* bb2[p^2] + 3*bb1[p]* (bb1[p^2]^2 +
	2*bb1[p^2]* bb2[p]^2 - bb2[p^2]^2), 
 bb2[p^5] -> -6*bb1[p]^2*bb1[p^2]* bb2[p] - 2*bb1[p]^3*bb2[p^2] +
	6*bb1[p]*(bb1[p^2] + bb2[p]^2)* bb2[p^2] + bb2[p]* (3*bb1[p^2]^2 + 2*bb1[p^2]* bb2[p]^2 -
	3*bb2[p^2]^2), 
 bb1[p^6] -> bb1[p]^6 + bb1[p^2]^3 - 3*bb1[p^2]^2*bb2[p]^2 - bb2[p]^6 - bb1[p]^4*(4*bb1[p^2] + 
     15*bb2[p]^2) + 16*bb1[p]^3* bb2[p]*bb2[p^2] - 4*bb1[p]* bb2[p]*(3*bb1[p^2] + 4*bb2[p]^2)*
    bb2[p^2] + 3*bb2[p]^2* bb2[p^2]^2 + 3*bb1[p]^2* (bb1[p^2]^2 + 8*bb1[p^2]* bb2[p]^2 + 5*bb2[p]^4 - 
     bb2[p^2]^2) - bb1[p^2]* (4*bb2[p]^4 + 3*bb2[p^2]^2), 
 bb2[p^6] -> 6*bb1[p]^5*bb2[p] - 4*bb1[p]^3*(4*bb1[p^2]*bb2[p] + 
     5*bb2[p]^3) - 4*bb1[p]^4* bb2[p^2] + 6*bb1[p]^2* (bb1[p^2] + 4*bb2[p]^2)*
    bb2[p^2] + 2*bb1[p]*bb2[p]* (3*bb1[p^2]^2 + 8*bb1[p^2]* bb2[p]^2 + 3*bb2[p]^4 - 
     3*bb2[p^2]^2) - bb2[p^2]* (-3*bb1[p^2]^2 + 6*bb1[p^2]* bb2[p]^2 + 4*bb2[p]^4 + bb2[p^2]^2),
 bb1[p^7]-> bb1[p]^7 - bb1[p]^5*(2*bb1[p^2] + 21*bb2[p]^2) + 10*bb1[p]^4* bb2[p]*bb2[p^2] + 4*bb1[p]^2*
    bb2[p]*(3*bb1[p^2] - 5*bb2[p]^2)* bb2[p^2] + 2*bb2[p]*bb2[p^2]* (-6*bb1[p^2]^2 - 2*bb1[p^2]*
      bb2[p]^2 + bb2[p]^4 + 2*bb2[p^2]^2) + bb1[p]^3* (-2*bb1[p^2]^2 + 20*bb1[p^2]*
      bb2[p]^2 + 35*bb2[p]^4 + 2*bb2[p^2]^2) + 
   bb1[p]*(4*bb1[p^2]^3 + 6*bb1[p^2]^2*bb2[p]^2 - 7*bb2[p]^6 - 6*bb2[p]^2* bb2[p^2]^2 - 2*bb1[p^2]*
      (5*bb2[p]^4 + 6*bb2[p^2]^2)), 
 bb2[p^7] -> 7*bb1[p]^6*bb2[p] - 5*bb1[p]^4*(2*bb1[p^2]*bb2[p] + 7*bb2[p]^3) - 2*bb1[p]^5*
    bb2[p^2] - 4*bb1[p]^3* (bb1[p^2] - 5*bb2[p]^2)* bb2[p^2] + 2*bb1[p]*bb2[p^2]*
    (6*bb1[p^2]^2 + 6*bb1[p^2]* bb2[p]^2 - 5*bb2[p]^4 - 2*bb2[p^2]^2) + bb1[p]^2*bb2[p]*
    (-6*bb1[p^2]^2 + 20*bb1[p^2]* bb2[p]^2 + 21*bb2[p]^4 + 6*bb2[p^2]^2) - 
   bb2[p]*(-4*bb1[p^2]^3 - 2*bb1[p^2]^2*bb2[p]^2 + bb2[p]^6 + 2*bb2[p]^2*
      bb2[p^2]^2 + 2*bb1[p^2]* (bb2[p]^4 + 6*bb2[p^2]^2)), 
 bb1[p^8] -> 3*bb1[p]^6*bb1[p^2] + bb1[p^2]^4 - 6*bb1[p^2]^3* bb2[p]^2 - 18*bb1[p]^5*bb2[p]*
    bb2[p^2] + 12*bb1[p]^3*bb2[p]* (6*bb1[p^2] + 5*bb2[p]^2)* bb2[p^2] + 9*bb2[p]^4*
    bb2[p^2]^2 + bb2[p^2]^4 - 6*bb1[p]*bb2[p]*bb2[p^2]* (6*bb1[p^2]^2 + 12*bb1[p^2]*
      bb2[p]^2 + 3*bb2[p]^4 - 2*bb2[p^2]^2) - 9*bb1[p]^4* (bb1[p^2]^2 + 5*bb1[p^2]*
      bb2[p]^2 - bb2[p^2]^2) - 3*bb1[p^2]^2*(3*bb2[p]^4 + 2*bb2[p^2]^2) - 3*bb1[p^2]*
    (bb2[p]^6 - 6*bb2[p]^2* bb2[p^2]^2) + 3*bb1[p]^2* (2*bb1[p^2]^3 + 18*bb1[p^2]^2*
      bb2[p]^2 - 18*bb2[p]^2* bb2[p^2]^2 + 3*bb1[p^2]* (5*bb2[p]^4 - 2*bb2[p^2]^2)), 
 bb2[p^8] -> 18*bb1[p]^5*bb1[p^2]* bb2[p] + 3*bb1[p]^6*bb2[p^2] - 9*bb1[p]^4*(2*bb1[p^2] + 
     5*bb2[p]^2)*bb2[p^2] - 12*bb1[p]^3*bb2[p]* (3*bb1[p^2]^2 + 5*bb1[p^2]* bb2[p]^2 - 3*bb2[p^2]^2) + 
   3*bb1[p]^2*bb2[p^2]* (6*bb1[p^2]^2 + 36*bb1[p^2]* bb2[p]^2 + 15*bb2[p]^4 - 
     2*bb2[p^2]^2) + 6*bb1[p]*bb2[p]* (2*bb1[p^2]^3 + 6*bb1[p^2]^2* bb2[p]^2 - 6*bb2[p]^2*
      bb2[p^2]^2 + 3*bb1[p^2]* (bb2[p]^4 - 2*bb2[p^2]^2)) - bb2[p^2]*(-4*bb1[p^2]^3 + 
     18*bb1[p^2]^2*bb2[p]^2 + 3*bb2[p]^2*(bb2[p]^4 - 2*bb2[p^2]^2) + 2*bb1[p^2]*
      (9*bb2[p]^4 + 2*bb2[p^2]^2)), 
 bb1[p^9] -> -bb1[p]^9 + 6*bb1[p]^7*(bb1[p^2] + 6*bb2[p]^2) - 42*bb1[p]^6* bb2[p]*bb2[p^2] + 30*bb1[p]^4*
    bb2[p]*(3*bb1[p^2] + 7*bb2[p]^2)* bb2[p^2] - 18*bb1[p]^2*bb2[p]^3* (10*bb1[p^2] + 7*bb2[p]^2)*
    bb2[p^2] + 6*bb1[p]^3*bb2[p]^2* (15*bb1[p^2]^2 + 35*bb1[p^2]* bb2[p]^2 + 14*bb2[p]^4 - 
     15*bb2[p^2]^2) - 9*bb1[p]^5* (bb1[p^2]^2 + 14*bb1[p^2]* bb2[p]^2 + 14*bb2[p]^4 - 
     bb2[p^2]^2) + bb1[p]* (5*bb1[p^2]^4 - 42*bb1[p^2]* bb2[p]^6 - 9*bb2[p]^8 + 
     45*bb2[p]^4*bb2[p^2]^2 + 5*bb2[p^2]^4 - 15*bb1[p^2]^2* (3*bb2[p]^4 + 2*bb2[p^2]^2)) + 
   2*bb2[p]*bb2[p^2]* (-10*bb1[p^2]^3 + 3*bb2[p]^6 + bb1[p^2]*(9*bb2[p]^4 + 10*bb2[p^2]^2)), 
 bb2[p^9] -> -9*bb1[p]^8*bb2[p] + 42*bb1[p]^6*bb2[p]*(bb1[p^2] + 2*bb2[p]^2) + 6*bb1[p]^7*
    bb2[p^2] - 18*bb1[p]^5* (bb1[p^2] + 7*bb2[p]^2)* bb2[p^2] + 30*bb1[p]^3*bb2[p]^2*
    (6*bb1[p^2] + 7*bb2[p]^2)* bb2[p^2] - 3*bb1[p]^4*bb2[p]* (15*bb1[p^2]^2 + 70*bb1[p^2]*
      bb2[p]^2 + 42*bb2[p]^4 - 15*bb2[p^2]^2) + 18*bb1[p]^2* bb2[p]^3*(5*bb1[p^2]^2 + 
     7*bb1[p^2]*bb2[p]^2 + 2*bb2[p]^4 - 5*bb2[p^2]^2) - 2*bb1[p]*bb2[p^2]* (-10*bb1[p^2]^3 + 21*bb2[p]^6 + 
     5*bb1[p^2]*(9*bb2[p]^4 + 2*bb2[p^2]^2)) - bb2[p]*(-5*bb1[p^2]^4 + 6*bb1[p^2]*bb2[p]^6 + 
     bb2[p]^8 - 9*bb2[p]^4* bb2[p^2]^2 - 5*bb2[p^2]^4 + bb1[p^2]^2*(9*bb2[p]^4 + 30*bb2[p^2]^2)), 
 bb1[p^10] -> -bb1[p]^10 + bb1[p^2]^5 - 10*bb1[p^2]^4* bb2[p]^2 + 3*bb1[p]^8* (bb1[p^2] + 15*bb2[p]^2) - 
   24*bb1[p]^7*bb2[p]*bb2[p^2] + 12*bb1[p]^5*bb2[p]*(-3*bb1[p^2] + 14*bb2[p]^2)*bb2[p^2] + 
   12*bb1[p]^3*bb2[p]*bb2[p^2]* (15*bb1[p^2]^2 + 10*bb1[p^2]* bb2[p]^2 - 14*bb2[p]^4 - 
     5*bb2[p^2]^2) + 3*bb1[p]^6* (bb1[p^2]^2 - 28*bb1[p^2]* bb2[p]^2 - 70*bb2[p]^4 - 
     bb2[p^2]^2) - 5*bb1[p^2]^3* (3*bb2[p]^4 + 2*bb2[p^2]^2) - 3*bb1[p^2]^2*(bb2[p]^6 - 
     20*bb2[p]^2*bb2[p^2]^2) + bb2[p]^2*(bb2[p]^8 + 3*bb2[p]^4* bb2[p^2]^2 - 10*bb2[p^2]^4) + 
   bb1[p^2]*(3*bb2[p]^8 + 45*bb2[p]^4*bb2[p^2]^2 + 5*bb2[p^2]^4) - 15*bb1[p]^4*
    (bb1[p^2]^3 + 3*bb1[p^2]^2* bb2[p]^2 - 14*bb2[p]^6 - 3*bb2[p]^2*bb2[p^2]^2 - 
     bb1[p^2]*(14*bb2[p]^4 + 3*bb2[p^2]^2)) + 4*bb1[p]*bb2[p]*bb2[p^2]* (-20*bb1[p^2]^3 - 45*bb1[p^2]^2*
      bb2[p]^2 + 6*bb2[p]^6 + 15*bb2[p]^2*bb2[p^2]^2 + bb1[p^2]*(-9*bb2[p]^4 + 
       20*bb2[p^2]^2)) + bb1[p]^2*(10*bb1[p^2]^4 + 90*bb1[p^2]^3*bb2[p]^2 + 15*bb1[p^2]^2*(3*bb2[p]^4 - 
       4*bb2[p^2]^2) - 6*bb1[p^2]* (14*bb2[p]^6 + 45*bb2[p]^2* bb2[p^2]^2) - 
     5*(9*bb2[p]^8 + 9*bb2[p]^4* bb2[p^2]^2 - 2*bb2[p^2]^4)), 
 bb2[p^10] -> -10*bb1[p]^9*bb2[p] + 24*bb1[p]^7*bb2[p]*(bb1[p^2] + 5*bb2[p]^2) + 3*bb1[p]^8*
    bb2[p^2] + 6*bb1[p]^6* (bb1[p^2] - 14*bb2[p]^2)* bb2[p^2] + 15*bb1[p]^4*bb2[p^2]*
    (-3*bb1[p^2]^2 - 6*bb1[p^2]* bb2[p]^2 + 14*bb2[p]^4 + bb2[p^2]^2) - 6*bb1[p]^5*bb2[p]*
    (-3*bb1[p^2]^2 + 28*bb1[p^2]* bb2[p]^2 + 42*bb2[p]^4 + 3*bb2[p^2]^2) + 2*bb1[p]^2*
    bb2[p^2]*(20*bb1[p^2]^3 + 135*bb1[p^2]^2*bb2[p]^2 - 42*bb2[p]^6 - 45*bb2[p]^2*
      bb2[p^2]^2 + 5*bb1[p^2]* (9*bb2[p]^4 - 4*bb2[p^2]^2)) + 12*bb1[p]^3*bb2[p]*
    (-5*bb1[p^2]^3 - 5*bb1[p^2]^2* bb2[p]^2 + 5*bb2[p]^2* (2*bb2[p]^4 + bb2[p^2]^2) + 
     bb1[p^2]*(14*bb2[p]^4 + 15*bb2[p^2]^2)) + 2*bb1[p]*bb2[p]*(10*bb1[p^2]^4 + 
     30*bb1[p^2]^3*bb2[p]^2 - 5*bb2[p]^8 - 9*bb2[p]^4* bb2[p^2]^2 + 10*bb2[p^2]^4 + 
     bb1[p^2]^2*(9*bb2[p]^4 - 60*bb2[p^2]^2) - 6*bb1[p^2]* (2*bb2[p]^6 + 15*bb2[p]^2*
        bb2[p^2]^2)) + bb2[p^2]*(5*bb1[p^2]^4 - 40*bb1[p^2]^3*bb2[p]^2 + 
     3*bb2[p]^8 + 15*bb2[p]^4* bb2[p^2]^2 + bb2[p^2]^4 - 5*bb1[p^2]^2*(9*bb2[p]^4 + 
       2*bb2[p^2]^2) + bb1[p^2]* (-6*bb2[p]^6 + 40*bb2[p]^2* bb2[p^2]^2))
}

(* general bad degree 1 euler factor *)
eulersubsbad[1,p_,""]:={bb1[p^2] -> bb1[p]^2 - bb2[p]^2, 
 bb2[p^2] -> 2*bb1[p]*bb2[p], 
 bb1[p^3] -> bb1[p]^3 - 3*bb1[p]*bb2[p]^2, 
 bb2[p^3] -> 3*bb1[p]^2*bb2[p] - bb2[p]^3, 
 bb1[p^4] -> bb1[p]^4 - 6*bb1[p]^2*bb2[p]^2 + bb2[p]^4, 
 bb2[p^4] -> 4*bb1[p]*bb2[p]*(bb1[p]^2 - bb2[p]^2), 
 bb1[p^5] -> bb1[p]^5 - 10*bb1[p]^3*bb2[p]^2 + 5*bb1[p]*bb2[p]^4, 
 bb2[p^5] -> 5*bb1[p]^4*bb2[p] - 10*bb1[p]^2*bb2[p]^3 + bb2[p]^5, 
 bb1[p^6] -> bb1[p]^6 - 15*bb1[p]^4*bb2[p]^2 + 15*bb1[p]^2*bb2[p]^4 - bb2[p]^6, 
 bb2[p^6] -> 6*bb1[p]^5*bb2[p] - 20*bb1[p]^3*bb2[p]^3 + 6*bb1[p]*bb2[p]^5, 
 bb1[p^7] -> bb1[p]^7 - 21*bb1[p]^5*bb2[p]^2 + 35*bb1[p]^3*bb2[p]^4 - 7*bb1[p]*bb2[p]^6, 
 bb2[p^7] -> 7*bb1[p]^6*bb2[p] - 35*bb1[p]^4*bb2[p]^3 + 21*bb1[p]^2*bb2[p]^5 - bb2[p]^7, 
 bb1[p^8] -> bb1[p]^8 - 28*bb1[p]^6*bb2[p]^2 + 70*bb1[p]^4*bb2[p]^4 - 28*bb1[p]^2*bb2[p]^6 + bb2[p]^8, 
 bb2[p^8] -> 8*bb1[p]*bb2[p]* (bb1[p]^6 - 7*bb1[p]^4*bb2[p]^2 + 7*bb1[p]^2*bb2[p]^4 - bb2[p]^6), 
 bb1[p^9] -> bb1[p]^9 - 36*bb1[p]^7*bb2[p]^2 + 126*bb1[p]^5*bb2[p]^4 - 84*bb1[p]^3*bb2[p]^6 + 
   9*bb1[p]*bb2[p]^8, 
 bb2[p^9] -> 9*bb1[p]^8*bb2[p] - 84*bb1[p]^6*bb2[p]^3 + 126*bb1[p]^4*bb2[p]^5 - 
   36*bb1[p]^2*bb2[p]^7 + bb2[p]^9, 
 bb1[p^10] -> bb1[p]^10 - 45*bb1[p]^8*bb2[p]^2 + 210*bb1[p]^6*bb2[p]^4 - 210*bb1[p]^4*bb2[p]^6 + 
   45*bb1[p]^2*bb2[p]^8 - bb2[p]^10, 
 bb2[p^10] -> 2*(5*bb1[p]^9*bb2[p] - 60*bb1[p]^7*bb2[p]^3 + 126*bb1[p]^5*bb2[p]^5 - 
    60*bb1[p]^3*bb2[p]^7 + 5*bb1[p]*bb2[p]^9)}

(* missing Euler factor *)
eulersubsbad[0,p_,""]:=Flaten[Table[{bb1[p^j] -> 0, bb2[p^j] -> 0},{j,1,10}]];

(* next step:  use  multsubs[sumlim] followed by substitutions for the good factors,
and then the bad factors *)

scalelist[numsolve_]:=Join[{0,1/100,1/100,1/100,1/10,1/10,1/10,1/10,1/4,1/4,1/4,1/4,1/10,1/10,1/10,1/10},
	Flatten[Table[{1,4,10},{j,1,numsolve/3}]],{100,100,100,100,100,100}];

scalelist[numsolve_]:=Join[{0,1/10000000,1/10000000,1/10000000,1/1000000,1/1000000,1/10000,1/10000},
	Flatten[Table[{1/1000,1/100,1/100},{j,1,numsolve/3}]],{1/10,1/10,1/10}];

scalelist[numsolve_]:={0,1/10000000,1,1,1};
	(* to speed up the zooming case *)

scalelist[numsolve_]:=Join[{0,1/10000000,1/100000,1/10000},
	Flatten[Table[{1/1000,1/100,1/10,1},{j,1,numsolve/5}]],Flatten[Table[{1/2,1/2,1,1,1,2,2,2,5,5,5,10,20,20},{j,1,numsolve}]]];

findsolone[eqns_, theunkns_, testunksIN_,numsolve_,eps_,{numcheck_,eps2_}] := Block[
{ct, foundsols={},residuals={},tmp(*,startvals1*)},
If[Length[eqns]<Length[testunksIN[[1]]],Print["not enough equations"];Return[{}]];
testunks={};
theeqns=Take[eqns,Length[testunksIN[[1]]]];
eqnPRECISION=Floor[Precision[theeqns/.Table[theunkns[[j]] -> 0, {j, 1, Length[theunkns]}]]];
For[jf=1,jf<=Length[testunksIN],++jf,
AppendTo[testunks,
Table[{testunksIN[[jf,j,1]],SetPrecision[testunksIN[[jf,j,2]],100]},{j,1,Length[testunksIN[[jf]]]}]]];

(* eqns : the equations
   theunkns     : the unknowns in the equations
   testunksIN     : starting values for (some of) the unknowns
   numsolve     : the number of times to solve the equations with each starting value
   eps  : threshold for considering the equations to be solved
   numcheck     : number of unknows we check to see if two solutions are the same
   eps2 : how close the numcheck unknowns have to be to be the same
*)

   thescalelist = Flatten[Union[{0},
               Table[1/10^j, {j, 0, 6}, {k, 1, 5}], {{1, 1, 1, 1, 3, 3, 3, 3, 3, 6,6,6,6,6,10,10,10,10,10,10}}]];
   For[nn=1,nn<=Length[testunks],++nn,
     ct=0;
     While[ct<Length[thescalelist],++ct;
        startvals1={};
        For[j=1,j<=Length[testunks[[nn]]],++j,
           AppendTo[startvals1,
             {testunks[[nn]][[j,1]],tmp=testunks[[nn]][[j,2]]+thescalelist[[ct]] RandomReal[{-1, 1}, WorkingPrecision -> 100]; tmp, tmp + 1/1000}]];
        For[j=Length[testunks[[nn]]]+1,j<=Length[theunkns],++j,
           AppendTo[startvals1,
         (*    {theunkns[[j]],tmp=thescalelist[[ct]] RandomReal[{-1, 1}, WorkingPrecision -> 100]; tmp, 101/100 tmp}]]; *)
             (* previous version was unused (because we pre-pad with 0s) but woudl have been wrong if tmp == 0 *)
             {theunkns[[j]],tmp=thescalelist[[ct]] RandomReal[{-1, 1}, WorkingPrecision -> 100]; tmp, 1/100 + tmp}]];
        test = FindRoot[theeqns, startvals1, Method -> "Secant",
                WorkingPrecision->(*-5+*)eqnPRECISION, MaxIterations->200];
        vec = theeqns /. test;
        If[Norm[vec]<eps,
              AppendTo[foundsols,test];
              AppendTo[residuals,Norm[vec]];
              Print["found a solution on iteration ",ct];
              Break[]  (* only want to find one solution *)
        ];
     ];
   ];
   foundsols (*  not caring about repeats =tossrepeats[foundsols,residuals,theunkns,numcheck,eps2] *)
];

findsolmult[eqns_, theunkns_, testunksIN_,numsolve_,eps_,{numcheck_,eps2_}] := Block[
{ct, foundsols={},residuals={},tmp(*,startvals1*)},
If[Length[eqns]<Length[testunksIN[[1]]],Print["not enough equations"];Return[{}]];
testunks={};
theeqns=Take[eqns,Length[testunksIN[[1]]]];
eqnPRECISION=Floor[Precision[theeqns/.Table[theunkns[[j]] -> 0, {j, 1, Length[theunkns]}]]];
For[jf=1,jf<=Length[testunksIN],++jf,
AppendTo[testunks,
Table[{testunksIN[[jf,j,1]],SetPrecision[testunksIN[[jf,j,2]],100]},{j,1,Length[testunksIN[[jf]]]}]]];

(* eqns : the equations
   theunkns     : the unknowns in the equations
   testunksIN     : starting values for (some of) the unknowns
   numsolve     : the number of times to solve the equations with each starting value
   eps  : threshold for considering the equations to be solved
   numcheck     : number of unknows we check to see if two solutions are the same
   eps2 : how close the numcheck unknowns have to be to be the same
*)
   thescalelist=scalelist[numsolve];
   For[nn=1,nn<=Length[testunks],++nn,
     ct=0;
     While[ct<Length[thescalelist],++ct;
        startvals1={};
        For[j=1,j<=Length[testunks[[nn]]],++j,
           originfactor = 1;  (* determines whether we use the input value, its negative, or 0 as the offset of the random part of the initial value *)
           If[ct>5, originfactor = 0];
           AppendTo[startvals1,
             {testunks[[nn]][[j,1]],tmp=originfactor testunks[[nn]][[j,2]]+thescalelist[[ct]] RandomReal[{-1, 1}, WorkingPrecision -> 100]; tmp, tmp + 1/1000}]];
        For[j=Length[testunks[[nn]]]+1,j<=Length[theunkns],++j,
           AppendTo[startvals1,
             {theunkns[[j]],tmp=thescalelist[[ct]] RandomReal[{-1, 1}, WorkingPrecision -> 100]; tmp, 101/100 tmp}]];
(*
Print["findroot with startvals1",N[startvals1]];
*)
        test = FindRoot[theeqns, startvals1, Method -> "Secant",
		WorkingPrecision->(*-5+*)eqnPRECISION];
        vec = theeqns /. test;
        If[Norm[vec]<eps,AppendTo[foundsols,test];AppendTo[residuals,Norm[vec]]];
     ];
   ];
   foundsols=tossrepeats[foundsols,residuals,theunkns,numcheck,eps2]
];

findclosest[target_,thelist_]:=thelist[[1]];

(* ???? *)
(*]*)

tossrepeats[sublists_,residuals_,allunknowns_,numunkn_,eps_]:=Block[{j,k,weededlist={},weededresiduals={}},
	checkunks=Take[allunknowns,numunkn];
	If[Length[sublists]<1, Return[{}]];
	weededlist={sublists[[1]] };
	weededresiduals={residuals[[1]] };
	For[j=2,j<=Length[sublists],++j,
	  foundmatch=0;
	  For[k=1,k<=Length[weededlist],++k,
		If[Norm[(checkunks/.sublists[[j]])-(checkunks/.weededlist[[k]])]<eps,
		  foundmatch=1;
		  If[residuals[[j]]<weededresiduals[[k]],
			weededresiduals[[k]]=residuals[[j]];
			weededlist[[k]]=sublists[[j]] ]];	
	  ];
	  If[foundmatch==0,
		AppendTo[weededresiduals,residuals[[j]]];
		AppendTo[weededlist,sublists[[j]]]];
	  ];
   weededlist
];
		

twodsecant[pt1_,pt2_,pt3_,pt4_]:=Block[{ans,data1,data2,x,y,plane1,plane2},  (* need the first segment to be horizontal? *)
data1={Flatten[{pt1[[1]],pt1[[2,1]]}], Flatten[{pt2[[1]],pt2[[2,1]]}],
        Flatten[{pt3[[1]],pt3[[2,1]]}], Flatten[{pt4[[1]],pt4[[2,1]]}]};
data2={Flatten[{pt1[[1]],pt1[[2,2]]}], Flatten[{pt2[[1]],pt2[[2,2]]}],
        Flatten[{pt3[[1]],pt3[[2,2]]}], Flatten[{pt4[[1]],pt4[[2,2]]}]};
plane1 = Fit[data1, {1, x, y}, {x, y}];
plane2 =  Fit[data2, {1, x, y}, {x, y}];
ans=Solve[{plane1==0,plane2==0},{x,y}];
({x,y}/.ans)[[1]]
];

nDsecant[pts_,indicators_]:=Block[{j,m,k,n,x,plane,data,vars,ans}, (* find common zero of indicators
				at k points in R^n *)
   k=Length[pts];
   n=Length[pts[[1]]];
   For[j=1,j<=k,++j,If[Length[pts[[j]]]!=n || Length[indicators[[j]]]!=n,Print["error, not all points in same dimension: ",j];
					Abort[]]];
   vars=Table[x[j],{j,1,n}];
   For[j=1,j<=n,++j,
	data[j]=Table[Append[pts[[m]],indicators[[m,j]]],{m,1,k}]
   ];
   For[j=1,j<=n,++j,
	plane[j]=Fit[data[j],Prepend[vars,1],vars]
   ];
   ans=Solve[Table[plane[j]==0,{j,1,n}],vars];
   If[Length[ans[[1]]] < Length[vars],
       Print["problem with equations in nDsecant. pts:", pts, "indicators", indicators, "the eqn",Table[plane[j]==0,{j,1,n}]];
       Table["",{j,1,n}]
   ,
       (vars/.ans)[[1]]
   ]
];
   

theunknowns[ep_,lim_]:=Block[{j,fi,ans={},deg,badps,badpunks,theprime,theexponent},
	ans={};
	deg=ep[[1,1]];  (* degree of the good local factors *)
	absdeg=Abs[deg];  (* degree could be negative, to indicate real coefficients  *)
	badps=ep[[2,1]];  (* bad primes *)
	badpunks=ep[[2,2]];  (* number of unknowns for each bad prime *)
	absbadpunks=Abs[ep[[2,2]]];  (* number is negative for real coefficients *)
  For[j=2,j<=lim,++j,
    fi=FactorInteger[j];
    If[Length[fi]==1,  (* ie, it is a prime power *)
	{theprime,theexponent}=fi[[1]];
	If[MemberQ[badps,theprime],  (* bad prime case *)
	   theposition=Position[badps,theprime][[1,1]];
	   If[theexponent<= (absbadpunks[[ theposition ]]),
		If[ (badpunks[[ theposition ]])>0,
			AppendTo[ans,{bb1[j],bb2[j]}], AppendTo[ans,bb1[j]]]
	     ],  (*else, good prime case *)
	   If[theexponent<=Floor[absdeg/2] || theexponent >=5, 
	     If[deg>0, AppendTo[ans,{bb1[j],bb2[j]}], AppendTo[ans,bb1[j]] ]
	     ]
	]  (* if which prime case *)
     ]  (* if we are at a prime power *)
    ];  (*for *)
  If[Length[ep[[2]]] > 2, (* case of known shape of bad factor *)
    If[Length[ep[[2,1]]] > 1 || ep[[2,2]] != {1}, Print["Error: unimplemented case"]];
  ];
Flatten[ans]
];
	

UNKNOWNSIGNtheunknowns[ep_,lim_]:=Block[{j,fi,ans={},deg,badps,badpunks,theprime,theexponent},
        ans={};
        deg=ep[[1,1]];  (* degree of the good local factors *)
        absdeg=Abs[deg];  (* degree could be negative, to indicate real coefficients  *)
        badps=ep[[2,1]];  (* bad primes *)
        badpunks=ep[[2,2]];  (* number of unknowns for each bad prime *)
        absbadpunks=Abs[ep[[2,2]]];  (* number is negative for real coefficients *)
  For[j=2,j<=lim,++j,
    fi=FactorInteger[j];
    If[Length[fi]==1,  (* ie, it is a prime power *)
        {theprime,theexponent}=fi[[1]];
        If[MemberQ[badps,theprime],  (* bad prime case *)
           theposition=Position[badps,theprime][[1,1]];
           If[theexponent<= (absbadpunks[[ theposition ]])|| theexponent >=11,
                If[ (badpunks[[ theposition ]])>0,
                        AppendTo[ans,{bb1[j],bb2[j]}], AppendTo[ans,bb1[j]]]
             ],  (*else, good prime case *)
           If[theexponent<=Floor[absdeg/2] || theexponent >=5,
             If[deg>0, AppendTo[ans,{bb1[j],bb2[j]}], AppendTo[ans,bb1[j]] ]
             ]
        ]  (* if which prime case *)
     ]  (* if we are at a prime power *)
    ];  (*for *)
ans = Flatten[ans];
  If[ep[[1,1]]==3 && Length[ep[[2,1]]] == 1,
    thebadprime = ep[[2,1,1]];
    If[ep[[2,2,1]] != ep[[1,1]]-1,
       Print["Error: only prime level case implemented"],
       (* so bb1[p] and bb2[p] are not unknowns, but theta and phi are *)
       ans = ans/.{bb1[thebadprime]->theta, bb2[thebadprime]->phi};
       ans = ans/.{bb1[thebadprime^2] -> Nothing, bb2[thebadprime^2] -> Nothing}
(* did not work
       ans = DeleteCases[ans, (# == bb1[thebadprime^2] || # == bb2[thebadprime^2]) &];
*)
    ]
  ];
  If[ep[[1,1]]==2 && Length[ep[[2,1]]] == 1,
    thebadprime = ep[[2,1,1]];
    If[ep[[2,2,1]] != ep[[1,1]]-1,
       Print["Error: only prime level case implemented"],
       (* so bb1[p] and bb2[p] are not unknowns, but theta and phi are *)
       ans = ans/.{bb1[thebadprime]->theta};
       ans = DeleteElement[ans,{bb2[thebadprime]}];
    ]
  ];
ans
];

fefromdata[dat_] := Block[{j},
   dat[[1, 2]]/.Table[XX[j] -> dat[[1, 1, j]], {j, 1, Length[dat[[1, 1]]]}]
(*
   If[Length[dat[[1,1]]] == 2,
      (dat[[1, 2]] /. {XX[1] -> dat[[1, 1, 1]], XX[2] -> dat[[1, 1, 2]]}),
      (dat[[1, 2]] /. {XX[1] -> dat[[1, 1, 1]], XX[2] -> dat[[1, 1, 2]], XX[3] -> dat[[1, 1, 3]]})]
*)
];

(* need to consolidate the Z, L, Lambda below *)
evaluateZfromAp[FE_, b_, s_, Ev_, gflag_, PRECIS_, ep_, ap_] :=
  Block[{numterms},
   numterms = If[ep[[1, 1]] > 0, NextPrime[Length[ap]] - 1, NextPrime[Length[ap]/2] - 1];
   thisEv = {Ev[[1]], numterms};
   rawZ = Z[FE, b, s, thisEv, gflag, PRECIS];
   rawZ/.substitutionsFromAp[ep, ap, numterms]
(*
   theseAnSubs = subsAn[anFromAp[ep, ap, numterms]];
   rawZ /. theseAnSubs
*)
];

evaluateLfromAp[FE_, b_, s_, Ev_, gflag_, PRECIS_, ep_, ap_] :=
  Block[{numterms},
   If[ep[[1, 1]] > 0, numterms = NextPrime[Length[ap]] - 1,
    numterms = NextPrime[Length[ap]/2] - 1];
   thisEv = {Ev[[1]], numterms};
   rawL = L[FE, b, s, thisEv, gflag, PRECIS];
   rawL/.substitutionsFromAp[ep, ap, numterms]
(*
   theseAnSubs = subsAn[anFromAp[ep, ap, numterms]];
   rawL /. theseAnSubs
*)
];

evaluateLambdafromAp[FE_, b_, s_, Ev_, gflag_, PRECIS_, ep_, ap_] :=
  Block[{numterms},
   If[ep[[1, 1]] > 0, numterms = NextPrime[Length[ap]] - 1,
    numterms = NextPrime[Length[ap]/2] - 1];
   thisEv = {Ev[[1]], numterms};
   rawLam = Lambda[FE, b, s, thisEv, gflag, PRECIS];
   rawLam/.substitutionsFromAp[ep, ap, numterms]
(*
   theseAnSubs = subsAn[anFromAp[ep, ap, numterms]];
   rawLam /. theseAnSubs
*)];

evaluateZfromApEXTRA[FE_, b_, s_, Ev_, gflag_, PRECIS_, ep_, ap_, extra_] := Block[{numterms},
   If[ep[[1, 1]] > 0, numterms = NextPrime[Length[ap]] - 1,
    numterms = NextPrime[Length[ap]/2] - 1];
   thisEv = {Ev[[1]], numterms + extra};
   rawZ = Z[FE, b, s, thisEv, gflag, PRECIS];
   rawZ/.substitutionsFromAp[ep, ap, numterms]
];

quickZ[ell_, b_, t_] := 
  evaluateZfromAp[fefromdata[ell], b, 1/2 + t I, {2, 0}, 1, 40, 
   ell[[1, 3]], ell[[1, 4]]];

quickZEXTRA[ell_, b_, t_, extra_] := Block[{},
  Expand[evaluateZfromApEXTRA[fefromdata[ell], b, 1/2 + t I, {2, 0}, 
    0, 40, ell[[1, 3]], ell[[1, 4]], extra]]];

boundtailsimple[obj_, lim_, deg_] := Block[{theerr},
  theerr = 0;
  For[j = 1, j <= lim, ++j, 
   theerr += 
    deg (Abs[Coefficient[obj, bb1[j]]] + Abs[Coefficient[obj, bb2[j]]])
   ];
  theerr
  ];

