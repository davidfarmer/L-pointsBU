
(*
4a:  new relativeerror0 by */100 instead of powers of 10.  10/7/19
4b: use refineandsave as the basis for testandsave, which is
    intended for finding new solutions, not for refining approximate
    solutions.
    Also, improve searchonce so that once it finds a solution at one point,
    use it as the starting value at the next point.   10/17/19
4c: At each test point, onlu use the first solution found
    This involves not sorting the results from findsolmult;
    Test whether the results are worth saving (based on star epsilon
    and clustering of most recent iteration).  Also save estimates of
    errors.  Eliminate foundsofar.
4d: Mode parameter: "poking" initially, then "zooming".
    If zooming and rel err is > 1, increase number of coefficients.
    If poking and the point does not move, also increase number of coefficients.
4e: Also save initguessIN to output file.
    Modify findgsandnumterms to force at least 4 coefficient unknowns  11/4/19
    In testandsave, prevent poking from moving more than wanderlimit=1/2 from the starting point. 11/6/19
    Added note about error on line 364.  Fix in next version
4f: Fixed error that put a wandering value in the wrong place
4g: Not save if not zooming; make restartct limit a parameter;
    and if zooming started on last step, then do one more (via NUMSTEPS+1). 12/10/19
4h: option to always save with categories: GOOD, BAD, and UGLY (i.e., Don't Know)
4i: when zooming, continue if relativeerror < 0.7 (i.e., 1/sqrt[2] improvement each iteration;
      otherwise, increase number of coefficients used if relativeerror < 0.7, not relativeerror < 1,
      as in version 4h  (1/27/20)
4j: add to output: no. of iterations
4k: changed itemtosave format (4/20/20)
4l: minor changes to output, add "warning" to filename when parameters have changed significantly
4m: saveBad flag to allow initial searching without creating lots of BAD files
4n: better tracking of wandering after zooming

4p: Have the epsilon for determining if a solution is valid, depend on the precision.
    findbadfactoreqns[EP]  for making the bad satake parameters the correct size.
    (need to also do the cases where we know the exact bad factor)

4pa: added "scaled cluster size:" in front of the output scaled cluster size (so that it
     is easy to see when perusing the output)  (4/21/23)
4q: When creating output data, don't take more coefficients then are available
    (Was causing an error, possibly only for BAD data)  (4/25/23)
    Initialize tmpY so that the error is calculated correctly when zooming
    finishes on the first iteration.  (4/26/23)
4r: offset star points (which appear to be a parameter but actually are
       set within the code).

4s: altered findstartingvalues, testandsave, and searchonce to use new FE notation.
    deleted many unused functions.
    lowered lowest zoom factor from 1/100 to 1/500
    rounded new starepsX to (2 digits)/100...001  at each step
    changed filenameX to filename.x for X = GOOD, BAD, UGLY

4t: changed file extension to ".m"
    added tossRepeats function.

    supplemented degree2 equations to include |eps| = 1

    Terminate search if wandering
    renamed "finddegree2equations" to "finddegree2andSIGNequations"

5a: major change in how FE signs are handled
      - the sign is always an unknown (even when actually known)
      - for every bad prime, eps_p is an unknown

    changed itemtosave format (added the unknowns as an entry)

    tossRepeats uses appropriate (custom) sorting
*)

debugging1 = False;  (* omits extra equaitons when degree is 2 *)

saveBad = True;  (* save runs that fail to converge.  Useful when refining, wasteful for initial search. *)

(* useWarnings should be set in the perl script *)
(* useWarnings = True;  put warning on filename when parameters have changed from initial values *)

testandsave[initguessIN_ (* the initial guess *),
        starpts_ (* configuration of nearby points to test *),
        starepsIN_ (* how far to move the starpts *),
        numtermsIN_ (* number of terms or truncation error *),
        gANDsList_ (* list of {g1, g2, s} to make equations or detectors, sorted *),
        FEin_ (* functional equation, for parts that are not being sesarched *),
        EvIN_  (* evaluation parameters *),
        coeffstartIN_  (* starting values of the coefficients *),
        gflag_ (* usually use gflag=1 for the Stefan version of g *),
        PRECISIN_ (* precision to use *),
        absflag_ (* 1 to use absolute errors in the truncation error, 0 for relative *),
        EP_ (* describes the shape of the Euler product *),
        type_,
        NUMSTEPS_,
        RESTARTMAX_,  (* new in 4g; formerly hard coded as 6 *)
        TARGETERR_
        ]:= Block[{PRECIS},
  Print["starting testandsave"];
  mode = "poking";
  PRECIS=PRECISIN;
  myct=1;
  restartct=0;
  starepsX=starepsIN;
  myguessX=initguessIN;
  numtermsX=numtermsIN;
  coeffstartX=coeffstartIN;
  FEinitial=FEnewtoold[FEin];
  (* tmpY = {}; *)
  (* in the case where the zooming finishes in one step, we need to compare
     with the input values, not the previous iteration.  So we make a fake
     tmpY (which will immediately be set to prevtmpY).  It does not have the
     correct length, but all that are needed are the 1st and 5th entries,
     which are the spectral parameters and the coefficients, respectively. *)
  tmpY = {initguessIN, 1, 1, {starepsIN, numtermsIN, PRECISIN}, coeffstartIN};
  prevtmpY = tmpY;

  wanderlimit = 1/2;  (* prevent poking from going too far *)

  (* main loop:  how many times to zoom *)
  While[myct<=NUMSTEPS + 1 && restartct <= RESTARTMAX && starepsX>TARGETERR,
    prevtmpY = tmpY;
    isWandering = False;

    Print["boxsize,truncation: ", N[{starepsX,numtermsX}]];

    tmpY=searchonce[myguessX,
            starpts (* configuration of nearby points to test *),
            starepsX (* how far to move the starpts *),
            numtermsX (* number of terms or truncation error *),
            gANDsList (* list of {g1, g2, s} to make equations or detectors, sorted *),
            FEin (* functional equation, for parts that are not being sesarched *),
            EvIN  (* evaluation parameters *),
            coeffstartX  (* starting values of the coefficients *),
            gflag (* usually use gflag=1 for the Stefan version of g *),
            PRECIS (* precision to use *),
            absflag (* 1 to use absolute errors in the truncation error, 0 for relative *),
            EP (* describes the shape of the Euler product *)];

    If[Length[tmpY]<4, restartct = restartct+1; Print["found nothing, try again"]; Continue[]];

    relativeerror=tmpY[[2]];

    For[jo=1, jo <= Length[initguessIN], ++jo,  (* check wandering in each coordinate separately *)
        If[Abs[tmpY[[1,jo]] - initguessIN[[jo]] ] > wanderlimit,
            Print["Wandering, so adjusting coordinate", jo];
            isWandering = True;
            If[mode == "zooming",
                Print["Wandering after zooming"];
                tmpY[[1,jo]] = initguessIN[[jo]] +  starepsX RandomReal[{-1,1}],
                tmpY[[1,jo]] = initguessIN[[jo]] + (wanderlimit/2) RandomReal[{-1,1}]
            ];
            If[mode == "zooming",
                tmpY[[5]] = coeffstartX,  (* use previous (while zooming) coefficients *)
                tmpY[[5]] = coeffstartIN   (* don't use the wandering coefficients *)
            ];
            relativeerror = 10 (* to prevent zooming, because we moved to a new random place *)
        ]
    ];

    If[isWandering, Print["wandering, to exiting this attempt"]; Break[]];
(*
    relativeerror0=10^(-1*Max[0,Floor[Log[10,1/relativeerror]]]);
*)
    relativeerror0 = Min[1, Ceiling[500 relativeerror]/500];

    Print["starepsX was", starepsX, "relativeerror is", relativeerror];
    starepsX = starepsX*relativeerror0;
    numtermsX = numtermsX*relativeerror0;
    Print["starepsX was", starepsX];
    starepsX = roundfraction[starepsX];
    Print["starepsX is", starepsX];
    numtermsX = roundfraction[numtermsX];

    coeffstartX=tmpY[[5]];  (* probably shodl be set before the wandering check, *)
    myguessX=tmpY[[1]];      (* and then these can be set directly if wandering *)

    If[mode == "zooming" && relativeerror > 0.7,  (* probably spinning, so use more coefficients *)
        Print["Seems to be spinning wheels after zooming.  Using more terms."];
        numtermsX = numtermsX/10];

    If[mode == "poking" && Length[prevtmpY] > 3 && Length[tmpY] > 3 && (* relativeerror > 1 && *)
              Abs[tmpY[[1,1]]-prevtmpY[[1,1]]] < 0.5 starepsX,
              If[Length[tmpY[[1]]] == 1 || Abs[tmpY[[1,2]]-prevtmpY[[1,2]]] < 0.5 starepsX,
                 Print["Seems to be spinning wheels before zooming.  Using more terms."];
                 numtermsX = numtermsX/10,
              Print["not really spinning"]]
     ];

    If[relativeerror < 1, mode = "zooming"; Print["in zooming mode"]];

    (* new in 4g:  skip last step if not zooming *)
    If[myct == NUMSTEPS && mode != "zooming",
       Break[]];

    relativeprecision=tmpY[[3]];
    If[relativeprecision < 10,
	Print["precision was", PRECIS];
        PRECIS=PRECIS+10,
        If[relativeprecision < 15, PRECIS = PRECIS+5]
    ];
    Print["Round ",myct," done, stareps: ", starepsX];
    Print[relativeerror, "   ",relativeerror0,"   "];
 (*   Print["RESULTS SO FAR: ", tmpY];  *)
    Print["RESULTS SO FAR: ", Take[tmpY,5]];
    ++myct;
    If[hasbeenfound[tmpY[[1]],0.01],
        Print["Has already found ",tmpY[[1]] ]; Return[]
    ];
  ];  (* While *)

  If[Length[tmpY]<4, Print["FAILED, ",restartct]; Return[]];

  Print["Time to save it."];
  myfilename=mydatadir<>"/"<>type<>"x";
  myfilename = myfilename<>numbertostring[myguessX[[1]],5];
  For[jh=2,jh<=Length[myguessX],++jh,
     myfilename = myfilename<>"_";
     myfilename = myfilename<>numbertostring[myguessX[[jh]],5]
  ];
(*  myfilename = myfilename<>"_";  *)
  If[starepsX <= TARGETERR,myfilename = myfilename<>".good"; Print["Saving as GOOD"],
    If[mode == "zooming", myfilename = myfilename<>".ugly"; Print["Saving as UGLY"],
      myfilename = myfilename<>".bad";
      If[saveBad, Print["Saving as BAD"],
        Print["Bad, but not saving"]; Return[tmpY]]
    ]
  ];

  If[useWarnings,
    If[Max[N[Abs[tmpY[[1]] - initguessIN]]] > 0.01,
      myfilename = myfilename <> "_eigen_warning"; Print["Eigenvalue may have wandered."]
    ];
    If[Max[Abs[N[Take[tmpY[[5]], 4] -  Take[coeffstartIN, 4]]]] > 0.2,
      myfilename = myfilename <> "_coeff_warning"; Print["Coefficients have changed."]
    ];
  ];

  numknowncoeffs = Length[prevtmpY[[5]]];
  Print["numknowncoeffs", numknowncoeffs, "which starts", prevtmpY[[5,1]],"compared to", Length[tmpY[[5]]], "which startts", tmpY[[5,1]]];

  itemtosave = {
    masterversion (* string designated in perl script to identify the symmetry type being investigated *),
    versionnumber (* version of the Mathematica code combination, set in perl script *),
    {type, versioneqns} (* string designated in perl script, version of equations, set in perl script *),
    {
      tmpY[[1]] (* final value of search parameter(s) *),
      FEin (* functional equation data, incluing XX[[n]] to be replaced by search parameters *),
      EP (* Euler product data *),
      tmpY[[5]] (* final values of coefficients *),
      tmpY[[8]] (* the unknowns *),
      masterversion
(*
      If[Length[tmpY[[5]]] >= numknowncoeffs,
         Take[tmpY[[5]], numknowncoeffs],
         tmpY[[5]] 
      ]
*)
    },
    {
      N[Abs[(tmpY[[1]]-prevtmpY[[1]])]] (* upper bound on uncertainty in final search parameters;
                                         determined by difference b/w last & next-to-last search parameters *),
      If[Length[tmpY[[5]]] >= numknowncoeffs,
         N[Abs[Take[tmpY[[5]], numknowncoeffs] - prevtmpY[[5]]]], (* upper bound on uncertainty in each coefficient *)
         N[Abs[   tmpY[[5]] - Take[ prevtmpY[[5]],Length[tmpY[[5]]] ]     ]]
      ]
    },
    {
      { (* intial search parameters *)
        {
          boxsizeIN (* initial size of search box *),
          truncationerrorIN (* selected size of truncation error used to set number of coefficients *)
        },
        {
          NUMSTEPS (* maximum number of iterations *),
          RESTARTMAX (* number of tries to find a solution *)
        },
        PRECISIN (* initial working precision *),
        TARGETERR (* goal for the final box size *)
      },
      { (* final search parameters *)
        {
          tmpY[[4,1]] (* final size of search box *),
          tmpY[[4,2]] (* final size of truncation error used to set number of coefficients *)
        },
        myct (* number of iterations that happened *),
        { tmpY[[4,3]], N[tmpY[[3]]] } (* {final working precision, remaining extra precision} *),
        "scaled cluster size:",
        N[tmpY[[2]]] (* scaled cluster size--scaled relative to serach box *)
      }
    },
    {
      N[Abs[tmpY[[1]] - initguessIN]] (* change in the search value, initial to final *),
      N[Take[tmpY[[5]], 4] -  Take[coeffstartIN, 4]] (* change in the first four coefficient unknowns, initial to final *)
    },
    {
      {sourcedatafile, candidatelistname, currentindex},  (* input file of candidates; name of list; and entry number *)
      initguessIN (* initial value of search parameter *),
      coeffstartIN (* inital values of coefficients *)
    }
  };  (* itemtosave *)


  Save[myfilename,itemtosave];

  Print["saved to ", myfilename];

  tmpY
];


numbertostring[num_, digits_] := Block[{ans},
  thenum = num;
  ans = "";
  If[thenum < 0, ans = ans <> "-"; thenum = Abs[thenum]];
  ans = ans <> IntegerString[Floor[thenum]];
  ans = ans <> ".";
  fracpart=IntegerString[Floor[10^digits FractionalPart[thenum]]];
  While[StringLength[fracpart]<digits, fracpart = "0"<>fracpart];
  ans = ans <> fracpart;
  ans]

(* we make a new searchzoom program called "searchonce" which is used to do
one step refinment of the functional equation parameters *)

(* this required making a new version of makeequations, called
makeequationsNEW  in gl3_4f.txt*)

(* the major change is that the "current guess", which is a list
of parameters, say {x_1,...,x_N}, are now interpreted as the
values of XX[1],...,XX[N] in the functional equation.  Thus,
the functional equation must be input with XX[1],...,XX[N]
as part of the format, to be substituted later. *)

(* 12 *)
searchonce[initguessIN_ (* the initial guess *),
        starpts_ (* configuration of nearby points to test *),
        starepsIN_ (* how far to move the starpts *),
        numtermsIN_ (* number of terms or truncation error *),
        gANDsList_ (* list of {g1, g2, s} to make equations or detectors, sorted *),
        FEin_ (* functional equation, for parts that are not being sesarched *),
        EvIN_  (* evaluation parameters *),
        coeffstart_  (* starting values of the coefficients *),
        gflag_ (* usually use gflag=1 for the Stefan version of g *),
        PRECIS_ (* precision to use *),
        absflag_ (* 1 to use absolute errors in the truncation error, 0 for relative *),
        EP_ (* describes the shape of the Euler product *)]:=
(* 11 *)
searchonce[initguessIN,
       (* starpts, omitted: automatic now --configuration of nearby points to test *)
        starepsIN,
        numtermsIN,
        gANDsList,
        FEin,
        EvIN,
        coeffstart,
        gflag,
        PRECIS,
        absflag,
        EP];

thestarpts[1] = {{-11/10}, {9/10}};

thestarpts[2] = {{-11/10, -21/20}, {9/10, -21/20}, {9/10, 19/20}, {-11/10, 19/20}};

(*
thestarpts[3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
      {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}};

thestarpts[3] = {{0,0,0},{1, -1, -1/2}, {-1, 1, -1/2}, {-1, -1, 1/2}, {1, 1, 1/2}};
*)
thestarpts[3] = {{9/10, -21/20, -31/30}, {-11/10, 19/20, -31/30}, {-11/10, -21/20, 29/30}, {9/10, 19/20, 29/30}};

thestarpts[4] = {{9/10, -21/20, -31/30,-41/40}, {-11/10, 19/20, -31/30,-42/41}, {-11/10, -21/20, 29/30, -16/15}, {-13/12, -14/13,-20/10, 12/13}, {9/10, 19/20, 29/30, 26/25}};

findstartingvalues[initguessIN_ (* the initial guess *),
	starepsIN_ (* how far to move the starpts *),
	numtermsIN_ (* number of terms or truncation error *),
	gANDsList_ (* list of {g1, g2, s} to make equations or detectors, sorted *),
	FEin_ (* functional equation, for parts that are not being sesarched *),
	EvIN_  (* evaluation parameters *),
	coeffstart_  (* starting values of the coefficients *),
	gflag_ (* usually use gflag=1 for the Stefan version of g *),
	PRECIS_ (* precision to use *),
	absflag_ (* 1 to use absolute errors in the truncation error, 0 for relative *),
	EP_ (* describes the shape of the Euler product *)]:=Block[
              {jz,kz,numterms,FE,numtermsTAB,Ev},
    starpts=thestarpts[Length[initguessIN]];
    DETECTPRECIS=PRECIS;
    detectedpoints={};
    hitpoints={};
    closepoints={};
    startcoeffs=coeffstart;
    starepssize=Floor[Log[10,starepsIN]];
    stareps=Ceiling[starepsIN/10^(starepssize-1)] 10^(starepssize-1);
    Print["Size of zooming neighborhood: ",stareps];
    initguess = Floor[10 initguessIN/stareps] stareps/10;
    Print["Initial guess: ",initguess, " which is approximately ",N[initguess,16]];
    Print["Initial coefficients: ", N[startcoeffs]];
    (* FE[[2]]=pairtotriple[initguess]; *)
    Print["Functional equation ", FEin];
    FE = FEnewtoold[FEin];
(*
    FE = (FEin/.Table[XX[j]->initguess[[j]],{j,1,Length[initguess]}]);
*)
    FE = (FE/.Table[XX[j]->initguess[[j]],{j,1,Length[initguess]}]);

    lev=FEin[[3]];
    nu=EvIN[[1]];
    istep=stepsizeRM[nu, PRECIS];
    ev={nu,istep};

    Print["finding starting values near ",N[initguess, PRECIS], " to ",stareps, ", level ",lev];

(* find how many terms are needed to make the equations *)

    {unknowns,gsCount,numterms,numdetectors} = findgsandnumterms[numtermsIN,FE,gANDsList,ev,gflag,PRECIS,
      absflag,EP,starpts];

    If[gsCount ==0 || gsCount<Length[unknowns]+numdetectors,Print["Error:  not enough detectors"];
Return[]];

(* now gsCount is the number of {g1,g2,s} pairs we need, numterms is the
number of terms we need in the Dirichlet series, and unknowns is the
list of unknowns.

We need to create the list of {g1,g2} and s's for making the equations,
and the same for the detectors.  We use 8 detectors.

*)

    degree2eqns= finddegree2eqns[EP,unknowns];
    signeqns= findsigneqns[FE,EP];
    signunknowns= findSIGNunknowns[FE, EP];

    unknowns = Flatten[{signunknowns, unknowns}];

(*
    badfactoreqns = findbadfactoreqns[EP];
    badfactorsubstitutions = findbadfactorsubstitutions[EP];
*)

    Ev={nu,numterms};

    eqngs={};
    detectgs={};
    detectmod=3;
    gs8=1+Floor[gsCount/8];
    For[j=1,j<=gsCount,++j,
      If[(Mod[j,gs8]==detectmod && Length[detectgs]<numdetectors) || Length[eqngs]>= Length[unknowns],
          AppendTo[detectgs,gANDsList[[gsCount-j+1]]],
          AppendTo[eqngs,gANDsList[[gsCount-j+1]]]
      ]
    ];

(* now turn detectgs and eqngs into separate {g1,g2} and s lists*)

    gtab=Table[{eqngs[[j,1]],eqngs[[j,2]]},{j,1,Length[eqngs]}];
    stab=Table[eqngs[[j,3]],{j,1,Length[eqngs]}];
    detectg=Table[{detectgs[[j,1]],detectgs[[j,2]]},{j,1,Length[detectgs]}];
    detects=Table[detectgs[[j,3]],{j,1,Length[detectgs]}];

    While[Length[unknowns]>Length[startcoeffs],
        AppendTo[startcoeffs,0]];  (* pad starting values with 0 if necessary*)

    detectpts={};

    For[starN=1,starN<=1,++starN,  (* startN indexes the starpts *)
        thept=initguess+stareps starpts[[starN]];
        AppendTo[detectpts,thept];

Print["testing ",starN,", :",thept];

        eq[starN] = makeequationsNEW[FEin, thept, gtab, stab, Ev,gflag,PRECIS];

        eq[starN] = Flatten[{degree2eqns, signeqns, eq[starN]}];

(*
        If[Length[degree2andSIGNeqns]>0 && Not[debugging1],
            Print["adding degree 2 equaitons"];
            eq[starN] = Flatten[{degree2andSIGNeqns,eq[starN]}];
        ];
*)

        eqsolv[starN] = converteqnsALL[EP, eq[starN], numterms, absflag];

tmpeqsolv = eqsolv[starN];

  (* should the badfactorsubstitutions be called from converteqnsALL? *)
  (*
        If[Length[badfactorsubstitutions]>0,
            Print["adding badfactorsubstitutions", badfactorsubstitutions];
            eqsolv[starN] = Expand[eqsolv[starN]/.badfactorsubstitutions];
        ];
   *)

        startvals={Table[{unknowns[[jz]],startcoeffs[[jz]]},{jz,1,Length[unknowns]}]};
        targeteps = 10.0^(-DETECTPRECIS/2);
(*
        targeteps = 10.0^(-DETECTPRECIS/3);
*)
        ans[starN]= findsolmult[eqsolv[starN], unknowns, startvals, 100,targeteps,{4,0.1}];
        Print["First up to 5 initial answers",If[Length[ans[starN]]>5,Take[ans[starN],5], ans[starN]]];

(*
        If[ans[starN]=={},
            Print["No solution at point ",starN,". Trying again."];
            ans[starN]= findsolmult[eqsolv[starN], unknowns, startvals,50,targeteps/10,{4,0.1}];
Print["all these starting values ",ans[starN]]
        ];
*)

        If[ans[starN]=={},
           Print["No solution at point ",starN,". Stopping."];
           Return[{}]
        ];
     ];
     unknowns/.ans[1]

];

(* 11 *)
searchonce[initguessIN_ (* the initial guess *),
	starepsIN_ (* how far to move the starpts *),
	numtermsIN_ (* number of terms or truncation error *),
	gANDsList_ (* list of {g1, g2, s} to make equations or detectors, sorted *),
	FEin_ (* functional equation, for parts that are not being sesarched *),
	EvIN_  (* evaluation parameters *),
	coeffstart_  (* starting values of the coefficients *),
	gflag_ (* usually use gflag=1 for the Stefan version of g *),
	PRECIS_ (* precision to use *),
	absflag_ (* 1 to use absolute errors in the truncation error, 0 for relative *),
	EP_ (* describes the shape of the Euler product *)]:=Block[{jz,kz,numterms,FE,numtermsTAB,Ev},
    starpts=thestarpts[Length[initguessIN]];
    DETECTPRECIS=PRECIS;
    detectedpoints={};
    hitpoints={};
    closepoints={};
    startcoeffs=coeffstart;
    starepssize=Floor[Log[10,starepsIN]];
    stareps=Ceiling[starepsIN/10^(starepssize-1)] 10^(starepssize-1);
    Print["Size of zooming neighborhood: ",stareps];
    initguess = Floor[10 initguessIN/stareps] stareps/10;
    Print["Initial guess: ",initguess, " which is approximately ",N[initguess,16]];
    Print["Initial coefficients: ", N[startcoeffs]];
    (* FE[[2]]=pairtotriple[initguess]; *)
    FE = FEnewtoold[FEin];
(*
    FE = (FEin/.Table[XX[j]->initguess[[j]],{j,1,Length[initguess]}]);
*)
    FE = (FE/.Table[XX[j]->initguess[[j]],{j,1,Length[initguess]}]);

    lev=FEin[[3]];
    nu=EvIN[[1]];
    istep=stepsizeRM[nu, PRECIS];
    ev={nu,istep};

    Print["testing near ",N[initguess, PRECIS], " to ",stareps, ", level ",lev];

(* find how many terms are needed to make the equations *)

    {unknowns,gsCount,numterms,numdetectors} = findgsandnumterms[numtermsIN,FE,gANDsList,ev,gflag,PRECIS,
      absflag,EP,starpts];

    If[gsCount ==0 || gsCount<Length[unknowns]+numdetectors,Print["Error:  not enough detectors"];
Return[]];

(* now gsCount is the number of {g1,g2,s} pairs we need, numterms is the
number of terms we need in the Dirichlet series, and unknowns is the
list of unknowns.

We need to create the list of {g1,g2} and s's for making the equations,
and the same for the detectors.  We use 8 detectors.

*)

    degree2eqns= finddegree2eqns[EP,unknowns];
    signeqns= findsigneqns[FE,EP];
    signunknowns= findSIGNunknowns[FE, EP];

    unknowns = Flatten[{signunknowns, unknowns}];

(*
    badfactoreqns = findbadfactoreqns[EP];
    badfactorsubstitutions = findbadfactorsubstitutions[EP];
*)

    Ev={nu,numterms};

    eqngs={};
    detectgs={};
    detectmod=3;
    gs8=1+Floor[gsCount/8];
    For[j=1,j<=gsCount,++j,
      If[(Mod[j,gs8]==detectmod && Length[detectgs]<numdetectors) || Length[eqngs]>= Length[unknowns],
          AppendTo[detectgs,gANDsList[[gsCount-j+1]]],
          AppendTo[eqngs,gANDsList[[gsCount-j+1]]]
      ]
    ];

(* now turn detectgs and eqngs into separate {g1,g2} and s lists*)

    gtab=Table[{eqngs[[j,1]],eqngs[[j,2]]},{j,1,Length[eqngs]}];
    stab=Table[eqngs[[j,3]],{j,1,Length[eqngs]}];
    detectg=Table[{detectgs[[j,1]],detectgs[[j,2]]},{j,1,Length[detectgs]}];
    detects=Table[detectgs[[j,3]],{j,1,Length[detectgs]}];

    While[Length[unknowns]>Length[startcoeffs],
        AppendTo[startcoeffs,0]];  (* pad starting values with 0 if necessary*)

    detectpts={};

    For[starN=1,starN<=Length[starpts],++starN,  (* startN indexes the starpts *)
        thept=initguess+stareps starpts[[starN]];
        AppendTo[detectpts,thept];

Print["testing ",starN,", :",thept];

(*
Print["eqution at",InputForm[{FEin, thept, gtab, stab, Ev,gflag,PRECIS}]];
*)

        eq[starN] = makeequationsNEW[FEin, thept, gtab, stab, Ev,gflag,PRECIS];

        eq[starN] = Flatten[{degree2eqns, signeqns, eq[starN]}];

(*
Print["is", InputForm[eq[starN]]];
*)
(*
        If[Length[degree2andSIGNeqns]>0 && Not[debugging1],
            Print["adding degree 2 andSIGNequaitons"];
            eq[starN] = Flatten[{degree2andSIGNeqns,eq[starN]}];
        ];
*)

        eqsolv[starN] = converteqnsALL[EP, eq[starN], numterms, absflag];

(*
        If[Length[badfactorsubstitutions]>0,
            Print["adding badfactorsubstitutions", badfactorsubstitutions];
            eqsolv[starN] = Expand[eqsolv[starN]/.badfactorsubstitutions];
        ];
*)

 (*       Print["checking on badfactorsubstitutions", N[eq[starN]]]; *)

(*
        Print["rechecking on badfactorsubstitutions", N[eqsolv[starN]]];
*)
        startvals={Table[{unknowns[[jz]],startcoeffs[[jz]]},{jz,1,Length[unknowns]}]};
        targeteps = 10.0^(-DETECTPRECIS/2);
        ans[starN]= findsolone[eqsolv[starN], unknowns, startvals,12,targeteps,{4,0.1}];
      (*  ans[starN]=Sort[ans[starN]];
       As per 4c *)
        Print[ans[starN]];

(*
        If[ans[starN]=={},
            Print["No solution at point ",starN,". Trying again."];
            targeteps = 10.0^(-DETECTPRECIS/2);
            ans[starN]= findsolone[eqsolv[starN], unknowns, startvals,5,targeteps,{4,0.1}];
Print[ans[starN]]
        ];
*)

        If[ans[starN]=={},
           Print["No solution at point ",starN,". Stopping."];
           Return[]
        ];

        detecteq[starN]=makeequationsNEW[FEin, thept, detectg, detects, Ev,gflag,DETECTPRECIS];
        detecteq[starN]=converteqnsALL[EP, detecteq[starN],numterms,1];

(*
        If[Length[badfactorsubstitutions]>0,
            Print["adding badfactorsubstitutions", badfactorsubstitutions];
            detecteq[starN] = Expand[detecteq[starN]/.badfactorsubstitutions];
        ];
*)

        (* we want to use the value we just found as the starting value
           for the next point *)
        tmpcoeffsubs = ans[starN][[1]];
        startcoeffs = Table[tmpcoeffsubs[[j, 2]], {j, 1, Length[tmpcoeffsubs]}]

    ]; (* For starN *)

    If[numdetectors==0,
        ans[1] (* just stop if only one point is being searched *)
      ,
        detectans=Table[ans[starN],{starN,1,Length[detectpts]}];
        thedetectors=Table[detecteq[starN],{starN,1,Length[detectpts]}];
           Print["about to check4pts", detectpts];
        detectANS=check4pts[detectpts,detectans,thedetectors];
        If[detectANS == "", Print["failure to detect"]; Return[""]];
        cp=closepts[stareps, detectANS[[1]]];
        themedian=Median[detectANS[[1]]];
           Print[cp,"  cp and median  ",{cp, themedian}];
        newcoeffs=multipleinterppairs[detectpts, detectans, themedian, 4];
           Print["newcoeffs", N[newcoeffs]];
           Print["detectANS", N[detectANS]];

        {themedian, cp, Precision[{themedian, cp, detectANS}],
            {starepsIN,numtermsIN,PRECIS},
            Transpose[newcoeffs[[1]]][[2]], detectANS[[2]], detectANS[[1]],
            unknowns}
    ]
];

printflag=0;


Ndigits[x_,n_]:=If[NumberQ[x],Floor[x]+N[1+x-Floor[x],n+1]-1,Print["Error in Ndigits"]; 1234];;

matchpair[xIN_,yIN_]:=Block[{len,diff,diffexp,diffpower},
len=Min[Length[xIN],Length[yIN]];
x=Take[xIN,len];
y=Take[yIN,len];
diff=Abs[x-y];
diffexp=Log[10,diff];
diffpower=Ceiling[diffexp];
(*
Print[N[Floor[y/10^diffpower]10^diffpower,1-diffpower]];
*)
If[NumberQ[x],
N[Floor[x/10^diffpower]10^diffpower,-diffpower],
(*
Table[N[Floor[x[[j]]/10^diffpower[[j]]]10^diffpower[[j]],-diffpower[[j]]],{j,1,Length[x]}]
*)
Table[Ndigits[x[[j]],-diffpower[[j]]],{j,1,Length[x]}]
]
];
(*
Floor[x/10^diffpower]10^diffpower];
*)

makeequationsR[FE_, glis_, svals_,Ev_,gflag_,PRECIS_] :=
Block[{v,w,j,k,sol,bvals,eqns,numeqns},
  FEtmp=FE;
  (*FEtmp[[2]]=eis; *)
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

boundtailR[obj_, deg_, limsum_] := Block[{}, subs0 = {};
  For[j = 1, j <= limsum, ++j, AppendTo[subs0, bb1[j] -> 0]; AppendTo[subs0,
bb2[j] -> 0]];
  For[j = 1, j <= Length[obj], ++j,
   ans1[j] = Re[(obj[[j, 2]] /. subs0)];
   errtot[j] = 0;
   For[k = 1, k <= limsum, ++k,
    errtot[j] +=
     Abs[Coefficient[obj[[j, 2]], bb1[k]] ] RamaBound[k, deg];
    errtot[j] +=
     Abs[Coefficient[obj[[j, 2]], bb2[k]] ] RamaBound[k, deg]
    ]];
  Table[{obj[[j, 1]], ans1[j] + errtot[j] err}, {j, 1, Length[obj]}]]

RamaBound[1, deg_] := 1

RamaBound[j_, deg_] := Block[{fi, m, k}, fi = FactorInteger[j];
  Product[m = fi[[k]];
   Binomial[m[[2]] + deg - 1, deg - 1], {k, 1, Length[fi]}]]

(* build in more restrictive vakues for the sharp Ramanujan bound for elliptic curves *)
boundtailEC[obj_, deg_, limsum_] := Block[{}, subs0 = {};
  For[j = 1, j <= limsum, ++j, AppendTo[subs0, bb1[j] -> 0]; AppendTo[subs0,
bb2[j] -> 0]];
  For[j = 1, j <= Length[obj], ++j,
   ans1[j] = Re[(obj[[j, 2]] /. subs0)];
   errtot[j] = 0;
   For[k = 1, k <= limsum, ++k,
    errtot[j] +=
     Abs[Coefficient[obj[[j, 2]], bb1[k]] ] ECBound[k, deg];
    errtot[j] +=
     Abs[Coefficient[obj[[j, 2]], bb2[k]] ] ECBound[k, deg]
    ]];
  Table[{obj[[j, 1]], ans1[j] + errtot[j] err}, {j, 1, Length[obj]}]]

boundtailECc[obj_, deg_, limsum_] := Block[{}, subs0 = {};
  For[j = 1, j <= limsum, ++j, AppendTo[subs0, bb1[j] -> 0]; AppendTo[subs0,
bb2[j] -> 0]];
  For[j = 1, j <= Length[obj], ++j,
   ans1[j] = (obj[[j, 2]] /. subs0);
   errtotR[j] = 0;
   errtotI[j] = 0;
   For[k = 1, k <= limsum, ++k,
    errtotR[j] +=
     Abs[Re[Coefficient[obj[[j, 2]], bb1[k]] ]] ECBound[k, deg];
    errtotI[j] +=
     Abs[Im[Coefficient[obj[[j, 2]], bb1[k]] ]] ECBound[k, deg];
    errtotR[j] +=
     Abs[Re[Coefficient[obj[[j, 2]], bb2[k]] ]] ECBound[k, deg];
    errtotI[j] +=
     Abs[Im[Coefficient[obj[[j, 2]], bb2[k]] ]] ECBound[k, deg];
    ]];
  Table[{obj[[j, 1]], ans1[j] + errtotR[j] err + I errtotI[j] err}, {j, 1, Length[obj]}]]






ECBound[1, deg_] := 1

ECBound[2,2] = Sqrt[2]
ECBound[3,2] = Sqrt[3]
ECBound[4,2] = 3
ECBound[5,2] = 4/Sqrt[5]
ECBound[6,2] = Sqrt[6]
ECBound[7,2] = 5/Sqrt[7]
ECBound[8,2] = 3/(2 Sqrt[2])
ECBound[9,2] = 2
ECBound[10,2] = 4 Sqrt[2]/Sqrt[5]
ECBound[j_, deg_] := Block[{fi, m, k}, fi = FactorInteger[j];
  Product[m = fi[[k]];
   Binomial[m[[2]] + deg - 1, deg - 1], {k, 1, Length[fi]}]]

coefficientsubs[EP_,solvedcoeff_,numterms_]:=Block[{(*aa,*)j,unknowns},
aa[1]=1;
deg=EP[[1,1]];
unknowns=theunknowns[EP,numterms];
charR=Re[EP[[1,-1,1]]];
charI=Im[EP[[1,-1,1]]];
unknownssubs=Table[unknowns[[j]]->solvedcoeff[[j]],{j,1,Length[unknowns]}];
For[n=2,n<=numterms,++n,
facn=FactorInteger[n];
(* the following is wrong: it assumes all primes are good *)
aa[n]=Product[p=facn[[j,1]];pn=facn[[j,1]]^facn[[j,2]];
   ((bb1[pn] + I bb2[pn])/.eulersubs[deg,p,charR,charI])/.unknownssubs,
        {j,1,Length[facn]}];
];
Table[a[n]->aa[n],{n,1,numterms}];
Table[{Re[aa[n]],Im[aa[n]]},{n,1,numterms}];
Flatten[Table[{bb1[n]->Re[aa[n]],bb2[n]->Im[aa[n]]},{n,1,numterms}]]
];



findgsandnumterms[numtermsIN_,FE_,gANDsList_,ev_,gflag_,PRECIS_,
      absflag_,EP_,starpts_]:=Block[{unknowns,gsCount,numterms,numdetectors},

  If[numtermsIN > 1, numterms1=numtermsIN,
    numterms1={findSumlim[
        FE, gANDsList[[1,1]], gANDsList[[1,3]], ev, Abs[numtermsIN],gflag,PRECIS,Sign[absflag]][[1]],
               findSumlim[
        FE, gANDsList[[1,2]], gANDsList[[1,3]], ev, Abs[numtermsIN],gflag,PRECIS,Sign[absflag]][[1]]};
  ];  (* If  *)

  Print["numterms1: ",numterms1];
  numterms=Max[numterms1];

  If[numterms < 3, numterms = 3];  (* later code assumes at least 4 unknowns.
                                      with numterms=3, the unknowns are bb1[2], bb2[2], bb1[3], bb2[3] *)

  unknowns=theunknowns[EP,numterms];
  Print["number of unknowns for first equation: ", Length[unknowns], " which are ", unknowns];

  (* now build up the list of {g1,g2,s} for solving and detecting *)
  (* keep adding pairs until there are enough, adjusting the target along the way*)

  gsCount=1;
  If[Length[starpts]==1,numdetectors=0,numdetectors=8];  (* usually 8 detectors *)
  If[numtermsIN > 1,
    numterms=numtermsIN; unknowns=theunknowns[EP,numterms];gsCount=Length[unknowns],
    While[(gsCount<Length[unknowns]+numdetectors && gsCount<Length[gANDsList]) || gsCount<8,
       ++gsCount;
       numterms2={findSumlim[
          FE, gANDsList[[gsCount,1]], gANDsList[[gsCount,3]], ev, Abs[numtermsIN],gflag,PRECIS,Sign[absflag]][[1]],
                  findSumlim[
          FE, gANDsList[[gsCount,2]], gANDsList[[gsCount,3]], ev, Abs[numtermsIN],gflag,PRECIS,Sign[absflag]][[1]]};
Print[gsCount,"   ",numterms2];
      If[Max[numterms2]> numterms,
         numterms=Max[numterms2];
         unknowns=theunknowns[EP,numterms];
         Print["Increasing numterms ", gsCount,"  ",numterms2, " num unknowns needed now: ",Length[unknowns]];
      ];  (* If Max *)
    ]; (* While *)
  ]; (* If numtermsIN > 1 *)

  If[gsCount<Length[unknowns]+numdetectors,
      Print["Error:  not enough detectors"];
      Return[{0,0,0,0}]
  ];  (* If *)

  {unknowns,gsCount,numterms,numdetectors}
];

findSIGNunknowns[feold_, ep_]:= Block[{signunknowns, j},
    signunknowns = {EpsilonR, EpsilonI};
    For[j=1, j<=Length[ep[[2,1]]], ++j,
        signunknowns = Flatten[{signunknowns, EPSILONpR[ ep[[2,1,j]] ], EPSILONpI[ep[[2,1,j]] ]}]
    ];
    signunknowns
];

finddegree2eqns[EP_,unknowns_]:=Block[{theconductor,degree2eqns,thecharacter},
  (* there are equations from degree 2 (char determines coeff sign),
     and also |eps| = 1 if the sign is unknown *)
  degree2eqns={};
  If[EP[[1,1]] == 2,
     (* degree 2 Euler product, so get a new equation when the character is nonzero *)
    thecharacter=EP[[1, 2]];
    theconductor=Length[thecharacter];
    For[jp=1,jp<=Length[unknowns],jp +=2,
     theindex=unknowns[[jp]][[-1]];
     If[thecharacter[[Mod[theindex,theconductor,1] ]] ==0,Continue[],
       thecharval= thecharacter[[Mod[theindex,theconductor,1] ]];
       If[thecharval == 1,AppendTo[degree2eqns,bb2[theindex]],
         If[thecharval == -1, AppendTo[degree2eqns,bb1[theindex]],
           AppendTo[degree2eqns,
             bb1[theindex] - (Re[thecharval] bb1[theindex] + Im[thecharval] bb2[theindex])]
           ]
         ]
       ]
    ]
  ];

  degree2eqns
];

findsigneqns[feold_,ep_] := Block[{},
    thesign = feold[[5]];
    If[Coefficient[thesign, EpsilonR] == 0,
        (* sign is known *)
        signeqns = {EpsilonR - Re[thesign], EpsilonI -  Im[thesign]}
      ,
        signeqns = {EpsilonR^2 + EpsilonI^2 - 1}
    ];

    newsigneqns = {};
    newEpApeqns = {};
    If[Length[ep[[2,1]]] > 0,
      fenew = FEoldtonew[feold];
      theinfinitysign = infinitysign[fenew];
      If[Length[ep[[2,1]]] == 1,
        thisprime = ep[[2,1,1]];
        newsigneqns = {Re[theinfinitysign] EPSILONpR[thisprime] - Im[theinfinitysign] EPSILONpI[thisprime] - EpsilonR,
	               Re[theinfinitysign] EPSILONpI[thisprime] + Im[theinfinitysign] EPSILONpR[thisprime] - EpsilonI};
        If[Length[ep[[2]]] > 2,
           newEpApeqns = ep[[2,3,1]];
        ]
        ,
        Print["error: multiple bad primes not implemented", ep, fenew]
       ]
   ];
   signeqns = Flatten[{newsigneqns, newEpApeqns, signeqns}];   
   signeqns
];

DONOTUSEfindbadfactoreqns[EP_]:= Block[{neweqns,thisp},
  neweqns = {};
  If[Length[EP] > 2,
    If[EP[[3,1]] == 3 && EP[[3,3]] == "IVa",
       thisp = EP[[3,2]];
       AppendTo[neweqns, bb1[thisp]^2 + bb2[thisp]^2 - bb1[1]],
       If[EP[[1,1]] == 3 && EP[[3,3]] == "IIIa",
         thisp = EP[[2,1,1]];
         Print["assigning specific value to a sub", thisp];
         Print["the values are", EP[[3,1]], "plus I times", EP[[3,2]]];
         AppendTo[neweqns, bb1[thisp] -  bb1[1]] EP[[3,1]];
         AppendTo[neweqns, bb2[thisp] - bb1[1] EP[[3,2]]]],
       Print["Error:  unimplemented bad factor"]
    ]
  ];

  neweqns
];

(* probably wrong, and needed to be rethought anyway   *)
DONTUSEfindbadfactorsubstitutions[EP_]:= Block[{neweqns,thisp},
neweqns = {};
If[Length[EP] > 2 && False,
    If[EP[[3,1]] == 3 && EP[[3,3]] == "IVa",
       thisp = EP[[3,2]];
       AppendTo[neweqns, bb1[thisp]^2 + bb2[thisp]^2 - bb1[1]],
       If[EP[[1,1]] == 3 && EP[[3,3]] == "IIIa",
         thisp = EP[[2,1,1]];
         Print["assigning specific value to a sub", thisp];
         Print["the values are", EP[[3,1]], "plus I times", EP[[3,2]]];
         AppendTo[neweqns, bb1[thisp] -  bb1[1]] EP[[3,1]];
         AppendTo[neweqns, bb2[thisp] - bb1[1] EP[[3,2]]]],
       Print["Error:  unimplemented bad factor"]
    ]
];
  If[EP[[1,1]]==3 && Length[EP[[2,1]]] == 1,
    thebadprime = EP[[2,1,1]];
    If[EP[[2,2,1]] != EP[[1,1]]-1,
       Print["Error: only prime level case implemented"],
       (* so bb1[p] and bb2[p] are not unknowns, but theta and phi are *)
       neweqns = {bb1[thebadprime] -> (Cos[theta]+Cos[phi]/Sqrt[thebadprime]), bb2[thebadprime] -> (Sin[theta]+Sin[phi]/Sqrt[thebadprime])};
       AppendTo[neweqns, bb1[thebadprime^2] -> (1/(2 thebadprime))(Cos[2 phi] + Cos[2 phi] + 2 Sqrt[thebadprime] Cos[phi] Cos[theta] + 2 thebadprime Cos[theta]^2 - 2 Sqrt[thebadprime] Sin[phi] Sin[theta] - 2 thebadprime Sin[theta]^2)];
       AppendTo[neweqns, bb2[thebadprime^2] -> (Cos[theta] Sin[phi])/Sqrt[thebadprime] + (Cos[phi] (2 Sin[phi] + Sqrt[thebadprime] Sin[theta]))/thebadprime + Sin[2 theta]];
    ]
Print["made findbadfactorsubstitutions",neweqns, "from", {EP[[2,1,1]],EP[[2,2,1]],EP[[1,1]]-1}];
  ];
  If[EP[[1,1]]==2 && Length[EP[[2,1]]] == 1,
    thebadprime = EP[[2,1,1]];
    If[EP[[2,2,1]] != EP[[1,1]]-1,
       Print["Error: only prime level case implemented"],
       (* so bb1[p] and bb2[p] are not unknowns, but theta and phi are *)
       neweqns = {bb1[thebadprime] -> (Cos[theta]/Sqrt[thebadprime]), bb2[thebadprime] -> (Sin[theta]/Sqrt[thebadprime])}
    ]
Print["made findbadfactorsubstitutions",neweqns, "from", {EP[[2,1,1]],EP[[2,2,1]],EP[[1,1]]-1}];
  ];
neweqns
];


tossRepeats[lis_] := tossRepeats[lis, 0.01];

tossRepeats[lis_, tolerance_] := Block[{slis, j, nn},
  slis = Sort[lis];
  If[NumberQ[slis[[1,1]]],
      slis = Sort[slis, #1[[1]] < #2[[1]]&]
    ,
      If[NumberQ[slis[[1,1,1]]],
        slis = Sort[slis, #1[[1,1]] < #2[[1,1]]&]
      ,
        slis = Sort[slis, #1[[1,1,1]] < #2[[1,1,1]]&]
      ]
  ];

  Print["Initial number of items: ", Length[slis]];
  counter = 0;
  fullcounter = 0;
  For[j = 1, j <= Length[slis] - 1, ++j,
   fullcounter += 1;
   If[NumberQ[slis[[1,1,1,1]]],  (* could be ldata or zdata *)
       thisitem = slis[[j]];
       nextitem = slis[[j + 1]]
     ,
       thisitem = slis[[j,1]];
       nextitem = slis[[j + 1,1]]
   ];
   (* check if the spectral parameters are close *)
   theyareclose = True;
   (* first check same funcitonal equation *)
   If[thisitem[[1, 2]] != nextitem[[1, 2]],
    theyareclose = False];
   If[Norm[thisitem[[1, 1]] - nextitem[[1, 1]]] > tolerance,
    theyareclose = False];
   For[nn = 1, nn <= 4, ++nn,
    If[Abs[thisitem[[1, 4, nn]] - nextitem[[1, 4, nn]]] > tolerance,
     theyareclose = False
     ]
    ];
   If[theyareclose,Print[fullcounter," ", N[thisitem[[1,1]]]," near ",N[nextitem[[1,1]]]," FE: ",InputForm[thisitem[[1,2]]]];
    counter += 1;
    If[thisitem[[2, 1, 1]] > nextitem[[2, 1, 1]],
       Print["deleting item ", j, " of ", Length[slis]];
       slis = Delete[slis, j];
       j -= 1
      ,
       Print["deleting item ", j + 1, " of ", Length[slis]];
       slis = Delete[slis, j + 1];
       j -= 1
      ,
       Print["deleting item ", j, " of ", Length[slis]];
       slis = Delete[slis, j];
       j -= 1
    ]
   ]
  ];
  Print["Omitted ", counter, " duplicates"];
  Print["Trimmed list length: ", Length[slis]];
  slis];

(*   for finding zeros and zero index after values have been computed  *)

tossRepeatsSimple[lis_, eps_] := Block[{slis},
  slis = Sort[lis];
  For[j = 1, j <= Length[slis] - 1, ++j,
   If[slis[[j + 1]] - slis[[j]] < eps,
    slis = Delete[slis, j];
    j -= 1]];
  slis]

trimList[lis_, lowerbound_, upperbound_] :=
 Cases[lis, s_ /; lowerbound < s < upperbound]

findallzeros[vallist_] := Block[{x(*,a,b,c,d*)},
  thisfunction = Interpolation[vallist];
  allroots = {};
  (* first find crossings *)
  For[j = 1, j <= Length[vallist] - 1, ++j,
   If[vallist[[j, 2]] vallist[[j + 1, 2]] <= 0,
    midpt = (vallist[[j, 1]] + vallist[[j + 1, 1]])/2;
    thisroot = x /. FindRoot[thisfunction[x], {x, midpt}];
    AppendTo[allroots, thisroot]
    ]
   ];
  For[j = 1, j <= Length[vallist] - 3, ++j,
   a = vallist[[j, 2]]; b = vallist[[j + 1, 2]];
   c = vallist[[j + 2, 2]]; d = vallist[[j + 3, 2]];
   If[((a - b) *(c - d) < 0 )(* local amx or min *)
     && ( Sign[a] == Sign[b] == Sign[c] == Sign[d])
     && ( a (b - a) < 0),
    thisroot = x /. FindRoot[thisfunction[x], {x, vallist[[j + 1, 1]]}];
    AppendTo[allroots, thisroot];
    thisroot = x /. FindRoot[thisfunction[x], {x, vallist[[j + 2, 1]]}];
    AppendTo[allroots, thisroot];
    ]
   ];
  allroots = trimList[allroots, vallist[[1, 1]], vallist[[-1, 1]]];
  checkedroots = {};
  For[j = 1, j <= Length[allroots], ++j,
   If[Abs[thisfunction[allroots[[j]]]] < 0.0001,
    AppendTo[checkedroots, allroots[[j]]]]
   ];
  checkedroots = tossRepeatsSimple[checkedroots, 0.00001];
  checkedroots
  ]

addzerodatatoZ[elem_] := Block[{i, j},
  thisitem = elem;
  thisfedata = thisitem[[1]];
  (*
  If[Abs[thisfedata[[1,1,2]]]<10^-8,Print["setting ", thisfedata[[1,1,
  2]], " to 0"];
  thisfedata[[1,1,2]]=0];
  *);
  thiseigs = thisfedata[[1, 1]];
  thisparamR = thisfedata[[1, 2, 1]];
  thisparamC = thisfedata[[1, 2, 2]];
  If[Length[thiseigs] == 1,
   thistrivzeroheights = Sort[{-thiseigs[[1]], thiseigs[[1]]}]
   ];
  If[Length[thiseigs] == 2 && Length[thisparamC] == 0,
   thistrivzeroheights = Sort[{-thiseigs[[1]], -thiseigs[[2]], thiseigs[[1]] + thiseigs[[2]]}]
   ];
  If[Length[thiseigs] == 2 && Length[thisparamC] == 1,
   thistrivzeroheights = Sort[{-thiseigs[[1]], -thiseigs[[ 2]], (thiseigs[[1]] + thiseigs[[2]])/2}]
   ];
  If[Length[thiseigs] == 3,
   thistrivzeroheights = Sort[{-thiseigs[[1]], -thiseigs[[2]], -thiseigs[[3]], thiseigs[[1]] + thiseigs[[2]] + thiseigs[[3]]}]
   ];
  thisvaluedata = thisitem[[2, 1]];
  thisvaluedata = Table[{thisvaluedata[[i, 1]], N[thisvaluedata[[i, 2]]]}, {i, 1, Length[thisvaluedata]}];
  thisprecisiondata = thisitem[[2, 2]];
  thisprecisiondata = Table[{thisprecisiondata[[i, 1]], {thisprecisiondata[[i, 2, 1]], N[thisprecisiondata[[i, 2, 2]]]}}, {i, 1, Length[thisprecisiondata]}];
  jt = 1; While[thisprecisiondata[[jt, 2, 2]] > 0.1, ++jt];
  thislowerlim = thisprecisiondata[[jt, 1]];
  jt = -1; While[thisprecisiondata[[jt, 2, 2]] > 0.1, --jt];
  thisupperlim = thisprecisiondata[[jt, 1]];
  thisvaluedata = Select[thisvaluedata, thislowerlim <= #[[1]] <= thisupperlim &];
  thesezeros = findallzeros[thisvaluedata];
  thiszerosig = Table[Length[
     Select[thesezeros, thistrivzeroheights[[j]] < # < thistrivzeroheights[[j + 1]] &]],
    {j, 1, Length[thistrivzeroheights] - 1}];
  thisspecialvalues = thisitem[[3]];
  theanswer = {thisfedata, {thisvaluedata, thisprecisiondata}, {thesezeros, thiszerosig}, thisspecialvalues};
  theanswer
  ]


