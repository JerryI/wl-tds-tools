BeginPackage["JerryI`TDSTools`nCPU`"]
Begin["`Private`"]

If[!(compiledQ // TrueQ),
  Print[Style["TDTools`nCPU >> Compiling to C", Italic] ];

  initialize = Compile[{{src, _Real, 1}, {meta, _Real, 1}}, Module[{out = Table[0.0, {5 meta[[5]]}]},
    With[{
      n0 = 1.0 + (0.029979 meta[[1]] / meta[[2]]),
      ft = 1.0 / meta[[2]]
    },
      out[[1;;10]] = {n0,0.0,src[[2]],src[[3]],0.0,n0,0.0,src[[5]],src[[6]], 0.0};


      Do[
        out[[5 i - 4]] = 1.0 + 0.159152 (src[[3 i - 0]] + meta[[4]]) ft / src[[3 i - 2]];
        out[[5 i - 3]] = - 0.159152 Log[src[[3 i - 1]] meta[[3]] ] ft / src[[3 i - 2]];
        out[[5 i - 2]] = src[[3 i - 1]];
        out[[5 i - 1]] = src[[3 i - 0]];
        out[[5 i - 0]] = 0.0;
      , {i, 2, meta[[5]]}];
    ];

    out
  ], "CompilationTarget" -> "C", "RuntimeOptions" -> "Speed"];

  solveNK = Compile[{{src, _Real, 1}, {dest, _Real, 1}, {meta, _Real, 1}, {iterations, _Integer}}, Module[{
    np = 0.0, kp = 0.0,
    logT = 0.0, ph = 0.0,
    n=0.0, k=0.0, denominator=0.0, arg=0.0, modulus=0.0, im=0.0, re=0.0, n2=0.0, fT=0.0,
    buffer = {0.0,0.0},

    ft = 1.0/meta[[2]],
    out = Table[0.0, {5 meta[[5]]}],
    n0 = 1.0 + (0.029979 meta[[1]] / meta[[2]])
  },

    out[[1;;10]] = dest[[1;;10]];

    Do[With[{
      freq = src[[3 i - 2]]
    },

    np = 0.0; kp = 0.0;
    logT = 0.0; ph = 0.0;
    n=0.0; k=0.0; denominator=0.0; arg=0.0; modulus=0.0; im=0.0; re=0.0; n2=0.0; fT=0.0;

      n = dest[[5 i - 4]];
      k = dest[[5 i - 3]];

      fT = 1.0 / (freq meta[[2]]);
      np = dest[[5 i - 4]];
      kp = dest[[5 i - 3]];
      logT = Log[dest[[5 i - 2]] meta[[3]] ];
      ph   = dest[[5 i - 1]] +  meta[[4]];

      Do[
            n2 = (1.0 + np);
            n2 = n2 n2;
            denominator = 1.0/(kp kp + n2);
            denominator = denominator denominator;

            re = denominator (np n2 + kp kp (2.0 + np));
            im = denominator (kp (kp kp + np np - 1.0));

            modulus = Sqrt[re re + im im];
            arg = Arg[I im + re];

            n = 1.0 + 0.159152  fT  (ph - arg);
            k = - 0.159152  fT  (logT - Log[4.0 modulus]);

            np = n; kp = k;  
      , {j, 1, iterations}];


      If[!NumericQ[k] || !NumericQ[n] || n < 1.0 || !RealValuedNumberQ[k] || !RealValuedNumberQ[n], 
        k = 0.0; n = n0;
        (* out[[5 i - 0]] = - 1.0; *)
      ];

      (* out[[4 i - 3]] = (buffer[[1]] + n)/2.0; *)
      (* out[[4 i - 2]] = (buffer[[2]] + k)/2.0; *)


      out[[5 i - 4]] = n;
      out[[5 i - 3]] = k;

      out[[5 i - 2]] = dest[[5 i - 2]];
      out[[5 i - 1]] = dest[[5 i - 1]];
      out[[5 i - 0]] = dest[[5 i - 0]];

    ], {i, 2, meta[[5]]}];

    out
  ], "CompilationTarget" -> "C", "RuntimeOptions" -> "Speed"];


  solveFP = Compile[{{src, _Real, 1}, {dest, _Real, 1}, {meta, _Real, 1}, {iterations, _Integer}}, Module[{
        var40, var41, var54, var57, var85, var90, var81, var82, var58, var93, var91, 
        var83, var84, var86, var87, var88, var89, var94, var80, var99, var98, 
        var92, abs, arg, fT=0.0, 
    CONSTF = 6.28331 * meta[[2]],
    out = Table[0.0, {5 meta[[5]]}]
  },

    out[[1;;10]] = dest[[1;;10]];

    Do[With[{
      freq = src[[3 i - 2]],
      t    = src[[3 i - 1]],
      ph   = src[[3 i - 0]],

      n = dest[[5 i - 4]],
      k = dest[[5 i - 3]]
    },


      (* Perform the optimized mathematical operations *)
      var40 = k^2;
      var41 = 1.0 + n;
      var54 = var41^2;
      var57 = var40 + var54;
      var85 = n^2;
      var90 = 2.0 * CONSTF * freq * n;
      var81 = -2.0 * CONSTF * freq * k;
      var82 = Exp[var81];
      var58 = 1.0 / (var57^2);
      var93 = -1.0 + var40 + var85;
      var91 = Cos[var90];
      var83 = -2.0 + k;
      var84 = var83 * k;
      var86 = -1.0 + var84 + var85;
      var87 = 2.0 + k;
      var88 = k * var87;
      var89 = -1.0 + var88 + var85;
      var94 = Sin[var90];
      var80 = var57^2;
      var99 = var40 + (-1.0 + n)^2;
      var98 = Exp[2.0 * CONSTF * freq * k];

      var92 = 4.0 * k * var93 * var94;

      abs = Sqrt[
        var58 * ((var99 / var98)^2 + var80 + var82 * (-2.0 * var86 * var89 * var91 + 2.0 * var92))
      ];

      arg = Arg[
        I(var82 * var58 * (4.0 * k * var93 * var91 + var86 * var89 * var94)) + 
        var82 * var58 * (var98 * var80 - var86 * var89 * var91 + var92)
      ];

      (* Handle non-finite values *)
      If[!NumericQ[abs] || !NumericQ[arg] || !FiniteQ[abs] || !FiniteQ[arg], 
        abs = 1.0; arg = 0.0;
        (* out[[5 i - 0]] = - 1.0; *)
      ];  

      (* Apply condition *)
      If[t abs > 1.2,
        out[[5 i - 2]] = 1.0;
        out[[5 i - 1]] = 0.0;
      ,
        out[[5 i - 2]] = t abs ;
        out[[5 i - 1]] = ph - arg ;
      ];

      out[[5 i - 0]] = dest[[5 i - 0]];

      out[[5 i - 4]] = dest[[5 i - 4]];
      out[[5 i - 3]] = dest[[5 i - 3]];

    ], {i, 2, meta[[5]]}];

    out

  ], "CompilationTarget" -> "C", "RuntimeOptions" -> "Speed"];


  movingAverage = Compile[{{src, _Real, 1}, {dest, _Real, 1}, {meta, _Real, 1}}, Module[{
    out = Table[0.0, {5 meta[[5]]}]
  },

    Do[With[{

    },


      If[i == meta[[5]],
        out[[5 i - 4]] = dest[[5 i - 4]];
        out[[5 i - 3]] = dest[[5 i - 3]];
        out[[5 i - 2]] = dest[[5 i - 2]];  
        out[[5 i - 1]] = dest[[5 i - 1]];
        out[[5 i - 0]] = dest[[5 i - 0]];
      ,
        out[[5 i - 4]] = (dest[[5 i - 4]] + dest[[5 (i+1) - 4]])/2.0;
        out[[5 i - 3]] = (dest[[5 i - 3]] + dest[[5 (i+1) - 3]])/2.0;  
        out[[5 i - 1]] = dest[[5 i - 1]];
        out[[5 i - 2]] = dest[[5 i - 2]];
        out[[5 i - 0]] = dest[[5 i - 0]];
      ];

    ], {i, 1, meta[[5]]}];

    out

  ], "CompilationTarget" -> "C", "RuntimeOptions" -> "Speed"];



  saveForDebug = Compile[{{dest, _Real, 1}, {meta, _Real, 1}}, Module[{
    out = dest
  },

    Do[With[{

    },
        out[[5 i - 0]] = dest[[5 i - 3]];

    ], {i, 1, meta[[5]]}];

    out

  ], "CompilationTarget" -> "C", "RuntimeOptions" -> "Speed"];

  clusterPhase = 
    Compile[
      {
        {input, _Real, 1}, 
        {low, _Integer, 0}, 
        {high, _Integer, 0}, 
        {PhaseTrashhold, _Real, 0}
      }, 
      Module[{phase = input, delta = 0.0},

        (* Adjust phase values based on PhaseTrashhold *)
        Do[
          If[phase[[i + 1]] - phase[[i]] > PhaseTrashhold,
            Do[phase[[j]] = phase[[j]] - 2 \[Pi];, {j, i + 1, Length[phase]}],

            If[phase[[i + 1]] - phase[[i]] < -PhaseTrashhold,
              Do[phase[[j]] = phase[[j]] + 2 \[Pi];, {j, i + 1, Length[phase]}]
            ]
          ],
          {i, low, high}
        ];

        (* Adjust phase before the 'low' index *)
        If[low > 1,
          delta = Differences[phase[[low ;; Floor[(3 high + low)/4] ]]  ] // Mean;

          Do[
            phase[[-i - 1]] = phase[[-i]] - delta/2,
            {i, -low, -2}
          ];
        ];

        (* Adjust phase after the 'high' index *)
        delta = Differences[phase[[Floor[(high + 3 low)/4] ;; high]]] // Mean;

        Do[
          phase[[i + 1]] = phase[[i]] + delta,
          {i, high, Length[phase] - 1}
        ];

        phase
    ], 
      "CompilationTarget" -> "C", 
      "RuntimeOptions" -> "Speed"
    ];

  Print[Style["TDTools`nCPU >> Done", Italic, Background->LightGreen] ];

  compiledQ = True;
];

End[]
EndPackage[]

{
  JerryI`TDSTools`nCPU`Private`initialize, 
  JerryI`TDSTools`nCPU`Private`solveNK, 
  JerryI`TDSTools`nCPU`Private`solveFP, 
  JerryI`TDSTools`nCPU`Private`movingAverage,
  JerryI`TDSTools`nCPU`Private`saveForDebug,
  JerryI`TDSTools`nCPU`Private`clusterPhase
}