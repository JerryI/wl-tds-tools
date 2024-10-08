BeginPackage["JerryI`TDSTools`nGPU`", {"OpenCLLink`"}]
Begin["`Private`"]

root = DirectoryName[$InputFileName];

If[!(compiledQ // TrueQ), Module[{},
  Print[Style["TDTools`nGPU >> Compiling to OpenCL", Italic] ];

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

  Print[Style["TDTools`nGPU >> Done", Italic, Background->LightGreen] ];

  compiledQ = True;
] ];

End[]
EndPackage[]

{
  JerryI`TDSTools`nGPU`Private`initialize, 
  JerryI`TDSTools`nGPU`Private`solveNK, 
  JerryI`TDSTools`nGPU`Private`solveFP, 
  JerryI`TDSTools`nGPU`Private`movingAverage,
  JerryI`TDSTools`nGPU`Private`saveForDebug,
  JerryI`TDSTools`nGPU`Private`clusterPhase
}