BeginPackage["JerryI`TDSTools`nGPU`", {"OpenCLLink`"}]
Begin["`Private`"]

root = DirectoryName[$InputFileName];
clFile = FileNameJoin[{root, "nGPU.cl"}];

If[!(compiledQ // TrueQ) && OpenCLQ[], Module[{},
  Print[Style["TDTools`nGPU >> Compiling to OpenCL", Italic] ];

  clRun = OpenCLFunctionLoad[File[clFile], "clRun", {
    {"Float", _, "InputOutput"}, 
    {"Float", _, "Input"}, 
    {"Float", _, "Input"},
    _Integer,
    _Integer,
    _Integer,
    _Integer,
    _Integer
  }, 256, "ShellOutputFunction"->Print];

  Print[Style["TDTools`nGPU >> Done", Italic, Background->LightGreen] ];



  (* compiledQ = True; *)
] ];

End[]
EndPackage[]

{
  JerryI`TDSTools`nGPU`Private`clRun, 
  JerryI`TDSTools`nGPU`Private`clusterPhase
}