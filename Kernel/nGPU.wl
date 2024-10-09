BeginPackage["JerryI`TDSTools`nGPU`", {"OpenCLLink`"}]
Begin["`Private`"]

root = DirectoryName[$InputFileName];
clFile = FileNameJoin[{root, "nGPU.cl"}];

If[!(compiledQ // TrueQ), Module[{},
  Print[Style["TDTools`nGPU >> Compiling to OpenCL", Italic] ];

  initialize = OpenCLFunctionLoad[File[clFile], "initialize", {
    {"Float", _, "InputOutput"}, 
    {"Float", _, "Input"}, 
    {"Float", _, "Input"}
  }, 256, "ShellOutputFunction"->Print];

  solveNK = OpenCLFunctionLoad[File[clFile], "solveNK", {
    {"Float", _, "InputOutput"}, 
    {"Float", _, "Input"}, 
    {"Float", _, "Input"},
    _Integer
  }, 256, "ShellOutputFunction"->Print];  

  movingAverage = OpenCLFunctionLoad[File[clFile], "movingAverage", {
    {"Float", _, "InputOutput"}, 
    {"Float", _, "Input"}
  }, 256, "ShellOutputFunction"->Print]; 

  saveForDebug = OpenCLFunctionLoad[File[clFile], "saveForDebug", {
    {"Float", _, "InputOutput"}, 
    {"Float", _, "Input"}
  }, 256, "ShellOutputFunction"->Print];       


  solveFP = OpenCLFunctionLoad[File[clFile], "solveFP", {
    {"Float", _, "InputOutput"}, 
    {"Float", _, "Input"}, 
    {"Float", _, "Input"},
    _Integer
  }, 256, "ShellOutputFunction"->Print]; 

  

  Print[Style["TDTools`nGPU >> Done", Italic, Background->LightGreen] ];



  (* compiledQ = True; *)
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