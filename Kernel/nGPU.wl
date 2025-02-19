BeginPackage["JerryI`TDSTools`nGPU`", {"OpenCLLink`"}]

run;
readyQ;
load;
unload;
get;

Begin["`Private`"]

root = DirectoryName[$InputFileName];
clFile = FileNameJoin[{root, "nGPU.cl"}];

readyQ = OpenCLQ;
load = OpenCLMemoryLoad;
unload = OpenCLMemoryUnload;
get = OpenCLMemoryGet;

If[OpenCLQ[], Module[{},
  Print[Style["TDTools`nGPU >> Compiling to OpenCL", Italic] ];

  run = OpenCLFunctionLoad[File[clFile], "clRun", {
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
