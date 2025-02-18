BeginPackage["JerryI`TDSTools`nCPULL`"]
Begin["`Private`"]


getLibraryLinkVersion[] := 
Which[
    $VersionNumber >= 14.1, 
        With[{n = LibraryVersionInformation[FindLibrary["demo"] ]["WolframLibraryVersion"]},
            If[!NumberQ[n], 7, n]
        ], 
    $VersionNumber >= 13.1, 
        7, 
    $VersionNumber >= 12.1, 
        6, 
    $VersionNumber >= 12.0, 
        5, 
    $VersionNumber >= 11.2, 
        4, 
    $VersionNumber >= 10.0, 
        3, 
    $VersionNumber >= 9.0, 
        2, 
    True, 
        1
]; 


$directory = DirectoryName[If[$InputFileName == "", 
        NotebookFileName[], 
        $InputFileName
    ], 2]

lib = FileNameJoin[{
        $directory, 
        "LibraryResources", 
        $SystemID <> "-v" <> ToString[getLibraryLinkVersion[] ],
        "nCPU" <> "." <> Internal`DynamicLibraryExtension[]
}]; 


cload = LibraryFunctionLoad[
 lib,  "load", 
 {  
   {Real, _, "Manual"}
 }, Integer];

cunload = LibraryFunctionLoad[
 lib,  "unload", 
 {  
   Integer
 }, Integer]

cget = LibraryFunctionLoad[
 lib,  "get", 
 {  
   Integer
 }, {Real, _, Automatic}]

 crun = LibraryFunctionLoad[
 lib,  "run", 
 {  
   Integer,
   Integer,
   Integer,

   Integer,
   Integer,
   Integer,
   Integer,
   Integer,

   Integer
 }, Integer]

load[dat_, _] :=  cload[ dat ] // llObject; 

unload[llObject[uid_] ] := cunload[uid]; 

get[llObject[uid_] ] := cget[uid]

run[llObject[d_], llObject[s_], llObject[m_], rest__] := crun[d,s,m,rest]


End[]
EndPackage[]

{
  JerryI`TDSTools`nCPULL`Private`run, 
  JerryI`TDSTools`nCPULL`Private`clusterPhase,
  JerryI`TDSTools`nCPULL`Private`clReadyQ,
  JerryI`TDSTools`nCPULL`Private`load,
  JerryI`TDSTools`nCPULL`Private`unload,
  JerryI`TDSTools`nCPULL`Private`get
}