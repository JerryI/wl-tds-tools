#!/usr/bin/env wolframscript

Get["CCompilerDriver`"]; 
Get["LibraryLink`"];

Echo[$CCompiler];
Echo[CCompilers[]];





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

build[lib_String, opts: OptionsPattern[CreateLibrary]] := 
Block[{$directory, $libSrc, $libDir, $linkerOptions}, 
    $directory = DirectoryName[$InputFileName, 2]; 

    $libSrc = File[FileNameJoin[{
        $directory, 
        "Kernel", 
        lib <> ".c"
    }]]; 

    $libDir = FileNameJoin[{
        $directory, 
        "LibraryResources", 
        $SystemID <> "-v" <> ToString[getLibraryLinkVersion[] ]
    }]; 


    If[!FileExistsQ[$libDir], CreateDirectory[]];

    CreateLibrary[$libSrc, lib, 
        "TargetDirectory" -> $libDir, 
        "Debug" -> True, 
        "ShellOutputFunction"->Print,
        "Compiler"->GenericCCompiler,
        "CompilerInstallation"->"/usr/bin",
        "CompilerName"->"gcc",
        opts
    ]
]; 

build["nCPU"];

Print["Finished! Yo"];
