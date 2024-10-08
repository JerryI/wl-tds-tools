BeginPackage["JerryI`TDSTools`Material`", {
  "JerryI`TDSTools`Utils`",
  "JerryI`TDSTools`Trace`",
  "JerryI`TDSTools`Transmission`"
}]

MaterialParameters::usage = "TransmissionObject[sam_TDTrace, ref_TDTrace, opts___] creates a transmission from sample and reference traces. See _TransmissionObject[\"Properties\"]"

Begin["`Private`"]

root = DirectoryName[$InputFileName];

{ 
  initialize, 
  solveNK, 
  solveFP, 
  movingAverage,
  saveForDebug,
  clusterPhase
} = Get[ FileNameJoin[{root, "nCPU.wl"}] ];

Options[MaterialParameters] = {"Target"->"CPU", "Tags"-><||>, "Model"->"Slab", "SolveNKEquations"->True, "MovingAverageFilter"->True, "FabryPerotCancellation"->True, "FabryPerotCycles"->8, "NKCycles"->30};
Options[materialParameters] = Options[MaterialParameters];

MaterialParameters[t_TransmissionObject, opts:OptionsPattern[]] := materialParameters[t, OptionValue["Model"], OptionValue["Target"], opts]

MaterialParameters[t:List[__TransmissionObject], opts:OptionsPattern[]] := materialParameters[t, OptionValue["Model"], OptionValue["Target"], opts]

MaterialParameters::unknown = "Unknown model or target";
MaterialParameters::unknownprop = "Unknown property `1`";

materialParameters[_, _, _, ___] := (Message[MaterialParameters::unknown]; $Failed)

MaterialParameters[n_Association][s_String] := If[!KeyExistsQ[n, s], 
  Message[MaterialParameters::unknownprop, s]; $Failed
,
  n[s]
]

MaterialParameters[n_Association]["Properties"] := {"Domain", "Transmission", "\[Alpha]", "Frequencies", "n", "k", "Properties", "Thickness", "Tags", "Gain", "PhaseShift", "Phase", "FrequencyDomainConfidenceInterval", "FDCI", "FDCI2", "FrequencyDomainConfidenceInterval2", "FPReduction"}

MaterialParameters[n_Association]["FrequencyDomainConfidenceInterval"] := MaterialParameters[n]["FDCI"]

MaterialParameters[n_Association]["FrequencyDomainConfidenceInterval2"] := MaterialParameters[n]["FDCI2"]

MaterialParameters[n_Association]["Transmission"] := QuantityArray[Transpose[{n["Frequencies"] // Normal, Normal[n["T"]]^2}], {1/"Centimeters", 1}]

MaterialParameters[n_Association]["Phase"] := QuantityArray[Transpose[{n["Frequencies"] // Normal, Normal[n["Phase"]]}], {1/"Centimeters", 1}]

MaterialParameters[n_Association]["n"] := QuantityArray[Transpose[{n["Frequencies"] // Normal, n["n"] // Normal}], {1/"Centimeters", 1}]
MaterialParameters[n_Association]["k"] := QuantityArray[Transpose[{n["Frequencies"] // Normal, n["k"] // Normal}], {1/"Centimeters", 1}]

MaterialParameters[n_Association]["\[Alpha]"] := QuantityArray[{#[[1]], (#[[2]] 4 \[Pi]  10^12 #[[1]])/(33.356 2.9979 10^10)} &/@ Transpose[{n["Frequencies"] // Normal, n["k"] // Normal}], {1/"Centimeters", 1/"Centimeters"}]

MaterialParameters[n_Association]["Frequencies"] := QuantityArray[n["Frequencies"] // Normal, 1/"Centimeters"]
MaterialParameters[n_Association]["Domain"] := Quantity[#, 1/"Centimeters"] &/@ MinMax[Normal @ n["Frequencies"]]


MaterialParameters[a_Association]["FDCI2"] := 
  With[{phase = MaterialParameters[a]["Phase"] // QuantityMagnitude}, 
    With[{r = Table[
        {phase[[q, 1]], 
          LinearModelFit[Take[phase, q], {1, x}, x]["RSquared"]}, 
        {q, Round[0.2 (phase // Length)], phase // Length, 25}]}, 
      With[{scape = {Drop[r[[All, 1]], -4], 
          MovingAverage[r[[All, 2]], 4] // Differences} // Transpose}, 
        Quantity[#, 1/"Centimeters"] & /@ {0.0, 
          scape[[FindPeaks[-scape[[All, 2]], 5][[All, 1]], 1]] // First}
      ]
    ]
  ]

MaterialParameters[a_Association]["FPReductionPlot"] := 
  With[{mFp = MaterialParameters[a]}, 
    With[{testRegion = 
      With[{range = findFDCIRanges[mFp]}, 
        {
            MovingAverage[
              Drop[
                Fourier[
                  (mFp[[1]]["k debug"] // Normal)[[range[[1]] ;; range[[2]]]]
                ] // Abs // dropHalf, 
                1
              ], 
              3
            ]
           , 
            MovingAverage[
              Drop[
                Fourier[
                  (mFp[[1]]["k"] // Normal)[[range[[1]] ;; range[[2]]]]
                ] // Abs // dropHalf, 
                1
              ], 
              3
            ]
        }
      ]}, 
      
      With[{r = # & /@ testRegion}, 
        ListLinePlot[r, PlotRange->Full, PlotLegends->{"FP", "no FP"}]
      ]
    ]
  ]

MaterialParameters[a_Association]["FPReduction"] := 
  With[{mFp = MaterialParameters[a]}, 
    With[{testRegion = 
      With[{range = findFDCIRanges[mFp]}, 
        {
          Interpolation[
            MovingAverage[
              Drop[
                Fourier[
                  (mFp[[1]]["k debug"] // Normal)[[range[[1]] ;; range[[2]]]]
                ] // Abs // dropHalf, 
                2
              ], 
              4
            ], 
            InterpolationOrder -> 1
          ], 
          Interpolation[
            MovingAverage[
              Drop[
                Fourier[
                  (mFp[[1]]["k"] // Normal)[[range[[1]] ;; range[[2]]]]
                ] // Abs // dropHalf, 
                2
              ], 
              4
            ], 
            InterpolationOrder -> 1
          ]
        }
      ]}, 
      
      With[{r = NIntegrate[
        #[x], {x, #["Domain"][[1, 1]], #["Domain"][[1, 2]]}] & /@ testRegion}, 
        r[[1]] / r[[2]]
      ]
    ]
  ]

   


MaterialParameters /: MakeBoxes[obj: MaterialParameters[a_Association], StandardForm] := With[{
  preview = With[{
    n = ArrayResample[Select[Normal[a["n"]]//dropHalf//keepMiddle, #>0.0 &], 150],
    k = ArrayResample[Select[Normal[a["k"]]//dropHalf//keepMiddle, #>0.0 &], 150]
  },
    With[{
      rangesK = MinMax[k],
      rangesN = MinMax[n]
    },
      ListLinePlot[{k, rangesK[[1]] + (rangesK[[2]]-rangesK[[1]])(#-rangesN[[1]])/(rangesN[[2]]-rangesN[[1]]) &/@ n}, PlotStyle->{ColorData[97][2], ColorData[97][3]}, PlotRange->Full,Axes -> None, ImagePadding->None]
    ]
  ]
},


    Module[{above, below},
        above = { 
          {BoxForm`SummaryItem[{"From ", Quantity[Part[a["Frequencies"],1], 1/"Centimeters"]//Round }]},
          {BoxForm`SummaryItem[{"To ", Quantity[Part[a["Frequencies"],-1], 1/"Centimeters"]//Round }]},
          {BoxForm`SummaryItem[{"Thickness ", a["Thickness"]}]},
          {BoxForm`SummaryItem[{"n ", Round[a["n0"], 0.01]}]},
          {BoxForm`SummaryItem[{"Gain ", a["Gain"]}]},
          {BoxForm`SummaryItem[{"PhaseShift ", 2Pi a["PhaseShift"]}]},
          
          If[Length[obj["Tags"]//Keys] > 0, {BoxForm`SummaryItem[{"Tags ", Style[#, Background->LightBlue]&/@obj["Tags"]//Keys}]}, Nothing]
        };

        BoxForm`ArrangeSummaryBox[
           MaterialParameters, (* head *)
           obj,      (* interpretation *)
           preview,    (* icon, use None if not needed *)
           (* above and below must be in a format suitable for Grid or Column *)
           above,    (* always shown content *)
           Null (* expandable content. Currently not supported!*)
        ]
    ]  
]

materialParameters[t_List, type_, "CPU", opts: OptionsPattern[]] := materialParameters[#, type, "CPU", opts] &/@ t
materialParameters[t_List, type_, "GPU", opts: OptionsPattern[]] := materialParameters[#, type, "GPU", opts] &/@ t

MaterialParameters::badgpu = "GPU acceleration is not available. Reason: ``";

materialParameters[t:TransmissionObject[a_], "Slab", "GPU", opts: OptionsPattern[]] := With[{},
  Message[MaterialParameters::badgpu, "not implemented"];
  materialParameters[t, "Slab", "CPU", opts]
]

materialParameters[TransmissionObject[a_], "Slab", "CPU", opts: OptionsPattern[]] := Module[{packed, src, dest}, With[{
 thickness = QuantityMagnitude[a["Thickness"], "Centimeters"] // N,
 \[Delta]t = QuantityMagnitude[a["\[Delta]t"], "Picoseconds"] // N,
 freqs = Normal[a["Frequencies"]],
 fpIterations = If[OptionValue["FabryPerotCancellation"] === True, 
   OptionValue["FabryPerotCycles"],
   0
 ]
},

  With[{
    src = Transpose[{freqs, Normal[a["T"]], Normal[a["Phase"]]}] // Flatten
  },
    (* dest = Transpose[{Array[0.0&, Length[freqs]], Array[0.0&, Length[freqs]]}] // Flatten; *)
    
    meta = {\[Delta]t, thickness, 1.0 a["Gain"], 2Pi Round[a["PhaseShift"]] // N, Length[freqs]};

    dest = initialize[src, meta];
    dest = solveNK[src, dest, meta, If[!OptionValue["SolveNKEquations"], 0, OptionValue["NKCycles"]]];
    dest = movingAverage[src, dest, meta];
    dest = saveForDebug[dest, meta];
    
    Do[
      dest = solveFP[src, dest, meta, 0]; 
      dest = solveNK[src, dest, meta, If[!OptionValue["SolveNKEquations"], 0, OptionValue["NKCycles"]]];
      dest = movingAverage[src, dest, meta]; 
    , {fpIterations}];

    dest = Partition[dest,5] // Transpose;

    MaterialParameters[<|
      "n0" -> 1.0 + (0.029979 meta[[1]] / meta[[2]]),
      "n" -> NumericArray[dest[[1]]],
      "k" -> NumericArray[dest[[2]]],
      "T" -> NumericArray[dest[[3]]],
      "k debug" -> NumericArray[dest[[5]]],
      "Phase" -> NumericArray[dest[[4]]],
      "Frequencies" -> NumericArray[a["Frequencies"]],
      "Thickness" -> a["Thickness"],
      "Gain" -> a["Gain"],
      "FDCI" -> a["FDCI"],
      "PhaseShift" -> a["PhaseShift"],
      "Date" -> Now,
      "Tags" -> Join[a["Tags"], OptionValue["Tags"]]
    |>]
  ]
]]

findFDCIRanges[mFp_MaterialParameters] := 
  With[{m = mFp["Frequencies"] // QuantityMagnitude}, 
    {
      FirstPosition[m, x_ /; (x > QuantityMagnitude[mFp["FDCI"][[1]]])] // First,
      Length[m] - FirstPosition[m // Reverse, x_ /; (x < QuantityMagnitude[mFp["FDCI"][[2]]])] // First
    }
  ];




End[]
EndPackage[]