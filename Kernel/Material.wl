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

{ 
  clRun,
  clclusterPhase,
  clReadyQ,
  clLoad,
  clUnload,
  clGet
} = Get[ FileNameJoin[{root, "nGPU.wl"}] ];

{ 
  llRun,
  llclusterPhase,
  llReadyQ,
  llLoad,
  llUnload,
  llGet
} = Get[ FileNameJoin[{root, "nCPULL.wl"}] ];

Options[MaterialParameters] = {"Target"->"CPU", "Tags"-><||>, "Model"->"Slab", "SolveNKEquations"->True, "MovingAverageFilter"->True, "FabryPerotCancellation"->True, "FabryPerotCycles"->8, "NKCycles"->30};
Options[materialParameters] = Join[Options[MaterialParameters], {"sharedSrcQ" -> False, "sharedSrc"->Null}];

MaterialParameters[t_TransmissionObject, opts:OptionsPattern[]] := validateOptions[opts] @  With[{r = materialParameters[{t}, OptionValue["Model"], OptionValue["Target"], opts]},
  If[Head[r] === List, ReleaseHold[r // First], r]
]

MaterialParameters[t:List[__TransmissionObject], opts:OptionsPattern[]] := validateOptions[opts] @  materialParameters[t, OptionValue["Model"], OptionValue["Target"], opts]

MaterialParameters::unknown = "Unknown model or target";
MaterialParameters::unknownprop = "Unknown property `1`";
MaterialParameters::invalidopts = "Invalid options provided"

materialParameters[_, _, _, ___] := (Message[MaterialParameters::unknown]; $Failed)

MaterialParameters[n_Association][s_String] := If[!KeyExistsQ[n, s], 
  Message[MaterialParameters::unknownprop, s]; $Failed
,
  n[s]
]

MaterialParameters /: Keys[t_MaterialParameters] :=  t["Properties"]


MaterialParameters[n_Association]["Properties"] := {"Domain", "Best Transmission", "Transmission", "Raw \[Alpha]", "Best Raw \[Alpha]", "\[Alpha]", "Best \[Alpha]", "Frequencies", "Best n", "n", "Best Raw k", "Raw k", "Best k", "k", "Properties", "Thickness", "Tags", "Gain", "PhaseShift", "Phase", "FrequencyDomainConfidenceInterval", "FDCI", "FDCI2", "FrequencyDomainConfidenceInterval2", "FPReduction"}

MaterialParameters[n_Association]["FrequencyDomainConfidenceInterval"] := MaterialParameters[n]["FDCI"]

MaterialParameters[n_Association]["FrequencyDomainConfidenceInterval2"] := MaterialParameters[n]["FDCI2"]

MaterialParameters[n_Association]["Transmission"] := QuantityArray[Transpose[{n["Frequencies"] // Normal, Normal[n["T"]]^2}], {1/"Centimeters", 1}]

MaterialParameters[n_Association]["Phase"] := QuantityArray[Transpose[{n["Frequencies"] // Normal, Normal[n["Phase"]]}], {1/"Centimeters", 1}]

MaterialParameters[a_Association]["Phase Features"] := With[{
  offset = offsetPhase[ a ],
  raw = Normal[ a["Phase"] ]
},
  QuantityArray[Transpose[{Normal @ a["Frequencies"], raw - offset }], {1/"Centimeters", 1}]
]


MaterialParameters[n_Association]["Best Transmission"] := SelectBestRange[n] @ QuantityArray[Transpose[{n["Frequencies"] // Normal, Normal[n["T"]]^2}], {1/"Centimeters", 1}]

MaterialParameters[n_Association]["Best Phase"] := SelectBestRange[n] @ QuantityArray[Transpose[{n["Frequencies"] // Normal, Normal[n["Phase"]]}], {1/"Centimeters", 1}]


MaterialParameters[a_Association]["Best Phase Features"] := With[{
  offset = offsetPhase[ a ],
  raw = Normal[ a["Phase"] ]
},
  SelectBestRange[n] @ QuantityArray[Transpose[{Normal @ a["Frequencies"], raw - offset }], {1/"Centimeters", 1}]
]


MaterialParameters[n_Association]["n"] := QuantityArray[Transpose[{n["Frequencies"] // Normal, n["n"] // Normal}], {1/"Centimeters", 1}]
MaterialParameters[n_Association]["k"] := QuantityArray[Transpose[{n["Frequencies"] // Normal, n["k"] // Normal}], {1/"Centimeters", 1}]

MaterialParameters[n_Association]["Raw k"] := QuantityArray[Transpose[{n["Frequencies"] // Normal, n["k debug"] // Normal}], {1/"Centimeters", 1}]
MaterialParameters[n_Association]["Best Raw k"] := SelectBestRange[n] @ QuantityArray[Transpose[{n["Frequencies"] // Normal, n["k debug"] // Normal}], {1/"Centimeters", 1}]
MaterialParameters[n_Association]["Raw Best k"] := SelectBestRange[n] @ QuantityArray[Transpose[{n["Frequencies"] // Normal, n["k debug"] // Normal}], {1/"Centimeters", 1}]

offsetPhase[t_] := With[{
  offset = 2 Pi (1/33.356) QuantityMagnitude[t["\[Delta]t"], "Picoseconds"],
  freqs = Normal[t["Frequencies"]]
},
  offset freqs 
]


SelectBestRange[n_][list_] := With[{ranges = findFDCIRanges[MaterialParameters[n] ]},
  list[[ranges[[1]] ;; ranges[[2]] ]]
]

MaterialParameters[n_Association]["Best n"] := SelectBestRange[n] @ QuantityArray[Transpose[{n["Frequencies"] // Normal, n["n"] // Normal}], {1/"Centimeters", 1}]
MaterialParameters[n_Association]["Best k"] := SelectBestRange[n] @ QuantityArray[Transpose[{n["Frequencies"] // Normal, n["k"] // Normal}], {1/"Centimeters", 1}]

MaterialParameters[n_Association]["Best \[Alpha]"] := SelectBestRange[n] @ QuantityArray[{#[[1]], (#[[2]] 4 \[Pi]  10^12 #[[1]])/(33.356 2.9979 10^10)} &/@ Transpose[{n["Frequencies"] // Normal, n["k"] // Normal}], {1/"Centimeters", 1/"Centimeters"}]
MaterialParameters[n_Association]["Best Raw \[Alpha]"] := SelectBestRange[n] @ QuantityArray[{#[[1]], (#[[2]] 4 \[Pi]  10^12 #[[1]])/(33.356 2.9979 10^10)} &/@ Transpose[{n["Frequencies"] // Normal, n["k debug"] // Normal}], {1/"Centimeters", 1/"Centimeters"}]
MaterialParameters[n_Association]["Raw Best \[Alpha]"] := MaterialParameters[n]["Best Raw \[Alpha]"]

MaterialParameters[n_Association]["\[Alpha]"] := QuantityArray[{#[[1]], (#[[2]] 4 \[Pi]  10^12 #[[1]])/(33.356 2.9979 10^10)} &/@ Transpose[{n["Frequencies"] // Normal, n["k"] // Normal}], {1/"Centimeters", 1/"Centimeters"}]

MaterialParameters[n_Association]["Raw \[Alpha]"] := QuantityArray[{#[[1]], (#[[2]] 4 \[Pi]  10^12 #[[1]])/(33.356 2.9979 10^10)} &/@ Transpose[{n["Frequencies"] // Normal, n["k debug"] // Normal}], {1/"Centimeters", 1/"Centimeters"}]


MaterialParameters[n_Association]["Frequencies"] := QuantityArray[n["Frequencies"] // Normal, 1/"Centimeters"]
MaterialParameters[n_Association]["Domain"] := Quantity[#, 1/"Centimeters"] &/@ MinMax[Normal @ n["Frequencies"]]


MaterialParameters[a_Association]["FDCI2"] := 
  With[{phase = QuantityMagnitude[MaterialParameters[a]["Phase"], {1/"Centimeters", 1}] }, 
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
                2
              ], 
              4
            ]
           , 
            MovingAverage[
              Drop[
                Fourier[
                  (mFp[[1]]["k"] // Normal)[[range[[1]] ;; range[[2]]]]
                ] // Abs // dropHalf, 
                2
              ], 
              4
            ]
        }
      ]}, 
      
      With[{r = # & /@ testRegion}, 
        ListLinePlot[r, PlotRange->Full, PlotLegends->{"FP", "no FP"}]
      ]
    ]
  ]

MaterialParameters[a_Association]["FPReduction"] := 
  With[{mFp = MaterialParameters[a], freqs = MaterialParameters[a][[1]]["Frequencies"]//Normal}, 
    With[{testRegion = 
      With[{range = findFDCIRanges[mFp]}, 
        {
          
            MovingAverage[
              Drop[
                Fourier[
                  (freqs[[range[[1]] ;; range[[2]]]]) (mFp[[1]]["k debug"] // Normal)[[range[[1]] ;; range[[2]]]]
                ] // Abs // dropHalf, 
                2
              ], 
              4
            ] // Interpolation, 
          
            MovingAverage[
              Drop[
                Fourier[
                  (freqs[[range[[1]] ;; range[[2]]]]) (mFp[[1]]["k"] // Normal)[[range[[1]] ;; range[[2]]]]
                ] // Abs // dropHalf, 
                2
              ], 
              4
            ] // Interpolation
        }
      ]}, 
      
      With[{r = NIntegrate[ #[x], {x, #["Domain"][[1,1]], #["Domain"][[1,2]]}] &/@ testRegion}, 
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

Options[validateOptions] = Options[MaterialParameters]
{"Target"->"CPU", "Tags"-><||>, "Model"->"Slab", "SolveNKEquations"->True, "MovingAverageFilter"->True, "FabryPerotCancellation"->True, "FabryPerotCycles"->8, "NKCycles"->30}

validateOptions[OptionsPattern[] ] := With[{},
  If[Or[
    !BooleanQ[OptionValue["SolveNKEquations"] ], 
    (!BooleanQ[OptionValue["MovingAverageFilter"] ] && !IntegerQ[OptionValue["MovingAverageFilter"] ]),
    !BooleanQ[OptionValue["FabryPerotCancellation"] ],
    !NumericQ[OptionValue["FabryPerotCycles"] ],
    !NumericQ[OptionValue["NKCycles"] ]
  ], 
    Message[MaterialParameters::invalidopts];
    $Failed &
  ,
    Identity
  ]
]


MaterialParameters::badgpu = "GPU acceleration is not available. Reason: ``";

materialParameters[l: List[__TransmissionObject], "Slab", "GPU", opts: OptionsPattern[] ] := Module[{dump = {}},
  If[!clReadyQ[], 
    Message[MaterialParameters::badgpu, "OpenCLQ[] returned False. Fallback to CPU"];
    materialParameters[l, "Slab", "CPU", opts]
  ,
    
    (* Group them by Phase (to avoid copying) *)
    With[{groups = SplitBy[l, Function[item,
      item[[1]]["Phase"]
    ] ]},
      With[{ results = Map[Function[group,
        With[{
          freqs = Normal[group[[1,1, "Frequencies"]] ]
        },
          With[{ (* take the first one only, cuz it will be the same for the rest *)
            src = clLoad[Transpose[{freqs, Normal[group[[1, 1, "T"]] ], Normal[group[[1, 1,"Phase"]] ]}] // Flatten, "Float"],
            dest = clLoad[Table[0., {Length[freqs] 5 Length[group] }], "Float"],
            meta = With[{rawMeta = generateMeta /@ group},
              clLoad[rawMeta // Flatten, "Float"]
            ]
          },
            AppendTo[dump, src];
            AppendTo[dump, dest];
            AppendTo[dump, meta // Last];

            materialParameters[group, dest, src, meta, {Length[freqs], Length[group]}, "Slab", "GPU", opts]
          ]
        ]
      ], groups]},

          With[{computed = ReleaseHold /@ results},
            clUnload /@ dump;
            ClearAll[dump];

            Flatten[computed, 1]
          ]

        
    ]
    ]
  ]
]

materialParameters[l: List[__TransmissionObject], "Slab", "CPU", opts: OptionsPattern[] ] := Module[{dump = {}},
    (* Group them by Phase (to avoid copying) *)
    With[{groups = SplitBy[l, Function[item,
      item[[1]]["Phase"]
    ] ]},
      With[{ results = Map[Function[group,
        With[{
          freqs = Normal[group[[1,1, "Frequencies"]] ]
        },
          With[{ (* take the first one only, cuz it will be the same for the rest *)
            src = llLoad[Transpose[{freqs, Normal[group[[1, 1, "T"]] ], Normal[group[[1, 1,"Phase"]] ]}] // Flatten, "Float"],
            dest = llLoad[Table[0., {Length[freqs] 5 Length[group] }], "Float"],
            meta = With[{rawMeta = generateMeta /@ group},
              llLoad[rawMeta // Flatten, "Float"]
            ]
          },
            AppendTo[dump, src];
            AppendTo[dump, dest];
            AppendTo[dump, meta // Last];

            materialParameters[group, dest, src, meta, {Length[freqs], Length[group]}, "Slab", "CPU", opts]
          ]
        ]
      ], groups]},

          With[{computed = ReleaseHold /@ results},
            llUnload /@ dump;
            ClearAll[dump];

            Flatten[computed, 1]
          ]

        
    ]
    ]
]




generateMeta[t: TransmissionObject[a_] ] := With[{
  thickness = QuantityMagnitude[a["Thickness"], "Centimeters"] // N,
  \[Delta]t = QuantityMagnitude[a["\[Delta]t"], "Picoseconds"] // N,
  size = a["Size"]
},
  {\[Delta]t, thickness, 1.0 a["Gain"], 2Pi Round[a["PhaseShift"] ] , size} // N
]


materialParameters[list: List[__TransmissionObject], destCl_, srcCl_, metaCl_, {itemSize_, groupSize_}, "Slab", "GPU", opts: OptionsPattern[] ] := Module[{
  array = list
}, With[{

 NKCycles = If[OptionValue["SolveNKEquations"], OptionValue["NKCycles"], 0], 
 FabryPerotCycles = If[OptionValue["FabryPerotCancellation"] === True, 
   OptionValue["FabryPerotCycles"],
   0
 ],

 MovingAverageFilter = If[MatchQ[OptionValue["MovingAverageFilter"], _Integer], 
  OptionValue["MovingAverageFilter"],
  If[OptionValue["MovingAverageFilter"], 1, 0]
 ]
},

  clRun[destCl, srcCl, metaCl, itemSize, groupSize, NKCycles, MovingAverageFilter,  FabryPerotCycles, groupSize 256];



  With[{
    dest = destCl,
    extraTags = OptionValue["Tags"],
    partition = itemSize
  }, Hold[
    With[{
      memory = (Transpose[ Partition[#, 5] ] &/@ Partition[clGet[destCl], 5 partition])
    },
      MapIndexed[Function[{item, itemSection}, With[{
        itemNumber = itemSection[[1]],
        dt = QuantityMagnitude[item[[1, "\[Delta]t"]], "Picoseconds"] // N,
        thickness = QuantityMagnitude[item[[1, "Thickness"]], "Centimeters"] // N
      },
        MaterialParameters[<|
            "n0" -> 1.0 + (0.029979 dt / thickness),
            "n" -> NumericArray[memory[[itemNumber, 1]]],
            "k" -> NumericArray[memory[[itemNumber, 2]]],
            "T" -> NumericArray[memory[[itemNumber, 3]]],
            "k debug" -> NumericArray[memory[[itemNumber, 5]]],
            "Phase" -> NumericArray[memory[[itemNumber, 4]]],
            "Frequencies" -> item[[1, "Frequencies"]],
            "Thickness" -> item[[1, "Thickness"]],
            "Gain" -> item[[1, "Gain"]],
            "FDCI" -> item[[1, "FDCI"]],
            "PhaseShift" -> item[[1, "PhaseShift"]],
            "Date" -> Now,
            "Tags" -> Join[item["Tags"], extraTags]
          |>]
        ]          
      ], array ] ]
    ]
  ] 
] ]


materialParameters[list: List[__TransmissionObject], destCl_, srcCl_, metaCl_, {itemSize_, groupSize_}, "Slab", "CPU", opts: OptionsPattern[] ] := Module[{
  array = list
}, With[{

 NKCycles = If[OptionValue["SolveNKEquations"], OptionValue["NKCycles"], 0], 
 FabryPerotCycles = If[OptionValue["FabryPerotCancellation"] === True, 
   OptionValue["FabryPerotCycles"],
   0
 ],

 MovingAverageFilter = If[MatchQ[OptionValue["MovingAverageFilter"], _Integer], 
  OptionValue["MovingAverageFilter"],
  If[OptionValue["MovingAverageFilter"], 1, 0]
 ]
},

  llRun[destCl, srcCl, metaCl, itemSize, groupSize, NKCycles, MovingAverageFilter,  FabryPerotCycles, groupSize 256];



  With[{
    dest = destCl,
    extraTags = OptionValue["Tags"],
    partition = itemSize
  }, Hold[
    With[{
      memory = (Transpose[ Partition[#, 5] ] &/@ Partition[llGet[destCl], 5 partition])
    },
      MapIndexed[Function[{item, itemSection}, With[{
        itemNumber = itemSection[[1]],
        dt = QuantityMagnitude[item[[1, "\[Delta]t"]], "Picoseconds"] // N,
        thickness = QuantityMagnitude[item[[1, "Thickness"]], "Centimeters"] // N
      },
        MaterialParameters[<|
            "n0" -> 1.0 + (0.029979 dt / thickness),
            "n" -> NumericArray[memory[[itemNumber, 1]]],
            "k" -> NumericArray[memory[[itemNumber, 2]]],
            "T" -> NumericArray[memory[[itemNumber, 3]]],
            "k debug" -> NumericArray[memory[[itemNumber, 5]]],
            "Phase" -> NumericArray[memory[[itemNumber, 4]]],
            "Frequencies" -> item[[1, "Frequencies"]],
            "Thickness" -> item[[1, "Thickness"]],
            "Gain" -> item[[1, "Gain"]],
            "FDCI" -> item[[1, "FDCI"]],
            "PhaseShift" -> item[[1, "PhaseShift"]],
            "Date" -> Now,
            "Tags" -> Join[item["Tags"], extraTags]
          |>]
        ]          
      ], array ] ]
    ]
  ] 
] ]

findFDCIRanges[mFp_MaterialParameters] := 
  With[{
    m = QuantityMagnitude[mFp["Frequencies"], 1/"Centimeters"] ,
    left = QuantityMagnitude[mFp["FDCI"][[1]], 1/"Centimeters"],
    right = QuantityMagnitude[mFp["FDCI"][[2]], 1/"Centimeters"]
  }, 
    {
      FirstPosition[m, x_ /; (x > left)] // First,
      Length[m] - FirstPosition[m // Reverse, x_ /; (x < right)] // First
    }
  ];



MaterialParameters /: Snippet[m_MaterialParameters] := With[{
  n = QuantityMagnitude[m["n"], {1/"Centimeters", 1}],
  k = QuantityMagnitude[m["\[Alpha]"], {1/"Centimeters", 1/"Centimeters"}],
  fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ m["FDCI"]
},

With[{
  nPlot = ListLinePlot[
    {#[[1]], Clip[#[[2]], {1,Infinity}]} &/@ Select[n, Function[item, item[[1]] > fdci[[1]] / E && item[[1]] < E fdci[[2]] ]], 
    PlotRange->Full, Frame->True, FrameLabel->{"wavenumber (1/cm)", ""}
  , PlotStyle->ColorData[97][2], Prolog->{
    Opacity[0.2], Green, Rectangle[{fdci[[1]], Min[n[[All,2]]]}, {fdci[[2]], Max[n[[All,2]]]}]
  }],
  
  kPlot = ListLinePlot[
    {#[[1]], Clip[#[[2]], {-20,Infinity}]} &/@ Select[k, Function[item, item[[1]] > fdci[[1]] / E && item[[1]] < E fdci[[2]] ]], 
    PlotRange->Full, Frame->True, FrameLabel->{"wavenumber (1/cm)", "absorption coefficient (1/cm)"}
  , Prolog->{
    Opacity[0.2], Green, Rectangle[{fdci[[1]], 0}, {fdci[[2]], Max[k[[All,2]]]}]
  }]
},
  Panel[Row[{kPlot, nPlot}]]
] ]


End[]
EndPackage[]