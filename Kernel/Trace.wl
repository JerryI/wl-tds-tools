BeginPackage["JerryI`TDSTools`Trace`", {
  "JerryI`TDSTools`Utils`"
}]

TDTrace::usage = "TDTrace[trace_QuantityArray] creates a time-trace object from table provided as trace. See _TDTrace[\"Properties\"]"

Begin["`Private`"]

tdValid[td: List[__List], OptionsPattern[]] := With[{
  q = QuantityMagnitude[OptionValue["Units"]]
},
  If[ListQ[td] && ListQ[td // First] && NumericQ[q],
    With[{a = ((Transpose[td] {q, 1}))},
      AllTrue[NumericQ][a[[1]]] && AllTrue[NumericQ][a[[2]]]
    ]
  ,
    Message[TDTrace::input];
    
    False
  ]
]

freqCut[test_] := 0.2008 (test[[{1,-1},1]]//Differences)/Length[test];
slowRate[x_] := Clip @ If[x > 0.8,  With[{s = Exp[-3(x-0.8)]},  s x + (1-s) Tanh[ x]], x]
transition[x_, y1_, y2_, y3_] := With[{j = If[x > 1.0, Exp[-3.0(x-1.0)], 1.0]}, With[{s = If[x <= 0.8, 1.0, Exp[-8.0(x-0.8)]]}, (1-s) y2 + s y1] j + y3 (1.0-j)];


pad["BaselineDecay", n_, test_List] := With[{int = Interpolation[test, InterpolationOrder->1], base = Interpolation[Transpose[{test[[All,1]], LowpassFilter[test[[All,2]], freqCut[test]]}], InterpolationOrder->1]},
  With[{min = test[[1,1]], max = test[[-1,1]], step = Abs[Min[Differences[test[[All,1]] ] ] ]},
    Table[{i, transition[
      (i-min)/(max - min),
      int[min + (max - min) slowRate[(i-min)/(max - min)]],
      base[min + (max - min) slowRate[(i-min)/(max - min)]],
      test[[1,2]]
    ]
      
    }, {i, min, max + n(max-min), step}]
  ]
]

pad["ZerosDecay", n_, test_List] := With[{int = Interpolation[test, InterpolationOrder->1], base = Interpolation[Transpose[{test[[All,1]], LowpassFilter[test[[All,2]], freqCut[test]]}], InterpolationOrder->1]},
  With[{min = test[[1,1]], max = test[[-1,1]], step = Abs[Min[Differences[test[[All,1]] ] ] ]},
    Table[{i, With[{x = (i-min)/(max - min)},
      If[x > 0.8, 
        If[i <= max, int[i], int[max]] Exp[-50.0(x - 0.8)^2] 
      , 
        int[i]
      ]
    ]
      
    }, {i, min, max + n(max-min), step}]
  ]
]

pad["Zeros", n_, test_List] := With[{int = Interpolation[test, InterpolationOrder->1], base = Interpolation[Transpose[{test[[All,1]], LowpassFilter[test[[All,2]], freqCut[test]]}], InterpolationOrder->1]},
  With[{min = test[[1,1]], max = test[[-1,1]], step = Abs[Min[Differences[test[[All,1]] ] ] ]},
    Table[{i, If[i < max,
      int[i]
    ,
      0
    ]
      
    }, {i, min, max + n(max-min), step}]
  ]
]

pad[test_List, None] := test;
pad[test_List, n_Integer] := pad["BaselineDecay", n, test]
pad[test_List, {name_String, n_Integer}] := pad[name, n, test]
pad[test_List, _] := test


ClearAll[TDTrace]

tdValid[_, OptionsPattern[]] := False

TDTrace::input = "Input data must be a QuantityArray or List[List[]..]";
TDTrace::unknownprop = "Unknown property `1`";

TDTrace[td_TDTrace, opts: OptionsPattern[]] := With[{
  units = QuantityMagnitude[OptionValue["Units"], "Picoseconds"],
  amp   = OptionValue["Gain"],
  offset = QuantityMagnitude[OptionValue["Offset"], "Picoseconds"],
  padding = OptionValue["PadZeros"]
},
  TDTrace[NumericArray[{#[[1]] units + offset, #[[2]] amp} &/@ pad[Normal[td[[1]]], padding] ] ]
]

TDTrace /: MakeBoxes[obj: TDTrace[n_NumericArray], StandardForm] := With[{
  preview = ListLinePlot[ArrayResample[n // Normal, 300], PlotRange->Full,Axes -> None, ImagePadding->None]
},
    Module[{above, below},
        above = { 
          {BoxForm`SummaryItem[{"Length ", Length[n]}]},
          {BoxForm`SummaryItem[{"From ", Quantity[Part[n,1,1], "Picoseconds"] }]},
          {BoxForm`SummaryItem[{"To ", Quantity[Part[n,-1,1], "Picoseconds"] }]}
        };

        BoxForm`ArrangeSummaryBox[
           TDTrace, (* head *)
           obj,      (* interpretation *)
           preview,    (* icon, use None if not needed *)
           (* above and below must be in a format suitable for Grid or Column *)
           above,    (* always shown content *)
           Null (* expandable content. Currently not supported!*)
        ]
    ]  
]

TDTrace[n_NumericArray][s_String]    := (Message[TDTrace::unknownprop, s]; $Failed)
TDTrace[n_NumericArray]["Trace"]    := QuantityArray[n // Normal, {"Picoseconds", 1}]
TDTrace[n_NumericArray]["Spectrum"] := QuantityArray[fourier2d[ QuantityMagnitude[TDTrace[n]["Trace"], {"Picoseconds", 1} ] ], {1/"Centimeters", 1}]


TDTrace[n_NumericArray]["FrequencyDomainConfidenceInterval"] := Module[{ model, x0, \[Sigma], A},
  With[{d = QuantityMagnitude[TDTrace[n]["PowerSpectrum"], {1/"Centimeters", 1}]}, 
    model = NonlinearModelFit[Drop[d,10], A Exp[-(*FB[*)(((*SpB[*)Power[(x0 - x)(*|*),(*|*)2](*]SpB*))(*,*)/(*,*)((*SpB[*)Power[\[Sigma](*|*),(*|*)2](*]SpB*)))(*]FB*)], {\[Sigma], A, x0}, x, ConfidenceLevel->0.5, Method->"NMinimize"];
    model = Association[model["BestFitParameters"]];
    Quantity[#, 1/"Centimeters"] &/@ {Clip[model[x0] - Abs[model[\[Sigma]]], {5.0, Infinity}], 1.2 Clip[model[x0] + Abs[model[\[Sigma]]], {30.0, Infinity}]}
  ]
]

TDTrace[n_NumericArray]["FDCI"] := TDTrace[n]["FrequencyDomainConfidenceInterval"]

TDTrace[n_NumericArray]["PowerSpectrum"] := QuantityArray[{#[[1]], Abs[#[[2]]]^2}&/@fourier2d[ QuantityMagnitude[TDTrace[n]["Trace"], {"Picoseconds", 1} ] ], {1/"Centimeters", 1}]

TDTrace[n_NumericArray]["Properties"] := {"Properties", "Spectrum", "PowerSpectrum", "Trace", "FrequencyDomainConfidenceInterval", "FDCI"}

TDTrace[q_QuantityArray, opts: OptionsPattern[]] := TDTrace[QuantityMagnitude[q, {"Picoseconds", 1}], opts]

TDTrace[td: List[__List], opts: OptionsPattern[] ] := With[{
  units = QuantityMagnitude[OptionValue["Units"], "Picoseconds"],
  amp   = OptionValue["Gain"],
  offset = QuantityMagnitude[OptionValue["Offset"], "Picoseconds"],
  padding = OptionValue["PadZeros"]
},

  TDTrace[NumericArray[{#[[1]] units + offset, #[[2]] amp} &/@ pad[td, padding] ] ]
  
] /; tdValid[td, opts]

Options[TDTrace] = {"Units"->Quantity[1, "Picoseconds"], "PadZeros"->None, "Gain"->1.0, "Offset"->Quantity[0, "Picoseconds"]};

Options[tdValid] = Options[TDTrace];

TDTrace /: Keys[t_TDTrace] :=  t["Properties"]

End[]
EndPackage[]