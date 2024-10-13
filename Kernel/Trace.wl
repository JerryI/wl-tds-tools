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

ClearAll[TDTrace]

tdValid[_, OptionsPattern[]] := False

TDTrace::input = "Input data must be a QuantityArray or List[List[]..]";
TDTrace::unknownprop = "Unknown property `1`";

TDTrace[td_TDTrace, opts: OptionsPattern[]] := With[{
  units = QuantityMagnitude[OptionValue["Units"], "Picoseconds"],
  amp   = OptionValue["Gain"],
  offset = QuantityMagnitude[OptionValue["Offset"], "Picoseconds"]
},
  TDTrace[NumericArray[{#[[1]] units + offset, #[[2]] amp} &/@ Normal[td[[1]]]]]
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
TDTrace[n_NumericArray]["Spectrum"] := QuantityArray[fourier2d[TDTrace[n]["Trace"] // QuantityMagnitude], {1/"Centimeters", 1}]


TDTrace[n_NumericArray]["FrequencyDomainConfidenceInterval"] := Module[{ model, x0, \[Sigma], A},
  With[{d = TDTrace[n]["PowerSpectrum"]//QuantityMagnitude}, 
    model = NonlinearModelFit[Drop[d,10], A Exp[-(*FB[*)(((*SpB[*)Power[(x0 - x)(*|*),(*|*)2](*]SpB*))(*,*)/(*,*)((*SpB[*)Power[\[Sigma](*|*),(*|*)2](*]SpB*)))(*]FB*)], {\[Sigma], A, x0}, x, ConfidenceLevel->0.5, Method->"NMinimize"];
    model = Association[model["BestFitParameters"]];
    Quantity[#, 1/"Centimeters"] &/@ {Clip[model[x0] - Abs[model[\[Sigma]]], {5.0, Infinity}], 1.2 Clip[model[x0] + Abs[model[\[Sigma]]], {30.0, Infinity}]}
  ]
]

TDTrace[n_NumericArray]["FDCI"] := TDTrace[n]["FrequencyDomainConfidenceInterval"]

TDTrace[n_NumericArray]["PowerSpectrum"] := QuantityArray[{#[[1]], Abs[#[[2]]]^2}&/@fourier2d[TDTrace[n]["Trace"] // QuantityMagnitude], {1/"Centimeters", 1}]

TDTrace[n_NumericArray]["Properties"] := {"Properties", "Spectrum", "PowerSpectrum", "Trace", "FrequencyDomainConfidenceInterval", "FDCI"}

TDTrace[q_QuantityArray, opts: OptionsPattern[]] := TDTrace[QuantityMagnitude[q, {"Picoseconds", 1}], opts]

TDTrace[td: List[__List], opts: OptionsPattern[]] := With[{
  units = QuantityMagnitude[OptionValue["Units"], "Picoseconds"],
  amp   = OptionValue["Gain"],
  offset = QuantityMagnitude[OptionValue["Offset"], "Picoseconds"]
},

  TDTrace[NumericArray[{#[[1]] units + offset, #[[2]] amp} &/@ td]]
  
] /; tdValid[td, opts]

Options[TDTrace] = {"Units"->Quantity[1, "Picoseconds"], "Gain"->1.0, "Offset"->Quantity[0, "Picoseconds"]};

Options[tdValid] = Options[TDTrace];

TDTrace /: Keys[t_TDTrace] :=  t["Properties"]

End[]
EndPackage[]