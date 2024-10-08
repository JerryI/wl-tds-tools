BeginPackage["JerryI`TDSTools`Transmission`", {
  "JerryI`TDSTools`Utils`",
  "JerryI`TDSTools`Trace`"
}]

TransmissionObject::usage = "TransmissionObject[sam_TDTrace, ref_TDTrace, opts___] creates a transmission from sample and reference traces. See _TransmissionObject[\"Properties\"]"
TransmissionUnwrap::usage = "TransmissionUnwrap[t_TransmissionObject, opts___] produces a new transsmission object with unwrapped phase"

Begin["`Private`"]

root = DirectoryName[$InputFileName];

{ 
  initalize, 
  solveNK, 
  solveFP, 
  movingAverage,
  saveForDebug,
  clusterPhase
} = Get[ FileNameJoin[{root, "nCPU.wl"}] ];

TransmissionObject::thickerr = "Thickness `1` is not valid";

TransmissionObject[sam_TDTrace, ref_TDTrace, opts: OptionsPattern[]] := Module[{}, With[{
  thickness = QuantityMagnitude[OptionValue["Thickness"], "Centimeters"],
  gain = OptionValue["Gain"]
},
  If[!NumericQ[thickness], Message[TransmissionObject::thickerr, OptionValue["Thickness"]]; Return[$Failed]];

  With[{
    minStep = Min @ ((Min @ Differences @ (QuantityMagnitude[#["Trace"], {"Picoseconds", 1}][[All,1]])) &/@ {sam, ref}),
    min = Min @ ((QuantityMagnitude[#["Trace"], {"Picoseconds", 1}][[1,1]]) &/@ {sam, ref}),
    max = Max @ ((QuantityMagnitude[#["Trace"], {"Picoseconds", 1}][[-1,1]]) &/@ {sam, ref}),
    t0Ref = MaximalBy[QuantityMagnitude[ref["Trace"], {"Picoseconds", 1}], Last]//First//First,
    t0Sam = MaximalBy[QuantityMagnitude[sam["Trace"], {"Picoseconds", 1}], Last]//First//First
  },
    With[{
      fsam = fourier2d[QuantityMagnitude[sam["Trace"], {"Picoseconds", 1}], {min, max, minStep}], 
      fref = fourier2d[QuantityMagnitude[ref["Trace"], {"Picoseconds", 1}], {min, max, minStep}]
    }, 
      With[{
        t = Abs[fsam[[All,2]] / fref[[All,2]]],
        freqs = fsam[[All,1]]
      },
        With[{
          pDiff = 2 Pi (1/33.356) freqs (t0Sam - t0Ref)
        },
        
        TransmissionObject[<|
          "\[Delta]t"->Quantity[(t0Sam-t0Ref), "Picoseconds"], 
          "T"->NumericArray[t], 
          "Thickness"->Quantity[thickness, "Centimeters"],
          "Gain"->gain,
          "FDCI"->sam["FDCI"],
          "PhaseShift"->0,
          "Date" ->Now,
          "Tags" -> OptionValue["Tags"],
          "Phase"->NumericArray[Arg[fsam[[All,2]]] - Arg[fref[[All,2]]] ], 
          "Frequencies"->NumericArray[freqs], (Sequence @@ Options[TransmissionObject]), 
          opts
        |>]
        ]
      ]
    ]
  ]
]]

TransmissionObject /: MakeBoxes[obj: TransmissionObject[a_Association], StandardForm] := With[{
  preview = ListLinePlot[ArrayResample[Select[obj["Transmission"]//QuantityMagnitude//dropHalf, #[[2]]<1.0 &], 300], PlotStyle->ColorData[97][2], PlotRange->Full,Axes -> None, ImagePadding->None],

  state = phaseState[QuantityMagnitude[obj["Phase"]]]
},
 

    Module[{above, below},
        above = { 
          {BoxForm`SummaryItem[{"From ", Quantity[Part[a["Frequencies"],1], 1/"Centimeters"]//Round }]},
          {BoxForm`SummaryItem[{"To ", Quantity[Part[a["Frequencies"],-1], 1/"Centimeters"]//Round }]},
          {BoxForm`SummaryItem[{"Thickness ", obj["Thickness"]}]},
          {BoxForm`SummaryItem[{"Gain ", obj["Gain"]}]},
          {BoxForm`SummaryItem[{"PhaseShift ", 2Pi obj["PhaseShift"]}]},
          {BoxForm`SummaryItem[{"Phase ", state}]},
          If[Length[obj["Tags"]//Keys] > 0, {BoxForm`SummaryItem[{"Tags ", Style[#, Background->LightBlue]&/@obj["Tags"]//Keys}]}, Nothing]
        };

        BoxForm`ArrangeSummaryBox[
           TransmissionObject, (* head *)
           obj,      (* interpretation *)
           preview,    (* icon, use None if not needed *)
           (* above and below must be in a format suitable for Grid or Column *)
           above,    (* always shown content *)
           Null (* expandable content. Currently not supported!*)
        ]
    ]  
]

TransmissionObject::unknownprop = "Unknown property `1`";

TransmissionObject[a_Association][prop_String] := If[!KeyExistsQ[a, prop],
  (Message[TransmissionObject::unknownprop, prop]; $Failed)
,
  a[prop]
]

TransmissionObject[a_Association]["Transmission"] := QuantityArray[Transpose[{Normal @ a["Frequencies"], (*SpB[*)Power[(a["Gain"] Normal @ a["T"])(*|*),(*|*)2](*]SpB*)}], {1/"Centimeters", 1}]

TransmissionObject /: Append[TransmissionObject[a_Association], props_Association] := TransmissionObject[Join[a, props]]
TransmissionObject /: Append[TransmissionObject[a_Association], prop_Rule] := TransmissionObject[Append[a, prop]]
TransmissionObject /: Append[TransmissionObject[a_Association], props_List] := TransmissionObject[Append[a, props]]

TransmissionObject[a_Association]["Frequencies"] := QuantityArray[Normal @ a["Frequencies"], 1/"Centimeters"]

TransmissionObject[a_Association]["Domain"] := Quantity[#, 1/"Centimeters"] &/@ MinMax[Normal @ a["Frequencies"]]


TransmissionObject[a_Association]["Phase"] := With[{
  shift = 2 Pi (1/33.356) QuantityMagnitude[a["\[Delta]t"], "Picoseconds"] Normal[a["Frequencies"]],
  off = a["PhaseShift"]
},
  QuantityArray[Transpose[{Normal @ a["Frequencies"], (Normal @ a["Phase"])}], {1/"Centimeters", 1}]
]

TransmissionObject[a_Association]["FrequencyDomainConfidenceInterval"] := TransmissionObject[a]["FDCI"]

TransmissionObject[a_Association]["FrequencyDomainConfidenceInterval2"] := TransmissionObject[a]["FDCI2"]


TransmissionObject[a_Association]["FDCI2"] := 
With[{phase = TransmissionObject[a]["Phase"] // QuantityMagnitude},
  With[{r = Table[{phase[[q, 1]], 
      LinearModelFit[Take[phase, q], {1, x}, x]["RSquared"]}, 
     {q, Round[0.2 (phase // Length)], phase // Length, 25}]},
    With[{scape = {Drop[r[[All, 1]], -4], 
       MovingAverage[r[[All, 2]], 4] // Differences} // Transpose},
     Quantity[#, 1/"Centimeters"] &/@ {0.0, scape[[FindPeaks[-scape[[All, 2]], 5][[All, 1]], 1]] // First}
     ]
    ]
  ]



phaseState[phase_List] := With[{},
  If[Fit[phase, {1, x}, x][[1]] > 10.0,
    Item[Style["Unwrapped", Bold], Background->LightGreen]
  ,
    Item[Style["Possibly Wrapped", Bold], Background->LightYellow]
  ]
]

TransmissionObject[a_Association]["Properties"] := Join[Options[TransmissionObject][[All,1]], {"Properties", "Frequencies", "Transmission", "Phase", "\[Delta]t", "Gain", "PhaseShift", "Thickness", "Domain", "FrequencyDomainConfidenceInterval", "FDCI", "FDCI2", "FrequencyDomainConfidenceInterval2"}]

Options[TransmissionObject] = {"Thickness"->Null, "Tags"-><||>, "Gain"->1.0, "PhaseShift"->0};


TransmissionUnwrap[t: TransmissionObject[a_], "Basic", OptionsPattern[]] := With[{
  offset = 2 Pi (1/33.356) QuantityMagnitude[t["\[Delta]t"], "Picoseconds"] ,
  th = OptionValue["PhaseThreshold"]//N,
  freqs = Normal[a["Frequencies"]],
  phaseShift = OptionValue["PhaseShift"]
},
  With[{unwrapped = With[{
    phase = Normal[a["Phase"]]
  },
    clusterPhase[phase - offset freqs, 1, Length[phase]-1, th] + offset freqs
  ]},
    Append[t, {"Phase"->NumericArray[unwrapped], "PhaseShift"->phaseShift}]
  ]
]

Options[TransmissionUnwrap] = {"PhaseThreshold"->5.6, "PhaseShift"->0};

End[]
EndPackage[]