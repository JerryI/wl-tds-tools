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
          "n0" -> 1.0 + (0.029979 (t0Sam-t0Ref) / thickness),
          "Gain"->gain,
          "FDCI"->sam["FDCI"],
          "Size" -> Length[freqs],
          "PhaseShift"->0,
          "Date" ->Now,
          "Tags" -> OptionValue["Tags"],
          "Phase"->NumericArray[   Arg[Exp[I (Arg[fsam[[All,2]]] - Arg[fref[[All,2]]]) ] Exp[ -I pDiff ] ] + pDiff ], 
          "Frequencies"->NumericArray[freqs], (Sequence @@ Options[TransmissionObject]), 
          opts
        |>]
        ]
      ]
    ]
  ]
]]

(* by Thies Heidecke see https://mathematica.stackexchange.com/questions/341/implementing-discrete-and-continuous-hilbert-transforms *)
HilbertSpectrum[0] := {};
HilbertSpectrum[n_Integer?Positive] := With[{nhalf = Quotient[n - 1, 2]},
  Join[{0}, ConstantArray[-I, nhalf], If[EvenQ[n], {0}, {}]
          , ConstantArray[ I, nhalf]
  ]
];

Hilbert[data_?VectorQ, padding_Integer?NonNegative] := Module[
  {fp = FourierParameters -> {1, -1}, n = Length[data], m, paddeddata},
  m = n + padding;
  paddeddata = PadRight[data, m];
  Re @ InverseFourier[ HilbertSpectrum[m] Fourier[paddeddata, fp], fp][[;;n]]
]

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
          {BoxForm`SummaryItem[{"n0 ", obj["n0"]}]},
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

TransmissionObject[a_Association]["Kramers-Kronig n"] := With[{
  n0 = a["n0"],
  thickness = QuantityMagnitude[a["Thickness"], "Centimeters"],
  freqs = Normal @ a["Frequencies"]
},
  With[{
    im\[Chi] =  (- (*FB[*)((1)(*,*)/(*,*)(#[[1]] thickness))(*]FB*) Log[#[[2]]]) &/@ Drop[QuantityMagnitude[TransmissionObject[a]["Transmission"], {1/"Centimeters", 1} ], 1]
  },
    QuantityArray[{freqs, Join[{n0}, (*SqB[*)Sqrt[n0^2 - Hilbert[im\[Chi], 0] ](*]SqB*)]} // Transpose, {1/"Centimeters", 1}]
  ]
]

TransmissionObject[a_Association]["Transmission"] := QuantityArray[Transpose[{Normal @ a["Frequencies"], (*SpB[*)Power[(a["Gain"] Normal @ a["T"])(*|*),(*|*)2](*]SpB*)}], {1/"Centimeters", 1}]

updateThicknessDependent[a_Association ] := With[{

},
  Join[a, <|
    "n0" -> 1.0 + (0.029979 QuantityMagnitude[a["\[Delta]t"], "Picoseconds" ] / QuantityMagnitude[a["Thickness"], "Centimeters"]),
    "Date" ->Now
  |>] 
]

TransmissionObject /: Append[TransmissionObject[a_Association], props_Association] := TransmissionObject[Join[a, props] // updateThicknessDependent ]
TransmissionObject /: Append[TransmissionObject[a_Association], prop_Rule] := TransmissionObject[Append[a, prop] // updateThicknessDependent ]
TransmissionObject /: Append[TransmissionObject[a_Association], props_List] := TransmissionObject[Append[a, props] // updateThicknessDependent ]

TransmissionObject[a_Association]["Frequencies"] := QuantityArray[Normal @ a["Frequencies"], 1/"Centimeters"]

TransmissionObject[a_Association]["Domain"] := Quantity[#, 1/"Centimeters"] &/@ MinMax[Normal @ a["Frequencies"]]


offsetPhase[t_] := With[{
  offset = 2 Pi (1/33.356) QuantityMagnitude[t["\[Delta]t"], "Picoseconds"],
  freqs = Normal[t["Frequencies"]]
},
  offset freqs 
]

TransmissionObject[a_Association]["Phase"] := With[{
  shift = 2 Pi (1/33.356) QuantityMagnitude[a["\[Delta]t"], "Picoseconds"] Normal[a["Frequencies"]],
  off = a["PhaseShift"]
},
  QuantityArray[Transpose[{Normal @ a["Frequencies"], (Normal @ a["Phase"])}], {1/"Centimeters", 1}]
]

TransmissionObject[a_Association]["Phase Features"] := With[{
  offset = offsetPhase[ a ],
  raw = Normal[ a["Phase"] ]
},
  QuantityArray[Transpose[{Normal @ a["Frequencies"], raw - offset }], {1/"Centimeters", 1}]
]


TransmissionObject[a_Association]["Approximated k"] := With[{
  thickness = QuantityMagnitude[a["Thickness"], "Centimeters"],
  gain = a["Gain"]
},
  QuantityArray[
    Join[{0, 0}, {#[[1]], - 0.159152 Log[#[[2]] gain ] / (#[[1]] thickness) } &/@ Drop[Transpose[Normal /@ {a["Frequencies"], a["T"]}], 1] ]
  , {1/"Centimeters", 1}]
]

TransmissionObject[a_Association]["Approximated \[Alpha]"] := With[{
  thickness = QuantityMagnitude[a["Thickness"], "Centimeters"],
  gain = a["Gain"]
},
  QuantityArray[
    Join[{0, 0}, {#[[1]], ((- 0.159152 Log[#[[2]] gain ] / (#[[1]] thickness)) 4 \[Pi]  10^12 #[[1]])/(33.356 2.9979 10^10)} &/@ Drop[Transpose[Normal /@ {a["Frequencies"], a["T"]}], 1] ]
  , {1/"Centimeters", 1/"Centimeters"}]
]


TransmissionObject[a_Association]["Approximated n"] := With[{
  shift = a["PhaseShift"],
  thickness = QuantityMagnitude[a["Thickness"], "Centimeters"],
  n0 = a["n0"]
},
  QuantityArray[
    Join[{0, n0}, {#[[1]], 1.0 + 0.159152 (#[[2]] + shift) / (#[[1]] thickness) } &/@ Drop[Transpose[Normal /@ {a["Frequencies"], a["Phase"]}],1] ]
  , {1/"Centimeters", 1}]
]



TransmissionObject[a_Association]["FrequencyDomainConfidenceInterval"] := TransmissionObject[a]["FDCI"]

TransmissionObject[a_Association]["FrequencyDomainConfidenceInterval2"] := TransmissionObject[a]["FDCI2"]


TransmissionObject[a_Association]["FDCI2"] := 
With[{phase = QuantityMagnitude[TransmissionObject[a]["Phase"], {1/"Centimeters", 1}]},
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

TransmissionObject /: Keys[t_TransmissionObject] :=  t["Properties"]

TransmissionObject[a_Association]["Properties"] := Join[Options[TransmissionObject][[All,1]], {"Properties", "n0", "Approximated n", "Approximated k", "Approximated \[Alpha]", "Kramers-Kronig n", "Frequencies", "Transmission", "Phase", "Phase Features", "\[Delta]t", "Gain", "PhaseShift", "Thickness", "Domain", "FrequencyDomainConfidenceInterval", "FDCI", "FDCI2", "FrequencyDomainConfidenceInterval2"}]

Options[TransmissionObject] = {"Thickness"->Null, "Tags"-><||>, "Gain"->1.0, "PhaseShift"->0};

TransmissionObject::invalidopts = "Invalid options provided"

Options[validateOptions] = Options[TransmissionObject] 
validateOptions[OptionsPattern[] ] := With[{},
  If[Or[
    !QuantityQ[OptionValue["Thickness"] ], 
    !NumericQ[OptionValue["Gain"] ],
    !NumericQ[OptionValue["PhaseShift"] ]
  ], 
    Message[TransmissionObject::invalidopts];
    $Failed &
  ,
    Identity
  ]
]

TransmissionUnwrap[t: TransmissionObject[a_], "Basic", OptionsPattern[]] := With[{
  offset = offsetPhase[a],
  th = OptionValue["PhaseThreshold"]//N,
  phaseShift = OptionValue["PhaseShift"]
},
  If[!NumericQ[OptionValue["PhaseThreshold"] ] || !NumericQ[OptionValue["PhaseShift"] ],
    Message[TransmissionObject::invalidopts];
    Return[$Failed];
  ];

  With[{unwrapped = With[{
    phase = Normal[a["Phase"] ]
  },
    clusterPhase[phase - offset, 1, Length[phase]-1, th] + offset
  ]},
    Append[t, {"Phase"->NumericArray[unwrapped], "PhaseShift"->phaseShift}]
  ]
]

Options[TransmissionUnwrap] = {"PhaseThreshold"->5.6, "PhaseShift"->0};

End[]
EndPackage[]