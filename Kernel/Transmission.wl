BeginPackage["JerryI`TDSTools`Transmission`", {
  "JerryI`TDSTools`Utils`",
  "JerryI`TDSTools`Trace`"
}]

TransmissionObject::usage = "TransmissionObject[sam_TDTrace, ref_TDTrace, opts___] creates a transmission from sample and reference traces. See _TransmissionObject[\"Properties\"]"
TransmissionUnwrap::usage = "TransmissionUnwrap[t_TransmissionObject, opts___] produces a new transsmission object with unwrapped phase"

Begin["`Private`"]

root = DirectoryName[$InputFileName];

clusterPhase = 
    Compile[
      {
        {input, _Real, 1}, 
        {low, _Integer, 0}, 
        {high, _Integer, 0}, 
        {PhaseTrashhold, _Real, 0}
      }, 
      Module[{phase = input, delta = 0.0},

        (* Adjust phase values based on PhaseTrashhold *)
        Do[
          If[phase[[i + 1]] - phase[[i]] > PhaseTrashhold,
            Do[phase[[j]] = phase[[j]] - 2 \[Pi];, {j, i + 1, Length[phase]}],

            If[phase[[i + 1]] - phase[[i]] < -PhaseTrashhold,
              Do[phase[[j]] = phase[[j]] + 2 \[Pi];, {j, i + 1, Length[phase]}]
            ]
          ],
          {i, low, high}
        ];

        (* Adjust phase before the 'low' index *)
        If[low > 1,
          delta = Differences[phase[[low ;; Floor[(3 high + low)/4] ]]  ] // Mean;

          Do[
            phase[[-i - 1]] = phase[[-i]] - delta/2,
            {i, -low, -2}
          ];
        ];

        (* Adjust phase after the 'high' index *)
        delta = Differences[phase[[Floor[(high + 3 low)/4] ;; high]]] // Mean;

        Do[
          phase[[i + 1]] = phase[[i]] + delta,
          {i, high, Length[phase] - 1}
        ];

        phase
    ], 
      "CompilationTarget" -> "C", 
      "RuntimeOptions" -> "Speed"
    ];

TransmissionObject::thickerr = "Thickness `1` is not valid";

(* :: Constructor :: *)

(* [TODO] Store Transmission as InterpolationFunction instead! *)

TransmissionObject /: Times[a_TransmissionObject, B_?NumberQ] := With[{
  TA = a[[1, "T"]]  // Normal,
  PhaseA = a[[1, "Phase"]]  // Normal
},

  Append[a, <|
          "T"->NumericArray[ TA Abs[B] ],
          "Phase"->NumericArray[ PhaseA + Arg[B] ]
  |>]

]

resampleSeries[{x_,y_}, {min_, max_, step_}] := With[{
  int = Interpolation[{x,y} // Transpose, InterpolationOrder->1]
},
  Table[{i, int[i]}, {i, min, max, step}]
]

TransmissionObject /: Plus[a_TransmissionObject, b_TransmissionObject] := With[{
  phaseShiftA = a[[1, "PhaseShift"]],
  sizeA = a[[1, "Size"]],
  nA = a[[1, "n0"]],
  gainA = a[[1, "Gain"]],
  TA = a[[1, "T"]] // Normal,
  PhaseA = a[[1, "Phase"]] // Normal,
  FreqsA = a[[1, "Frequencies"]] // Normal,
  delta = QuantityMagnitude[a[[1, "\[Delta]t"]], "Picoseconds"],

  phaseShiftB = b[[1, "PhaseShift"]],
  sizeB = b[[1, "Size"]],
  nB = b[[1, "n0"]],
  gainB = b[[1, "Gain"]],
  TB = b[[1, "T"]] // Normal,
  PhaseB = b[[1, "Phase"]] // Normal,
  FreqsB = b[[1, "Frequencies"]] // Normal
},

  (* resample and apply gain / phase shifts *)
  With[{
    step = Min[Differences[FreqsA]//Abs//Min, Differences[FreqsA]//Abs//Min],
    min = Min[FreqsA//Min, FreqsB//Min],
    max = Max[FreqsA//Max, FreqsB//Max]
  },
    With[{
      resampledTA = {#[[1]], #[[2]] gainA} &/@ resampleSeries[{FreqsA, TA}, {min, max, step}],
      resampledTB = {#[[1]], #[[2]] gainB} &/@ resampleSeries[{FreqsB, TB}, {min, max, step}],
      resampledPA = {#[[1]], #[[2]] + 2Pi phaseShiftA} &/@ resampleSeries[{FreqsA, PhaseA}, {min, max, step}],
      resampledPB = {#[[1]], #[[2]] + 2Pi phaseShiftB} &/@ resampleSeries[{FreqsB, PhaseB}, {min, max, step}]

      
    },
      With[{
        result =  resampledTA[[All,2]] Exp[I resampledPA[[All,2]] ] +  resampledTB[[All,2]] Exp[I resampledPB[[All,2]] ],
        pDiff = 2 Pi (1/33.356) resampledTA[[All,1]] delta
      
      },
        Append[a, <|
          "T"->NumericArray[Abs[result] ],
          "Frequencies" -> NumericArray[resampledTA[[All,1]]],
          "Gain" -> 1,
          "PhaseShift"->0,
          "Size" -> Length[resampledTA],
          "Phase"->NumericArray[ Arg[result Exp[- I pDiff] ] + pDiff ]
        |>]
      ]
    ]
  ]

]

TransmissionObject[sam_TDTrace, ref_TDTrace, opts: OptionsPattern[] ] := Module[{}, With[{
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
] ]

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


(* :: Snippet :: *)
TransmissionObject /: Snippet[m_TransmissionObject] := With[{
  phase = QuantityMagnitude[m["Phase Features"], {1/"Centimeters", 1}],
  transm = QuantityMagnitude[m["Transmission"], {1/"Centimeters", 1}],
  fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ m["FDCI"]
},

With[{
  phasePlot = ListLinePlot[
    Select[phase, Function[item, item[[1]] < E fdci[[2]] ]], 
    PlotRange->Full, Frame->True, FrameLabel->{"wavenumber (1/cm)", "Radians"}
  , PlotStyle->ColorData[97][2], Prolog->{
    Opacity[0.2], Green, Rectangle[{fdci[[1]], Min[phase[[All,2]]]}, {fdci[[2]], Max[phase[[All,2]]]}]
  }],
  
  transmPlot = ListLinePlot[
    {#[[1]], Clip[#[[2]], {0,1}]} &/@ Select[transm, Function[item, item[[1]] < E fdci[[2]] ]], 
    PlotRange->Full, Frame->True, FrameLabel->{"wavenumber (1/cm)", "T"}
  , Prolog->{
    Opacity[0.2], Green, Rectangle[{fdci[[1]], -0.1}, {fdci[[2]], 1.1}]
  }]
},
  Panel[Row[{transmPlot, phasePlot}]]
] ]


(* :: Calculated properties :: *)

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

(* :: Transition methods :: *)

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

(* :: Normal properties :: *)

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


(* :: Approximated properties :: *)

TransmissionObject[a_Association]["Approximated k"] := With[{
  thickness = QuantityMagnitude[a["Thickness"], "Centimeters"],
  gain = a["Gain"]
},
  QuantityArray[
    Join[{{0, 0}}, {#[[1]], - 0.159152 Log[#[[2]] gain ] / (#[[1]] thickness) } &/@ Drop[Transpose[Normal /@ {a["Frequencies"], a["T"]}], 1] ]
  , {1/"Centimeters", 1}]
]

TransmissionObject[a_Association]["Approximated \[Alpha]"] := With[{
  thickness = QuantityMagnitude[a["Thickness"], "Centimeters"],
  gain = a["Gain"],
  n0 = Clip[a["n0"], {0, Infinity}]
},
{
  freshnel = (1/thickness) Log[16 n0 n0 / ((n0 + 1)^4)]
},
  QuantityArray[
    Join[{{0, 0}}, {#[[1]], freshnel + ((- 0.159152 Log[#[[2]] gain ] / (#[[1]] thickness)) 4 \[Pi]  10^12 #[[1]])/(33.356 2.9979 10^10)} &/@ Drop[Transpose[Normal /@ {a["Frequencies"], a["T"]}], 1] ]
  , {1/"Centimeters", 1/"Centimeters"}]
]


TransmissionObject[a_Association]["Approximated n"] := With[{
  shift = a["PhaseShift"],
  thickness = QuantityMagnitude[a["Thickness"], "Centimeters"],
  n0 = a["n0"]
},
  QuantityArray[
    Join[{{0, n0}}, {#[[1]], 1.0 + 0.159152 (#[[2]] + shift) / (#[[1]] thickness) } &/@ Drop[Transpose[Normal /@ {a["Frequencies"], a["Phase"]}],1] ]
  , {1/"Centimeters", 1}]
]

(* :: FDCI Functions :: *)

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

(* :: internal :: *)

phaseState[phase_List] := With[{},
  If[Fit[phase, {1, x}, x][[1]] > 10.0,
    Item[Style["Unwrapped", Bold], Background->LightGreen]
  ,
    Item[Style["Possibly Wrapped", Bold], Background->LightYellow]
  ]
]

(* :: Properties list :: *)

TransmissionObject /: Keys[t_TransmissionObject] :=  t["Properties"]

TransmissionObject[a_Association]["Properties"] := Join[Options[TransmissionObject][[All,1]], {"Properties", "n0", "Approximated n", "Approximated k", "Approximated \[Alpha]", "Kramers-Kronig n", "Frequencies", "Transmission", "Phase", "Phase Features", "\[Delta]t", "Gain", "PhaseShift", "Thickness", "Domain", "FrequencyDomainConfidenceInterval", "FDCI", "FDCI2", "FrequencyDomainConfidenceInterval2"}]

Options[TransmissionObject] = {"Thickness"->Null, "Tags"-><||>, "Gain"->1.0, "PhaseShift"->0};


(* :: Options validator :: *)

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

(* :: Automatic phase unwrapping :: *)

TransmissionUnwrap[t: TransmissionObject[a_], "Basic" | Automatic, OptionsPattern[]] := With[{
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


(* :: Semi-automatic phase unwrapping :: *)

TransmissionUnwrap[t: TransmissionObject[a_], "Held" | "Hold" | "Manual", OptionsPattern[]] := With[{
  th = OptionValue["PhaseThreshold"]//N,
  phaseShift = OptionValue["PhaseShift"]
},
  applyDefferedBranching[t, th, phaseShift]
]


splitPhase[tr_, th_:3.16] := With[{
  phase = QuantityMagnitude[tr["Phase Features"], {1/"Centimeters", 1}]
},
  Split[phase, Not[Abs[(#2[[2]]-#1[[2]])] > th] &]
]

autoAdjust[joints_] := Map[Function[item, {item[[1]], item[[1]]}], joints];
recombinePhase[sausages_, joints_] := With[{
  accumulated = Accumulate[joints[[All,2]]] 
}, 
  Flatten[Join[{sausages[[1,All,2]]}, Table[
    Map[Function[p, p[[2]] + 2Pi accumulated[[i-1]]], sausages[[i]]]
  , {i, 2, Length[sausages]}]], 1]
]

applyDefferedBranching[tr_, th_:5.7, shift_:0] := Module[{sausages, joints, offset, freqs},
  offset = 2 Pi (1/33.356) QuantityMagnitude[tr["\[Delta]t"], "Picoseconds"];
  freqs = Normal[tr[[1]]["Frequencies"]];
  sausages = splitPhase[tr, th];
  joints = Table[{-Sign[sausages[[i+1, 1, 2]] - sausages[[i, -1, 2]]], 0}, {i, Length[sausages]-1}] // autoAdjust;


  With[{sausages = sausages, joints = joints, off = offset freqs},
    transmissionPhaseRecombine[tr, shift, off][{sausages, joints}] // Hold 
  ]
]

transmissionPhaseRecombine[tr_, shift_, off_][{sausages_, joints_}] := (Append[tr, {"Phase"->NumericArray[recombinePhase[sausages, joints] + off], "PhaseShift"->shift}]);


End[]
EndPackage[]
