BeginPackage["JerryI`TDSTools`Wizard`", {
  "JerryI`TDSTools`Trace`",
  "JerryI`TDSTools`Transmission`",
  "JerryI`TDSTools`Material`",

  "Notebook`Kernel`Inputs`",
  "JerryI`Misc`Events`",
  "JerryI`Misc`Events`Promise`",
  "JerryI`Misc`Language`",
  "JerryI`Misc`Async`",
  "JerryI`Misc`WLJS`Transport`",
  "Notebook`Editor`Boxes`",
  "Notebook`EditorUtils`"
}]

TDSWizard::usage = "TDSWizard"

Begin["`Private`"]


TDSWizard[TransmissionUnwrap][t_TransmissionObject] := TDSWizard[TransmissionUnwrap, Method->Automatic][t]
TDSWizard[TransmissionUnwrap][{t__TransmissionObject}] := TDSWizard[TransmissionUnwrap, Method->Automatic][{t}]

TDSWizard[TransmissionUnwrap][p_Promise] := TDSWizard[TransmissionUnwrap, Method->Automatic][p]

TDSWizard[TransmissionUnwrap, opts__][t_TransmissionObject] := With[{
  parent = EvaluationCell[],
  fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ t["FDCI"],
  aopts = KeyDrop[Association[opts], Method] // Normal,
  cpts  = Association[opts],
  promise = Promise[]
},
  If[!KeyExistsQ[cpts, Method],
    automaicUnwrap[parent, Null, fdci, aopts, t, Function[object,
      CellPrint[ToString[object, StandardForm], "After" -> parent, "Type" -> "Output"];
      EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> object, "Length" -> 1|>];
    ]]
  ,
    If[MatchQ[cpts[Method], "Manual" | "Held" | "Hold"],
      assistedUnwrap[parent, Null, fdci, aopts, t, Function[object,
        CellPrint[ToString[object, StandardForm], "After" -> parent, "Type" -> "Output"];
        EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> object, "Length" -> 1|>];
      ]]
    ,
      automaicUnwrap[parent, Null, fdci, aopts, t, Function[object,
        CellPrint[ToString[object, StandardForm], "After" -> parent, "Type" -> "Output"];
        EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> object, "Length" -> 1|>];
      ]]
    ]
  ];
  
  promise
]

TDSWizard[TransmissionUnwrap, opts__][p_Promise] := With[{
  parent = EvaluationCell[],
  fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ t["FDCI"],
  aopts = KeyDrop[Association[opts], Method] // Normal,
  cpts  = Association[opts],
  promise = Promise[]
},

  Then[p, Function[input,
    If[input["Length"] == 1,
    
      With[{
        t = input["Object"],
        fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ (input["Object"]["FDCI"]),
        parent = input["Cell"]
      },
        If[!KeyExistsQ[cpts, Method],
          automaicUnwrap[parent, Null, fdci, aopts, t, Function[object,
            CellPrint[ToString[object, StandardForm], "After" -> parent, "Type" -> "Output"];
            EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> object, "Length" -> 1|>];
          ]]
        ,
          If[MatchQ[cpts[Method], "Manual" | "Held" | "Hold"],
            assistedUnwrap[parent, Null, fdci, aopts, t, Function[object,
              CellPrint[ToString[object, StandardForm], "After" -> parent, "Type" -> "Output"];
              EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> object, "Length" -> 1|>];
            ]]
          ,
            automaicUnwrap[parent, Null, fdci, aopts, t, Function[object,
              CellPrint[ToString[object, StandardForm], "After" -> parent, "Type" -> "Output"];
              EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> object, "Length" -> 1|>];
            ]]
          ]
        ];
      ];
      
    ,

      LeakyModule[{
        index = 1, 
        checkStack, 
        results
      },
        With[{
          stack = input["Object"],
          resultsName = ToString[results],
          parent = input["Cell"]
        },

          checkStack := If[index > Length[stack],
            CellPrint[resultsName, "After" -> parent, "Type" -> "Output"];
            EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> results, "Length" -> Length[stack]|>];
            ClearAll[index]; ClearAll[checkStack];
          ,
            With[{
              fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ stack[[index]]["FDCI"]
            },


              If[!KeyExistsQ[cpts, Method],
                automaicUnwrap[parent, Null, fdci, aopts, stack[[index]], Function[object,
                    AppendTo[results, object];
                    checkStack;
                ]]
              ,
                If[MatchQ[cpts[Method], "Manual" | "Held" | "Hold"],
                  assistedUnwrap[parent, Null, fdci, aopts, stack[[index]], Function[object,
                    AppendTo[results, object];
                    checkStack;
                  ]]
                ,
                  automaicUnwrap[parent, Null, fdci, aopts, stack[[index]], Function[object,
                    AppendTo[results, object];
                    checkStack;
                  ]]
                ]
              ];

            
            ];
            
            index++;
          ];
        
        
          checkStack;

        ]
      ]
    ]
  ]];  
  
  promise
]

TDSWizard[TransmissionUnwrap, opts__][{t__TransmissionObject}] := LeakyModule[{results, index=1, checkStack}, With[{
  parent = EvaluationCell[],
  fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ t["FDCI"],
  aopts = KeyDrop[Association[opts], Method] // Normal,
  cpts  = Association[opts],
  stack = List[t],
  promise = Promise[],
  resultsName = ToString[results]
},

  results = {};

  checkStack := If[index > Length[stack],
    CellPrint[resultsName, "After" -> parent, "Type" -> "Output"];
    EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> results, "Length" -> Length[stack]|>];
    ClearAll[index]; ClearAll[checkStack];
  ,
    With[{
      fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ stack[[index]]["FDCI"]
    },


      If[!KeyExistsQ[cpts, Method],
        automaicUnwrap[parent, Null, fdci, aopts, stack[[index]], Function[object,
            AppendTo[results, object];
            checkStack;
        ]]
      ,
        If[MatchQ[cpts[Method], "Manual" | "Held" | "Hold"],
          assistedUnwrap[parent, Null, fdci, aopts, stack[[index]], Function[object,
            AppendTo[results, object];
            checkStack;
          ]]
        ,
          automaicUnwrap[parent, Null, fdci, aopts, stack[[index]], Function[object,
            AppendTo[results, object];
            checkStack;
          ]]
        ]
      ];

    
    ];
    
    index++;
  ];


  checkStack;  


  
  promise
] ]

capturePhaseTransform[Hold[callback_[{parts_, joints_}]]] := Module[{
      myJoints = joints
},
      (* modify points *)

      {callback, parts, myJoints}
];

ClearAll[automaicUnwrap];

toBoxes[expr_] := ToString[expr, StandardForm]
automaicUnwrap[parent_, _, fdci_, aopts_, t_, cbk_] := 
 Module[{Global`phase = QuantityMagnitude[
     TransmissionUnwrap[t, Automatic, Sequence @@ aopts]["Phase Features"], 
     {1/"Centimeters", 1}
   ],
   object = t,
   cell
   },
   
   With[{xMinMax = Global`phase[[All, 1]] // MinMax,
         yMinMax = Global`phase[[All, 2]] // MinMax,
         initial = "PhaseThreshold" /. aopts /. Options[TransmissionUnwrap]
       },
      
     cell = CellPrint[
       toBoxes @ Panel[
         Column[{
           Row[{
             EventHandler[
               InputRange[3.14, 6.2, 0.025, initial], 
               Function[value, 
                 Global`phase = QuantityMagnitude[
                   (object = TransmissionUnwrap[t, Automatic, "PhaseThreshold" -> value, Sequence @@ aopts])["Phase Features"], 
                   {1/"Centimeters", 1}
                 ];       
               ]
             ],
             
             EventHandler[
               InputButton["Proceed"], 
               Function[Null,
                 cell // Delete;
                 cbk[object];
                 
               ]
             ]
           }],
           
           Graphics[{
             {Green, Opacity[0.5], 
              Rectangle[{fdci[[1]], yMinMax[[1]]}, {E fdci[[2]], yMinMax[[2]]}]},
             ColorData[97][4], 
             Line[Global`phase // Offload]
           }, PlotRange -> {xMinMax, yMinMax}, Frame -> True, FrameLabel -> {"wavenumber (1/cm)", "Radians"}]
         }], 
         Style["Phase unwrapping", 10]
       ], 
       "After" -> parent, "Type" -> "Output"
     ];
   ]
]

assistedUnwrap[parent_, _, fdci_, aopts_, t_, cbk_] := 
 Module[{Global`phase = QuantityMagnitude[
     TransmissionUnwrap[t, Automatic, Sequence @@ aopts]["Phase Features"], 
     {1/"Centimeters", 1}
   ],
   object = t,
   cell, recombine, parts, joints, Global`jointsTable, makeGrid, 
   lastUpdated = 0, debounce = Null, updateEvent = CreateUUID[]
   },
   
   {recombine, parts, joints} = 
    TransmissionUnwrap[t, "Held", Sequence @@ aopts] // capturePhaseTransform;
    
   makeGrid := With[{grid = MapIndexed[
       Function[{item, index}, 
         {Style[Round[parts[[index[[1]], -1, 1]], 1] // ToString, 10], item[[2]]}
       ], 
       Take[joints, Min[20, Length[joints]]]
     ]},
     Pane[TableForm[grid // Transpose, TableHeadings -> {{"f", "d"}, None}], ImageSize -> 340]
   ];
   
   Global`jointsTable = "";
   Global`jointsTable = toBoxes @ makeGrid;

   With[{xMinMax = Global`phase[[All, 1]] // MinMax,
         yMinMax = Global`phase[[All, 2]] // MinMax,
         initial = "PhaseThreshold" /. aopts /. Options[TransmissionUnwrap]
       },
     
     EventHandler[
       updateEvent, 
       Function[value,
         With[{int = ToExpression[value][[2, 2 ;;]] // Quiet},
           If[ListQ[int],
             With[{length = Length[int]},
               If[AllTrue[int, IntegerQ],
                 joints[[1 ;; length, 2]] = int;
                 object = recombine[{parts, joints}];
                 Global`phase = QuantityMagnitude[object["Phase Features"], {1/"Centimeters", 1}];
               ];
             ]
           ]
         ]
       ]
     ];
     
     cell = CellPrint[
       toBoxes @ Panel[
         Column[{
           Row[{
             EventHandler[
               InputRange[3.14, 6.2, 0.025, initial], 
               Function[value,
                 {recombine, parts, joints} = TransmissionUnwrap[t, "Held", "PhaseThreshold" -> value, Sequence @@ aopts] // capturePhaseTransform;
                 object = recombine[{parts, joints}];
                 Global`phase = QuantityMagnitude[object["Phase Features"], {1/"Centimeters", 1}];
                 
                 lastUpdated = AbsoluteTime[];
                 
                 If[debounce =!= Null, Return[]];
                 
                 debounce = SetInterval[
                   If[AbsoluteTime[] - lastUpdated > 0.8,
                     Global`jointsTable = "";
                     Global`jointsTable = toBoxes @ makeGrid;
                     TaskRemove[debounce];
                     debounce = Null;
                   ];
                   , 300
                 ];
               ]
             ],
             
             EventHandler[
               InputButton["Proceed"], 
               Function[Null,
                 cell // Delete;
                 cbk[object];
               ]
             ]
           }],
           
           Graphics[{
             {Green, Opacity[0.3], 
              Rectangle[{fdci[[1]], yMinMax[[1]]}, {E fdci[[2]], yMinMax[[2]]}]},
             ColorData[97][4], 
             Line[Global`phase // Offload]
           }, PlotRange -> {xMinMax, yMinMax}, Frame -> True, FrameLabel -> {"wavenumber (1/cm)", "Radians"}],
           
           EditorView[Global`jointsTable // Offload, "Event" -> updateEvent]
         }], 
         Style["Assisted phase unwrapping", 10]
       ], 
       "After" -> parent, "Type" -> "Output"
     ];
   ]
]




TDSWizard[MaterialParameters][t_TransmissionObject] := TDSWizard[MaterialParameters, Method->Automatic][t]
TDSWizard[MaterialParameters][p_Promise] := TDSWizard[MaterialParameters, Method->Automatic][p]
TDSWizard[MaterialParameters][{t__TransmissionObject}] := TDSWizard[MaterialParameters, Method->Automatic][{t}]


TDSWizard[MaterialParameters, opts__][p_Promise] := With[{
  aopts = KeyDrop[Association[opts], Method] // Normal,
  cpts  = Association[opts],
  promise = Promise[]
},

  Then[p, Function[input,
    If[input["Length"] == 1,
    
      With[{
        t = input["Object"],
        fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ (input["Object"]["FDCI"]),
        parent = input["Cell"]
      },
        makeWidgetMaterials[t, fdci, opts, parent, Function[result, 
          CellPrint[ToString[result["Object"], StandardForm], "After" -> parent, "Type" -> "Output"];
          EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> result["Object"], "Length" -> 1|>];
        ], <||>];
      ];
      
    ,

      LeakyModule[{
        index = 1, 
        checkStack, 
        results
      },
        With[{
          stack = input["Object"],
          resultsName = ToString[results],
          parent = input["Cell"]
        },

          checkStack := If[index > Length[stack],
            CellPrint[resultsName, "After" -> parent, "Type" -> "Output"];
            EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> results, "Length" -> Length[stack]|>];
            ClearAll[index]; ClearAll[checkStack];
          ,
            With[{
              fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ stack[[index]]["FDCI"]
            },
        
              makeWidgetMaterials[stack[[index]], fdci, opts, parent, Function[result, 
                AppendTo[results, result["Object"]];
                checkStack;
              ], <||>];
            
            ];
            
            index++;
          ];
        
        
          checkStack;

        ]
      ]
    ]
  ]];
  
  
  promise
]


TDSWizard[MaterialParameters, opts__][t_TransmissionObject] := With[{
  parent = EvaluationCell[],
  fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ t["FDCI"],
  aopts = KeyDrop[Association[opts], Method] // Normal,
  cpts  = Association[opts],
  promise = Promise[]
},

  makeWidgetMaterials[t, fdci, opts, parent, Function[result, 
    CellPrint[ToString[result["Object"], StandardForm], "After" -> parent, "Type" -> "Output"];
    EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> result["Object"], "Length" -> 1|>];
  ], <||>];
  
  
  promise
]

TDSWizard[MaterialParameters, opts__][{t__TransmissionObject}] := LeakyModule[{
  index = 1, 
  checkStack, 
  results
}, With[{
  parent = EvaluationCell[],
  aopts = KeyDrop[Association[opts], Method] // Normal,
  cpts  = Association[opts],
  promise = Promise[],
  stack = List[t],
  resultsName = ToString[results]
},

  results = {};

  checkStack := If[index > Length[stack],
    CellPrint[resultsName, "After" -> parent, "Type" -> "Output"];
    EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> results, "Length" -> Length[stack]|>];
    ClearAll[index]; ClearAll[checkStack];
  ,
    With[{
      fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ stack[[index]]["FDCI"]
    },

      makeWidgetMaterials[stack[[index]], fdci, opts, parent, Function[result, 
        AppendTo[results, result["Object"]];
        checkStack;
      ], <||>];
    
    ];
    
    index++;
  ];


  checkStack;
  
  
  promise
]]


makeWidgetMaterials[t_TransmissionObject, fdci_, opts_, parent_, cbk_, meta_] := Module[{cell},
  cell = CellPrint[toBoxes @ With[{
    initialPhase = If[StringQ[#], t["PhaseShift"], #] &@ ("PhaseShift" /. List[opts]),
    initialThickness = QuantityMagnitude[If[QuantityQ[#], #, t["Thickness"]] &@ ("Thickness" /. List[opts]), "Millimeters"],
    initialGain = If[StringQ[#], t["Gain"], #] &@ ("Gain" /. List[opts])
  },
  
  Module[{
    material,
    group,
    Global`alpha,
    Global`rawAlpha,
    Global`n,
    minmaxA,
    minmaxN
  },

    group = EventHandler[InputGroup[<|
        "Done" -> InputButton["Proceed"],
        "Phase" -> InputRange[initialPhase-5, initialPhase+5, 1, initialPhase, "Label"->"Phase shift"],
        "Thickness" -> With[{h = initialThickness}, InputRange[0.7 h, 1.3 h, 0.01 h, h, "Label"->"Thickness"]],
        "Gain" -> InputRange[0.3, 1.5, 0.05, initialGain, "Label"->"Gain"]
    |>], Function[assoc,

      If[assoc["Done"],
        cell // Delete;
        cbk[<|"Cell" -> parent, "Object" -> material, "Length" -> 1|>];
        Return[];
      ];
    
      material = MaterialParameters[
        Append[t, {"PhaseShift"->assoc["Phase"], "Gain"->assoc["Gain"], "Thickness"->Quantity[assoc["Thickness"], "Millimeters"]}]
        , "FabryPerotCancellation"->True, "Target"->"GPU"];
      
      Global`alpha = Select[QuantityMagnitude[material["\[Alpha]"], {1/"Centimeters", 1/"Centimeters"}], Function[item, item[[1]] > (*FB[*)((fdci[[1]])(*,*)/(*,*)(E))(*]FB*) && item[[1]] <  E fdci[[2]]] ];

      Global`rawAlpha = Select[QuantityMagnitude[material["Raw \[Alpha]"], {1/"Centimeters", 1/"Centimeters"}], Function[item, item[[1]] > (*FB[*)((fdci[[1]])(*,*)/(*,*)(E))(*]FB*) && item[[1]] <  E fdci[[2]]] ];

      Global`n = Select[QuantityMagnitude[material["n"], {1/"Centimeters", 1}], Function[item, item[[1]] > (*FB[*)((fdci[[1]])(*,*)/(*,*)(E))(*]FB*) && item[[1]] <  E fdci[[2]]] ];
      
    ]];

    EventFire[group];

    minmaxA = MinMax[Global`alpha[[All,2]]];
    minmaxN = MinMax[Select[Global`n, Function[item, item[[1]] > fdci[[1]] && item[[1]] <  (*SpB[*)Power[E(*|*),(*|*)1/2](*]SpB*) fdci[[2]]] ][[All,2]]];
    minmaxf = MinMax[Global`n[[All,1]]];
  
    Panel[Row[{
      group,

      Graphics[{
        {Opacity[0.2], Green, Rectangle[{fdci[[1]], minmaxA[[1]]}, {fdci[[2]], 2 minmaxA[[2]]}]},
        ColorData[97][1], Line[Global`rawAlpha // Offload], 
        ColorData[97][2], Line[Global`alpha // Offload]
      }, Frame->True, PlotRange->{minmaxf, minmaxA}, FrameLabel->{"wavenumber (1/cm)", "absorption (1/cm)"}],

      Graphics[{
        {Opacity[0.2], Green, Rectangle[{fdci[[1]], minmaxN[[1]]}, {fdci[[2]], 2 minmaxN[[2]]}]},
        ColorData[97][4], Line[Global`n // Offload]
      }, Frame->True, PlotRange->{minmaxf, minmaxN}, FrameLabel->{"wavenumber (1/cm)", "refractive index"}]

    }], Style["Material parameters", 10]]
  ]], "Type"->"Output", "After"->parent];
]




TDSWizard::frontend = "TDSWizard is not supported on this platform"

If[!TrueQ[Internal`Kernel`WLJSQ],
    Message[TDSWizard::frontend];
    ClearAll[TDSWizard];
];


End[]
EndPackage[]