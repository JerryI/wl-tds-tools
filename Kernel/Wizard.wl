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
  aopts = KeyDrop[Association[opts], {Method, "InheritParameters"}] // Normal,
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
  aopts = KeyDrop[Association[opts], {Method, "InheritParameters"}] // Normal,
  cpts  = Association[opts],
  phaseTh = If[NumericQ[#], #, 5.7] &@ ("PhaseThreshold" /. List[opts]),
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

      Module[{
        index = 1, 
        checkStack, 
        results,
        PhaseThreshold = phaseTh
      },

        ClearAttributes[index, Temporal];
ClearAttributes[checkStack, Temporal];
ClearAttributes[results, Temporal];
ClearAttributes[PhaseThreshold, Temporal];
AppendTo[dump, Hold[{index, checkStack, results, PhaseThreshold}]];



        With[{
          stack = input["Object"],
          resultsName = ToString[results],
          parent = input["Cell"]
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
                automaicUnwrap[parent, PhaseThreshold, fdci, aopts, stack[[index]], Function[{object, ph},
                    PhaseThreshold = ph;
                    AppendTo[results, object];
                    checkStack;
                ]]
              ,
                If[MatchQ[cpts[Method], "Manual" | "Held" | "Hold"],
                  assistedUnwrap[parent, PhaseThreshold, fdci, aopts, stack[[index]], Function[{object, ph},
                     PhaseThreshold = ph;
                    AppendTo[results, object];
                    checkStack;
                  ]]
                ,
                  automaicUnwrap[parent, PhaseThreshold, fdci, aopts, stack[[index]], Function[{object, ph},
                   PhaseThreshold = ph;
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

TDSWizard[TransmissionUnwrap, opts__][{t__TransmissionObject}] := Module[{results, index=1, checkStack, PhaseThreshold}, With[{
  parent = EvaluationCell[],
  aopts = KeyDrop[Association[opts], {Method, "InheritParameters"}] // Normal,
  cpts  = Association[opts],
  stack = List[t],
  promise = Promise[],
    phaseTh = If[NumericQ[#], #, 5.7] &@ ("PhaseThreshold" /. List[opts]),
  resultsName = ToString[results]
},



  ClearAttributes[index, Temporal];
ClearAttributes[checkStack, Temporal];
ClearAttributes[results, Temporal];
ClearAttributes[PhaseThreshold, Temporal];
AppendTo[dump, Hold[{index, checkStack, results, PhaseThreshold}]];

PhaseThreshold = phaseTh;

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
        automaicUnwrap[parent, PhaseThreshold, fdci, aopts, stack[[index]], Function[{object, th},
            AppendTo[results, object];
            PhaseThreshold = th;
            checkStack;
        ]]
      ,
        If[MatchQ[cpts[Method], "Manual" | "Held" | "Hold"],
          assistedUnwrap[parent, PhaseThreshold, fdci, aopts, stack[[index]], Function[{object, th},
            AppendTo[results, object];
            PhaseThreshold = th;
            checkStack;
          ]]
        ,
          automaicUnwrap[parent, PhaseThreshold, fdci, aopts, stack[[index]], Function[{object, th},
            AppendTo[results, object];
            PhaseThreshold = th;
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
automaicUnwrap[parent_, initialPhase_, fdci_, aopts_, t_, cbk_] := 
 Module[{Global`phase = QuantityMagnitude[
     TransmissionUnwrap[t, Automatic, Sequence @@ aopts]["Phase Features"], 
     {1/"Centimeters", 1}
   ],
   object = t,
   ph = initialPhase,
   cell
   },
   
   With[{xMinMax = Global`phase[[All, 1]] // MinMax,
         yMinMax = Global`phase[[All, 2]] // MinMax,
         initial = If[initialPhase =!= Null, initialPhase, "PhaseThreshold" /. aopts /. Options[TransmissionUnwrap] ],
         kramer = QuantityMagnitude[t["Kramers-Kronig n"], {1/"Centimeters", 1}] // Re
       },
      
     cell = CellPrint[
       toBoxes @ Panel[
         Column[{
           Row[{
             EventHandler[
               InputRange[3.14, 6.2, 0.025, initial], 
               Function[value, 
                 ph = value;
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
                 cbk[object, ph];
                 
               ]
             ]
           }],
           
           Graphics[{
             {Green, Opacity[0.5], 
              Rectangle[{fdci[[1]], yMinMax[[1]]}, {E fdci[[2]], yMinMax[[2]]}]},
             ColorData[97][4], 
             Line[Global`phase // Offload],
             ColorData[97][7], Opacity[0.3], Line[kramer]
           }, PlotRange -> {xMinMax, yMinMax}, Frame -> True, FrameLabel -> {"wavenumber (1/cm)", "Radians"}]
         }], 
         Style["Phase unwrapping", 10]
       ], 
       "After" -> parent, "Type" -> "Output"
     ];
   ]
]

assistedUnwrap[parent_, initialPhase_, fdci_, aopts_, t_, cbk_] := 
 Module[{Global`phase = QuantityMagnitude[
     TransmissionUnwrap[t, Automatic, Sequence @@ aopts]["Phase Features"], 
     {1/"Centimeters", 1}
   ],
   object = t,
   ph = initialPhase,
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
         initial = If[initialPhase =!= Null, initialPhase, "PhaseThreshold" /. aopts /. Options[TransmissionUnwrap] ],
         kramer = QuantityMagnitude[t["Kramers-Kronig n"], {1/"Centimeters", 1}] // Re
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
                 ph = value;
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
                 cbk[object, ph];
               ]
             ]
           }],
           
           Graphics[{
             {Green, Opacity[0.3], 
              Rectangle[{fdci[[1]], yMinMax[[1]]}, {E fdci[[2]], yMinMax[[2]]}]},
             ColorData[97][4], 
             Line[Global`phase // Offload],
             ColorData[97][7], Opacity[0.3], Line[kramer]
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
  aopts = KeyDrop[Association[opts], {Method, "InheritParameters"}] // Normal,
  cpts  = Association[opts],
  promise = Promise[],
  inherit = TrueQ["InheritParameters" /. List[opts] ]
},

  Then[p, Function[input,
    If[input["Length"] == 1,
    
      With[{
        t = input["Object"],
        fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ (input["Object"]["FDCI"]),
        parent = input["Cell"]
      },
        makeWidgetMaterials[t, Null, fdci, opts, parent, Function[result, 
          CellPrint[ToString[result["Object"], StandardForm], "After" -> parent, "Type" -> "Output"];
          EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> result["Object"], "Length" -> 1|>];
        ], <||>];
      ];
      
    ,

      Module[{
        index = 1, 
        checkStack, 
        results,
        parameters = Null
      },
        
        ClearAttributes[index, Temporal];
ClearAttributes[checkStack, Temporal];
ClearAttributes[results, Temporal];
ClearAttributes[parameters, Temporal];
AppendTo[dump, Hold[{index, checkStack, results, parameters}]];



        With[{
          stack = input["Object"],
          resultsName = ToString[results],
          parent = input["Cell"]
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
        
              makeWidgetMaterials[stack[[index]], parameters, fdci, opts, parent, Function[{result, params}, 
                If[inherit, parameters = params ];
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
  aopts = KeyDrop[Association[opts], {Method, "InheritParameters"}] // Normal,
  cpts  = Association[opts],
  promise = Promise[]
},

  makeWidgetMaterials[t, Null, fdci, opts, parent, Function[{result, parameters}, 
    CellPrint[ToString[result["Object"], StandardForm], "After" -> parent, "Type" -> "Output"];
    EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> result["Object"], "Length" -> 1|>];
  ], <||>];
  
  
  promise
]

TDSWizard[MaterialParameters, opts__][{t__TransmissionObject}] := Module[{
  index = 1, 
  checkStack, 
  results, 
  parameters = Null
}, With[{
  parent = EvaluationCell[],
  aopts = KeyDrop[Association[opts], {Method, "InheritParameters"}] // Normal,
  cpts  = Association[opts],
  inherit = TrueQ["InheritParameters" /. List[opts] ],
  promise = Promise[],
  stack = List[t],
  resultsName = ToString[results]
},

  ClearAttributes[index, Temporal];
ClearAttributes[checkStack, Temporal];
ClearAttributes[results, Temporal];
ClearAttributes[parameters, Temporal];
AppendTo[dump, Hold[{index, checkStack, results, parameters}]];

  results = {};

  checkStack := If[index > Length[stack],
    CellPrint[resultsName, "After" -> parent, "Type" -> "Output"];
    EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> results, "Length" -> Length[stack]|>];
    ClearAll[index]; ClearAll[checkStack];
  ,
    With[{
      fdci = QuantityMagnitude[#, 1/"Centimeters"] &/@ stack[[index]]["FDCI"]
    },

      makeWidgetMaterials[stack[[index]], parameters, fdci, opts, parent, Function[{result, params}, 
        If[inherit, parameters = params];
        AppendTo[results, result["Object"]];
        checkStack;
      ], <||>];
    
    ];
    
    index++;
  ];


  checkStack;
  
  
  promise
]]


makeWidgetMaterials[t_TransmissionObject, params_, fdci_, opts_, parent_, cbk_, meta_] := Module[{cell},
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
    minmaxN,
    parameters
  },

    parameters = If[params === Null, <|"initialPhase" -> initialPhase, "initialThickness" -> initialThickness, "initialGain" -> initialGain|>, params];

    group = EventHandler[InputGroup[<|
        "Done" -> InputButton["Proceed"],
        "Phase" -> InputRange[parameters["initialPhase"]-5, parameters["initialPhase"]+5, 1, parameters["initialPhase"], "Label"->"Phase shift"],
        "Thickness" -> With[{h = parameters["initialThickness"]}, InputRange[0.7 h, 1.3 h, 0.0025 h, h, "Label"->"Thickness"]],
        "Gain" -> InputRange[0.3, 1.5, 0.05, parameters["initialGain"], "Label"->"Gain"]
    |>], Function[assoc,

      If[assoc["Done"],
        cell // Delete;
        cbk[<|"Cell" -> parent, "Object" -> material, "Length" -> 1|>, parameters];
        Return[];
      ];
    
      material = MaterialParameters[
        Append[t, {"PhaseShift"->assoc["Phase"], "Gain"->assoc["Gain"], "Thickness"->Quantity[assoc["Thickness"], "Millimeters"]}]
        , "FabryPerotCancellation"->True, "Target"->"GPU"];

      parameters = <|"initialPhase" -> assoc["Phase"], "initialThickness" -> assoc["Thickness"], "initialGain" -> assoc["Gain"]|>;
      
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


TDSWizard::samref = "Input must be in a form of association with two keys \"Sample\" and \"Reference\"";
TDSWizard::samrefprovided = "Provided sample of reference is not TDTrace"

TDSWizard[TransmissionObject, opts__][_] := (Message[TDSWizard::samref]; $Failed)
TDSWizard[TransmissionObject, opts__][__] := (Message[TDSWizard::samref]; $Failed)

TDSWizard[TransmissionObject][any__] := TDSWizard[TransmissionObject, Method->Automatic, "InheritParameters"->False][any]

TDSWizard[TransmissionObject, opts__][sam_TDTrace, ref_TDTrace] := TDSWizard[TransmissionObject, opts][<|"Sample"->sam, "Reference"->ref|>]

TDSWizard[TransmissionObject, opts__][{sam_TDTrace, ref_TDTrace}] := TDSWizard[TransmissionObject, opts][<|"Sample"->sam, "Reference"->ref|>]

TDSWizard[TransmissionObject, opts__][{list__List}] := TDSWizard[TransmissionObject, opts][Map[(<|"Sample"->#[[1]], "Reference"->#[[2]]|>) &, {list}]]

TDSWizard[TransmissionObject, opts__][a_Association] := Module[{}, With[{
  sample = a["Sample"],
  reference = a["Reference"],
  parent = EvaluationCell[],
  aopts = KeyDrop[Association[opts], {Method, "InheritParameters"}] // Normal,
  promise = Promise[],
  thickness = QuantityMagnitude[If[QuantityQ[#], #, Quantity[1, "Millimeters"]] &@ ("Thickness" /. {opts}), "Millimeters"]
},
  If[!MatchQ[sample, _TDTrace] || !MatchQ[reference, _TDTrace],
    Message[TDSWizard::samrefprovided]; Return[$Failed];
  ];


  makeTDSTransmissionWidget[sample, reference, thickness, aopts, parent, Function[object,
      CellPrint[ToString[object, StandardForm], "After"->parent, "Type"->"Output"];
      EventFire[promise, Resolve, <|"Cell"->parent, "Length"->1, "Object"->object|>];
  ]];

  promise
]]

toBoxes[expr_] := ToString[expr, StandardForm]

makeTDSTransmissionWidget[sample_, reference_, thickness_, aopts_, parent_, cbk_] := Module[{
  object,
  setThickness,
  cell
},

  AppendTo[dump, Hold[{object, setThickness, cell}]];

  setThickness = thickness;

  object = TransmissionObject[sample, reference, "Thickness"->Quantity[setThickness, "Millimeters"], aopts];

  cell = CellPrint[toBoxes @ Panel[Column[{
    EventHandler[InputText[thickness//ToString, "Label"->"Thickness in mm"], Function[value,
      With[{t = ToExpression[value]//Quiet}, 
        If[NumericQ[t],
          setThickness = t
        ]
      ]
    ]],
  
    EventHandler[InputButton["Proceed"], Function[Null,
      Delete[cell];
      object = TransmissionObject[sample, reference, "Thickness"->Quantity[setThickness, "Millimeters"], aopts];
      cbk[object];
    ]],

    object // Snippet
  }], Style["Initial settings", 10]], "After"->parent, "Type"->"Output"];

]

dump = {};


TDSWizard[TransmissionObject, opts__][p_Promise] := With[{
  aopts = KeyDrop[Association[opts], {Method, "InheritParameters"}] // Normal,
  cpts  = Association[opts],
  promise = Promise[],
  thickness = QuantityMagnitude[If[QuantityQ[#], #, Quantity[1, "Millimeters"]] &@ ("Thickness" /. {opts}), "Millimeters"]
},

  Then[p, Function[input,
    If[input["Length"] == 1,
    
      With[{
        t = input["Object"],
        parent = input["Cell"]
      },

        makeTDSTransmissionWidget[t["Sample"], t["Reference"], thickness, aopts, parent, Function[object,
          CellPrint[ToString[object, StandardForm], "After"->parent, "Type"->"Output"];
          EventFire[promise, Resolve, <|"Cell"->parent, "Length"->1, "Object"->object|>];
        ]];
      ];
      
    ,

      Module[{
        index = 1, 
        checkStack, 
        results
      },

        ClearAttributes[index, Temporal];
ClearAttributes[checkStack, Temporal];
ClearAttributes[results, Temporal];
AppendTo[dump, Hold[{index, checkStack, results}]];

        With[{
          stack = input["Object"],
          resultsName = ToString[results],
          parent = input["Cell"]
        },

            results={};

          checkStack := If[index > Length[stack],
            CellPrint[resultsName, "After" -> parent, "Type" -> "Output"];
            EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> results, "Length" -> Length[stack]|>];
            ClearAll[index]; ClearAll[checkStack];
          ,
            With[{

            },
      
              If[index > 1 && TrueQ["InheritParameters" /. List[opts]],
                AppendTo[results, TransmissionObject[stack[[index]]["Sample"], stack[[index]]["Reference"], "Thickness"->(results[[1]]["Thickness"]), aopts]];
                index++;
                checkStack;
                Return[];
              ,
                makeTDSTransmissionWidget[stack[[index]]["Sample"], stack[[index]]["Reference"], thickness, aopts, parent, Function[object,
                  AppendTo[results, object];
                  checkStack;
                ]];
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


TDSWizard[TransmissionObject, opts__][{a__Association}] := Module[{
  index,
  checkStack,
  results
}, With[{
  stack = List[a],
  promise = Promise[],
  parent = EvaluationCell[],
  aopts = KeyDrop[Association[opts], {Method, "InheritParameters"}] // Normal,
  resultsName = ToString[results],
  thickness = QuantityMagnitude[If[QuantityQ[#], #, Quantity[1, "Millimeters"]] &@ ("Thickness" /. {opts}), "Millimeters"]
},

  
  ClearAttributes[index, Temporal];
ClearAttributes[checkStack, Temporal];
ClearAttributes[results, Temporal];
AppendTo[dump, Hold[{index, checkStack, results}]]; (* a problem with Garbage collector *)



  index = 1;
  
  If[!AllTrue[stack, Function[item, MatchQ[item["Sample"], _TDTrace] && MatchQ[item["Reference"], _TDTrace]]],
    Message[TDSWizard::samrefprovided];
    Return[$Failed];
  ];

  results = {};

  

  checkStack := ( If[index > Length[stack],
    CellPrint[resultsName, "After" -> parent, "Type" -> "Output"];
    EventFire[promise, Resolve, <|"Cell" -> parent, "Object" -> results, "Length" -> Length[stack]|>];
    ClearAll[index]; ClearAll[checkStack];
  ,
    With[{
      sample = stack[[index]]["Sample"],
      reference = stack[[index]]["Reference"]
    },


      If[index > 1 && TrueQ["InheritParameters" /. List[opts]],
        AppendTo[results, TransmissionObject[sample, reference, "Thickness"->(results[[1]]["Thickness"]), aopts]];
        index++;
        checkStack;
        Return[];
      ,
        makeTDSTransmissionWidget[sample, reference, thickness, aopts, parent, Function[object,
              AppendTo[results, object];
              checkStack;
        ]];
      ]

    
    ];
    
    index++;
  ]);


  checkStack;  

  promise

]]


TDSWizard::frontend = "TDSWizard is not supported on this platform"

If[!TrueQ[Internal`Kernel`WLJSQ],
    Message[TDSWizard::frontend];
    ClearAll[TDSWizard];
];


End[]
EndPackage[]