BeginPackage["JerryI`TDSTools`Utils`"]

dropHalf;
skip10;
keepMiddle;
fourier2d;

Begin["`Private`"]

dropHalf[list_List] := Drop[list, -Round[Length[list]/2]]

skip10[list_List] := Drop[list, Round[Length[list]/9]]

keepMiddle[list_List] := Drop[Drop[list, -Round[Length[list]/4]], Round[Length[list]/5]]

fourier2d[sample_NumericArray] := fourier2d[sample // Normal]
fourier2d[sample_List] := fourier2d[sample, Differences[sample[[All,1]]] // Min]
fourier2d[sample_List, step_] := With[{int = Interpolation[sample]},
  Transpose[{
    Table[33.356 i, {i, 0, (*FB[*)((1)(*,*)/(*,*)(step))(*]FB*) - (*FB[*)((1)(*,*)/(*,*)(step Length[sample]))(*]FB*), (*FB[*)((1)(*,*)/(*,*)(step Length[sample]))(*]FB*)}],
    Fourier[Table[int[i], {i, sample[[1,1]], sample[[-1,1]], step}]]
  }] // dropHalf
]

fourier2d[sample_List, {min_, max_, step_}] := With[{int = Interpolation[sample], len = (max-min)/step},
  Transpose[{
    Table[33.356 ( (*FB[*)((1)(*,*)/(*,*)(step len step))(*]FB*) (i - min)),  {i, min, max, step}],
    Fourier[Table[int[i], {i, min, max, step}]]
  }] // dropHalf
]
End[]
EndPackage[]