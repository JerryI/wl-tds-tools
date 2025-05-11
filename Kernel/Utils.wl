BeginPackage["JerryI`TDSTools`Utils`"]

dropHalf;
skip10;
keepMiddle;
fourier2d;
inverseFourier2d;

Begin["`Private`"]

dropHalf[list_List] := Drop[list, -Round[Length[list]/2]]

skip10[list_List] := Drop[list, Round[Length[list]/9]]

keepMiddle[list_List] := Drop[Drop[list, -Round[Length[list]/4]], Round[Length[list]/5]]

fourier2d[sample_NumericArray] := fourier2d[sample // Normal]
fourier2d[sample_List] := fourier2d[sample, Differences[sample[[All,1]]] // Min]
fourier2d[sample_List, step_] := With[{int = Interpolation[sample], resampledLength =  Length[Table[0, {i, sample[[1,1]], sample[[-1,1]], step}]]},
  Transpose[{
    Table[33.356 i, {i, 0, (*FB[*)((1)(*,*)/(*,*)(step))(*]FB*) - (*FB[*)((1)(*,*)/(*,*)(step resampledLength))(*]FB*), (*FB[*)((1)(*,*)/(*,*)(step resampledLength))(*]FB*)}],
    Fourier[Table[int[i], {i, sample[[1,1]], sample[[-1,1]], step}]]
  }] // dropHalf
]

fourier2d[sample_List, {min_, max_, step_}] := With[{int = Interpolation[sample], len = Length[Table[0,  {i, min, max, step}]]},
  Transpose[{
    Table[33.356 ( (*FB[*)((1)(*,*)/(*,*)(step len step))(*]FB*) (i - min)),  {i, min, max, step}],
    Fourier[Table[int[i], {i, min, max, step}]]
  }] // dropHalf
]


inverseFourier2d[sample_NumericArray] := inverseFourier2d[sample // Normal]
inverseFourier2d[image_List] := With[
  {
    step = Differences[image[[All, 1]]] // Abs // Min,
    min = image[[All,1]] // Min,
    max = image[[All,1]] // Max,
    int = Interpolation[image, InterpolationOrder->1]
  },
    inverseFourier2d[Table[{x,  int[x]}, {x, min, max, step}], {min, max, step}]
]

inverseFourier2d[image_List, {min_, max_, step_}] := Module[
  {
    freqs, fourierData, fullFourier, reconstructed, timeDomain
  },
  
  {freqs, fourierData} = Transpose[image];
  
  fullFourier = Join[fourierData,   {0.}, Drop[Conjugate[Reverse[fourierData]], -1]];
  reconstructed = InverseFourier[fullFourier];

  timeDomain = Transpose@{
    Table[33.356 i, {i, 0, (*FB[*)((1)(*,*)/(*,*)(step))(*]FB*) - (*FB[*)((1)(*,*)/(*,*)(step Length[reconstructed]))(*]FB*), (*FB[*)((1)(*,*)/(*,*)(step Length[reconstructed]))(*]FB*)}],
    Re[reconstructed]
  };
  
  timeDomain
]

End[]
EndPackage[]
