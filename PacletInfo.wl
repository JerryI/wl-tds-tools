(* ::Package:: *)

PacletObject[
  <|
    "Name" -> "JerryI/TDSTools",
    "Description" -> "Tools for working with time-domain THz spectroscopy (transmission)",
    "Creator" -> "Kirill Vasin",
    "License" -> "MIT",
    "PublisherID" -> "JerryI",
    "Version" -> "0.0.5",
    "WolframVersion" -> "13+",
    "PrimaryContext" -> "JerryI`TDSTools`",
    "Extensions" -> {
      {
        "Kernel",
        "Root" -> "Kernel",
        "Context" -> {
          {"JerryI`TDSTools`Material`", "Material.wl"}, 
          {"JerryI`TDSTools`Transmission`", "Transmission.wl"}, 
          {"JerryI`TDSTools`Trace`", "Trace.wl"}, 
          {"JerryI`TDSTools`Wizard`", "Wizard.wl"},  
          {"JerryI`TDSTools`Utils`", "Utils.wl"}                  
        },
        "Symbols" -> {}
      },
 
      {
        "Asset",
        "Assets" -> {}
      }
    }
  |>
]