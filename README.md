# CW_ICA
This is a robust algorithm to automatically determine the optimal number of ICs from mixed signals.

This repository contains three code files:

The defined function "CWICA_fastinfo" is used for automatically determining optimal number of ICs.
The input arguments of this function are 
1) Number of repetitions: replicate (numerical value)
2) Mixed signals: X (length of singal by number of mixed signals)
3) Maximum number of ICs: maxic
4) Choice of ICA method: algorithm (fastICA or Infomax),
5) Generate signal-correlation plot: plt=T
The output of the function is the detected optimal number of ICs


Simulated EEG Data are obtained by linear combination of simulated EEG components, which are generated based on characteristics of each component.
The related code can be found in "Simulated Data Generation.R"



All simulated experiment-related code can be found in the file "Plot of simulation result.R", including:
  1. Coef Comparion test
  2. Accuracy test
  3. Robustness test

