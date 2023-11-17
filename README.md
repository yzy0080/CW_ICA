# CW_ICA
This is a robust algorithm to automatically determine the optimal number of ICs from mixed signals.

Users are required to install RStudio and run the defined CW_ICA function.
Additionally, there are some packages that are also required to install before conducting CW_ICA:
1)fastICA
2)steadyICA
3)JADE
4)ggpubr
5)ggplot2

This repository contains three R code files:

The defined function "CWICA_fastinfo" is used for automatically determining optimal number of ICs.
The input arguments of this function are 
1) Number of repetitions: replicate (numerical value)
2) Mixed signals: X (length of singal by number of mixed signals)
3) Maximum number of ICs: maxic
4) Choice of ICA method: algorithm (fastICA, Infomax, JADE),
5) Generate signal-correlation plot: plt=T
The output of the function is the detected optimal number of ICs


Simulated EEG Data are obtained by linear combination of simulated EEG components, which are generated based on characteristics of each component.
The related code can be found in "Simulated Data Generation.R"



All simulated experiment-related code can be found in the file "Plot of simulation result.R", including:
  1. Coef Comparion test
  2. Accuracy test
  3. Robustness test

