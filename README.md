# CW_ICA
This is a robust algorithm to automatically determine the optimal number of ICs from mixed signals.

Users are required to install RStudio and run the defined CW_ICA function.
Additionally, there are some packages that are also required to install before conducting CW_ICA:
1)fastICA_1.2-3
2)steadyICA_1.0
3)JADE_2.0-4
4)ggpubr_0.6.0
5)ggplot2_3.4.3
6)MASS_7.3-60         


This repository contains four R code files:

The tutorial on the CW_ICA method. The simulation results in the manuscript can be reproduced using this tutorial.


The defined function "CWICA" is used for automatically determining optimal number of ICs.
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

