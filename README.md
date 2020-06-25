# SepsiCool
This project contains software related to the manuscript titled "Robust Adaptive Incorporation of Historical Control Data in a Randomized Trial of External Cooling to Treat Septic Shock" by Thomas A. Murray, Peter F. Thall, Frederique Schortgen, Sarah Zohar and Sandrine Katsahian

This repository contains items (i)-(v) describe herein which facilitate reproducing Table 1 and Figures 2-4 of Murray et al. (2018)

(i) Historical-Control-Data.csv 
	-Contains the 761 historical control observations from the trial reported by Asfar et al. (2014) that will be adaptively incorporated into the new trial
	-y is the time to death/censoring, delta is the indicator for death

(ii) Results folder
	-Contains the simulation results that we obtained and used to generate Table 1 and Figures 2-4

(iii) functions.R 
	-Contains the relevant R functions

(iv) simulation.R
	-Contains R code for carrying out our computer simulation study described in Section 6 and reproducing results already contained in the 'Results' folder
	-Also contains R code for analyzing the contents in the Results folder to generate Table 1 and Figures 2-4 

Work-flow: 
	Keep all these contents in the same working directory 
	To reproduce the contents of the Results folder, execute the "Computer Simulation" section in simulation.R as described within
	To reproduce Tables 1 and Figures 2-4 based on the contents of the Results folder, execute the "Analyze Results" section of simulation.R

R-package dependencies:
dplyr, survival, BayesLogit, jsonlite, gsDesign, mgcv
