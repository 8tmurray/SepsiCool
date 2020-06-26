This project contains software related to the manuscript titled "Robust Adaptive Incorporation of Historical Control Data in a Randomized Trial of External Cooling to Treat Septic Shock" by Thomas A. Murray, Peter F. Thall, Frederique Schortgen, Sarah Zohar and Sandrine Katsahian

(i) Historical-Control-Data.csv 

	-761 historical control observations from the trial reported by Asfar et al. (2014) that will be adaptively incorporated into the new trial
	-y is the time to death/censoring, delta is the indicator for death

(ii) Results folder

	-Simulation results used to generate Tables and Figures

(iii) functions.R 

	-R functions definitions

(iv) simulation.R

	-R code for carrying out computer simulation study described in Section 6 (Results of this are provided in the Results folder so that these do not need to be re-run)
	-R code for analyzing the contents in the Results folder to generate Tables and Figures

Work-flow: 

	Keep all these contents in the same working directory 
	To reproduce the contents of the Results folder, execute the "Computer Simulation" section in simulation.R as described within
	To reproduce Tables and Figures, execute the "Analyze Results" section of simulation.R

R-package dependencies:
dplyr, survival, BayesLogit, jsonlite, gsDesign, mgcv
