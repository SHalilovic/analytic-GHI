# analytic-GHI
An analytical approach for estimating the global horizontal from the global tilted irradiance

Implementation of the function calc.GHI(), which analytically computes the global horizontal (Gh) from the global tilted irradiance (Gc) using the new approaches introduced in the paper 'An analytical approach for estimating the global horizontal from the global tilted irradiance'.

Installation:
Download all the files from the repository and place them into one folder.

The files are the following:
- 'Input_example.rds': an exemplary input data, which contains information about timestamps, Gh and Gc for a one day with a 1h resolution.
- 'example.R': R-script with a simple usage of the function calc.GHI(). Here, the input data from 'Input_example.rds' are first read and then prepared for the computation (estimation) of Gh from Gc, which is done using calc.GHI() function. Thereafter, the results are plotted using simple time chart. 
- 'analytical_approach.R': R-script with the implementation of calc.GHI() function. All the information about this function (e.g. inputs, output, etc.) are provided at the beginning of 'analytical_approach.R'.

Usage:
1) Install the R-packages: "maptools", "insol" and "ggplot2".
2) Run the script 'analytical_approach.R', this will load the function calc.GHI().
3) Run the script 'example.R'.
4) In pop-up window select 'Input_example.rds'.
