# Old

## Correlation_time.py

Contains a function correlation_time that calculates the correlation time for a given timeseries. It takes:
- timeseries: the timeseries to calculate the correlation time for
- dt: the time step between each data point in the timeseries. Default is 5 GMc^-3
- angle: if the timeseries is of an angular quantity. Should be in radians. Default is False
  Note that the angular correlation time may not be formally correct if the distribution wraps around 
  the 2pi boundary, only in the small angle limit. I hope to update to be formally correct soon, 
  but in the meantime, it may be easier to convert anguaar quantities to the imaginary plane and 
  calculate a correlation time for the x and y components separately.
- plot: if True, will plot the timeseries and the autocorrelation function. Default is False

Credit: Nicholas Conroy.

## Pizza_plotting.py

Pizza plot code. This contains a series of functions to produce pizza plots. 

The first function, constraint_plot, is for Sgr A* Paper 5. The other functions are alternate versions.
It takes pass_table as its primary imput. Pass table is a 4d array. 2 (mad or sane) x 5 (spin) x n (inclination =5 values?) x n (rhigh = 4 values).
0 corresponds to a pass, 1 and -1 correspond to failing (although there's different options... I believe 1 and -1 are used for mutli=False keyword). 
When multi=False, 1 Failure corresponds to too high vaue, and -1 corresponds to too low. Multi=True ignores this and does something else. There are a few different failure states. 

Credit: Michi Baubock and others. I've been told the current pizza plot code is finnicky. Edit with caution.


## load_hdf5.py

Standard hdf5 loading code. Takes hdf5 files of standard ipole format, and produces 
(1) Iall. An array of shape [x, y, t] of all stokes I images for a model. 
(2) Various useful quantities, e.g. the field of view, flux scale, pixel size, time.
(3) a cylinder plot, q, which shows brightness temperature along the ring as a function of PA and time. 

You may need to change a few quantities for your needs, including
-the for loop values, if youre looking for specific models or another library (v3 vs v5)
-the list_directory, where the files are stored
-A few other changes, e.g. uncommenting out the cylinder plot section, or setting a save_directory or plot_directory in your home folder.

Credit: Based on cylinder.py. Nicholas Conroy, Charles Gammie, Michi Baubock. 
