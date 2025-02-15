'''  Calculates correlation time for a given timeseries. It takes:
- timeseries: the timeseries to calculate the correlation time for
- dt: the time step between each data point in the timeseries. Default is 5 GMc^-3
- angle: if the timeseries is of an angular quantity. Should be in radians. Default is False
  Note that the angular correlation time may not be formally correct if the distribution wraps around 
  the 2pi boundary, only in the small angle limit. I hope to update to be formally correct soon, 
  but in the meantime, it may be easier to convert anguaar quantities to the imaginary plane and 
  calculate a correlation time for the x and y components separately.
- plot: if True, will plot the timeseries and the autocorrelation function. Default is False

Credit: Nicholas Conroy.

'''

from scipy.stats import circmean
import numpy.fft as fft

def correlation_time(timeseries, dt=5, angle=False, plot=False): 
    ## to calculate the correlation time, the mean must be 0
    if angle: ## if timeseries is of an angular quantity. Should be in radians. 
        timeseries_zeroed = timeseries - circmean(timeseries)
    else: ## if timeseries is of a non-angular quantity
        timeseries_zeroed = timeseries - np.mean(timeseries)

    ## calculate the autocorrelation
    qk = fft.fft(timeseries_zeroed)   ## take fourier transform 
    Pk = np.absolute(qk)**2           ## take square of the magnitude of the FT
    acf = np.real( fft.ifft(Pk) )     ## reverse transform, take real component which drops ~0 imaginary component
    acf = acf/acf[0]                  ## normalize
    index_0 = int(len(acf)/2)
    racf = np.roll(acf, index_0)      ## roll, so that peak is at center

    threshold = 1/np.e                ## find when autocorrelation drops to 1/e
    peak_index = np.argmax(racf)
    correlation_index = np.argmax(racf[peak_index:] <= threshold) 

    ## do linear interpolation to get precise correlation time
    correlation_time_interp = np.interp(threshold, [racf[peak_index + correlation_index], racf[peak_index + correlation_index - 1]], [correlation_index, correlation_index -1])

    if plot:
        t_axis = np.arange(-dt*len(racf)/2, dt*len(racf)/2, dt)
        
        fig, axs = plt.subplots(2, 1)
        axs[0].scatter( t_axis + dt*len(racf)/2, timeseries_zeroed)
        axs[0].set_xlabel("Time [GMc^-3]")
        axs[0].set_ylabel("f")

        axs[1].scatter( t_axis, racf)
        axs[1].axhline(y=threshold, color='r', linestyle='--', label='1/e')
        axs[1].axvline(x=correlation_time_interp*dt, color='g', linestyle='--', label='Correlation Time')
        axs[1].axvline(x=peak_index*dt - dt*len(racf)/2, color='grey', linestyle='--', label='Delta t = 0')
        axs[1].set_xlabel("Delta Time [GMc^-3]")
        axs[1].set_ylabel("RACF")
        plt.show()
    
    ## correlation_time_interp is in units of indices. multiply by dt to get in physical units. standard is dt=5 GMc^-3
    return correlation_time_interp * dt 