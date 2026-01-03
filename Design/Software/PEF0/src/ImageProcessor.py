import numpy as np
from scipy.signal import correlate
import rawpy

#This class contains all the methods that are needed to transform the photo into a spectrum

class ProcessThis:
    
    #function that loads the photo, demosaices it and adds the green, red and blue components. Function with side effects
    @staticmethod
    def get_signal(pathToPhoto):
        raw = rawpy.imread(pathToPhoto) # Load the image
        rgb = raw.postprocess(demosaic_algorithm=0,gamma=(1,1), no_auto_bright=True, output_bps=16) # demosaicing with 1,1 gamma compression (exponential 1 and slope one to ensure linearity. Max bit depth (16 bit))
        Signal = rgb[:,:,0]+rgb[:,:,1]+rgb[:,:,2] # compose signal from red, blue and green components
        return(Signal)
    
    #Function used in intermediate steps. Gets the pixel shift between the first and the nth spectra by the cross correlation
    @staticmethod
    def get_pixel_shift(A, B):
        # regularize datasets by subtracting mean and dividing by s.d.
        A = A - A.mean() # center the first spectrum
        A = A / A.std() # scale the first spectrum
        B = B - B.mean() # center the nth spectrum
        B = B / B.std() # Scale the nth spectrum
        xcorr = correlate(A, B) # get cross correlation
        # delta time array to match xcorr
        nsamples = A.size 
        dt = np.arange(1 - nsamples, nsamples) #generate indices
        recovered_time_shift = dt[xcorr.argmax()] 
        return (recovered_time_shift)
    
    # internal function used to align all spectra
    @staticmethod
    def align_spectra(F):
        size = len(F[:, 0])
        Out = F[:, :]  # select first row
        A = F[0, :]  # select first row
        for i in range(1, size):  # and then iterate through the rest, comparing against the first
            B = F[i, :]
            shift = ProcessThis.get_pixel_shift(B, A)
            # calculates the lag between both rows
            Out[i, :] = np.roll(F[i, :],-shift)  # then shifts the row to match both spectra, replacing the spectrum in the original data
        return (Out)
    
    #Wraper function that takes the signal and outputs the spectra ready to use. 
    # The Baseline correction is by default turned off because up to now there is not substantial evidence about the relevance of the effect, but yes about hte effect on the noise
    @staticmethod
    def get_spectrum(Signal,corrected=False):
        SignalSP = Signal[0:450,:]
        Aligned = ProcessThis.align_spectra(SignalSP)
        SpectrumFull = np.mean(Aligned[:,0:2560],axis=0)
        test = np.reshape(SpectrumFull,(512,5))
        Spectrum = np.mean(test,axis=1)
        if corrected:
            Spectrum = Spectrum - ProcessThis.get_baseline(Signal)
        return(Spectrum)
    
    #Get the noise of the uniluminated pixels to calculate the signal to noise ratio

    @staticmethod
    def get_noise(Signal):
        SignalNoise = Signal[1500:1900,0:2560]
        test = np.reshape(SignalNoise,(204800,5)) # Group the spectra into the 5 contiguous pixels cells
        Spectrum = np.mean(test,axis=1) # Average those pixels
        Spectrum = np.reshape(Spectrum,(400,512)) # rearrange the data into the original format (each column has information about waveletgths and rows are replicate measurements)
        Noise = np.std(Spectrum,axis=0) # get the noise as standard deviation of each column
        return(Noise)

    #Get the baseline of signal from the uniluminated pixels
    @staticmethod
    def get_baseline(Signal):
        BaseLine = np.mean(Signal[:,0:2560],axis=0)
        BaseLine = np.reshape(BaseLine,(512,5))
        BaseLine = np.mean(BaseLine,axis=1)
        return(BaseLine)
   # Calculate the signal to noise ratio
    @staticmethod
    def get_SNR(Signal):
        SignalSNR = ProcessThis.get_spectrum(Signal)
        NoiseSNR = ProcessThis.get_noise(Signal)
        SNR = SignalSNR/NoiseSNR
        return(SNR)
