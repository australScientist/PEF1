import numpy as np
from time import *
from gpiozero import Button, LED, PWMLED
from signal import pause
import os
import subprocess
from SpectralData import Spectrum
from ImageProcessor import ProcessThis
import datetime

class Pump:
    pumpPin = 20
    pumpTime = 10
    time = pumpTime
    pmp = LED(pumpPin)
    pmp.off()

class Source:
    
    def __init__ (self, name, pin, exposureTime, gain = 1):
        self.name = name
        self.pin = pin
        self.exposureTime = exposureTime
        self.gain = gain
        self.LightSource = PWMLED(self.pin)
        self.LightSource.value = 0
    def camera(self):
        #takes a raw picture with the parameters passed and those present at the config file
        command = 'libcamera-still -r -o ' + 'measurement.dng' + ' --shutter ' + str(self.exposureTime) + ' --gain '+ str(self.gain) + ' --awbgains 1,1 --immediate --nopreview'
        #some cameras freeze, halting the program. When this happens the process remains open, preventing future usage of the camera. To prevent this a timeout is implemented. In case of error or time limit exceeded the process is killed
        try:
            # Wait 15 seconds, otherwise rise error
            subprocess.run(command, capture_output=True, timeout=15,shell=True)
            return(True,1)
        except subprocess.CalledProcessError as e:
            # code that handles the CalledProcessError exception
            subprocess.run('kill $(fuser /dev/video0)',shell=True)
            return(False,"Error: lib-raw command failed with return code" + e.returncode)
        except subprocess.TimeoutExpired as e:
            # code that handles the TimeoutExpired exception
            subprocess.run('kill $(fuser /dev/video0)',shell=True)
            return (False,"Error: lib-raw command exceeded the specified timeout of" + str(e.timeout) + "seconds")

    def measure(self,save=True,pumpOn=True,Pump = Pump(),ledVal = 1, Tag = "NA"):
        if pumpOn:
            Pump.on()
            sleep(Pump.time)
            Pump.off()
        self.LightSource.value = ledVal
        name = "measurement.dng"
        if os.path.exists(name):
            os.remove(name)
        i = 0
        while i < 10:
            Measurement = self.camera()
            if Measurement[0]:
                break
            print(Measurement[1])
        if i == 9:
            print("Fatal error in camera. Ending the execution")
            exit()
        self.LightSource.value = 0
        Signal = ProcessThis.get_signal("./measurement.dng")
        newSp = Spectrum(spectrum=ProcessThis.get_spectrum(Signal=Signal),SNR=ProcessThis.get_SNR(Signal=Signal),Noise=ProcessThis.get_noise(Signal=Signal),ExposureTime= self.exposureTime ,LED= self.name,BaseLine=ProcessThis.get_baseline(Signal=Signal),tag=Tag)
        if save:
            newSp.save_results()
        return(newSp)

    def just_max_signal_please(self,ExTime,pumpOn=True,Pump = Pump(),ledVal = 1):
        if pumpOn:
            Pump.on()
            sleep(Pump.time)
            Pump.off()
        self.LightSource.value = ledVal
        name = "measurement.dng"
        command = 'libcamera-still -r -o ' + name + ' --shutter ' + str(ExTime) + ' --gain '+ str(self.gain) + ' --awbgains 1,1 --immediate --nopreview'
        os.system(command)
        self.LightSource.value = 0
        Signal = ProcessThis.get_signal("./measurement.dng")
        Spectrum = ProcessThis.get_spectrum(Signal)
        maxValue = max(Spectrum)
        return(maxValue)

