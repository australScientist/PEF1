import os
from dataclasses import dataclass
import json
import numpy as np
from datetime import datetime

#this class is ment to contain and structure the data resulting from each measurement toghether with the method to save the results to a file
@dataclass
class Spectrum:
    def __init__(self, spectrum,SNR,Noise,BaseLine,ExposureTime,LED,tag = "NA"):
        Spectrum.check_existence()
        self.ID = Spectrum.get_id()
        self.date = datetime.now()
        self.date = self.date.strftime("%d/%m/%Y %H:%M:%S")
        self.LED = LED
        self.spectrum = spectrum
        self.SNR = SNR
        self.noise = Noise
        self.BaseLine = BaseLine
        self.ExposureTime = ExposureTime
        self.tag = tag
    
    @staticmethod
    def check_existence():
        maxSizeFile = 10000000 # maximum size of results file. Otherwise it grows too large and makes saving files slow
        if  os.path.getsize('./Results/results.json') > maxSizeFile:
            date = datetime.now()
            date = date.strftime("%d-%m-%Y_%H:%M:%S")
            newFileName = './Results/results_'+ date +'.json'
            os.rename('./Results/results.json',newFileName)
        if not  os.path.isfile('./Results/results.json'):
            holder = []
            with open('./Results/results.json',"w") as File:
                json.dump(holder,File)
                return(None)
    @staticmethod
    def get_id():
        with open("./Results/results.json", "r") as file:
            results = json.load(file)
        if len(results) == 0:
            return(1)
        lastMeasurement = results[-1]
        identification = lastMeasurement["ID"] + 1
        return(identification)
    
    def formatSpectra(self):
        readyToSave = {
                    "ID" :self.ID,
                    "date" : self.date,
                    "LED" : self.LED,
                    "spectrum" : self.spectrum.tolist(),
                    "SNR" : self.SNR.tolist(),
                    "noise" : self.noise.tolist(),
                    "BaseLine" : self.BaseLine.tolist(),
                    "ExposureTime" : self.ExposureTime,
                    "tag" : self.tag
                }
        return(readyToSave)
    
    def save_results(self):
        with open("./Results/results.json", "r") as file:
            results = json.load(file)
            results.append(self.formatSpectra())            
        
        with open("./Results/results.json", "w") as file:
            json.dump(results,file,indent=4)
