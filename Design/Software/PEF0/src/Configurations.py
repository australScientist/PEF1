import json
from Measurement import Source
import math


class LEDConfig:

    @staticmethod
    def read_config():
        with open('./Config/config.json', "r") as configFile:
            return(json.load(configFile))
    
    @staticmethod
    def write_config(newConfig): #,path='./Config/config.json'):
        with open('./Config/config.json',"w") as File:
            json.dump(newConfig,File)
            

# finds exposure time that maximizes the signal of the most luminous sample. Least concentrated for transmitance/absorbance measurements or most concentrated for fluorescence measurements.

def resetLEDConfig():
    configFile = {
    "395":{"expTime":10000,"gain":1,"pin":10},
    "440":{"expTime":10000,"gain":1,"pin":9},
    "470":{"expTime":10000,"gain":1,"pin":11},
    "530":{"expTime":10000,"gain":1,"pin":5},
    "600":{"expTime":10000,"gain":1,"pin":6},
    "660":{"expTime":10000,"gain":1,"pin":13},
    "280":{"expTime":10000,"gain":1,"pin":19},
    "WA":{"expTime":10000,"gain":1,"pin":26},
    "365":{"expTime":10000,"gain":1,"pin":21}
    }
    LEDConfig.write_config(newConfig= configFile)

def InitializeLEDsFromFile():
    LEDs={}
    try: 
        config = LEDConfig.read_config()
    except:
        print("Create a config file")
    for key in config:
        LEDs[key] = Source(name=key,pin=config[key]["pin"],exposureTime=config[key]["expTime"],gain=config[key]["gain"])
    print(type(LEDs))
    return(LEDs)

def saveNewConfigsToFile(LEDsHolder):
    configFile = {
    "395":{"expTime":LEDsHolder["395"].exposureTime,"gain":LEDsHolder["395"].gain,"pin":10},
    "440":{"expTime":LEDsHolder["440"].exposureTime,"gain":LEDsHolder["440"].gain,"pin":9},
    "470":{"expTime":LEDsHolder["470"].exposureTime,"gain":LEDsHolder["470"].gain,"pin":11},
    "530":{"expTime":LEDsHolder["530"].exposureTime,"gain":LEDsHolder["530"].gain,"pin":5},
    "600":{"expTime":LEDsHolder["600"].exposureTime,"gain":LEDsHolder["600"].gain,"pin":6},
    "660":{"expTime":LEDsHolder["660"].exposureTime,"gain":LEDsHolder["660"].gain,"pin":13},
    "280":{"expTime":LEDsHolder["280"].exposureTime,"gain":LEDsHolder["280"].gain,"pin":19},
    "WA":{"expTime":LEDsHolder["WA"].exposureTime,"gain":LEDsHolder["WA"].gain,"pin":26},
    "365":{"expTime":LEDsHolder["365"].exposureTime,"gain":LEDsHolder["365"].gain,"pin":21}
    }
    LEDConfig.write_config(newConfig= configFile)

class Calibrate:
    @staticmethod
    def expTime(LightSource = Source,Time = 10000):
        maxValue = LightSource.just_max_signal_please(ExTime=Time)
        print(type(maxValue))
        nIter = 0
        while maxValue <31900 or maxValue > 32100:
            nIter += 1
            if nIter == 100:
                break
            try:
                q = 32000/maxValue
            except:
                q = 32000
            step = 14000*math.log(q)
            Time = Time+step
            maxValue = LightSource.just_max_signal_please(ExTime=Time)
            print(maxValue)
            print(Time)
        LightSource.exposureTime = Time
