from Measurement import Source
from Configurations import *
import json
from Measurement import Source

lSources = InitializeLEDsFromFile()

values = [0,0.01,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1]
board = input("input board number:\n")
carcase = input("input carcase number:\n")
cell = input("input cell number:\n")
camera = input("input camera number:\n")
diff = input("input diffraction grating number:\n")

for X in values:
    TAG = board + "-" + carcase + "-" + cell + "-" + camera+ "-" + diff + "-" + str(X)
    print(TAG)
    lSources["WA"].measure(pumpOn=False,ledVal = X, Tag = TAG)
    
print("DONE")
    
