from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as sg
import matplotlib, threading
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from ImageProcessor import ProcessThis as ImP
from itertools import compress
import time
import json
from gpiozero import LED

from Configurations import *
from Measurement import Source, Pump


import RPi.GPIO as GPIO
GPIO.setmode(GPIO.BCM)
GPIO.setwarnings(False)

## FUNCTIONS
def fig_maker(window,data,SampleID):
    for key in data:
        plt.plot(data[key].spectrum, label=key+' nm')
    plt.title(SampleID)
    plt.ylabel = "Light Intensity (AU)"
    plt.xlabel = "Pixel"
    plt.legend()
    window.write_event_value('-THREAD-', 'done.')
    time.sleep(1)
    return plt.gcf()


def draw_figure(canvas, figure, loc=(0, 0)):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg


def delete_fig_agg(fig_agg):
    fig_agg.get_tk_widget().forget()
    plt.close('all')


def preview(lightSource):
    Spectra = {}
    global fig
    global fig_agg
    global window
    global values
    TAG = values['ID'] + '-' + lightSource
    Spectra[lightSource] = lSources[lightSource].measure(pumpOn=False,ledVal = 1, Tag = TAG,save=False)
    #Signal = ImP.get_signal('../measurement.dng')
    #Spectra[lightSource] = ImP.get_spectrum(Signal=Signal)*random.uniform(a=0.5,b=1.5)
    if fig_agg is not None:
        delete_fig_agg(fig_agg)
    fig = fig_maker(window,data=Spectra,SampleID=values['ID'])
    fig_agg = draw_figure(window['canvas'].TKCanvas, fig)

def updateTime(lightSource):
    global window
    global values
    global lSources
    window['time.'+lightSource].update(values['newTime.'+lightSource])
    #change value in lightSource
    lSources[lightSource].exposureTime = values['newTime.'+lightSource]

def calibrateTime(lightSource):
    global window
    global values
    global lSources
    exTime = Calibrate.expTime(lSources[lightSource])
    #exTime = random.uniform(10000,2000000)
    window['time.'+lightSource].update(exTime)
    #change value in lightSource
    lSources[lightSource].exposureTime = exTime


#### INITIALIZE EVERYTHING ######

LEDs = ['WA','280','365','395','440','470','530','600','660']
#PUMP = LED(20)
GPIO.setup(20,GPIO.OUT)
lSources = InitializeLEDsFromFile()
expTimeWA = lSources['WA'].exposureTime
expTime280 = lSources['280'].exposureTime
expTime365 = lSources['365'].exposureTime
expTime395 = lSources['395'].exposureTime
expTime440 = lSources['440'].exposureTime
expTime470 = lSources['470'].exposureTime
expTime530 = lSources['530'].exposureTime
expTime600 = lSources['600'].exposureTime
expTime660 = lSources['660'].exposureTime

#expTimeWA = random.uniform(10000,2000000)
#expTime280 = random.uniform(10000,2000000)
#expTime365 = random.uniform(10000,2000000)
#expTime395 = random.uniform(10000,2000000)
#expTime440 = random.uniform(10000,2000000)
#expTime470 = random.uniform(10000,2000000)
#expTime530 = random.uniform(10000,2000000)
#expTime600 = random.uniform(10000,2000000)
#expTime660 = random.uniform(10000,2000000)


if __name__ == '__main__':

    layout = [
        [sg.Text("Sample ID"), sg.Input(key='ID'),sg.Button("Flush",key='Flush'), sg.Text('Set Time:',key='FlushMessage'),sg.Input(key='FlushTime',s=3)],
        [sg.Button('Preview',key='pview.WA'),sg.Checkbox("WA",key="on.WA"), sg.Text('Time'),sg.Input(key='newTime.WA',s=7),sg.Button('Update',key='update.WA'),sg.Text(expTimeWA,key='time.WA'),sg.Button('Calibrate',key='calib.WA')],
        [sg.Button('Preview',key='pview.280'),sg.Checkbox("280",key="on.280"), sg.Text('Time'),sg.Input(key='newTime.280',s=7),sg.Button('Update',key='update.280'),sg.Text(expTime280,key='time.280'),sg.Button('Calibrate',key="calib.280"),sg.Button('Preview',key='pview.365'),sg.Checkbox("365",key="on.365"), sg.Text('Time'),sg.Input(key='newTime.365',s=7),sg.Button('Update',key='update.365',s=7),sg.Text(expTime365,key='time.365'),sg.Button('Calibrate',key="calib.365")],
        #[],
        [sg.Button('Preview',key='pview.395'),sg.Checkbox("395",key="on.395"), sg.Text('Time'),sg.Input(key='newTime.395',s=7),sg.Button('Update',key='update.395'),sg.Text(expTime395,key='time.395'),sg.Button('Calibrate',key="calib.395"),sg.Button('Preview',key='pview.440'),sg.Checkbox("440",key="on.440"), sg.Text('Time'),sg.Input(key='newTime.440',s=7),sg.Button('Update',key='update.440'),sg.Text(expTime440,key='time.440'),sg.Button('Calibrate',key="calib.440")],
        #[],
        [sg.Button('Preview',key='pview.470'),sg.Checkbox("470",key="on.470"), sg.Text('Time'),sg.Input(key='newTime.470',s=7),sg.Button('Update',key='update.470'),sg.Text(expTime470,key='time.470'),sg.Button('Calibrate',key="calib.470"),sg.Button('Preview',key='pview.530'),sg.Checkbox("530",key="on.530"), sg.Text('Time'),sg.Input(key='newTime.530',s=7),sg.Button('Update',key='update.530'),sg.Text(expTime530,key='time.530'),sg.Button('Calibrate',key="calib.530")],
        #[],
        [sg.Button('Preview',key='pview.600'),sg.Checkbox("600",key="on.600"), sg.Text('Time'),sg.Input(key='newTime.600',s=7),sg.Button('Update',key='update.600'),sg.Text(expTime600,key='time.600'),sg.Button('Calibrate',key="calib.600"),sg.Button('Preview',key='pview.660'),sg.Checkbox("660",key="on.660"), sg.Text('Time'),sg.Input(key='newTime.660',s=7),sg.Button('Update',key='update.660'),sg.Text(expTime660,key='time.660'),sg.Button('Calibrate',key="calib.660")],
        #[],
        [sg.Canvas(size=(240,120), key='canvas'),[sg.Button('Measure'), sg.Button('Stop', key="-STOP-"), sg.Button('Exit', key="-EXIT-")]],
        
    ]
    

    # create the form and show it without the plot
    window = sg.Window('Ecophysiometer 1.0',
                        layout, finalize=True, resizable=True)

    fig_agg = None
    while True:
        event, values = window.read()
        if event is None:  # if user closes window
            break
########################### CONFIGURE EXPOSURE TIME #############################
        if 'update' in event:
            updateTime(lightSource=event.split('.')[1])   

############################## CALIBRATE ###############################################
        if 'calib' in event:
            calibrateTime(lightSource=event.split('.')[1])
            
############################### PREVIEW #############################################
        if 'pview' in event:
            preview(lightSource=event.split('.')[1])
        
############################## Measure ###############################################
        if event == "Measure":
            Active = [
                values["on.WA"],
                values["on.280"],
                values["on.365"],
                values["on.395"],
                values["on.440"],
                values["on.470"],
                values["on.530"],
                values["on.600"],
                values["on.660"]
            ]
            if sum(Active) == 0:
                next
            LEDsActive = list(compress(LEDs,Active))
            Spectra = {}
            for lightSource in LEDsActive:
                TAG = values['ID'] + '-' + lightSource
                Spectra[lightSource] = lSources[lightSource].measure(pumpOn=False,ledVal = 1, Tag = TAG)
                #Signal = ImP.get_signal('../measurement.dng')
                #Spectra[lightSource] = ImP.get_spectrum(Signal=Signal)*random.uniform(a=0.5,b=1.5)

            if fig_agg is not None:
                delete_fig_agg(fig_agg)
            fig = fig_maker(window,data=Spectra,SampleID = values['ID'])
            fig_agg = draw_figure(window['canvas'].TKCanvas, fig)
###################################################################################
        if event == "Flush":
            window['FlushMessage'].update('Flushing')
            #PUMP.on()
            GPIO.output(20,True)
            time.sleep(int(values['FlushTime']))
            GPIO.output(20,False)
            #PUMP.off()
            window['FlushMessage'].update('Set Time:')
        if event == "-EXIT-":
            break
        if event == 'Show':
            # Update the "output" text element to be the value of "input" element
            window['-OUTPUT-'].update(values['-IN-'])
            print(values['-IN-'])
    
    window.close() 
