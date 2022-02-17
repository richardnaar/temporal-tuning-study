﻿#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2021.2.3),
    on February 17, 2022, at 20:47
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

from __future__ import absolute_import, division

from psychopy import locale_setup
from psychopy import prefs
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

from psychopy.hardware import keyboard



# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '2021.2.3'
expName = 'UR_ssvep_19_stairs_CALIB_bayes'  # from the Builder filename that created this script
expInfo = {'sound': '1', 'participant': '999', 'medium gabor': '0', 'high': '15', 'rotate positions': '1', 'flickering cue': '1', 'low': '4', 'eeg': '0', 'looming fixation': '1', 'intro': '1', 'jumper': '1'}
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='C:\\Users\\Richard Naar\\Documents\\dok\\ssvep\\Visit to York\\PsychoPy\\uncertainty reduction ssvep 13-05-19\\UR_ssvep_2022_Tartu_lastrun.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# Setup the Window
win = visual.Window(
    size=[1920, 1080], fullscr=True, screen=0, 
    winType='pyglet', allowGUI=False, allowStencil=False,
    monitor='kodus', color=[0,0,0], colorSpace='rgb',
    blendMode='avg', useFBO=True, 
    units='pix')
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# Setup eyetracking
ioDevice = ioConfig = ioSession = ioServer = eyetracker = None

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard()

# Initialize components for Routine "main_intro"
main_introClock = core.Clock()
np.random.seed(int(expInfo[ 'participant'])) # this needs to be tested
if int(expInfo[ 'eeg' ]):
    #if the library isn't in the path already
    #import sys
    #sys.path.append('/Users/tyrion/Documents/HaoTing/nbackmindwandering/src/pypixxlib-3.4.4448')
    #sys.path.append('/Users/tyrion/Documents/HaoTing/nbackmindwandering/src/pypixxlib-2.0.3917')

    from pypixxlib.viewpixx import VIEWPixx3D
# enable ViewPixx3D pixel mode (demands firmware revision 31 or higher)
    my_device = VIEWPixx3D()
    my_device.dout.enablePixelMode()
    my_device.updateRegisterCache()

# these values are used to set the RGB triggers
rgbValsRed = [0,2,4,8,16,32,64] # values taken by the green. Each one corresponding primarily to a specific routine in the task (e.g. ITI and 1; CUE and 4 ...)
rgbValsGreen = [0, 1,2,4,8,16,32,128] # values taken by the red. Each one corresponding primarily to a condition in the task (e.g. CUE_HIGH and 4; CUE_LOW and 16 ...)
rgbValsBlue = [0,1] # Values taken by blue. 
# these values correspond to the trigger outputs for previously defined inputs
trigValOutGreen = [57, 25,17,41,89,105,121,81]
trigValOutRed = [57,59,60,62,61,63]
trigValOutBlue = [57, 185]
# to get more information about the triggering consult ViewPixxEEG manual

# from the ITI routine

# set values for the staircase (to measure the point of subjective 
# equality (average of last 10 trials are used to set the individual
# delta values for each participant):

stairesN = 65 # in order to have equal number of repetitions in each condition keep it divisible by 4

# srate == screen refresh
srate = int(win.getActualFrameRate())# get the actual monitor fRate

# save the selected frequencies to a variable
high_freq = float(expInfo[ 'high' ])
low_freq = float(expInfo[ 'low' ])

# SET ROUTINE DURATIONS FOR ITI, CUE, PRED, STIM, FB

cueDur = srate # 1 s

# PERDICTION and ITI duration (based on screen refresh rate and 
# least common multiplier

# ITI duration
if (120-srate) < 2:
    fix_durs = [60, 120, 180] # from 0.5 to 1.5 s
    pred_durs = [240,300,360,420,480, 600,720,840,960] 
    lcm_frex = 60
elif (60-srate) < 2:
    fix_durs = [30, 60, 90]
    pred_durs =  [120,150,180,210,240,300,360,420,480]
    lcm_frex = 30
else:
    print('Check the refresh rate!')


# predictionDur = int(2*srate)

# STIM duration 
stim_duration = srate+lcm_frex

# feedback duration
fbDur=srate # 1 s


# select random orientation for target gratings (different orientation for each)
orients =  randint(1,360)
annulusOri = orients + 90

# set etalon opacity 
#opa1 = 0.4
# equal opacity for the second grating at the beginning of the experiment
#opa2 = 0.4
# default opacities 
#opas_default = [opa1, opa2] # opas_default will remain the same throughout the experiment
#opas = [opa1, opa2] # one of the elements in opas will change each trial (the value of 
# 'level' from the staircase will be added to it)

# set etalon opacity 
cont1 = 0.4
# equal opacity for the second grating at the beginning of the experiment
cont2 = 0.4
# default opacities 
conts_default = [cont1, cont1] # opas_default will remain the same throughout the experiment
conts = [cont1, cont2] # one of the elements in opas will change each trial (the value of 
# 'level' from the staircase will be added to it)

# size of the grating in degrees
gabor_size = 2
# gratings midpoint distance from the centre (radius of the circle in degrees)
gabor_dist = 5
# default sizes
gsizes_default = [(gabor_size, gabor_size), (gabor_size, gabor_size)]

# number of gratings in each arrey
distractorN = 2
contAnnulus = 0.5/2 # 0.7/2 #0.5 # opacity of the annulus is hard coded to 0.4
opaTarget = 0.625/2 # 0.875/2#0.625


# build an arrey of gratings. this information will be sent
# to the graphics card for better timing

# these are the target gratings drawn in the STIM routine
# colors, orientations, positions and opacities will be changing on each iteration.
# everything else remains as it is defined here.
# color variable will be used to give the gratings a desired flickering rate on each trial.
gabors = visual.ElementArrayStim(win, units='deg', fieldPos=(0.0, 0.0), 
fieldSize=(20, 20), fieldShape='circle', nElements=2, 
sizes=gsizes_default, xys=[(-gabor_dist,0),(gabor_dist,0)], colors=([0.5, 0.5, 0.5]) , 
colorSpace='rgb', opacities=opaTarget, oris=orients, sfs=2.0, contrs=conts_default, 
phases=0, elementTex='sin',elementMask='circle', texRes=128, 
interpolate=True, name=None, autoLog=None, maskParams=None)

distractors = visual.ElementArrayStim(win, units='deg', fieldPos=(0.0, 0.0), 
fieldSize=(20, 20), fieldShape='circle', nElements=2, 
sizes=gsizes_default, xys=[(-gabor_dist,0),(gabor_dist,0)], colors=([0.5, 0.5, 0.5]) , 
colorSpace='rgb', opacities=opaTarget, oris=orients, sfs=2.0, contrs=conts_default, 
phases=0, elementTex='sin',elementMask='circle', texRes=128, 
interpolate=True, name=None, autoLog=None, maskParams=None)

# these are the large gratings drawn in all the routines
# colors will be changing on each iteration.
# everything else remains as it is defined here.
# color variable will be used to give the two overlapping gratings different flickering rates.
# note: the opacities of the lower and top gratings are 1 and 0.5 respectively.
# LARGE 
gabors_large = visual.ElementArrayStim(win, units='deg', fieldPos=(0.0, 0.0), 
fieldSize=(30, 30), fieldShape='circle', nElements=distractorN, 
sizes=20, xys=[(0,0), (0,0)], colors=([0.5, 0.5, 0.5]) , 
colorSpace='rgb', opacities=[1,0.5], oris=annulusOri, sfs=2.0, contrs=contAnnulus, 
phases=0, elementTex='sin',elementMask='circle', texRes=128,
interpolate=True, name=None, autoLog=None, maskParams=None)


# SMALL
gabors_small = visual.ElementArrayStim(win, units='deg', fieldPos=(0.0, 0.0), 
fieldSize=(30, 30), fieldShape='circle', nElements=distractorN, 
sizes=7.5, xys=[(0,0), (0,0)], colors=([0.5, 0.5, 0.5]) , 
colorSpace='rgb', opacities=[1, 0.5], oris=annulusOri, sfs=2.0, contrs=contAnnulus, 
phases=0, elementTex='sin',elementMask='circle', texRes=128, 
interpolate=True, name=None, autoLog=None, maskParams=None)

# circle used to cover two larger gabors
circle = visual.Circle(
win=win,
units="deg",
radius=6.5,
fillColor=[0, 0, 0],
lineColor=None,
edges=99
)

# gabor(s) used as fill-in's
# contrast of the first grating will be set to zero if HIGH/LOW CUE
gaborsMedium = visual.ElementArrayStim(win, units='deg', fieldPos=(0.0, 0.0), 
fieldSize=(30, 30), fieldShape='circle', nElements=2, 
sizes=13, xys=[(0,0), (0,0)], colors=([0.5, 0.5, 0.5]) , 
colorSpace='rgb', opacities=[1, 0.5], oris=orients, sfs=2.0, contrs=contAnnulus, 
phases=0, elementTex='sin',elementMask='circle', texRes=128, 
interpolate=True, name=None, autoLog=None, maskParams=None)

# thos is used to present the RGB trigger on screen
# fillColor value will be changed based on the trial information
trigger = visual.Rect(win, units = 'norm', width=0.01, height=0.01, 
autoLog=None, fillColor=[0, 0, 0], fillColorSpace='rgb255', lineColor=None, pos = [-1,1],
opacity = 1)

# a function to draw the large gabors and fill-in on screen
def drawGabors(expInfo,t,frameN, low_freq, high_freq, annulusOri, 
                 cueText, fliCue,t1,t2):
    
    
    if expInfo['flickering cue'] == '1':
        flickering = True
    else:
        flickering = False  # binary
    
    
    if int(expInfo[ 'jumper' ]):
        # change the stimulus colour based on the desired frequency
        # where t is time, trial_freq is the frequency
        sine_high = sin(2*pi*(high_freq)*t)
        stimCol_high = [sine_high,sine_high,sine_high]# [x * sine_high for x in col]
        sine_low = sin(2*pi*(low_freq)*t)
        stimCol_low = [sine_low,sine_low,sine_low]#[x * sine_low for x in col]
        # update gratings
        gabors_large.colors=[stimCol_low, stimCol_high]
        gabors_small.colors=[stimCol_low, stimCol_high]
        # set positions
        gabors_large.oris=annulusOri
        gabors_small.oris=annulusOri
        # set contrast
        gabors_large.contrs=contAnnulus
        gabors_small.contrs=contAnnulus
        gabors_large.draw()
        if flickering and fliCue == 1 and frameN >= t1 and frameN <= t2:
            if cueText == 'fast':
                gaborsMedium.colors = stimCol_high
                gaborsMedium.contrs = [0, contAnnulus]
            elif cueText == 'slow':
                gaborsMedium.colors = stimCol_low
                gaborsMedium.contrs = [0, contAnnulus]
            elif cueText == '?':
                gaborsMedium.colors = [stimCol_low, stimCol_high]
                gaborsMedium.contrs = [contAnnulus, contAnnulus]
#            # TRIGGER: start prediction (show on screen)
#            triggers(expInfo, cueText, 2, 1)
            gaborsMedium.oris = annulusOri+90
            gaborsMedium.draw()
        elif int(expInfo[ 'medium gabor' ]):
            gaborsMedium.oris = annulusOri
            gaborsMedium.contrs = [contAnnulus, 0]
            gaborsMedium.draw()
        else:
            circle.draw()
        gabors_small.draw()

# this is a function used to draw the triggers
def triggers(expInfo, cueText,routine, resp, corr):
    if int(expInfo[ 'eeg' ]):
        my_device.updateRegisterCache()
    if not resp:
        if cueText == 'fast':
            msgTrigger = str(trigValOutRed[routine-1])
            trigRGB = [rgbValsRed[routine-1],0,0]
        elif cueText == 'slow':
            msgTrigger = str(trigValOutRed[routine+2])
            trigRGB = [rgbValsRed[routine+2],0,0]
        elif cueText == '?' and frex == 'low':
            msgTrigger = str(trigValOutGreen[routine])
            trigRGB = [0,rgbValsGreen[routine],0]
        elif cueText == '?' and frex == 'high':
            msgTrigger = str(trigValOutGreen[routine+3])
            trigRGB = [0,rgbValsGreen[routine+3],0]
    else:
        if corr:
            msgTrigger = str(trigValOutGreen[7])
            trigRGB = [0,rgbValsGreen[7],0]
        else:
            msgTrigger = str(trigValOutGreen[7])
            trigRGB = [4,rgbValsGreen[1],0]#28
    trigNum(win, msgTrigger)
    trigger.fillColor = trigRGB

# the body of the trigNum function can be commented in to present
# the trigger numbers in the uppermost left corner. these are the values
# you should see in the EEG recording
def trigNum(win, msgTrigger):
#    trigNum = visual.TextStim(win, msgTrigger, pos=(-0.9, 0.9), units='norm',color=[1,1,1], colorSpace='rgb')
#    trigNum.draw()
    pass
    
key_resp_intro = keyboard.Keyboard()
text_main_intro = visual.TextStim(win=win, name='text_main_intro',
    text='',
    font='Arial',
    units='deg', pos=[0,0], height=1, wrapWidth=25, ori=0, 
    color='black', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-2.0);
pageNr = visual.TextStim(win=win, name='pageNr',
    text='',
    font='Arial',
    units='deg', pos=(0, -10.5), height=1, wrapWidth=None, ori=0, 
    color='black', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-3.0);

# Initialize components for Routine "ITI"
ITIClock = core.Clock()
trialCounter = 0
ellips_fix = visual.ShapeStim(
    win=win, name='ellips_fix', vertices=99,units='deg', 
    size=(3, 1.5),
    ori=0, pos=(0, 0),
    lineWidth=0,     colorSpace='rgb',  lineColor=[0.8,0.8,0.8], fillColor='white',
    opacity=1, depth=-1.0, interpolate=True)

# Initialize components for Routine "cue_s"
cue_sClock = core.Clock()
import math
white_ellips_cfi_2 = visual.ShapeStim(
    win=win, name='white_ellips_cfi_2', vertices=99,units='deg', 
    size=(3,1.5),
    ori=0, pos=(0, 0),
    lineWidth=0,     colorSpace='rgb',  lineColor=[0.8,0.8,0.8], fillColor=[0.1,0.1,0.1],
    opacity=1, depth=-1.0, interpolate=True)
cue_txt_cfi_2 = visual.TextStim(win=win, name='cue_txt_cfi_2',
    text='',
    font='Arial',
    units='deg', pos=(0, 0), height=1, wrapWidth=None, ori=0, 
    color=[0.75,0.75,0.75], colorSpace='rgb', opacity=1.0, 
    languageStyle='LTR',
    depth=-2.0);

# Initialize components for Routine "pred_event_s"
pred_event_sClock = core.Clock()
# just to ensure that the cumulative duration of the prediction time would be 
# equal for all the conditions


preDurN = math.ceil(float(stairesN)/len(pred_durs))
pred_durs_low = np.repeat(pred_durs,preDurN)
pred_durs_high = np.repeat(pred_durs,preDurN)
pred_durs_low_rnd = np.repeat(pred_durs,preDurN)
pred_durs_high_rnd = np.repeat(pred_durs,preDurN) 

pred_durs_low = pred_durs_low[0:int(stairesN)-5]
pred_durs_high = pred_durs_high[0:int(stairesN)-5]
pred_durs_low_rnd = pred_durs_low_rnd[0:int(stairesN)-5]
pred_durs_high_rnd = pred_durs_high_rnd[0:int(stairesN)-5]

pred_durs_low = pred_durs_low.tolist()
pred_durs_high = pred_durs_high.tolist()
pred_durs_low_rnd = pred_durs_low_rnd.tolist()
pred_durs_high_rnd = pred_durs_high_rnd.tolist()

shuffle(pred_durs_low)
shuffle(pred_durs_high)
shuffle(pred_durs_low_rnd)
shuffle(pred_durs_high_rnd)

training_low = pred_durs[0:5]
shuffle(training_low)
training_high = pred_durs[0:5]
shuffle(training_high)
training_low_rnd = pred_durs[0:5]
shuffle(training_low_rnd)
training_high_rnd = pred_durs[0:5]
shuffle(training_high_rnd)

training_low.extend(pred_durs_low)
training_high.extend(pred_durs_high)
training_low_rnd.extend(pred_durs_low_rnd)
training_high_rnd.extend(pred_durs_high_rnd)

pred_durs_low = training_low
pred_durs_high = training_high
pred_durs_low_rnd = training_low_rnd
pred_durs_high_rnd = training_high_rnd

# delta_no = np.matlib.repmat(delta,1,nTrials)

count_fast = 0
count_slow = 0
count_slow_rnd = 0
count_fast_rnd = 0
fix_prediction_2 = visual.ShapeStim(
    win=win, name='fix_prediction_2', vertices=50,units='deg', 
    size=[1.0, 1.0],
    ori=0, pos=(0, 0),
    lineWidth=0,     colorSpace='rgb',  lineColor=[1,1,1], fillColor=[0.25,0.25,0.25],
    opacity=1, depth=-1.0, interpolate=True)

# Initialize components for Routine "stim_s"
stim_sClock = core.Clock()
# all the values taken by 'level' will be appended to this variable
contr = []
correctResps = []
resp_stim_staires = keyboard.Keyboard()
fix_stim_2 = visual.ShapeStim(
    win=win, name='fix_stim_2', vertices=99,units='deg', 
    size=[1.0, 1.0],
    ori=0, pos=(0, 0),
    lineWidth=0,     colorSpace='rgb',  lineColor=[1,1,1], fillColor=[0.25,0.25,0.25],
    opacity=1, depth=-2.0, interpolate=True)

# Initialize components for Routine "pause_s"
pause_sClock = core.Clock()
pauseCounter = 0

pauseMsg = ' '

pauseMsgText = visual.TextStim(win, pauseMsg,
pos=(0, -5.5), units='deg',color=[1,1,1], height = 1, 
bold=False, colorSpace='rgb', italic = True, wrapWidth=35)
text_pause_2 = visual.TextStim(win=win, name='text_pause_2',
    text='',
    font='Arial',
    units='deg', pos=(0, 5), height=1, wrapWidth=35, ori=0, 
    color='black', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);
any_key_end_pause_2 = keyboard.Keyboard()

# Initialize components for Routine "bye"
byeClock = core.Clock()
text = visual.TextStim(win=win, name='text',
    text='',
    font='Arial',
    units='norm', pos=[0,0], height=0.075, wrapWidth=None, ori=0, 
    color='black', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);
key_resp_2 = keyboard.Keyboard()

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# ------Prepare to start Routine "main_intro"-------
continueRoutine = True
# update component parameters for each repeat
from psychopy import sound

clock = core.Clock()
# these are just for presenting instroctions and are not commented at the moment
# same appliest to the whole introduction routine (also "each frame" and "end routine") 

introTextSize = 1
introVol = 0
introCounter = 0
introMSG = []
counter = 0
musicStart = 99
contIntro = conts_default
contIntroBlip = [cont1, 1]
shuffle(contIntroBlip)
stimCol=0
posYtext = 0
introFrex = 4
distrFrex = 15
line = '-'
sine = 0
counterGo = 0
introTxtCol = 'black'
fill_in = 0
keyDown = 0

introTxt = 'Kui Sa pole kindel, millises võres muutus toimus, siis lähtu kõhutundes. Ilmselt märkasid, et võred ei seisa paigal, vaid vilguvad. Võredel on kaks kiirust. Vajuta klaviatuuril \"a\" või \"k\" tähte, et ekraanile esitatud võrede kiirust vahetada.'

introTxt2 = 'Sellel lehel näed nelja väiksemat võret, mis paiknevad suuremate võrede sees olevas hallikas alas. Samamoodi kuvatakse võresid ka katse ajal.'

introTxt2_1 = 'Vastamiseks on 3 sekundit. \n\n\
Vastuse registreerumisest annab märku helitoon. Õige vastuse korral mängib programm kõrget ja vale vastuse korral madalat heli. Kõrge ja madala heli võrdlemiseks vajuta kordamööda tagasi ja edasi nooleklahve klaviatuuril.'

pressLeft = visual.TextStim(win, "Eelmine", 
pos=(-0.7, -0.7), units='norm',color=[-1,-1,-1], height = 0.05, 
bold=True,colorSpace='rgb')

pressRight = visual.TextStim(win, "Järgmine",
pos=(+0.7, -0.7), units='norm',color=[-1,-1,-1], height = 0.05, 
bold=True, colorSpace='rgb')

introText = visual.TextStim(win, introTxt,
pos=(0, -5.5), units='deg',color=[-1,-1,-1], height = 1, 
bold=False, colorSpace='rgb',wrapWidth=25)

introText2 = visual.TextStim(win, introTxt2,
pos=(-17, 0), units='deg',color=[-1,-1,-1], height = 1, 
bold=False, colorSpace='rgb',wrapWidth=10)

introText2_1 = visual.TextStim(win, introTxt2_1,
pos=(17, 0), units='deg',color=[-1,-1,-1], height = 1, 
bold=False, colorSpace='rgb',wrapWidth=10)

go = visual.TextStim(win, '',
pos=(0, 0), units='deg',color=introTxtCol, height = 0.6, 
bold=False,  colorSpace='rgb')

# circle used to cover two larger gabors
introDot = visual.Circle(
win=win,
units="deg",
radius=.15,
fillColor=[0.25, 0.25, 0.25],
lineColor=None,
edges=99
)

introCue = visual.Polygon(win, edges=99, units = 'deg' , pos = (0,0),
size = (3, 1.5), fillColor = (0.1,0.1,0.1), lineColor = None)


# rotate gratings around the center if 'rotate positions' set 1
if int(expInfo[ 'rotate positions' ]) == 1:
    center = [0,0] # relative center
    pos1 = [0,5] # coordinates for the first point
    pos2 = [0,-5] # coordinates for the second point
    pos1m = [center, pos1] # this will be multiplied with the rotation matrix
    pos2m = [center, pos2] # this will be multiplied with the rotation matrix
    angle = randint(45,135) # take random integer from 45 to 135
    theta = np.radians(angle) # convert to radians
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c,-s), (s, c))) # rotation matrix
    pos1 = np.matmul(pos1m, R) # multiply first vector with rotation matrix
    pos2 = np.matmul(pos2m, R) # multiply second vector with rotation matrix
#    gabors.xys = [pos1[1], pos2[1]] # extract x and y and update the gratings
#    if angle < 55 | angle > 125:
#        theta2 = [theta + pi/2.75, theta - pi/2.75] # convert to radians
#        shuffle(theta2)
#        theta2 = theta2[0]
#    else:
#        theta2 = theta + pi/2 # convert to radians
    if angle < 90:
        angle2 = randint(30,60)+angle
    else:
        angle2 = randint(-60,-30)+angle
#    theta2 = theta + pi/3 * multi[0]
    theta2 = np.radians(angle2)
    c2, s2 = np.cos(theta2), np.sin(theta2)
    R2 = np.array(((c2,-s2), (s2, c2))) # rotation matrix
    pos12 = np.matmul(pos1m, R2) # multiply first vector with rotation matrix
    pos22 = np.matmul(pos2m, R2) # multiply second vector with rotation matrix
    distractors.xys = [pos12[1], pos22[1]] # [(0, 5), (0, -5)] # 
key_resp_intro.keys = []
key_resp_intro.rt = []
_key_resp_intro_allKeys = []
# keep track of which components have finished
main_introComponents = [key_resp_intro, text_main_intro, pageNr]
for thisComponent in main_introComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
main_introClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "main_intro"-------
while continueRoutine:
    # get current time
    t = main_introClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=main_introClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    keys_resp_i = event.getKeys()
    
    if int(expInfo[ 'intro' ]):
        posYtext = 0
        if int(introCounter) == 0:
            introMSG = 'Tere tulemast! Täname, et otsustasid enda aega ja energiat meie uurimistöösse investeerida. Selleks hetkeks oled suure tõenäosusega juba läbi teinud üsna elektroodide paigaldamise protseduuri.\n\
            \nAitäh kannatlikkuse eest!'
        elif int(introCounter) == 1:
            gabors.xys=[(-gabor_dist,-1),(gabor_dist,-1)]
            introMSG = 'Katses esitatakse ekraanile neli väikest võret ja üks neist võredest muudab lühikeseks ajahetkeks kontrasti. Sinu ülesandeks on vasaku või parema nooleklahviga märku anda, kas muutus toimus vastavalt vasakul või paremal pool ekraani. Kui vajutad klaviatuuril alla noolt, siis ilmuvad ekraanile katses kasutatavad võred. Üles ja alla noolt kordamööda vajutades näed ka, kuidas kontrasti muutus katses välja nägema hakkab. Pane tähele, et katse jooksul võib muutus olla vaevu märgatav ja seetõttu on hea, kui prooviksid tähelepanu kogu katse vältel katse juures hoida.'
            posYtext = 7
            if 'down' in keys_resp_i: 
                keyDown = 1
                gabors.colors=stimCol
                gabors.draw()
            elif keyDown:
                gabors.colors=stimCol
                gabors.draw()
            introText.draw()
    # target gratings will be presented if 'up' or 'down' pressed on keyboard
        elif int(introCounter) == 2:
            introMSG = ''
            introText2.draw()
            introText2_1.draw()
            drawGabors(expInfo,t, frameN,low_freq, high_freq, annulusOri, 'cueText', fill_in,frameN,frameN+srate)
            gabors.colors=stimCol
            gabors.xys = [pos1[1], pos2[1]]
            distractors.colors=stimCol2
            distractors.draw()
            gabors.draw()
            introDot.draw()
        elif int(introCounter) == 3:
            introMSG  = 'Me ei ütle sulle kummal ekraanipoolel muutus toimub, sest see teeks katse liiga lihtsaks ja läheks ka vastuollu meie uurimiseesmärkidega. Selle asemele vihjame hoopis, kas muutus toimub kiiresti või aeglaselt vilkuvate võrede hulgas. Kui tahaksid meenutada, kuidas kiiresti ja aeglaselt vilkuvad võred välja nägidki, siis liigu tagasi noolega esimesele slaidile, et kiirusi uue pilguga vaadata.\n\
            \nKatses on kolme tüüpi vihjeid. Kui vihjena kuvatakse sõna „kiire“, siis tähendab see, et muutus toimub kiiresti vilkuvate võrede hulgas. Sõna „aeglane“ omakorda aga seda, et muutus toimub aeglaselt vilkuvate võrede hulgas. Kui ekraanile ilmub küsimärk, siis tähendab see, et muutus võib toimuda nii kiiresti, kui aeglaselt vilkuvate võrede hulgas.'
        elif int(introCounter) == 4:
            introMSG = 'Olemegi peaaegu valmis alustama...\n\n\
            \nVeel vaid viimane palve. Palun väldi suuremaid liigutusi. Liigutused jätavad EEG signaali palju müra ja see raskendab saadud andmete analüüsi. Päris soolasambaks muutuda me ei palu ja mõned liigutused on katses ka möödapääsmatud (nt klaviatuuriklahvi vajutamine), kuid palume võimalusel asendi kohendamise või teiste suuremate liigutustega oodata pausini. Katses on kokku 9 pausi. Võta endale pausi ajal täpselt nii palju aega, et jaksaksid igas järgmises katseosas rahulikult ülesandele keskenduda.'
        elif int(introCounter) == 5:
            introMSG  = 'Pöördu palun katse läbiviija poole, et mõõtmine käivitada.\n\
            \nInstruktsioon oli küllalt pikk ja on loomulik, kui end veel mõnes katsega seonduvas aspektis ebakindlalt tunned. Aruta seda küsimust julgelt katse läbiviijaga, sest see tõstab potentsiaalselt nii tulevaste andmete kui katsekogemuse kvaliteeti.'
    else:
        introMSG = 'Siit edasi liikudes algab eksperiment.\n\
        \n Alustamiseks vajuta palun: \"g\"'
    
    if 'g' in keys_resp_i and (introCounter == 0 or introCounter == 5):
        introMSG = ''
        if not counterGo:
            ready = clock.getTime() 
            counterGo = counterGo + 1
            introTextSize = 0.6
            introTxtCol = 'red'
            if int(expInfo[ 'sound' ]):
                efl = sound.Sound(value="Efl", secs=0.1) 
                g = sound.Sound(value="G", secs=0.1) 
                g.play()
                efl.play()
                continueRoutine =  False
    elif 'left' in keys_resp_i and introCounter > 0:
        introCounter -= 1
        introFrex = low_freq
        if int(expInfo[ 'sound' ]):
            c = sound.Sound(value="C", secs=0.1) 
            c.play()
    elif 'right' in keys_resp_i and introCounter < 5:
        introCounter +=  1
        introFrex = low_freq
        core.wait(0.2)
        if int(expInfo[ 'sound' ]):
            g = sound.Sound(value="G", secs=0.1)
            g.play()
    elif 'p' in keys_resp_i:
        introVol=0
    elif key_resp_intro.keys == 'up' and (int(introCounter) == 2 or int(introCounter) == 1):
        gabors.colors=stimCol
        gabors.contrs=contIntroBlip
        gabors.draw()
    elif key_resp_intro.keys == 'down' and (int(introCounter) == 2 or int(introCounter) == 1):
        gabors.colors=stimCol
        gabors.contrs=contIntro
        gabors.draw()
    elif (int(introCounter) == 2 or int(introCounter) == 1) and 'k' in keys_resp_i:
        introFrex = high_freq
        gabors.colors=stimCol
    elif (int(introCounter) == 2 or int(introCounter) == 1) and 'a' in keys_resp_i:
        introFrex = low_freq
        gabors.colors=stimCol
    elif 'c' in keys_resp_i and int(introCounter) == 2:
        fill_in = 1
        drawGabors(expInfo,t, frameN,low_freq, high_freq, annulusOri, 'cueText', fill_in,frameN,frameN+srate)
    
    if clock.getTime()  > 0.5 and counterGo < 2:
        pressLeft.draw()
        pressRight.draw()
    
    # change the stimulus colour based on the desired frequency
    # where t is time, trial_freq is the frequency
    sine = sin(2*pi*(introFrex)*t)
    stimCol = [sine,sine,sine]
    
    sine2 = sin(2*pi*(distrFrex)*t)
    stimCol2 = [sine2,sine2,sine2]
    
    # *key_resp_intro* updates
    waitOnFlip = False
    if key_resp_intro.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        key_resp_intro.frameNStart = frameN  # exact frame index
        key_resp_intro.tStart = t  # local t and not account for scr refresh
        key_resp_intro.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(key_resp_intro, 'tStartRefresh')  # time at next scr refresh
        key_resp_intro.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(key_resp_intro.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(key_resp_intro.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if key_resp_intro.status == STARTED and not waitOnFlip:
        theseKeys = key_resp_intro.getKeys(keyList=['left', 'right', 'p', 'up', 'down', 't', 'k', 'a', 'g', 'c'], waitRelease=False)
        _key_resp_intro_allKeys.extend(theseKeys)
        if len(_key_resp_intro_allKeys):
            key_resp_intro.keys = _key_resp_intro_allKeys[-1].name  # just the last key pressed
            key_resp_intro.rt = _key_resp_intro_allKeys[-1].rt
    
    # *text_main_intro* updates
    if text_main_intro.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text_main_intro.frameNStart = frameN  # exact frame index
        text_main_intro.tStart = t  # local t and not account for scr refresh
        text_main_intro.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text_main_intro, 'tStartRefresh')  # time at next scr refresh
        text_main_intro.setAutoDraw(True)
    if text_main_intro.status == STARTED:  # only update if drawing
        text_main_intro.setPos((0, posYtext), log=False)
        text_main_intro.setText(introMSG, log=False)
    
    # *pageNr* updates
    if pageNr.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        pageNr.frameNStart = frameN  # exact frame index
        pageNr.tStart = t  # local t and not account for scr refresh
        pageNr.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(pageNr, 'tStartRefresh')  # time at next scr refresh
        pageNr.setAutoDraw(True)
    if pageNr.status == STARTED:  # only update if drawing
        pageNr.setText(line + str(introCounter) + line, log=False)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in main_introComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "main_intro"-------
for thisComponent in main_introComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
gabors.contrs=conts_default

expStartTime = core.getTime()
thisExp.addData('Exp time', expStartTime)
thisExp.addData('text_main_intro.started', text_main_intro.tStartRefresh)
thisExp.addData('text_main_intro.stopped', text_main_intro.tStopRefresh)
thisExp.addData('pageNr.started', pageNr.tStartRefresh)
thisExp.addData('pageNr.stopped', pageNr.tStopRefresh)
# the Routine "main_intro" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of trials etc
conditions = data.importConditions('conditions_quest.xlsx')
trials_2 = data.MultiStairHandler(stairType='QUEST', name='trials_2',
    nTrials=stairesN,
    conditions=conditions,
    method='random',
    originPath=-1)
thisExp.addLoop(trials_2)  # add the loop to the experiment
# initialise values for first condition
level = trials_2._nextIntensity  # initialise some vals
condition = trials_2.currentStaircase.condition

for level, condition in trials_2:
    currentLoop = trials_2
    # abbreviate parameter names if possible (e.g. rgb=condition.rgb)
    for paramName in condition:
        exec(paramName + '= condition[paramName]')
    
    # ------Prepare to start Routine "ITI"-------
    continueRoutine = True
    # update component parameters for each repeat
    # TRIGGER: start ITI
    if int(expInfo[ 'eeg' ]):
        my_device.updateRegisterCache()
    
    # RGB values are hard coded, because the value of the ITI
    # trigger do not depend on the trial information (its always 8)
    trigRGB = [0, 0, rgbValsBlue[1]]
    trigger.fillColor = trigRGB
    msgTrigger = str(trigValOutBlue[1])
    
    # set fixation duration randomly 0.5,1 or 1.5 s
    shuffle(fix_durs)
    fixation_duration = int(fix_durs[0])
    
    #give a random orientation to target gratings
    orients = randint(1,360) # randint(1,360)
    annulusOri = orients + 90# 0 # this variable can be used to change the background 
    # grating orientations on each trial
    fliCue = 0 # if this is set zero no fill-in will be presented
    
    trialCounter += 1
    
    # this bit is just for timing test
    lastT = t 
    #print('ITI')
    
    sumOfCount = 0
    frameDurSum = 0
    # keep track of which components have finished
    ITIComponents = [ellips_fix]
    for thisComponent in ITIComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    ITIClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "ITI"-------
    while continueRoutine:
        # get current time
        t = ITIClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=ITIClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        # this bit is just for timing test
        frameDuration = t - lastT
        lastT = t
        sumOfCount += 1
        frameDurSum +=  frameDuration
        
        # draw/send the trigger
        if int(expInfo[ 'eeg' ]):
            trigger.draw()
        
        # draw the background gabors and fill-in (drawGabors() gets defined in ITI routine)
        drawGabors(expInfo,t, frameN,low_freq, high_freq, annulusOri, cueText, fliCue,0,fixation_duration)
        
        # this calls the trigNum function (see ITI code_fix: Begin Experiment)
        trigNum(win, msgTrigger)
        
        # *ellips_fix* updates
        if ellips_fix.status == NOT_STARTED and frameN >= 0.0:
            # keep track of start time/frame for later
            ellips_fix.frameNStart = frameN  # exact frame index
            ellips_fix.tStart = t  # local t and not account for scr refresh
            ellips_fix.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(ellips_fix, 'tStartRefresh')  # time at next scr refresh
            ellips_fix.setAutoDraw(True)
        if ellips_fix.status == STARTED:
            if frameN >= fixation_duration:
                # keep track of stop time/frame for later
                ellips_fix.tStop = t  # not accounting for scr refresh
                ellips_fix.frameNStop = frameN  # exact frame index
                win.timeOnFlip(ellips_fix, 'tStopRefresh')  # time at next scr refresh
                ellips_fix.setAutoDraw(False)
        if ellips_fix.status == STARTED:  # only update if drawing
            ellips_fix.setFillColor([0.1,0.1,0.1], log=False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in ITIComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "ITI"-------
    for thisComponent in ITIComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # this bit is just for timing test
    refR = float(sumOfCount)/float(frameDurSum)
    #print(refR)
    thisExp.addData('cue refR', refR) # write average srate to the file
    trials_2.addOtherData('ellips_fix.started', ellips_fix.tStartRefresh)
    trials_2.addOtherData('ellips_fix.stopped', ellips_fix.tStopRefresh)
    # the Routine "ITI" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # ------Prepare to start Routine "cue_s"-------
    continueRoutine = True
    # update component parameters for each repeat
    
    # shuffle prediction time durations
    # shuffle(pred_durs)
    pred_duration = pred_durs[0]
    # prepare the frequency information for this trial (based on the values defined
    # in the previous loop)
    if frex == 'high':
        trial_freq = high_freq
        dist_freq = low_freq
    elif frex == 'low':
        trial_freq = low_freq
        dist_freq = high_freq
    
    # write into the file
    thisExp.addData('cue', cueText) # cue name
    
    # TRIGGER: start cue
    triggers(expInfo, cueText, 1,0,0)
    
    # this value will be read by the drawGabors() function and enables fill-in cue
    # presentation if set to 1
    fliCue = 1
    
    # this bit is just for timing test
    lastT = t 
    #print('cue')
    
    sumOfCount = 0
    frameDurSum = 0 
    cue_txt_cfi_2.setText(cueText)
    # keep track of which components have finished
    cue_sComponents = [white_ellips_cfi_2, cue_txt_cfi_2]
    for thisComponent in cue_sComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    cue_sClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "cue_s"-------
    while continueRoutine:
        # get current time
        t = cue_sClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=cue_sClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        # this bit is just for timing test
        frameDuration = t - lastT
        lastT = t
        sumOfCount += 1
        frameDurSum +=  frameDuration
        
        # TRIGGER: start cue (show on screen)
        triggers(expInfo, cueText, 1,0,0)
        
        # draw/send the trigger
        if int(expInfo[ 'eeg' ]):
            trigger.draw()
        
        # draw the background gabors and fill-in (drawGabors() gets defined in ITI routine)
        drawGabors(expInfo,t,frameN, low_freq, high_freq, annulusOri, cueText, fliCue,0,cueDur) # t1
        
        #cueText = frameN
        
        # *white_ellips_cfi_2* updates
        if white_ellips_cfi_2.status == NOT_STARTED and frameN >= 0:
            # keep track of start time/frame for later
            white_ellips_cfi_2.frameNStart = frameN  # exact frame index
            white_ellips_cfi_2.tStart = t  # local t and not account for scr refresh
            white_ellips_cfi_2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(white_ellips_cfi_2, 'tStartRefresh')  # time at next scr refresh
            white_ellips_cfi_2.setAutoDraw(True)
        if white_ellips_cfi_2.status == STARTED:
            if frameN >= cueDur:
                # keep track of stop time/frame for later
                white_ellips_cfi_2.tStop = t  # not accounting for scr refresh
                white_ellips_cfi_2.frameNStop = frameN  # exact frame index
                win.timeOnFlip(white_ellips_cfi_2, 'tStopRefresh')  # time at next scr refresh
                white_ellips_cfi_2.setAutoDraw(False)
        
        # *cue_txt_cfi_2* updates
        if cue_txt_cfi_2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
            # keep track of start time/frame for later
            cue_txt_cfi_2.frameNStart = frameN  # exact frame index
            cue_txt_cfi_2.tStart = t  # local t and not account for scr refresh
            cue_txt_cfi_2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(cue_txt_cfi_2, 'tStartRefresh')  # time at next scr refresh
            cue_txt_cfi_2.setAutoDraw(True)
        if cue_txt_cfi_2.status == STARTED:
            if frameN >= cueDur:
                # keep track of stop time/frame for later
                cue_txt_cfi_2.tStop = t  # not accounting for scr refresh
                cue_txt_cfi_2.frameNStop = frameN  # exact frame index
                win.timeOnFlip(cue_txt_cfi_2, 'tStopRefresh')  # time at next scr refresh
                cue_txt_cfi_2.setAutoDraw(False)
        if cue_txt_cfi_2.status == STARTED:  # only update if drawing
            cue_txt_cfi_2.setOpacity(1, log=False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in cue_sComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "cue_s"-------
    for thisComponent in cue_sComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # this bit is just for timing test
    refR = float(sumOfCount)/float(frameDurSum)
    #print(refR)
    thisExp.addData('cue refR', refR) # write average srate to the file
    
    trials_2.addOtherData('white_ellips_cfi_2.started', white_ellips_cfi_2.tStartRefresh)
    trials_2.addOtherData('white_ellips_cfi_2.stopped', white_ellips_cfi_2.tStopRefresh)
    trials_2.addOtherData('cue_txt_cfi_2.started', cue_txt_cfi_2.tStartRefresh)
    trials_2.addOtherData('cue_txt_cfi_2.stopped', cue_txt_cfi_2.tStopRefresh)
    # the Routine "cue_s" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # ------Prepare to start Routine "pred_event_s"-------
    continueRoutine = True
    # update component parameters for each repeat
    # TRIGGER: start prediction
    triggers(expInfo, cueText, 2,0,0)
    
    # this value will be read by the drawGabors() function and enables fill-in cue
    # presentation if set to 1
    fliCue = 0
    
    # take a random duration from the durations in the pred_durs
    # shuffle(pred_durs)
    # pred_duration = int(pred_durs[0])
    
    # this bit is just for timing test
    lastT = t 
    #print('predict')
    
    sumOfCount = 0
    frameDurSum = 0 
    
    
    # this loop will ensure that the cumulative duration of the prediction time would be 
    # equal for all the conditions
    
    if label == 'high':
        pred_duration = pred_durs_high[count_fast]
        count_fast += 1
    elif label == 'low':
        pred_duration = pred_durs_low[count_slow]
        count_slow += 1
    elif label == 'high50':
        pred_duration = pred_durs_high_rnd[count_fast_rnd]
        count_fast_rnd += 1
    elif label == 'low50':
        pred_duration = pred_durs_low_rnd[count_slow_rnd]
        count_slow_rnd += 1
    else:
        print('Something went wrong with prediction durations')
    
    thisExp.addData('prediction duration', pred_duration) 
    # keep track of which components have finished
    pred_event_sComponents = [fix_prediction_2]
    for thisComponent in pred_event_sComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    pred_event_sClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "pred_event_s"-------
    while continueRoutine:
        # get current time
        t = pred_event_sClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=pred_event_sClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        # this bit is just for timing test
        frameDuration = t - lastT
        lastT = t
        sumOfCount += 1
        frameDurSum +=  frameDuration
        
        # TRIGGER: start prediction (show on screen)
        triggers(expInfo, cueText, 2,0,0)
        
        # draw/send the trigger
        if int(expInfo[ 'eeg' ]):
            trigger.draw()
        
        # draw the prediction time fixation
        if int(expInfo[ 'looming fixation' ]):
            fixSize = (t/8, t/8)
        else:
            fixSize = (0.3, 0.3)
        
        # draw the background gabors and fill-in (drawGabors() gets defined in ITI routine)
        drawGabors(expInfo,t, frameN, low_freq, high_freq, annulusOri, cueText, fliCue,0,pred_duration) # predictionDur
        
        
        # *fix_prediction_2* updates
        if fix_prediction_2.status == NOT_STARTED and frameN >= 0:
            # keep track of start time/frame for later
            fix_prediction_2.frameNStart = frameN  # exact frame index
            fix_prediction_2.tStart = t  # local t and not account for scr refresh
            fix_prediction_2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(fix_prediction_2, 'tStartRefresh')  # time at next scr refresh
            fix_prediction_2.setAutoDraw(True)
        if fix_prediction_2.status == STARTED:
            if frameN >= pred_duration:
                # keep track of stop time/frame for later
                fix_prediction_2.tStop = t  # not accounting for scr refresh
                fix_prediction_2.frameNStop = frameN  # exact frame index
                win.timeOnFlip(fix_prediction_2, 'tStopRefresh')  # time at next scr refresh
                fix_prediction_2.setAutoDraw(False)
        if fix_prediction_2.status == STARTED:  # only update if drawing
            fix_prediction_2.setSize(fixSize, log=False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in pred_event_sComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "pred_event_s"-------
    for thisComponent in pred_event_sComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # this bit is just for timing test
    refR = float(sumOfCount)/float(frameDurSum)
    #print(refR)
    thisExp.addData('prediction refR', refR) # write average srate to the file
    trials_2.addOtherData('fix_prediction_2.started', fix_prediction_2.tStartRefresh)
    trials_2.addOtherData('fix_prediction_2.stopped', fix_prediction_2.tStopRefresh)
    # the Routine "pred_event_s" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # ------Prepare to start Routine "stim_s"-------
    continueRoutine = True
    # update component parameters for each repeat
    # TRIGGER: start stim
    triggers(expInfo, cueText, 3,0,0)
    
    # this sets the trial duration. stim_duration time will be extended
    # by lcm_frex and until the maximal duration has reached (~3 s). 
    stim_duration  = srate + lcm_frex 
    stim_duration2 = srate + lcm_frex 
    
    # write orientations to the file
    thisExp.addData('orientation', orients)
    
    multi = [1, -1]
    shuffle(multi)
    
    # rotate gratings around the center if 'rotate positions' set 1
    if int(expInfo[ 'rotate positions' ]) == 1:
        center = [0,0] # relative center
        pos1 = [0,5] # coordinates for the first point
        pos2 = [0,-5] # coordinates for the second point
        pos1m = [center, pos1] # this will be multiplied with the rotation matrix
        pos2m = [center, pos2] # this will be multiplied with the rotation matrix
        angle = randint(45,135) # take random integer from 45 to 135
        theta = np.radians(angle) # convert to radians
        c, s = np.cos(theta), np.sin(theta)
        R = np.array(((c,-s), (s, c))) # rotation matrix
        pos1 = np.matmul(pos1m, R) # multiply first vector with rotation matrix
        pos2 = np.matmul(pos2m, R) # multiply second vector with rotation matrix
        gabors.xys = [pos1[1], pos2[1]] # extract x and y and update the gratings
    #    if angle < 55 | angle > 125:
    #        theta2 = [theta + pi/2.75, theta - pi/2.75] # convert to radians
    #        shuffle(theta2)
    #        theta2 = theta2[0]
    #    else:
    #        theta2 = theta + pi/2 # convert to radians
        if angle < 90:
            angle2 = randint(30,60)+angle
        else:
            angle2 = randint(-60,-30)+angle
    #    theta2 = theta + pi/3 * multi[0]
        theta2 = np.radians(angle2)
        c2, s2 = np.cos(theta2), np.sin(theta2)
        R2 = np.array(((c2,-s2), (s2, c2))) # rotation matrix
        pos12 = np.matmul(pos1m, R2) # multiply first vector with rotation matrix
        pos22 = np.matmul(pos2m, R2) # multiply second vector with rotation matrix
        distractors.xys = [pos12[1], pos22[1]] # [(0, 5), (0, -5)] # 
    
    # opacities
    delta = level
    cont2 = cont1 + delta # opa1 is defined in ITI code component (see Begin Experiment) 
    conts = [cont1, cont2]
    shuffle(conts)
    
    # save the angle of the  target position to the file
    if conts[0] == cont1:
        targetAngle = 180+angle
    else:
        targetAngle = angle
    
    thisExp.addData('target angle', targetAngle)
    thisExp.addData('trial_freq', trial_freq)
    
    # make the target gabors look the same at the beginning
    gabors.contrs= conts_default 
    
    # get information about the correct key
    if conts[0] < conts[1]:
        corrKey = 'left'
        thisExp.addData('corrKey', 'left')
    else:
        corrKey = 'right'
        thisExp.addData('corrKey', 'right')
    
    # changeT defines the blip duration in screen frames from sixt of a
    # screen refresh rate to whole (0.2-1s from the start of the presentation)
    changeT = randint(srate/6,srate)
    
    # updating the orientations and opacities in target gabors
    gabors.oris = orients
    # distractors
    distractors.oris = orients
    
    
    # this value will be read by the drawGabors() function and enables fill-in cue
    # presentation if set to 1
    fliCue = 0
    
    # this bit is just for timing test
    lastT = t
    #print('stim')
    
    sumOfCount = 0
    frameDurSum = 0 
    soundNotPlayed = 1
    resp_stim_staires.keys = []
    resp_stim_staires.rt = []
    _resp_stim_staires_allKeys = []
    # keep track of which components have finished
    stim_sComponents = [resp_stim_staires, fix_stim_2]
    for thisComponent in stim_sComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    stim_sClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "stim_s"-------
    while continueRoutine:
        # get current time
        t = stim_sClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=stim_sClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        # this bit is just for timing test
        frameDuration = t - lastT
        lastT = t
        sumOfCount += 1
        frameDurSum +=  frameDuration
        
        # change the trigger RGB based on if the response is recorded or not
        # and if it was correct or not
        if not resp_stim_staires.keys:
        # TRIGGER: start stim
            triggers(expInfo, cueText, 3,0,0)
        else:
        # TRIGGER: response
            if resp_stim_staires.corr:
                triggers(expInfo, cueText, 3,1,1)
                if soundNotPlayed and int(expInfo[ 'sound' ]):
                    g = sound.Sound(value="G", secs=0.1)
                    g.play()
                    soundNotPlayed = 0
            elif not resp_stim_staires.corr:
                triggers(expInfo, cueText, 3,1,0)
                if soundNotPlayed and int(expInfo[ 'sound' ]):
                    c = sound.Sound(value="C", secs=0.1) 
                    c.play()
                    soundNotPlayed = 0
        
        
        # draw/send the trigger
        if int(expInfo[ 'eeg' ]):
            trigger.draw()
        
        # keep looming for the duration of stim_duration2
        if int(expInfo[ 'looming fixation' ]):
            if t < stim_duration2/srate:
                fixSizeStim = (fixSize[0]+t/16, fixSize[1]+t/16)
        else:
            fixSizeStim = (0.3, 0.3)
        
        
        # change the stimulus colour based on the desired frequency
        # where t is time, trial_freq is the frequency
        sine = sin(2*pi*(trial_freq)*t)
        stimCol = [sine,sine,sine]
        
        # distractors
        sine_d = sin(2*pi*(dist_freq)*t)
        stimCol_d = [sine_d,sine_d,sine_d]
        
        
        # extend the duration by lcm until the response or the maximal duration has reached (~3 s) 
        if not resp_stim_staires.keys and frameN >= stim_duration and stim_duration < (lcm_frex*6):
            stim_duration  = stim_duration+lcm_frex
        
        # draw the background gabors and fill-in (drawGabors() gets defined in ITI routine)
        drawGabors(expInfo,t,frameN, low_freq, high_freq, annulusOri,cueText, fliCue,0,stim_duration)
        
        # update gratings
        if frameN < stim_duration2:
            gabors.colors=stimCol
            distractors.colors=stimCol_d
            distractors.contrs=conts_default
        ##### make the blip happen
            if changeT <= frameN <= changeT+srate*0.2: #
                gabors.contrs=conts
            else:
                gabors.contrs=conts_default
            gabors.draw()
            distractors.draw()
        
        
        # *resp_stim_staires* updates
        waitOnFlip = False
        if resp_stim_staires.status == NOT_STARTED and frameN >= 0:
            # keep track of start time/frame for later
            resp_stim_staires.frameNStart = frameN  # exact frame index
            resp_stim_staires.tStart = t  # local t and not account for scr refresh
            resp_stim_staires.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(resp_stim_staires, 'tStartRefresh')  # time at next scr refresh
            resp_stim_staires.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(resp_stim_staires.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(resp_stim_staires.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if resp_stim_staires.status == STARTED:
            if frameN >= stim_duration:
                # keep track of stop time/frame for later
                resp_stim_staires.tStop = t  # not accounting for scr refresh
                resp_stim_staires.frameNStop = frameN  # exact frame index
                win.timeOnFlip(resp_stim_staires, 'tStopRefresh')  # time at next scr refresh
                resp_stim_staires.status = FINISHED
        if resp_stim_staires.status == STARTED and not waitOnFlip:
            theseKeys = resp_stim_staires.getKeys(keyList=['left', 'right'], waitRelease=False)
            _resp_stim_staires_allKeys.extend(theseKeys)
            if len(_resp_stim_staires_allKeys):
                resp_stim_staires.keys = _resp_stim_staires_allKeys[0].name  # just the first key pressed
                resp_stim_staires.rt = _resp_stim_staires_allKeys[0].rt
                # was this correct?
                if (resp_stim_staires.keys == str(corrKey)) or (resp_stim_staires.keys == corrKey):
                    resp_stim_staires.corr = 1
                else:
                    resp_stim_staires.corr = 0
        
        # *fix_stim_2* updates
        if fix_stim_2.status == NOT_STARTED and frameN >= 0:
            # keep track of start time/frame for later
            fix_stim_2.frameNStart = frameN  # exact frame index
            fix_stim_2.tStart = t  # local t and not account for scr refresh
            fix_stim_2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(fix_stim_2, 'tStartRefresh')  # time at next scr refresh
            fix_stim_2.setAutoDraw(True)
        if fix_stim_2.status == STARTED:
            if frameN >= stim_duration:
                # keep track of stop time/frame for later
                fix_stim_2.tStop = t  # not accounting for scr refresh
                fix_stim_2.frameNStop = frameN  # exact frame index
                win.timeOnFlip(fix_stim_2, 'tStopRefresh')  # time at next scr refresh
                fix_stim_2.setAutoDraw(False)
        if fix_stim_2.status == STARTED:  # only update if drawing
            fix_stim_2.setSize(fixSizeStim, log=False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in stim_sComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "stim_s"-------
    for thisComponent in stim_sComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # this bit is just for timing test
    refR = float(sumOfCount)/float(frameDurSum)
    #print(refR)
    thisExp.addData('stim refR', refR) # write average srate to the file
    
    # convert responses to 0 and 1 values, and write to the file
    
    if resp_stim_staires.keys == 'right':
        whichResponse = 1
    elif resp_stim_staires.keys == 'left':
        whichResponse = 0
    else:
        whichResponse = 'No response'
    
    # write to the file (changeT == time at which the blip happened)
    thisExp.addData('whichResponse', whichResponse)
    thisExp.addData('changeT', changeT)
    
    # save the RT taking blip time into accountt
    if not resp_stim_staires.keys:
        thisExp.addData('trueRT', 'no response')
    else:
        thisExp.addData('trueRT', resp_stim_staires.rt-changeT/srate)
    
    # keep all the levels used in the staircase
    contr.append(level)
    correctResps.append(resp_stim_staires.corr)
    
    # check responses
    if resp_stim_staires.keys in ['', [], None]:  # No response was made
        resp_stim_staires.keys = None
        # was no response the correct answer?!
        if str(corrKey).lower() == 'none':
           resp_stim_staires.corr = 1;  # correct non-response
        else:
           resp_stim_staires.corr = 0;  # failed to respond (incorrectly)
    # store data for trials_2 (MultiStairHandler)
    trials_2.addResponse(resp_stim_staires.corr)
    trials_2.addOtherData('resp_stim_staires.rt', resp_stim_staires.rt)
    trials_2.addOtherData('resp_stim_staires.started', resp_stim_staires.tStartRefresh)
    trials_2.addOtherData('resp_stim_staires.stopped', resp_stim_staires.tStopRefresh)
    trials_2.addOtherData('fix_stim_2.started', fix_stim_2.tStartRefresh)
    trials_2.addOtherData('fix_stim_2.stopped', fix_stim_2.tStopRefresh)
    # the Routine "stim_s" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # ------Prepare to start Routine "pause_s"-------
    continueRoutine = True
    # update component parameters for each repeat
    # trigger 27 is pause
    trigRGB = [rgbValsRed[1], rgbValsGreen[1], rgbValsBlue[0]]
    msgTrigger = '?'
    trigger.fillColor = trigRGB
    
    if trialCounter > 32:
    # only pause for a rest on every 30th trial: 
        if  trialCounter %  30  !=  0: # 
            continueRoutine =  False # so don’t run the pause routine this time.
        else:
            pauseCounter = pauseCounter + 1
    else:
        continueRoutine =  False
    
    
    pauseTxt = 'See on paus. Jätkamiseks vajuta palun nooleklahvi...\n\
    \nÕigete vastuste protsent: '
    
    pauseTxt2 = '\n\nHinnanguliselt on katse lõpuni jäänud veel umbes '
    
    soundNotPlayed = 1
    currentTime = core.getTime()
    timePassed = currentTime - expStartTime
    
    trialRatio = stairesN*4/(trialCounter) # thisN+1
    expDuration = trialRatio * timePassed
    timeLeft = round((expDuration - timePassed)/60)
    
    if pauseCounter == 1:
        expDuration += 3
    text_pause_2.setText(pauseTxt + str(int(np.mean(correctResps)*100)) + str(' %') + pauseTxt2 + str(timeLeft) + str(' minutit'))
    any_key_end_pause_2.keys = []
    any_key_end_pause_2.rt = []
    _any_key_end_pause_2_allKeys = []
    # keep track of which components have finished
    pause_sComponents = [text_pause_2, any_key_end_pause_2]
    for thisComponent in pause_sComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    pause_sClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "pause_s"-------
    while continueRoutine:
        # get current time
        t = pause_sClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=pause_sClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        # draw/send the trigger
        if int(expInfo[ 'eeg' ]):
            trigger.draw()
        
        if int(pauseCounter) == 10:
            pauseMsg = 'Did you know?\n\n\
        \nThe first recording of a human electroencephalography is often attributed to a German psychiatrist Hans Berger, around 1920s. He is also well known for his discovery of the alpha wave or sometimes called \"Berger wave\". What is less known, is how he got involved with EEG in the first place. To be continued...'
            pauseMsgText.setText(pauseMsg)
            pauseMsgText.draw()
        elif int(pauseCounter) == 20:
            pauseMsg = 'Did you know?\n\n\
            \nThe story goes that Hans Berger got interested in EEG because he believed to have experienced spontaneous telepathy with her sister. He reasoned that if it was truly telepathy, recording the electrical activity of the brain might help to bring light to the mechanisms of it. Maybe this belief is to be blamed or something else, but his recordings of the alpha wave were left relatively unnoticed by other investigators for another decade or so. To be continued...'
            pauseMsgText.setText(pauseMsg)
            pauseMsgText.draw()
        elif int(pauseCounter) == 30:
            pauseMsg = 'Did you know?\n\n\
            \nHans Berger\'s work was later picked up by English investigator, Edgar Douglas Adrian, who was awarded for the Nobel prize with another investigator Sir Charles Sherrington, for their discoveries regarding the functions of neurons. One of the things that they investigated was all-or-nothing principle in the nerve cells.'
            pauseMsgText.setText(pauseMsg)
            pauseMsgText.draw()
        elif int(pauseCounter) == 40:
            pauseMsg = 'Did you know?\n\n\
            \nThe same guy, Edgar Douglas Adrian, who replicated Hans Berger\'s findings, was also involved in another matter close to the hearts of the researchers involved in the current experiment. He described something that in today\'s terms would be called \"steady-state visual evoked potentials\", which is a fancy term used to describe brains frequency specific response to a visual flicker.'
            pauseMsgText.setText(pauseMsg)
            pauseMsgText.draw()
        elif int(pauseCounter) == 50:
            pauseMsg = 'Did you know?\n\n\
            \nResearchers are using the Fourier transformation to analyse brain responses to the flickering stimuli, but also in more standard filtering. This method dates back to 19th century French mathematician, Jean-Baptiste Joseph Fourier. He came up with an idea that any signal could be represented as a collection of sinusoidal waves without losing any information. Or another way to think of it is that Fourier transformasion tells you what frequencies are present in your signal and in what proportions.'
            pauseMsgText.setText(pauseMsg)
            pauseMsgText.draw()
        elif int(pauseCounter) == 60:
            pauseMsg = 'Did you know?\n\n\
            \nFourier transformation, that at the time of Jean-Baptiste, attracted attention by few mathematics enthusiasts, is now widely used method for analysing signals across the disciplines, sprinkling from genetics to quantum mechanics, but has also found it\'s way into practical applications and gadgets that we use every day.'
            pauseMsgText.setText(pauseMsg)
            pauseMsgText.draw()
        elif int(pauseCounter) == 70:
            pauseMsg = 'Did you know?\n\n\
            \nThere are a huge number of nerve cells in the brain. According to best of our estimates, this number is around 86 billion (give or take 8 billion). It is often said, to bring that number into perspective, that this is about the same as the number of the stars in the Milky Way galaxy. Latest estimates taught us, that there might be up to 400 billion stars in our galaxy. So lets try to update the popular statement. There are up to as many stars in the Milky Way galaxy as there are nerve cells in the brains of you and your three friends combined. Sounds catchy, doesn\'t it? If you did the math, you found that there is a little miscalculation. I\'d take it as a compliment :)'
            pauseMsgText.setText(pauseMsg)
            pauseMsgText.draw()
        elif int(pauseCounter) == 80:
            pauseMsg = 'Did you know?\n\n\
            \nWelcome to the department of psychology at the University of York! The department holds it\'s roots in the year 1974 when Founding Head of Department Peter Venables appointed the first lecturers. He helped to pioneer the application of physiological measures to psychological problems, paving the way into the top class research in biological basis of psychology and behaviour. '
            pauseMsgText.setText(pauseMsg)
            pauseMsgText.draw()
        
        
        
        # *text_pause_2* updates
        if text_pause_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_pause_2.frameNStart = frameN  # exact frame index
            text_pause_2.tStart = t  # local t and not account for scr refresh
            text_pause_2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_pause_2, 'tStartRefresh')  # time at next scr refresh
            text_pause_2.setAutoDraw(True)
        
        # *any_key_end_pause_2* updates
        waitOnFlip = False
        if any_key_end_pause_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            any_key_end_pause_2.frameNStart = frameN  # exact frame index
            any_key_end_pause_2.tStart = t  # local t and not account for scr refresh
            any_key_end_pause_2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(any_key_end_pause_2, 'tStartRefresh')  # time at next scr refresh
            any_key_end_pause_2.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(any_key_end_pause_2.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(any_key_end_pause_2.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if any_key_end_pause_2.status == STARTED and not waitOnFlip:
            theseKeys = any_key_end_pause_2.getKeys(keyList=None, waitRelease=False)
            _any_key_end_pause_2_allKeys.extend(theseKeys)
            if len(_any_key_end_pause_2_allKeys):
                any_key_end_pause_2.keys = _any_key_end_pause_2_allKeys[-1].name  # just the last key pressed
                any_key_end_pause_2.rt = _any_key_end_pause_2_allKeys[-1].rt
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in pause_sComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "pause_s"-------
    for thisComponent in pause_sComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    trials_2.addOtherData('text_pause_2.started', text_pause_2.tStartRefresh)
    trials_2.addOtherData('text_pause_2.stopped', text_pause_2.tStopRefresh)
    # the Routine "pause_s" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    thisExp.nextEntry()
    
# all staircases completed


# ------Prepare to start Routine "bye"-------
continueRoutine = True
# update component parameters for each repeat
# trigger ??? is bye
trigRGB = [rgbValsRed[3], rgbValsGreen[3], 0]
msgTrigger = '???'
trigger.fillColor = trigRGB

outroVol = 0

outroTxt = 'Thank you for participation! \n\
\nYou got a remarkable '

outroTxt2 = ' % of correct trials.\n\nCongratulations!\n\
\nThis work is a collaboration between the psychology departments of the University of York and the University of Tartu.\n\
\nYork, 2019'


expEndTime = core.getTime()
thisExp.addData('Exp time', expEndTime)
text.setText(outroTxt + str(int(np.mean(correctResps)*100)) + outroTxt2)
key_resp_2.keys = []
key_resp_2.rt = []
_key_resp_2_allKeys = []
# keep track of which components have finished
byeComponents = [text, key_resp_2]
for thisComponent in byeComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
byeClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "bye"-------
while continueRoutine:
    # get current time
    t = byeClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=byeClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    # draw/send the trigger
    if int(expInfo[ 'eeg' ]):
        trigger.draw()
    
    # *text* updates
    if text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text.frameNStart = frameN  # exact frame index
        text.tStart = t  # local t and not account for scr refresh
        text.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text, 'tStartRefresh')  # time at next scr refresh
        text.setAutoDraw(True)
    if text.status == STARTED:  # only update if drawing
        text.setPos((0, cos((t-12)/10)+tan((t-12)/10)), log=False)
    
    # *key_resp_2* updates
    waitOnFlip = False
    if key_resp_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        key_resp_2.frameNStart = frameN  # exact frame index
        key_resp_2.tStart = t  # local t and not account for scr refresh
        key_resp_2.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(key_resp_2, 'tStartRefresh')  # time at next scr refresh
        key_resp_2.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(key_resp_2.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(key_resp_2.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if key_resp_2.status == STARTED and not waitOnFlip:
        theseKeys = key_resp_2.getKeys(keyList=None, waitRelease=False)
        _key_resp_2_allKeys.extend(theseKeys)
        if len(_key_resp_2_allKeys):
            key_resp_2.keys = _key_resp_2_allKeys[-1].name  # just the last key pressed
            key_resp_2.rt = _key_resp_2_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in byeComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "bye"-------
for thisComponent in byeComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('text.started', text.tStartRefresh)
thisExp.addData('text.stopped', text.tStopRefresh)
# the Routine "bye" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv', delim='auto')
thisExp.saveAsPickle(filename)
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
