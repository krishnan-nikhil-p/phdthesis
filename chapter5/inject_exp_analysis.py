###prequisites
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageSequence
import glob
import matplotlib.patches as patches
import cv2 as cv

# load the images
dhagIm = Image.open('inject_experiment/dHag_tl_bkg.tif')
beadsIm = Image.open('inject_experiment/beads_tl.tif')
wtIm = Image.open('inject_experiment/wt_tl_bkg.tif')

coords = [[280,80],[220,210],[100,80],[330,170],[140,340],[50, 600]] # coordinates for each box
hT = [8,20,33] 

boxMeansClip = np.zeros((3,6,308)) # establish array to save data
boxDataClip = np.zeros((2,5,308,30,30))
for i, dPage,bPage, wPage in zip(np.arange(308),ImageSequence.Iterator( dhagIm), 
                          ImageSequence.Iterator( beadsIm), ImageSequence.Iterator( wtIm)):
    

    for j, page in enumerate([dPage,bPage, wPage]): #iterate through cannels
        for k, coord in enumerate(coords):
            boxMeansClip[j,k,i] = np.mean((np.clip(np.array(page),0,hT[j]) )[coord[1]:coord[1]+30,coord[0]:coord[0]+30]) #save the mean signal to array