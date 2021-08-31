#!/usr/bin/env python
# coding: utf-8

# In[15]:


#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import scipy as sp
import random
from scipy.ndimage import distance_transform_bf
import matplotlib.pyplot as plt
import multiprocessing as mp
import sys
import os
import glob

# Download dataset from https://github.com/lorinanthony/SECT
dirNameSECT = "./SECT/" # change the SECT directory as necessary
# ID for patients
patName = os.listdir(dirNameSECT+"Data/MITKSegmentations")
# folder for patients
dirName = glob.glob(dirNameSECT+"Data/MITKSegmentations/*/")
# Directory to save results
dirNameSEDT = os.path.join(dirNameSECT+"data/brain/SEDT-2-rearranged/")
os.makedirs(dirNameSEDT,exist_ok=True)

# Define function that can be run in parallel
def computedist(iter):
    # compute distnace
    for i in range(0,len(dirName)):
        subfolder = 'baseline/Segmentations/enh/*.png'
        subdirName = dirName[i] + subfolder
        fileName = glob.glob(subdirName)
        np.random.seed(1000+iter) # seed
    
        for j in range(0,len(fileName)):
            idf = plt.imread(fileName[j])
            rate=np.shape(idf)[0]/256 # rescale image
            idf2 = np.reshape(idf, (-1, 1)) 
            np.random.shuffle(idf2) # rearrange pixels
            idf3 = np.reshape(idf2,idf.shape)
            
            if (np.sum(idf)>=100): # images with more than 100 tumor pixels
                distimgn=distance_transform_bf(idf3,metric='euclidean')/rate
                distimgp=distance_transform_bf(1-idf3,metric='euclidean')/rate
                distimgp = distimgp.astype(np.float64)
                distimgn = distimgn.astype(np.float64)   
                distimg=distimgp-distimgn
                per_disimg=np.ravel(distimg)

                # save as np array
                per_disimg_fin=np.array(per_disimg.flatten())
                info=np.array([2,idf.shape[1],idf.shape[0]])

                base=os.path.basename(fileName[j])
                filename = os.path.splitext(base)[0]

                # write txt file (for computing persistent homology)
                f= open(dirNameSEDT + patName[i] + "_" + filename+ "_" + str(iter) + ".txt","w+")
                for ll in range(0,len(info)):
                    f.write("%d\n" % (info[ll]))
                for mm in range(0,len(per_disimg_fin)):  
                    f.write("%f\n" % (per_disimg_fin[mm]))
                f.close()  

# Only run on main thread
if __name__ == '__main__':
    jobs = []
    
    # Launch processes
    for i in range(50):
        p = mp.Process(target = computedist, args = (i,))
        jobs.append(p)
        p.start()

