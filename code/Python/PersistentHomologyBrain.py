#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scipy as sp
import random
import math
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import distance_transform_bf
import sys
import os
import glob
import gudhi as gd # for installation of gudhi, see http://gudhi.gforge.inria.fr/


# In[2]:


# convert list to array
def pdarray(pd):
    pd_array=np.zeros((len(pd),3))
    for i in range(0,len(pd)):
        pd_array[i,0]=np.asarray(pd[i][0])
        pd_array[i,1]=np.asarray(pd[i][1][0])
        pd_array[i,2]=np.asarray(pd[i][1][1])
    return pd_array;


# In[3]:


# Download dataset from https://github.com/lorinanthony/SECT
dirNameSECT = "./SECT/" # change the SECT directory as necessary
# ID for patients
patName = os.listdir(dirNameSECT+"Data/MITKSegmentations")
# folder for patients
dirName = glob.glob(dirNameSECT+"Data/MITKSegmentations/*/")

# Directory to save results
dirNameSEDT = os.path.join(dirNameSECT+"data/brain/SEDT-2/")
os.makedirs(dirNameSEDT,exist_ok=True)
dirNamePD = os.path.join(dirNameSECT+"data/brain/PersistenceDiagram/") 
os.makedirs(dirNamePD,exist_ok=True)


# In[4]:


# compute distnace
for i in range(0,len(dirName)):
    subfolder = 'baseline/Segmentations/enh/*.png'
    subdirName = dirName[i] + subfolder
    fileName = glob.glob(subdirName)
    for j in range(0,len(fileName)):
        idf = plt.imread(fileName[j])
        rate=np.shape(idf)[0]/256
        if (np.sum(idf)>=100):
            distimgn=distance_transform_bf(idf,metric='euclidean')/rate
            distimgp=distance_transform_bf(1-idf,metric='euclidean')/rate
            distimgp = distimgp.astype(np.float64)
            distimgn = distimgn.astype(np.float64)   
            distimg=distimgp-distimgn
            per_disimg=np.ravel(distimg)
        
            # save as np array
            per_disimg_fin=np.array(per_disimg.flatten())
            info=np.array([2,idf.shape[1],idf.shape[0]])
        
            base=os.path.basename(fileName[j])
            filename = os.path.splitext(base)[0]
        
            # write txt file
            f= open(dirNameSEDT + patName[i] + "_" + filename+ ".txt","w+")
            for ll in range(0,len(info)):
                f.write("%d\n" % (info[ll]))
            for mm in range(0,len(per_disimg_fin)):  
                f.write("%f\n" % (per_disimg_fin[mm]))
            f.close()  


# In[5]:


# SEDT-2 file paths
dirPath_sedt = os.path.join(dirNameSEDT,"*.txt")
filePath_sedt = glob.glob(dirPath_sedt)

# compute persistent homology using gudhi
for i in range(0,len(filePath_sedt)):
    # compute PH
    md_cubical_complex = gd.CubicalComplex(perseus_file=filePath_sedt[i])
    # result
    md_cc_diag=md_cubical_complex.persistence()
    
    pd_array=pdarray(md_cc_diag)
    
    # filename
    base=os.path.basename(filePath_sedt[i])
    filename = os.path.splitext(base)[0]
    
    # write txt file
    f= dirNamePD + filename + "_pd.txt"
    np.savetxt(f,pd_array,fmt='%1.6f')    

