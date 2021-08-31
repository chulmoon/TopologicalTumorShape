#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import scipy as sp
import math
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import distance_transform_bf
import sys
import os
import pyreadr
import glob
import gudhi as gd # for installation of gudhi, see http://gudhi.gforge.inria.fr/


# In[3]:


# convert list to array
def pdarray(pd):
    pd_array=np.zeros((len(pd),3))
    for i in range(0,len(pd)):
        pd_array[i,0]=np.asarray(pd[i][0])
        pd_array[i,1]=np.asarray(pd[i][1][0])
        pd_array[i,2]=np.asarray(pd[i][1][1])
    return pd_array;


# In[4]:


# base file directory
dirName = "./TopologicalTumorShape/" # change as necessary
# read lung data image files
dirPath = os.path.join(dirName, "data/lung/image/*.Rdata")
filePath = glob.glob(dirPath)
# Directory to save results
dirNameSEDT = os.path.join(dirName+"data/lung/SEDT-3/")
os.makedirs(dirNameSEDT,exist_ok=True)
dirNamePD = os.path.join(dirName+"data/lung/PersistenceDiagram/") 
os.makedirs(dirNamePD,exist_ok=True)


# In[ ]:


for i in range(0,len(filePath)):   
    imgdat = pyreadr.read_r(filePath[i])
    idf = imgdat["dataset"]
    
    # how images are coded
    ## normal regions: 0
    ## tumor regions: 1 
    ## empty: 2    
    dfnormal=idf.copy()
    dftumor=idf.copy()
    
    matdf = idf.values
    
    # for normal, normal -1, otherwise: 0
    matnormal = dfnormal.values
    # first make tumor or empty regions 1 and plug in 1-matnormal
    ## tumor or empty regions
    matnormal[matnormal>1]=1

    # for tumor, tumor: 1, otherwise:0
    mattumor = dftumor.values
    ## make empty regions 0 (only tumors are 1)
    mattumor[mattumor>1]=0
    
    # compute distance
    ## negative: tumor
    distimgn=distance_transform_bf(mattumor,metric='euclidean')
    ## positive: normal
    distimgp=distance_transform_bf(1-matnormal,metric='euclidean')
    
    distimgp = distimgp.astype(np.float64)
    distimgn = distimgn.astype(np.float64)
    
    ## aggregated distance
    distimg=distimgp-distimgn
    per_disimg=np.ravel(distimg)
    
    # replace empty cells with inf
    per_matdf=np.ravel(matdf)
    per_disimg[per_matdf>1]=np.inf
    
    # filename
    base=os.path.basename(filePath[i])
    filename = os.path.splitext(base)[0]
    
    # save as np array
    per_disimg_fin=np.array(per_disimg.flatten())
    info=np.array([2,idf.shape[1],idf.shape[0]])

    # write txt file (for persistent homology computation)
    f= open(dirNameSEDT + filename + ".txt","w+")
    for ll in range(0,len(info)):
        f.write("%d\n" % (info[ll]))
    for mm in range(0,len(per_disimg_fin)):  
        f.write("%f\n" % (per_disimg_fin[mm]))
    f.close()


# In[ ]:


# SEDT-3 file paths
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
    f= dirNamePD + filename+ "_pd.txt"
    np.savetxt(f,pd_array,fmt='%1.6f')

