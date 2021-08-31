#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import scipy as sp
import random
import multiprocessing as mp
import sys
import os
import glob
import gudhi as gd

# base file directory
dirName = "./TopologicalTumorShape/" # change as necessary
# read lung data image files
dirPath = os.path.join(dirName, "data/lung/image/*.Rdata")
filePath = glob.glob(dirPath)

# Directory to save results
dirNameSEDT = os.path.join(dirName+"data/lung/SEDT-3-rearranged/")
dirNamePD = os.path.join(dirName+"data/lung/PersistenceDiagram-rearranged/") 
os.makedirs(dirNamePD,exist_ok=True)

dirPath = os.path.join(dirNameSEDT+"*.txt") # change the directory of the SEDT transformed images as needed
filePath = glob.glob(dirPath)

# convert list to array
def pdarray(pd):
    pd_array=np.zeros((len(pd),3))
    for i in range(0,len(pd)):
        pd_array[i,0]=np.asarray(pd[i][0]) # dimension
        pd_array[i,1]=np.asarray(pd[i][1][0]) # birth
        pd_array[i,2]=np.asarray(pd[i][1][1]) # death
    return pd_array;

def computepd(i):
    # compute persistent homology using cubical complex
    md_cubical_complex = gd.CubicalComplex(perseus_file=filePath[i])
    # result
    md_cc_diag=md_cubical_complex.persistence()
    # convert to array
    pd_array=pdarray(md_cc_diag)
    # filename
    base=os.path.basename(filePath[i])
    filename = os.path.splitext(base)[0]
    # write txt file
    f= dirNamePD + filename+ "_pd.txt" # change the output directory as needed
    np.savetxt(f,pd_array,fmt='%1.6f')    
      
# Only run on main thread
if __name__ == '__main__':
    jobs = []
    
    # Launch processes
    for i in range(len(filePath)):
        p = mp.Process(target = computepd, args = (i,))
        jobs.append(p)
        p.start()

