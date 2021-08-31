#!/usr/bin/env python
# coding: utf-8

# In[3]:


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

# Download dataset from https://github.com/lorinanthony/SECT
dirNameSECT = "./SECT/" # change the SECT directory as necessary
# ID for patients
patName = os.listdir(dirNameSECT+"Data/MITKSegmentations")
# folder for patients
dirName = glob.glob(dirNameSECT+"Data/MITKSegmentations/*/")
# Directory that saved SEDT-2 results
dirNameSEDT = os.path.join(dirNameSECT+"data/brain/SEDT-2-rearranged/")
# Directory to save persistent homology results
dirNamePD = os.path.join(dirNameSECT+"data/brain/PersistenceDiagram-rearranged/") 
os.makedirs(dirNamePD,exist_ok=True)

dirPath = os.path.join(dirNameSEDT+'*.txt') 
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
    f= dirNamePD + filename+ "_pd.txt"
    np.savetxt(f,pd_array,fmt='%1.6f')    
      
# Only run on main thread
if __name__ == '__main__':
    jobs = []
    
    # Launch processes
    for i in range(len(filePath)):
        p = mp.Process(target = computepd, args = (i,))
        jobs.append(p)
        p.start()

