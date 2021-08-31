#!/usr/bin/env python
# coding: utf-8

# In[6]:


#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import scipy as sp
import random
from scipy.ndimage import distance_transform_bf
import multiprocessing as mp
import sys
import os
import pyreadr
import glob

# base file directory
dirName = "./TopologicalTumorShape/" # change as necessary
# read lung data image files
dirPath = os.path.join(dirName,"data/lung/image/*.Rdata")
filePath = glob.glob(dirPath)
# Directory to save results
dirNameSEDT = os.path.join(dirName+"data/lung/SEDT-3-rearranged/")
os.makedirs(dirName,exist_ok=True)

# Define function that can be run in parallel
def computedist(iter):
    np.random.seed(iter) # seed

    for i in range(0,len(filePath)):
        imgdat = pyreadr.read_r(filePath[i])
        idf = imgdat["dataset"]
        idf2 = np.reshape(idf.values, (-1, 1))
        np.random.shuffle(idf2) # rearrange pixels
        idfrdm = np.reshape(idf2,idf.shape)
        idfrdm_pd = pd.DataFrame(idfrdm)

        # normal cell: 0
        # tumor cell: 1
        # empty: 2
        dfnormal=idfrdm_pd.copy()
        dftumor=idfrdm_pd.copy()

        matdf = idfrdm_pd.values

        # for normal, normal -1, otherwise: 0
        matnormal = dfnormal.values
        # first make tumor or empty regions 1 and plug in 1-matnormal
        ## tumor or empty regions
        matnormal[matnormal>1]=1


        # for tumor, tumor: 1, otherwise:0
        mattumor = dftumor.values
        ## make empty regions 0 (only tumors are 1)
        mattumor[mattumor>1]=0

        # current version
        # negative: tumor
        distimgn=distance_transform_bf(mattumor,metric='euclidean')
        # positive: normal
        distimgp=distance_transform_bf(1-matnormal,metric='euclidean')

        distimgp = distimgp.astype(np.float64)
        distimgn = distimgn.astype(np.float64)

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

        # write txt file (for computing persistent homology)
        f= open(dirNameSEDT + filename + "_" + str(iter) + ".txt","w+")
        for ll in range(0,len(info)):
            f.write("%d\n" % (info[ll]))
        for mm in range(0,len(per_disimg_fin)):  
            f.write("%f\n" % (per_disimg_fin[mm]))
        f.close()

# Only run on main thread
if __name__ == '__main__':
    jobs = []

    # Launch processes
    for pind in range(50):
        p = mp.Process(target = computedist, args = (pind,))
        jobs.append(p)
        p.start()

