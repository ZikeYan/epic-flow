#! /usr/bin/env python2
#-*- coding: UTF-8 -*-
# Requirements: Numpy as PIL/Pillow
import os
import os.path
import sys
import time
import numpy as np
import subprocess

TAG_FLOAT = 202021.25
TAG_CHAR = 'PIEH'


time1 = time.time()
filepath_img = sys.argv[1]
filepath_edge = sys.argv[2]
filepath_stereo = sys.argv[3]
filepath_flow = sys.argv[4]
filepath_output = sys.argv[5]

if not os.path.exists(filepath_img)  or not os.path.exists(filepath_edge) or not os.path.exists(filepath_stereo)  or not os.path.exists(filepath_flow):
    print "The path does not exist"

else:
    dirlist_img = []
    # ----------get the folder list
    nameFolds = os.listdir(filepath_img)
    for name in nameFolds:
        dir = os.path.join(filepath_img, name)
        if os.path.isdir(dir):
            dirlist_img.append(name)
    # ----------traverse all the folders
    dirlist_img.sort()
    for l in dirlist_img:
        imgFile = []
        #----------traverse all the files
        for a in os.listdir(filepath_img + '/' + l):
            if a.endswith(".png"):
                imgFile.append(os.path.splitext(a)[0])
        imgFile.sort()
        img1 = imgFile[0:len(imgFile)-1]
        for i in img1:
            imgFile_1 = filepath_img + "/" + l + "/" + i + ".png"
            imgFile_2 = filepath_img + "/" + l + "/" + imgFile[img1.index(i)+1] + ".png"
            edgeFile = filepath_edge + "/" + l + "/" + i
            stereoFile = filepath_stereo + "/" + l + "/" + i + ".flo"
            flowFile = filepath_flow + "/" + l + "/" + i + ".flo"
            outputFile = filepath_output + "/" + l + "/" + i + ".flo"
            #print os.path.exists(imgFile_1), os.path.exists(imgFile_2), os.path.exists(edgeFile), os.path.exists(matchFile), os.path.exists(outputFile)
            ######!!!!!! when shell = True, subprocess.Popen() needs a char file as 'python server.py'. when shell = False, it needs a list as ['python', 'server.py']
            tmp = " ".join(["./epicflow", imgFile_1, imgFile_2, edgeFile, stereoFile, flowFile, outputFile])
            os.system(tmp)
        print "calculated folder:", l
    print '%.2f' %(time.time()-time1), "seconds wall time"
