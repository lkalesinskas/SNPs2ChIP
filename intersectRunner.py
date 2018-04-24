from os import listdir
from os.path import isfile, join
import os
folderPath = "/home/larry/Desktop/ChipChip/Blastoid_5_filtered_19/"
keyFile = "/home/larry/Desktop/ChipChip/chr19.bed"
outputDirectory = "/home/larry/Desktop/ChipChip/FinalBed/"
files = [f for f in listdir(folderPath) if isfile(join(folderPath, f))]
for fileName in files:
    fullPath = folderPath + fileName
    command = "bedtools intersect -a " + fullPath + " -b " + keyFile + " -wb > " + outputDirectory + "final_" + fileName
    os.system(command)