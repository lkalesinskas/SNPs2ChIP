from os import listdir
from os.path import isfile, join
import numpy as np
folderPath = "C:\\Users\\Larry\\PycharmProjects\\ChipChip\\Blastoid_5\\"
outputFolder = "C:\\Users\\Larry\\PycharmProjects\\ChipChip\\Blastoid_5_filtered_score\\"
onlyfiles = [f for f in listdir(folderPath) if isfile(join(folderPath, f))]

def chromosomeFilter(folderPath, files, outputFolder, chromosome):
    count = 0
    for fileName in files:
        fullPath = folderPath + fileName
        outputFileName = outputFolder + "filtered_chr_" + fileName
        f = open(fullPath, "r")
        outputFile = open(outputFileName, "w")
        for line in f:
            arraySplit = line.split("\t")
            if(arraySplit[0] == chromosome):
                outputFile.write(line)
        outputFile.close()
        f.close()
        print(str(count) + "/" + str(len(files)))
        count+=1

def scoreFilter(folderPath, files, outputFolder, topNumber):
    count = 0
    for fileName in files:
        fullPath = folderPath + fileName
        outputFileName = outputFolder + "filtered_score_" + fileName
        outputFile = open(outputFileName, "w")
        data = np.genfromtxt(fullPath, dtype=None)
        sortedData = sorted(data, key=lambda x: x[4], reverse=True)
        finalSet = sortedData[0:topNumber]
        for item in finalSet:
            for element in item:
                outputFile.write(str(element) + "\t")
            outputFile.write("\n")
        outputFile.close()
        print(str(count) + "/" + str(len(files)))
        count+=1

scoreFilter(folderPath, onlyfiles, outputFolder, 1000)