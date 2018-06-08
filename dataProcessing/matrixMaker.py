from os import listdir
from os.path import isfile, join
import numpy as np
folderPath = "C:\\Users\\Larry\\PycharmProjects\\ChipChip\\FinalBed\\"
files = [f for f in listdir(folderPath) if isfile(join(folderPath, f))]
#Number of files
sampleNumber = len(files)
#by number of possible bins
binNumber = 58617

ourMatrix = np.zeros((sampleNumber, binNumber), dtype=bool)
print(ourMatrix.shape)

count = 0
for fileName in files:
    fullPath = folderPath + fileName
    f = open(fullPath, "r")
    for line in f:
        split = line.split()
        binToAdd = int(split[len(split)-1])
        ourMatrix[count][binToAdd] = True
    count += 1
    f.close()

ourMatrix = ourMatrix.astype(int)
print(ourMatrix[0][60:80])

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

pca = PCA(n_components=20)
pca.fit(ourMatrix)
print(pca.components_)
print(pca.explained_variance_)

outputFile = open("output.txt", "w")
for i in range(0, 58617):
    outputFile.write(str(i) + "," + str(pca.components_[0][i]) + "\n")
outputFile.close()




pca = PCA().fit(ourMatrix)
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('number of components')
plt.ylabel('cumulative explained variance');
plt.show()