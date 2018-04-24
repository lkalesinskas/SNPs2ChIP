output = open("chr19.bed", "w")
count = 0
binSize = 1000
number = 0
while count < 58617616:
    output.write("chr19" + "\t" + str(count) + "\t" + str(count + binSize) + "\t" + str(number) + "\n")
    number += 1
    count += binSize
output.close()
print(number)