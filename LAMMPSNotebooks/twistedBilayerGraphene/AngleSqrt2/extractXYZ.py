f = open("visualizeInit.xyz", "r")
g = open("angle1.08Init.xyz", "w")
h = open("forHeight1.08Init.xyz", "w")
hhh = open("differentIndexPerLayerInit.xyz", "w")
#f = open("visualize.xyz", "r")
#g = open("angle1.08.xyz", "w")
#h = open("forHeight1.08.xyz", "w")
#hhh = open("differentIndexPerLayer.xyz", "w")

g.write("11908" + "\n")
#g.write("\n")
line = "Lattice=\" 134.222253 0.0 0.0 67.1111265 116.23988085 0.0 0.0 0.0 0.0 35.00000 # typea \" " 
g.write(line + "\n")

shift = 19.15-3.3

for line in f:
    a = line.split()
    if (a[1]=="1"):
       g.write("C1   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       h.write("1   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       if (float(a[4])-shift)<3.3:
            hhh.write(a[0] + "  " + "1   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
    elif (a[1]=="2"):
       g.write("C2   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       h.write("2   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       if (float(a[4])-shift)>0.0:
           hhh.write(a[0] + "  " + "2   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")

f.close()
g.close()
h.close()
