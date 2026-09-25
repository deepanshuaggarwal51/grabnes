#f = open("visualizeInit.xyz", "r")
#g = open("angle1.08Init.xyz", "w")
#h = open("forHeightAngle1.08Init.xyz", "w")
#hh = open("removeOutside2Init.xyz", "w")
f = open("visualize.xyz", "r")
g = open("angle1.08.xyz", "w")
h = open("forHeightAngle1.08.xyz", "w")
#hh = open("removeOutside2.xyz", "w")
hh = open("removeOutside.xyz", "w")

g.write("16746" + "\n")
hh.write("16746" + "\n")
#g.write("\n")
line = "Lattice=\" 57.5345566 0.0 0.0 28.76727829 49.82638759 0.0 0.0 0.0 0.0 35.00000 # typea \" " 
g.write(line + "\n")

shift = 19.15-3.3

for line in f:
    a = line.split()
    if (a[1]=="1"):
       g.write("C   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       h.write("1   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       if (float(a[4])-shift)<3.3:
            hh.write("1   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
    elif (a[1]=="2"):
       g.write("C   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       h.write("2   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       if (float(a[4])-shift)>0.0:
           hh.write("2   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")

f.close()
g.close()
h.close()
hh.close()
