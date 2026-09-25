#f = open("visualizeInit.xyz", "r")
#g = open("angle1.5Init.xyz", "w")
#h = open("forHeightAngle1.5Init.xyz", "w")
#hh = open("removeOutside2Init.xyz", "w")
#hhh = open("visualizeInsideInit.xyz", "w")
f = open("visualize.xyz", "r")
g = open("angle1.5.xyz", "w")
h = open("forHeightAngle1.5.xyz", "w")
hh = open("removeOutside2.xyz", "w")
hhh = open("visualizeInside.xyz", "w")

g.write("11096" + "\n")
hh.write("11096" + "\n")
#g.write("\n")
line = "Lattice=\" 91.61642429 0.0 0.0 45.80821215 79.34215084 0.0 0.0 0.0 0.0 35.00000 # typea \" " 
g.write(line + "\n")

shift = 19.15-3.3

for line in f:
    a = line.split()
    if (a[1]=="1"):
       g.write("C   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       h.write("1   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       if (float(a[4])-shift)<3.3:
            hh.write("C   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
            hhh.write(line)
    elif (a[1]=="2"):
       g.write("C   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       h.write("2   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       if (float(a[4])-shift)>0.0:
           hh.write("C   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
           hhh.write(line)

f.close()
g.close()
h.close()
hh.close()
