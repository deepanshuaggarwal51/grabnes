f = open("visualizeInit.xyz", "r")
g = open("angleSqrt2Init.xyz", "w")
h = open("forHeightAngleSqrt2Init.xyz", "w")
hh = open("removeOutside2Init.xyz", "w")
#f = open("visualize.xyz", "r")
#g = open("angleSqrt2.xyz", "w")
#h = open("forHeightAngleSqrt2.xyz", "w")
##hh = open("removeOutside2.xyz", "w")
#hh = open("removeOutside.xyz", "w")

g.write("9942" + "\n")
hh.write("9942" + "\n")
#g.write("\n")
line = "Lattice=\" 100.137412 0.0 0.0 50.0687058 86.72154231 0.0 0.0 0.0 0.0 35.00000 # typea \" " 
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
