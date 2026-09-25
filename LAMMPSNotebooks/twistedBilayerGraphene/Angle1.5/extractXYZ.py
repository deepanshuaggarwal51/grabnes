f = open("visualizeInit.xyz", "r")
g = open("angle1.5Init.xyz", "w")
h = open("forHeight1.5Init.xyz", "w")
#f = open("visualize.xyz", "r")
#g = open("angle1.5.xyz", "w")
#h = open("forHeight1.5.xyz", "w")

g.write("11908" + "\n")
#g.write("\n")
line = "Lattice=\" 187.49331828 0.0 0.0 93.74665914 162.37397667 0.0 0.0 0.0 0.0 35.00000 # typea \" " 
g.write(line + "\n")

shift = 19.15

for line in f:
    a = line.split()
    if (a[1]=="1"):
       g.write("C1   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       h.write("1   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
    elif (a[1]=="2"):
       g.write("C2   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")
       h.write("2   " + a[2] + "\t" + a[3] + "\t" + str(float(a[4])-shift) + "\n")

f.close()
g.close()
h.close()
