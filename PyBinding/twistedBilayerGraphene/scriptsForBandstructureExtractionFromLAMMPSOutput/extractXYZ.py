#f = open("visualizeInit.xyz", "r")
#g = open("angle6.01Init.xyz", "w")
f = open("visualize.xyz", "r")
g = open("angle6.01.xyz", "w")

g.write("27208" + "\n")
g.write("\n")

for line in f:
    a = line.split()
    if (a[1]=="1"):
       g.write("C   " + a[2] + "  " + a[3] + "  " + a[4] + "\n")
    elif (a[1]=="2"):
       g.write("B   " + a[2] + "  " + a[3] + "  " + a[4] + "\n")
    elif (a[1]=="3"):
       g.write("N   " + a[2] + "  " + a[3] + "  " + a[4] + "\n")

f.close()
g.close()
