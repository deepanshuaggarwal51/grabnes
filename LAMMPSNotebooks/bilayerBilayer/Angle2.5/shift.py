f = open("BLBL.cart", "r")
g = open("angle2.5Init", "w")



shift = 19.15-3.3

for line in f:
    a = line.split()
    if len(a)>0:
        if a[0]=="C" or a[0]=="C1" or a[0]=="C2":
                g.write(a[0] + "  " + a[1] + "   " + a[2] + "   " + str(float(a[3])-shift) + "\n")
