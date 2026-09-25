f = open("inputInit.xyz", "r")
g = open("inputInitShifted.xyz", "w")

for line in f:
    a = line.split()
    if len(a)>0:
        if a[0]=="C":
                g.write(a[0] + "  " + a[1] + "   " + a[2] + "   " + str(float(a[3])-22.107050) + "\n") 
        else:
            g.write(line)
    else:
	g.write(line)
