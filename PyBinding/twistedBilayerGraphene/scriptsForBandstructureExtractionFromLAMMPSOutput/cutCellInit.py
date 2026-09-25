def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

#f = open("input.frac", "r")
#f2 = open("input.cart", "r")
f = open("inputInit.frac", "r")
f2 = open("inputInit.cart", "r")

g = open("unitCellInit.xyz", "w")


#cutoff = 0.986
#cutoff2 = -0.048
cutoff = 1.0
cutoff2 = 0.0
#cutoff = 1.5
#cutoff2 = 0.5
counter = 0

#shift = 22.107050
#shift = 0.0
#shift = 19.2818
shift = 19.15

maxZ = -1000

for line, line2 in zip(f,f2):
    a = line.split()
    b = line2.split()
    if len(a)>0:
        if a[0]=="C1" or a[0]=="C2" or a[0]=="C" or a[0]=="B" or a[0]=="N":
            x = float(a[1])
            y = float(a[2])
            z = float(a[3])
            if (x < cutoff and y < cutoff and z < cutoff and x > cutoff2 and y>cutoff2 and z>cutoff2):
                #g.write(line2)
                g.write(b[0] + "  " + b[1] + "   " + b[2] + "   " + str(float(b[3])-shift) + "\n") 
                #g.write("C  " + b[1] + "   " + b[2] + "   " + str(float(b[3])-shift) + "\n") 
                counter = counter+1
                if float(b[3])>maxZ:
                    maxZ = float(b[3])

print(maxZ)

print("number of atoms in unit cell: ", counter)
f.close()
f2.close()
g.close()

line = "Lattice=\" 18.450000000 6.391267480 0.00000000 3.690000000 19.173802440 0.000000000  0.0  0.0  35.000000 # type \" "
line_prepender("unitCellInit.xyz", line)
line_prepender("unitCellInit.xyz", "    " + str(counter))

