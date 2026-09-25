Steps to get the flat structure in LAMMPS (assuming that is the most stable configuration) as well as get symmetric displacements to get there.

Note that for single layer graphene, the relaxed lattice constant is 2.46019.

- Start with this value to build your commensurate supercell. Relax the system using `fix 1 all box/relax iso 0.0 vmax 0.001` to allow for the simulation to change size as well. Indeed, the optimal cell size changes with the ratio of AA and AB stacking in the moire cell. Be careful. If you put the two layers too far from each other, they might not interact and you will not get the optimal cell size that leads to AA/AB reconstruction. Yet, too close might lead to unnecessary instabilities.
- Extract the relaxed cell size parameters from your dump.minimization file. Use the following code snippet to obtain the new lx, ly, lz, xy, etc.

```
(xlo_bound, xhi_bound, xy) = (-2.9933985680883461e+02, 5.9865741656547743e+02, -2.9933242445800494e+02)
(ylo_bound, yhi_bound, xz) = (-6.4366046283187948e-03, 5.1845253090935398e+02, 0.0000000000000000e+00)
(zlo_bound, zhi_bound, yz) = (-4.3453150571666401e-04, 3.5000434531505725e+01, 0.0000000000000000e+00)


xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
ylo = ylo_bound - min(0.0,yz)
yhi = yhi_bound - max(0.0,yz)
zlo = zlo_bound
zhi = zhi_bound

#print("abc:", a, b, c)
lx = xhi-xlo
ly = yhi-ylo
lz = zhi-zlo

a = lx
b = np.sqrt(ly**2 + xy**2)
c = np.sqrt(lz**2 + xz**2 + yz**2)
alpha = np.arccos((xy*xz + ly*yz)/(b*c))
beta = np.arccos(xz/c)
gamma = np.arccos(xy/b)


print("------------------ Parameters to be used by LAMMPS ------------------")
print("lx = ", lx)
print("ly = ", ly)
print("lz = ", lz)
print("xy = ", xy)
print("xz = ", xz)
print("yz = ", yz)
print("\n")

print("------------------ Conventional triclinic cell notation ------------------")
print("a= ", a)
print("b= ", b)
print("c= ", c)
print("alpha: ", alpha*rad)
print("beta: ", beta*rad)
print("gamma: ", gamma*rad)
print("\n")
```
- The new lattice constant will be different than the one for ideal single layer graphene. Obtain this value using the following code snippet.

```
a2=  598.6962214082658 # this is the value you got for the supercell built using lattice vector 2.46019

ratio = a2/a 
print(ratio)

print(2.46019/ratio)
```

- Use this new simulation cell to do a new relaxation (allow the cell the change) with the atoms located at the non-relaxed positions for this simulation box (don't use the output from the previous relaxation... if you do, your system will be flat, but the displacements are not symmetric).
- Do at least one continuation run (still allow the cell to change) based on the output of previous calculation. Check you have reached the same total energy minimum as the one in the first simulation. This continuation run should be flat and symmetric, as well as give the lowest energy.
