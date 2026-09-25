How to set up the triclinic cell:
- Obtain your commensurate unit cell as you did before, including the positions of the atoms in reduced coordinates
- Input the unit cell vectors in one of the cells given in twistedBilayerGrapheneReport.ipynb
- Use the rotated unit cell vectors to set up your cell in the lammps input file (a has to be in the positive x direction)
- Use the obtained lx, ly and lz to set up your simulation box, as well as the xy, xz, yz

Here, my for lx = 134.22225299857, ly =  116.2398808502863, lz =  40.40453015592109, with an additional shift of 0.1 to avoid issues with atoms being exactly on the boundary of the simulation box (adding 0.1 to xlo, xhi and ylo, yhi):

    region mybox prism 0.0 134.22225299857 0.0 116.2398808502863 0.0 40.40453015592109 -67.1111264989935  0.0 0.0 units box
    create_box 2 mybox
    create_atoms 2 region mybox &
    basis 1 1 &
    basis 2 1 &
    ...
    
You might have to change your tilt factors if their value is not in the interval \[-lenght/2, lenght/2\].
Also note the special case were one would be exactly at the border of the interval (https://lammps.sandia.gov/threads/msg29769.html). I could avoid this problem taking advantage of machine precision. But this fix might be needed if we ever do a simulation where the box changes during the simulation.

Now, if all went well, the atoms in the output file will be the same as the ones you put in the input file, thus avoiding the problems with cutting out atoms from a rectangular cell.

WARNING: PyBinding needs to have the atoms organized per layer. However, LAMMPS messes up the order due to its domain decomposition method. This little code snippet can be used to reorder the atoms read from the (disordered) .xyz file

    x_coord, y_coord, z_coord, atom_type, l1, l2, l3 = load_ovito_lattice(name)

    zind = np.argsort(z_coord)[::-1] # [::-1] is used to reverse the order.
    x_coord = x_coord[zind]
    y_coord = y_coord[zind]
    z_coord = z_coord[zind]
    atom_type = atom_type[zind]









Extended XYZ format (usually the second line is ignored. In this PyBinding code, we use it to define the unit cell vectors).

    148
    Lattice=" 12.297500 8.519958 0.000000 -1.229750 14.909925 0.000000 0.000000 0.000000 10.000000 # type  "
    C1 1.229750 1.419989 0.000000
    C1 -0.000001 3.549979 0.000000
    C1 2.459500 3.549979 0.000000
    C1 4.919000 3.549979 0.000000
    C1 1.229749 5.679969 0.000000
    ...
    C1 x y z
    

Obsolete:

How to cut the supercell data out of the data dumped by LAMMPS (should be combined into one single script)
- Create a file visualize.dat and visualizeInit.dat with the final and initial positions, respectively, extracted from the dump.mimimization file (better to this on the server to avoid saving too much data locally).
- Run extractXYZ.py on both of these files.
- Copy the angleXXX.xyz data into a template input.cart (or inputInit.cart) file (also modify the unit cell vector information and number of atoms)
- Run fracToCart.py with input.cart as argument input
- Run cutCell.py or cutCellInit.py and check if number of output atoms matches with the number of expected atoms. If not, tincker with the cutoff variables where the difference between cutoff and cutoff2 should be around 1. 1.2-0.2=1.0 or 1.0-0.0=1.0 are all ok combinations. Sometimes, after relaxation, you need to play with the decimals. A better way should be found for this, but if the relaxation is done perfectly and very robustly, one should expect it too work.
- Read the unitCell.xyz or unitCellInit.xyz into the ipynb pybinding code.

Solved (by using the triclinic cell):
- properly cutting out the hexagonal supercell from the rectangular cell from LAMMPS. This works for the unrelaxed system, but fails after relaxation. Especially when considering GG instead of GBN. The trick for GBN is to use box -2 2 instead of box 0 4. Let's check if that also works for GG.
- solve the relaxation issues I am getting for twisted GG. I am getting large deformations that don't seem to necesssarily have the periodicity of the supercell (_box_ in lammps).
