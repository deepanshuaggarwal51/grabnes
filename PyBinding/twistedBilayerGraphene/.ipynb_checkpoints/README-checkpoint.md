There are two issues that I need to solve:
- properly cutting out the hexagonal supercell from the rectangular cell from LAMMPS. This works for the unrelaxed system, but fails after relaxation. Especially when considering GG instead of GBN. The trick for GBN is to use box -2 2 instead of box 0 4. Let's check if that also works for GG.
- solve the relaxation issues I am getting for twisted GG. I am getting large deformations that don't seem to necesssarily have the periodicity of the supercell (_box_ in lammps).






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
