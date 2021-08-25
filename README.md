# pynics
## Python script for computing Nuclear Independent Chemical Shifts and buildup functions from Castep data

This script was created to carry out the calculations in the "The Lorentz Sphere visualised" paper. 

#### How to install

Numpy and ASE are requirements. Install with pip:

    pip install ./ --user

in the pynics folder.

#### How to use

1. Run a Castep NMR calculation on the system of interest and store the currents in a binary file - this is done by using the 
`MAGRES_WRITE_RESPONSE: TRUE` parameter in the `.param` file;
2. Create a `.nicslist` file, formatted like the `.cell` file, with a single `NICS_POINTS_FRAC` block, listing the
fractional coordinates of all points for which the NICS buildup functions are desired, or a `NICS_POINTS_ABS` block, which can use cartesian coordinates.
3. Run the command line tool:

    `lorentz_buildup_nics <seedname>`
    
with the seedname of the system, in the same folder as the Castep calculation.
Use `lorentz_buildup_nics -h` for a list of additional command line parameters.
