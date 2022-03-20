# pynics
## Python script for computing Nuclear Independent Chemical Shifts and buildup functions from Castep data

This script was created to carry out the calculations in the "The Lorentz Sphere visualised" paper. 

#### How to install

Numpy and ASE are requirements. Install with pip:

    pip install ./ --user

in the pynics folder.

#### How to use (lorentz_buildup_nics)

1. Run a Castep NMR calculation on the system of interest and store the currents in a binary file - this is done by using the 
`MAGRES_WRITE_RESPONSE: TRUE` parameter in the `.param` file;
2. Create a `.nicslist` file, formatted like the `.cell` file, with a single `NICS_POINTS_FRAC` block, listing the
fractional coordinates of all points for which the NICS buildup functions are desired, or a `NICS_POINTS_ABS` block, which can use cartesian coordinates.
3. Run the command line tool:

    `lorentz_buildup_nics <seedname>`
    
with the seedname of the system, in the same folder as the Castep calculation.
Use `lorentz_buildup_nics -h` for a list of additional command line parameters.

#### How to use (nicsanalyse)

1. A single unit cell of the system of interest should be geometry optimized using Castep.
2. This unit cell should then be propagated to a supercell such that there is a separation $>10\\A$ between a central molecule and the next cell.
3. The `splitcell` utility can then be run on this, with the following options

```
python splitcell.py 	-struct <structure file> \
				-onemol <onemol output file> \
				-nomol <nomol output file> \
				-supercell <supercell output file> \
				-nicslist <file to output nicslist> 
```
4. The output cell files should then have Castep GIPAW calculations done performed on them, following the same steps as above for Pynics (e.g., pass the `MAGRES_WRITE_RESPONSE:True` parameter). The `.nicslist` should be the output nicslist of `splitcell.py`.
5. On the resulting current files, run the `nicsanalyse` utility as

```
python nicsanalyse.py -nicslist <path to nicslist> \
                      -nomol_current <path to nomol current> \
                      -nomol_magres <path to nomol magres> \
                      -onemol_current <path to onemol current> \
                      -onemol_magres <path to onemol magres> \
                      -supercell_current <path to supercell current> \ 
                      -supercell_magres <path to supercell magres> \
                      -output <output file> 
                      ```
