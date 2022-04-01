# pynics
## Python scripts for computing Nuclear Independent Chemical Shifts and buildup functions from Castep data

The chemical shift at which a nuclear environment appears in solid-state NMR is dependent not only on the local intramolecular interactions affecting that environment, but also on interactions with longer range features such as ring currents and hydrogen bonding. Utilising DFT GIPAW calculations, it is possible to partition the contributions of these different features to chemical shift. Further, it is possible to observe the spatial effects of these contributions, and the distance at which these contributions arise.

See the below 'how-to's for how pynics can perform such analysis and see the CCP-NC website for the full protocol: https://www.ccpnc.ac.uk/



### How to install

pynics is available on the Python Package Index. You can install the latest stable release by using `pip`:

    pip install pynics --user

It's recommended to use a virtual environment to reduce the chance of conflicts with the rest of your environment. With conda you could do:

 1. Create conda virtual environment   
    `conda create --name pynicsenv`
 2. Activate the conda environment   
    `conda activate pynicsenv`
 3. Make sure pip is installed   
    `conda install pip`
 4. pip-install pynics with it's dependencies.   
    `pip install pynics`

You could also use virtualenv to manage your virtual environments. 

Alternatively, you can get the latest version (not guaranteed to be stable) from github:

    git clone https://github.com/CCP-NC/pynics.git
    pip install ./pynics --user
 
This approach should work even on machines for which one does not possess admin privileges (such as HPC clusters), as long as Python and `pip` are present.
Numpy, ASE and soprano are requirements. Install with pip:

    pip install ./ --user

in the pynics folder.

### How to use (lorentz_buildup_nics)

1. Run a Castep NMR calculation on the system of interest and store the currents in a binary file - this is done by using the 
`MAGRES_WRITE_RESPONSE: TRUE` parameter in the `.param` file;
2. Create a `.nicslist` file, formatted like the `.cell` file, with a single `NICS_POINTS_FRAC` block, listing the
fractional coordinates of all points for which the NICS buildup functions are desired, or a `NICS_POINTS_ABS` block, which can use cartesian coordinates.
3. Run the command line tool:

    `lorentz_buildup_nics <seedname>`
    
with the seedname of the system, in the same folder as the Castep calculation.
Use `lorentz_buildup_nics -h` for a list of additional command line parameters.

### How to use (nicsanalyse)

1. A single unit cell of the system of interest should be geometry optimized using Castep.
2. This unit cell should then be propagated to a supercell such that there is a separation > 10 Å between a central molecule and the next cell.
3. The `splitcell` utility can then be run on this, with the following options

```
splitcell 	--struct <structure file> \
		--nicslist <file to output nicslist>
```

4. The output cell files should then have Castep GIPAW calculations done performed on them, following the same steps as above for lorentz_buildup_nics (e.g., pass the `MAGRES_WRITE_RESPONSE:True` parameter). The `.nicslist` should be the output nicslist of `splitcell.py`.
5. On the resulting current files, run the `nicsanalyse` utility as

```
nicsanalyse        --nicslist <path to nicslist> \
		   --nomol_magres <path to nomol magres> \
		   --onemol_magres <path to onemol magres> \
		   --supercell_magres <path to supercell magres> \
		   --output <output file> 
```

### Cite
This library was originally created to carry out the calculations in the "The Lorentz Sphere visualised" paper. If you use this code in your work, please cite:

S. Sturniolo and J. R. Yates , "The Lorentz sphere visualised", J. Chem. Phys. 150, 094103 (2019) https://doi.org/10.1063/1.5080298

```
@article{doi:10.1063/1.5080298,
	author = {Sturniolo,S.  and Yates,J. R. },
	title = {The Lorentz sphere visualised},
	journal = {The Journal of Chemical Physics},
	volume = {150},
	number = {9},
	pages = {094103},
	year = {2019},
	doi = {10.1063/1.5080298},
	URL = {https://doi.org/10.1063/1.5080298}
}
```
