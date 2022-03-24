from setuptools import setup, find_packages

"""
pynics - Computing Nuclear Independent Chemical Shifts and buildup functions using python
by Simone Sturniolo and Ben Tatman 
"""

long_description = """
Python scripts for computing Nuclear Independent Chemical Shifts and buildup functions from Castep data """

if __name__ == "__main__":

    setup(
        name="pynics",
        version='0.1.1',
        description="Computing Nuclear Independent Chemical Shifts",
        long_description=long_description,
        url="https://www.ccpnc.ac.uk/software/",
        author="Simone Sturniolo",
        author_email="simone.sturniolo@stfc.ac.uk",
        license="MIT",
        classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            "Development Status :: 4 - Beta",
            # Indicate who your project is intended for
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Scientific/Engineering :: Physics",
            # Pick your license as you wish (should match "license" above)
            "License :: OSI Approved :: MIT License",
            # Specify the Python versions you support here. In particular,
            # ensure that you indicate whether you support Python 2, Python 3
            # or both.
            "Programming Language :: Python :: 3",
        ],
        keywords=["crystallography", "ccpnc", "computational chemistry", "nmr", "nics"],
        packages=find_packages(),
        entry_points={
          'console_scripts': [
              'lorentz_buildup_nics = pynics.__main__:nics_buildup',
              'nicsanalyse = pynics.nicsanalyse:nics_analyse',
              'splitcell = pynics.splitcell:split_cell'
          ]
        },
        # Requirements
        install_requires=[
          'numpy',
          'ase',
          'soprano>=0.8.10'
        ],
        python_requires=">=3.6.*",
    )