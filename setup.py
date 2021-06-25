from setuptools import setup, find_packages

setup(name='pynics',
      version='0.1.0',
      packages=find_packages(),
      install_requires=[
          'numpy',
          'ase'
      ],      
      entry_points={
          'console_scripts': [
              'lorentz_buildup_nics = pynics.__main__:nics_buildup'
          ]
      },
      include_package_data=True,
)
