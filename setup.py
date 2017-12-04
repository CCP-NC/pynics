from setuptools import setup

setup(name='pynics',
      version='0.1.0',
      packages=['pynics'],
      entry_points={
          'console_scripts': [
              'lorentz_buildup_nics = pynics.__main__:nics_buildup'
          ]
      },
      )
