from setuptools import setup

setup(name='PV-DER',
      version=open("pvder/_version.py").readlines()[-1].split()[-1].strip("\"'"),
      packages=['pvder',],
      description='Utility for simulating PV-DER',
      author = 'Siby Jose Plathottam',
      author_email='sibyjackgrove@gmail.com',
      license= 'LICENSE.txt',
      install_requires=['scipy>=1.0.0','numpy>=1.15.1','matplotlib==2.0.2'],#And any other dependencies required	  
      )