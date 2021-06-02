from setuptools import setup

setup(name='pvder',
      version=open("pvder/_version.py").readlines()[-1].split()[-1].strip("\"'"),
      packages=['pvder',],
      description='Utility for simulating PV-DER',
      url ='https://github.com/tdcosim/SolarPV-DER-simulation-tool',
      author = 'Siby Jose Plathottam',
      author_email='sibyjackgrove@gmail.com',
      license= 'LICENSE.txt',
      classifiers=[
        'License :: OSI Approved :: BSD-3 License',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
      ],
      install_requires=['scipy>=1.0.0','numpy>=1.15.1','matplotlib>=2.0.2'],#And any other dependencies required
      extras_require={"docs": ['sphinx-rtd-theme','nbsphinx','nbsphinx-link']}
      )
