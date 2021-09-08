import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(name='pvder',
      version=open("pvder/_version.py").readlines()[-1].split()[-1].strip("\"'"),
      packages=['pvder',],
      include_package_data=True,
      description='Utility for simulating PV-DER',
      long_description=README,
      long_description_content_type="text/markdown",
      url ='https://github.com/tdcosim/SolarPV-DER-simulation-tool',
      author = 'Siby Jose Plathottam',
      author_email='sibyjackgrove@gmail.com',
      license= 'LICENSE.txt',
      classifiers=[
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
      ],
      install_requires=['scipy>=1.0.0','numpy>=1.15.1','matplotlib>=2.0.2','numba>=0.53.0'],#And any other dependencies required
      extras_require={"docs": ['sphinx-rtd-theme','nbsphinx','nbsphinx-link']}
      )
