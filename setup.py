from setuptools import find_packages, setup

from wcx2cytosure.__version__ import __version__

setup(
	name='wcx2cytosure',
	version=__version__,
	author='Denise Ravinale',
	author_email='denise.ravinale@regionstockholm.se',
	url='https://github.com/denrav99/wcx2cytosure.git',
	description='Convert WCX output with structural variations to CytoSure format',
	license='MIT',
	packages=find_packages(exclude=["tests/", "dist/", "build/"]),
	install_requires=['lxml','pandas', 'cyvcf2'],
	entry_points={'console_scripts': ['wcx2cytosure=wcx2cytosure.wcx2cytosure:main']},
	classifiers=[
		"Development Status :: 3 - Alpha",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
