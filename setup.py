#!/usr/bin/env python
from setuptools import setup

setup(  name='GeneList2COVID19',
		version='1.0.0',
		description='To tell whether a list of input genes are significantly connected to the SARS-CoV-2 proteins',
		author='Jun Ding',
		author_email='jund@andrew.cmu.edu',
		url="https://github.com/phoenixding/genelist2covid19",
		license='MIT',
		packages=['GeneList2COVID19'],
		entry_points={'console_scripts':['GeneList2COVID19=GeneList2COVID19.GeneList2COVID19:main']},
		install_requires=['scipy>=1.4.1','numpy>=1.18.4','matplotlib>=3.2.1','seaborn>=0.10.1','networkx>=2.4','pandas>=1.0.3'],
		classifiers=[
			'License :: OSI Approved :: MIT License',
			'Programming Language :: Python :: 3',
		],
		)
		
