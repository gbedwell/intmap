from setuptools import setup, find_packages

setup(
	name='intmap',
	version='1.0.0',
	author='Greg Bedwell',
	author_email='gregoryjbedwell@gmail.com',
	packages=find_packages(),
	python_requires='>=3.10',
	scripts=[
    	'intmap/intmap',
     	'intmap_demux/intmap_demux',
      	'intmap_multi/intmap_multi', 
       	'intmap_se/intmap_se'
    ],
	url='',
	license='LICENSE.txt',
	description='A CLI app for mapping retroviral integration sites from NGS data.',
	include_package_data = True,
	install_requires=[
		'biopython',
		'joblib',
		'pysam',
		'RapidFuzz',
		'numpy',
		'faiss-cpu',
		'datasketch',
		'pybloom-live'
	],
 	setup_requires=[
        "pytest-runner",
    ],
    tests_require=[
        "pytest",
    ],
    extras_require={
        'test': [
            'pytest',
            'pytest-cov',
            ],
	}
)