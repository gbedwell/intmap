from setuptools import setup

setup(
	name='intmap',
	version='1.0.0',
	author='Greg Bedwell',
	author_email='gregoryjbedwell@gmail.com',
	packages=['intmap','intmap-demux','intmap-multi'],
	python_requires='>=3.10',
	scripts=['intmap/intmap','intmap-demux/intmap-demux','intmap-multi/intmap-multi'],
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
        'datasketch'
        ]
    )