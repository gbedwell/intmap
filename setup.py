from setuptools import setup, find_packages

setup(
    name='intmap',
    version='2.0.0',
    author='Greg Bedwell',
    author_email='gregoryjbedwell@gmail.com',
    packages=find_packages(),
    python_requires='>=3.11,<3.12',
    scripts=[
        'intmap/intmap',
        'intmap_demux/intmap_demux',
        'intmap_multi/intmap_multi'
    ],
    url='',
    license='LICENSE.txt',
    description='A CLI app for mapping genomic integration sites from NGS data.',
    include_package_data = True,
    install_requires=[
        'pybloom-live==4.0.0',
        'ruptures==1.1.10',
        'pyranges==0.1.4',
        'mmh3==5.2.1'
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