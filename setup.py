import os
from setuptools import setup, find_packages

setup(
    name='geo_subsampler',
    packages=find_packages(),
    include_package_data=True,
    package_data={'geo_subsampler': [os.path.join('..', 'README.md'), os.path.join('..', 'LICENCE')]},
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
    ],
    version='0.1',
    description='Subsampling of rooted phylogenetic trees using phylogenetic diversity and location proportions.',
    author='Anna Zhukova',
    author_email='anna.zhukova@pasteur.fr',
    url='https://github.com/evolbioinfo/geo_subsampler',
    keywords=['phylogenetics', 'subsampling', 'phylogenetic diversity', 'phylogeography'],
    requires=['pandas', 'numpy', 'ete3'],
    entry_points={
        'console_scripts': [
            'geo_subsample = geo_subsampler.phylogenetic_diversity:main',
        ]
    },
)
