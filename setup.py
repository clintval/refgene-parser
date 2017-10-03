import setuptools

from packagename.version import Version


setuptools.setup(
    name='refgene_parser',
    version=Version('1.0.0').number,
    description='Python parser to RefGene files',
    long_description=open('README.md').read().strip(),
    author='clintval',
    author_email='valentine.clint@gmail.com',
    url='https://github.com/clintval/refgene-parser',
    py_modules=['refgene_parser'],
    install_requires=[

    ],
    license='MIT License',
    zip_safe=False,
    keywords='refgene bioinformatics',
)
