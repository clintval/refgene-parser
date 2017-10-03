import setuptools

setuptools.setup(
    name='refgene_parser',
    packages=['refgene_parser'],
    version='0.0.1',
    description='Quickly parse genes and exons from a RefGene file',
    long_description=open('README.md').read().strip(),
    author='clintval',
    author_email='valentine.clint@gmail.com',
    url='https://github.com/clintval/refgene-parser',
    download_url='https://github.com/clintval/refgene-parser/archive/0.0.1.tar.gz',  # noqa
    py_modules=['refgene_parser'],
    install_requires=[],
    license='MIT License',
    zip_safe=True,
    keywords='refgene parser bioinformatics',
    classifiers=[],
)
