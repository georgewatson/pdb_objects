import setuptools

with open("README.md", 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='pypdb',
    version='0.1.0',
    author='George Watson',
    author_email='george@georgewatson.me',
    description='Object-oriented processing of Protein Databank (PDB) files',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/georgewatson/pypdb',
    license='MIT',
    packages=setuptools.find_packages(),
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    include_package_data=True,
    zip_safe=False)
