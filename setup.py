from setuptools import setup

setup(
    name='molecule_rxn',
    version='0.0.1',    
    description='A molecule process Python package',
    url='https://github.com/cchang373/molecule_rxn',
    author='Chaoyi Chang',
    author_email='cchang373@gatech.edu',
    #license='',
    packages=['molecule_rxn'],
    install_requires=['networkx',
                      'imolecule',
                      'ase',
                      'bitarray'
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        #'License :: OSI Approved :: BSD License',  
        'Programming Language :: Python :: 3.6',
    ],
)
