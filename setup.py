from setuptools import setup

setup(
    name='HOSE_code_generator',
    url='https://github.com/Ratsemaat/HOSE_code_generator',

    # Needed to actually package something
    packages=['hosegen'],

    # Needed for dependencies
    install_requires=['rdkit'],

    # *strongly* suggested for sharing
    version='0.1',

    # The license can be anything you like
    description='HOSE code generator producing stereo enhanced HOSE code',
    # long_description=open('README.txt').read(),
    license_files = ('LICENSE.txt',),
)