from setuptools import setup

setup(
    name='HOSE_code_generator',
    url='https://github.com/Ratsemaat/HOSE_code_generator',

    # Needed to actually package something
    packages=['hose_generator'],

    # Needed for dependencies
    install_requires=['rdkit'],

    # *strongly* suggested for sharing
    version='0.1',

    # The license can be anything you like
    description='An example of a python package from pre-existing code',
    # long_description=open('README.txt').read(),
)