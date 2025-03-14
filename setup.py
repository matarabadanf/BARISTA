from setuptools import setup, find_packages

setup(
    name="barista",
    version="0.1.0",
    packages=find_packages(),
    #packages=find_packages(include=['barista', 'barista.*']),
    install_requires=[
        'numpy',
        'pandas',  
        'matplotlib',
        'plotly',
        'ase'
    ],
)
