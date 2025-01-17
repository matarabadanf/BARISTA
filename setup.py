from setuptools import setup, find_packages

setup(
    packages=find_packages(include=['barista', 'barista.*']),
    install_requires=[
        'numpy',
        'pandas',  
        'matplotlib'
    ],
)
