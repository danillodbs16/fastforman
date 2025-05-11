from setuptools import setup, Extension

# Define your Cython extension modules here
extensions = [
    Extension("fastforman", ["fastforman.pyx"]),
    # Add more extensions if you have multiple Cython modules
]

setup(
    name='fastforman',
    author='Danillo Barros de Souza',
    author_email='danillo.dbs16@gmail.com',
    description='fastforman - An efficient algorithm for computing Forman-Ricci curvature from point cloud data.',
    url='https://github.com/danillodbs16/fastforman',
    version='0.5.0',
    ext_modules=extensions,
    install_requires=[
    'numpy',
    'cython',
    'gudhi',
    ],
)
