from setuptools import setup, Extension

# Define your Cython extension modules here
extensions = [
    Extension("fastforman", ["fastforman.pyx"]),
    # Add more extensions if you have multiple Cython modules
]

setup(
    name='fastforman',
    version='0.1',
    ext_modules=extensions,
)
