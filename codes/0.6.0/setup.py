from setuptools import setup, Extension
from Cython.Build import cythonize
import os

# Helper to collect .pyx files and construct Extension objects
def find_pyx_modules(package_dir):
    extensions = []
    for root, dirs, files in os.walk(package_dir):
        for file in files:
            if file.endswith(".pyx"):
                module_path = os.path.join(root, file)
                module_name = module_path.replace(os.path.sep, ".").replace(".pyx", "")
                extensions.append(Extension(module_name, [module_path]))
    return extensions

extensions = find_pyx_modules("fastforman")

setup(
    name="fastforman",
    version="0.6.0",
    ext_modules=cythonize(extensions, language_level="3"),
    packages=["fastforman"],
    install_requires=[
        "cython",
        "numpy",
        "scikit-learn",
        "scipy",
        "networkx"
    ],
    zip_safe=False,
)
