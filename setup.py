import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="crystal-genome-util",
    version="0.2.0",
    description=(
        "Common functions used by Crystal Genome test drivers"
    ),
    author=["ilia Nikiforov", "Ellad B. Tadmor", "Daniel S. Karls"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    license="LGPL",
    install_requires=["numpy >= 1.13.1", "ase >= 3.19.0b1", "spglib >= 2.1.0", "kim-edn >= 1.3.1", "kim-property >= 2.6.0"],
    python_requires=">=3.8",
    include_package_data=True
)
