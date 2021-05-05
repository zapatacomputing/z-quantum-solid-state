import setuptools
import os

readme_path = os.path.join("README.md")
with open(readme_path, "r") as f:
    long_description = f.read()

setuptools.setup(
    name="z-quantum-solid-state",
    version="0.1.0",
    author="Zapata Computing, Inc.",
    author_email="info@zapatacomputing.com",
    description="Library for solid state calculations for Orquestra.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zapatacomputing/z-quantum-solid-state ",
    packages=setuptools.find_namespace_packages(include=["zquantum.*"]),
    package_dir={"": "src/python"},
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ),
    install_requires=[
        "z-quantum-core",
    ],
)
