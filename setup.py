import os

import setuptools

dev_requires = [
    "pytest>=3.7.1",
    "pytest-cov>=2.5.1",
    "tox>=3.2.1",
    "flake8>=3.7.9",
    "black>=19.3b0",
    "pre_commit>=2.10.1",
]

extras_require = {
    "dev": dev_requires,
}


def _this_path():
    return os.path.abspath(os.path.dirname(__file__))


def _read_readme():
    with open(os.path.join(_this_path(), "README.md")) as f:
        return f.read()


setuptools.setup(
    name="z-quantum-solid-state",
    version="0.1.0",
    author="Zapata Computing, Inc.",
    author_email="info@zapatacomputing.com",
    description="Library for solid state calculations for Orquestra.",
    license="Apache-2.0",
    long_description=_read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/zapatacomputing/z-quantum-core",
    packages=setuptools.find_namespace_packages(
        include=["zquantum.*"], where="src/python"
    ),
    package_dir={"": "src/python"},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Operating System :: OS Independent",
    ],
    install_requires=["z-quantum-core"],
    extras_require=extras_require,
)
