# File: setup.py
# Date: 1-Oct-2019
#
# Updates:
#
#
import re

from setuptools import find_packages
from setuptools import setup

packages = []
thisPackage = "rcsb.ccmodels"

with open("rcsb/ccmodels/search/__init__.py", "r") as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE).group(1)

if not version:
    raise RuntimeError("Cannot find version information")

setup(
    name=thisPackage,
    version=version,
    description="RCSB Python Chemical Component Model Utilities",
    long_description="See:  README.md",
    author="John Westbrook",
    author_email="john.westbrook@rcsb.org",
    url="https://github.com/rcsb/py-rcsb_ccmodels",
    #
    license="Apache 2.0",
    classifiers=(
        "Development Status :: 3 - Alpha",
        # 'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
    ),
    entry_points={
        "console_scripts": [
            "cc_models_cli=rcsb.ccmodels.search.ChemCompModelExec:main",
        ]
    },
    #  The following is somewhat flakey --
    dependency_links=["https://pypi.anaconda.org/OpenEye/simple#egg=OpenEye-toolkits-2020.2.0"],
    install_requires=[
        "mmcif >= 0.61",
        "rcsb.utils.io >= 1.01",
        "rcsb.utils.multiproc >= 0.18",
        "rcsb.utils.chem >= 0.70",
        "rcsb.utils.chemref >= 0.66",
        "OpenEye-toolkits>=2020.2.0",
        "wrapt_timeout_decorator >= 1.3.1",
    ],
    packages=find_packages(exclude=["rcsb.mock-data", "rcsb.utils.tests-chem", "rcsb.utils.tests-*", "tests.*"]),
    package_data={
        # If any package contains *.md or *.rst ...  files, include them:
        "": ["*.md", "*.rst", "*.txt", "*.cfg"]
    },
    #
    test_suite="rcsb.ccmodels.tests",
    tests_require=["tox"],
    #
    # Not configured ...
    extras_require={"dev": ["check-manifest"], "test": ["coverage"]},
    # Added for
    command_options={"build_sphinx": {"project": ("setup.py", thisPackage), "version": ("setup.py", version), "release": ("setup.py", version)}},
    # This setting for namespace package support -
    zip_safe=False,
)
