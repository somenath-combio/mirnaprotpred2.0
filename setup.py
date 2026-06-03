from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="mirnaprotpred2",
    version="2.0.0",
    author="Somenath Dutta, Sudipta Sardar",
    author_email="somenath@pusan.ac.kr, sudipta@pusan.ac.kr",
    description="Two-stage viral miRNA-CTS prediction: thermodynamic scan + ML scoring",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/somenath-combio/mirnaprotpred2",
    packages=find_packages(),
    package_data={
        "mirnaprotpred2.SeqFinder": ["data/*.pkl", "data/*.csv"],
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "SeqFinder2=mirnaprotpred2.SeqFinder.seqfinder:cli",
            "validator2=mirnaprotpred2.validator.validator:cli",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.21",
        "scikit-learn>=1.0",
        "xgboost>=1.6",
        "scipy>=1.7",
        "biopython>=1.79",
        "pyfiglet>=1.0.4",
        "openpyxl>=3.1.5",
        "lightgbm>=3.3.0",
        "catboost>=1.0.0",
    ],
)
