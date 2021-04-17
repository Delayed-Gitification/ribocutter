from setuptools import setup

# This call to setup() does all the work
setup(
    name="ribocutter",
    version="0.0.1",
    description="Design oligos to produce sgRNAs for abundant sequences in a fastq file",
    #long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/Delayed-Gitification/ribocutter.git",
    author="Oscar Wilkins",
    author_email="oscar.wilkins@crick.ac.uk",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=["ribocutter"],
    include_package_data=True,
    install_requires=["pandas", "dnaio"],
    entry_points={
        "console_scripts": [
            "ribocutter=ribocutter.__main__:main",
        ]
    },
)
