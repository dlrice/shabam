
from setuptools import setup

setup(
    name = "shabam",
    version = "1.0.0",
    author = "Daniel Rice",
    author_email = "daniel.rice@sanger.ac.uk",
    description = ("Easy sequence alignment plots"),
    license = "GPL",
    packages=["shabam"],
    install_requires=['pysam >= 0.9.0',
                      'cairocffi >= 0.7.0',
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
