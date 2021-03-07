from setuptools import setup, find_packages

setup(
    name="cygnus",
    version="0.1",
    author="Bernie Lee",
    author_email="bernadette.lee.cy@gmail.com",
    description=("cygnus: boop boop"),
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=[
        "numpy>=1.17.0",
        "matplotlib",
        "pandas",
        "seaborn",
    ]

)
