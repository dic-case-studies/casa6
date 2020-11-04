import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="casatestutils", # Replace with your own username
    version="0.1.10",
    author="A. Wells",
    author_email="awells@nrao.edu",
    description="Tools for use with casatest and testing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/browse",
    download_url="https://casa.nrao.edu/download/",
    packages=setuptools.find_packages(),
    package_data={"casatestutils": ["component_to_test_map.json"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[ 'scipy', 'numpy' ]

)
