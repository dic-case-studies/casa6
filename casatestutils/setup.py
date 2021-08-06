import setuptools
import argparse
import sys
from subprocess import Popen, PIPE


with open("README.md", "r") as fh:
    long_description = fh.read()

parser=argparse.ArgumentParser()
parser.add_argument('--version', help='version')
parser.add_argument('bdist_wheel', help='bdist_wheel')
args=parser.parse_args()

# Remove the "non-standard" arguments from sys.argv so as not to confuse dist_tools
for arg in sys.argv:
    if (arg.startswith("--version")):
        sys.argv.remove(arg)

pyversion = float(sys.version_info[0]) + float(sys.version_info[1]) / 10.0

if pyversion < 3:
    str_encode = str
    str_decode = str
    def pipe_decode(output):
        return output
else:
    def str_encode(s):
        return bytes(s,sys.getdefaultencoding())
    def str_decode(bs):
        return bs.decode(sys.getdefaultencoding(),"strict")
    def pipe_decode(output):
        if isinstance(output,bytes) or isinstance(output,bytearray):
            return str_decode(output)
        elif isinstance(output,tuple):
            return (str_decode(output[0]),str_decode(output[1]))
        else:
            return ("","")


def compute_version( ):
    if (args.version != None ):
        print (args.version.split("."))
        (major, minor, patch, feature) = args.version.split(".")
        return(int(major), int(minor), int(patch), int(feature),"","","")
    else:
        proc = Popen( [ "./version" ], stdout=PIPE, stderr=PIPE )
        out,err = pipe_decode(proc.communicate( ))
        print(out)
        devbranchtag = out.split(" ")[0].strip()
        print(devbranchtag)
        releasetag = out.split(" ")[1].strip()
        dirty=""
        if (len(out.split(" ")) == 3):
            print("Latest commit doesn't have a tag. Adding -dirty flag to version string.")
            dirty="+" + out.split(" ")[2].strip() # "+" denotes local version identifier as described in PEP440
        print(releasetag)
        devbranchversion = ""
        devbranchrevision = ""
        if (devbranchtag != releasetag):
            devbranchrevision = devbranchtag.split("-")[-1]
            if (devbranchtag.startswith("CAS-")):
                devbranchversion=devbranchtag.split("-")[1]
            else:
                devbranchversion=100
            devbranchrevision = devbranchtag.split("-")[-1]
        else:
            isDevBranch = False
        (major, minor, patch, feature) = releasetag.split(".")
        #print(major, minor, patch, feature, devbranchversion, devbranchrevision, dirty)
        return(int(major), int(minor), int(patch), int(feature), devbranchversion, devbranchrevision, dirty)

(casatestutils_major,casatestutils_minor,casatestutils_patch,casatestutils_feature, devbranchversion, devbranchrevision, dirty) = compute_version( )
print(casatestutils_major, casatestutils_minor, casatestutils_patch, casatestutils_feature, devbranchversion, devbranchrevision, dirty)
casatestutils_version = '%d.%d.%d.%d%s' % (casatestutils_major,casatestutils_minor,casatestutils_patch,casatestutils_feature,dirty)
if devbranchversion !="":
    casatestutils_version = '%d.%d.%d.%da%s.dev%s%s' % (casatestutils_major,casatestutils_minor,casatestutils_patch,casatestutils_feature,devbranchversion,devbranchrevision,dirty)
print(casatestutils_version)

# Copy runtest to casatestutils during setup.py so we can call 'casatestutils.runtest'
import shutil, os
runtest_path = os.path.join(os.getcwd(),"casatestutils","runtest.py")
shutil.copyfile("runtest.py", runtest_path)

setuptools.setup(
    name="casatestutils", # Replace with your own username
    version=casatestutils_version,
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
    install_requires=[ ] #'scipy','numpy', 'six' 
)

# Delete Copy of File to avoid accidental commit to this location
if os.path.exists(runtest_path):
    os.remove(runtest_path)

