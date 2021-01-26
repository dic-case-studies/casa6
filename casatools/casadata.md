
## CASA Runtime Data

CASA uses binary tables which are updated regularly. These tables contain information about [the earth's orientation](https://www.iers.org/IERS/EN/DataProducts/data.html), the path of ephemera and other support data. 

A subset of this runtime data is updated regularly because things like the earth's orientation changes as time progresses. After a couple of months, this data is sufficiently out of date to begin affecting the images which CASA produces.

If you came here after getting an error like:
```
>>> import casatools
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/export/hypnos/tmp/casa6-home-data/casa6/casatools/build/lib.linux-x86_64-3.6/casatools/__init__.py", line 153, in <module>
    raise ImportError('measures data is not available, visit http://go.nrao.edu/casadata-info for more information')
ImportError: measures data is not available, visit http://go.nrao.edu/casadata-info for more information
>>> 
```
the easiest solution is to run:
```
-bash-4.2$ python3 -m casatools --update-user-data
```
This will install CASA's runtime data in ~/.casa/data, and the next time you import casatools you should not encounter this error.

## Installing Runtime Data

There are three options for installing CASA's runtime data:

  * use an copy of data maintained by your site or institution
  * fetch your own copy using the casatools package
  * install a copy using [pip](https://pypi.org/project/pip/)

### Use Existing Data

If you already have a CASA data repository that contains CASA's runtime data, you can point casatools to this data by setting ```rundata``` in CASA's config.py file, ```~/.casa/config.py```. For example:
```
rundata = '/home/casa/data/trunk'
```
All [ALMA Regional Centers (ARC)](https://www.almaobservatory.org/en/about-alma-at-first-glance/global-collaboration/) and all NRAO sites have a copy of the runtime data used by casatools. For users at these locations, using the local copy of the runtime data is the best option because it is kept up to date by each location.

### user data installation

As described above, you can install your own copy of CASA's runtime data using casatools with:
```
-bash-4.2$ python3 -m casatools --update-user-data
```
This will install the runtime data into ```~/.casa/data``` (or a user specified alternative ```rcdir```). If a copy of the runtime data already exists there, it will be updated. If it does not already exist, all of the runtime data will be downloaded (approximately 830MB).

### casadata package

The casadata package can also be used to supply the runtime data for casatools. It can be installed from CASA's internal PyPI server. This can be done like:
```
-bash-4.2$ pip3 install --index-url=https://go.nrao.edu/pypi casadata
```
This will make the ```casadata``` package available within the current python environment.
