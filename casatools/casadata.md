
## CASA Runtime Data

CASA uses binary tables which are updated regularly. These tables contain information about [the earth's orientation](https://www.iers.org/IERS/EN/DataProducts/data.html), the path of ephemera and other support data. 

A subset of this runtime data is updated regularly because things like the earth's orientation changes as time progresses. After a couple of months, this data is sufficiently out of date to begin affecting the images which CASA produces.

## Initial Setup

There are two paths to getting CASA's runtime data. One is by installing a Python package that is only available from NRAO's PyPI repository. It is not possible to put the data on standard Python package servers because of the data size and the frequency of the updates.

### Use Existing Data

If you already have a CASA data repository that contains CASA's runtime data, you can point casatools to this data by setting ```rundata``` in CASA's config.py file, ```~/.casa/config.py```. For example:
```
rundata = '/home/casa/data/trunk'
```

### user data installation

You can install your own copy of CASA's runtime data using casatools with:
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
