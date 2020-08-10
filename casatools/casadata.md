
## CASA Runtime Data

CASA uses binary tables which are updated regularly. These tables contain information about [the earth's orientation](https://www.iers.org/IERS/EN/DataProducts/data.html), the path of ephemera and other support data. 

A subset of this runtime data is updated regularly because things like the earth's orientation changes as time progresses. After a couple of months, this data is sufficiently out of date to begin affecting the images which CASA produces.

## Initial Setup

There are two paths to getting CASA's runtime data. One is by installing a Python package that is only available from NRAO's PyPI repository. Because the packaged data is 336M and is updated every week, it is not possible to put the data on standard Python package servers because of the data size and the frequency of the updates.

### Use Existing Data

If you already have a CASA data repository that contains CASA's runtime data, you can point casatools to this data by setting ```rundata``` in CASA's config directory. For example:
```
rundata = '/home/casa/data/trunk'
```

### user data installation

You can install the data in CASA's runtime data using casatools with:
```
-bash-4.2$ python3 -m casatools --update-user-data
```
This will install the runtime data into ```~/.casa/data``` (or wherever ```rcdir``` specifies) if it does not already exist. If the directory already exists, it will update the runtime data to the most recent version.

### casadata package

The casadata package can also be used to supply the runtime data for casatools. It can be installed from CASA's internal PyPI server.

