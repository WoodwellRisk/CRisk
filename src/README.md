This directory is for routines that are shared by the whole repository.
Each directory that needs to be 'found' by the package installation should contain a `__init__.py` file.
Functions that are within files here can be imported into scripts using something like:

```
from CRISK.track_analysis import windspeed_analysis
```
