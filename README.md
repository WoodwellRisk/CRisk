## Extreme Coastal Risk System (TO BE RENAMED)

### Install Repository

The easiest way to set everything up for this package is to use the `install.sh` script.
This script will install pip, install dependencies for the `tc_analysis` scripts and install
climada to your active python environment.

Steps:

* Clone this repository somewhere onto your system using `git clone`.
* Create a new Python environment using `conda` or `mamba` and activate it.
* Change directory into the repository directory (where you cloned it to).
* Run `./install.sh`. You may need to change permissions on the file using `chmod` first.

This scripts will also create `input` and `output` directories in the repository directory. You can keep data wherever you want
by using the arguments in the scripts, however they are setup to look for the above directories by
default.

To install the package for development (if you want to modify the code) you will need to follow the climada installation instructions on the
climada website. Then you can install this package using `pip install -e .`
