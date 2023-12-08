# EpithelialPolarisation_Org3D

## Model of crypt formation on 3D intestinal organoids with 4 different cell types 

Before looking at this, you may wish to look at some of the https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials for Chaste.

### Getting the code and installing dependencies 

Before running this model you will need to install Chaste following these instructionds: https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/InstallGuide
The easiest way to do this is using an Ubuntu machine (or an Ubuntu virtual machine) as discussed on https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/UbuntuPackage
Note that Chaste is only fully supported on Linux/Unix systems, so users of Windows or Mac OS X may need to follow the virtual machine route.
For manual installation of each dependency, on any version of Linux, see DeveloperInstallGuide. Alternatively, Chaste can also be run using Docker, and the dockerfile can be accessed at https://github.com/Chaste/chaste-docker, where you'll find instructions for installation and usage.

Then you will need to add this project to the projects folder in the Chaste directory:
```sh
git clone https://github.com/slmontes/EpithelialPolarisation_Org3D.git
```

### Documentation 
There are two folders - `src` and `test`.
 1. The `src` folder contains the classes necesary to run the simulation. These define the aditional forces and boundary conditions not in the core chaste code.
 1. The `test` folder contains:
  - TestIntestinalOrganoid_3D_4CT.hpp this file can be run to generate the in-silico results
  
### Running tests
You can then run tests and simulations from the terminal with:
```sh
cd <Chaste_build path>
make -j4 TestIntestinalOrganoid_3D_4CT
ctest -j4 -V -R TestIntestinalOrganoid_3D_4CT
```
----
> Note: the current code was developed with release version 2021, but will work on release version 2019 and 2018.1. It will not work with with release version 3.4 or under.

For further information on using Chaste, see https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides.
You may also wish to look at some of the https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials.