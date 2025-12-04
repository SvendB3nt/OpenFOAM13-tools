# OpenFOAM13-tools
Various convenience tools for OpenFOAM.

BCs for turbulence model of k-epsilon and k-omega types, in pipes.
The following BCs are based on an extension of the inletOutlet BC, that allows for automatic calulation 
of the hydraulic diameter of the patch for simple shapes.

- turbulentPipeIntensityKineticEnergyInlet
- turbulentPipeMixingLengthDissipationRateInlet
- turbulentPipeMixingLengthFrequencyInlet

To download and compile, assuming that openFOAM is installed from the debian package:
```bash
source . /opt/openfoam13/etc/bashrc \
&& cd $WM_PROJECT_USER_DIR \
&& git clone https://github.com/SvendB3nt/OpenFOAM13-tools \
&& cd src \
&& ./Allwmake
&& cd ..
```
