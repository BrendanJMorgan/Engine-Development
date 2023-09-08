# Engine-Development
LOx-Ethanol Gas Generator Cycle Rocket Engine Analysis

This folder centers around designing a gas generator cycle rocket engine. By default the propellants are liquid oxygen and ethanol, though it should be fairly easy to change the code to other common propellants. The solidworks CAD reflects design parameters outputted by the MATLAB code. 

To use, open the main.m file and run. Engine target parameters can be adjusted right from this file. For example, target thrust, chamber pressure, OF ratio, and some basic dimensions for both the main combustion chamber (CC) and gas generator (GG) can be changed here.

The MATLAB code requires CoolProp integration to run. CoolProp is a python library with a MATLAB wrapper. The easiest installation method I know of can be followed from this video: https://www.youtube.com/watch?v=XvR10Fjph7U
The code also reuires NASA's CEAM to run, though this is already included in the folder so no further action should be reuired. NASA CEAM is approved for both domestic and international release.

Further documentation for running the MATLAB code can be found in Documentation.pptx, found in the MATLAB folder.
