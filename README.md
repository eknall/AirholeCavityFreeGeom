Authors: Mike Burek, Bart Machielse, Michelle Chalupnik, Erik Knall 
ReadMe Updated: 07/09/19

The master branch of the AirholeCavityFreeGeom repository contains simulation scripts for optimizing an overcoupled airhole cavity in diamond. The optimizeNanobeamCavity.m employs a figure of merit which is rightsidedness*guidedness*Qscat/modeVolume. That is, rightsidedness = Qx2/Qx1 and guidedness = Qscat/Qtot.

In order to run the optimizer on the harvard odyssey cluster, it is necessary to use the wrapper optimizeWrapper which simply calls optimizeNanobeamCavity. This is because the bash script which is used to launch jobs on the cluster can only call matlab executable scripts; to the best of my knowledge, it can't launch files that contain functions.

There are forks of this repository which use different figures of merit. This can be for optimizing a different type of cavity (eg. a different geometry restriction) or for optimizing the same type of cavity but exploring the impact of changing one of the optimization params.
