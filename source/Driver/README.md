Driver   ([back to main page](/source/CFS_Library_Documentation.md))
===

In CFS++ we differ between the different analysis types, which are implemented as *drivers*:

* StaticDriver: performs a static analysis
* TransientDriver: performs a transient analysis
* HarmonictDriver: performs a harmonic analysis
* EigenFrequencyDriver: performs an eigen-frequency analysis

>

Thereby, for each analysis (driver) it is defined how a single step (e.g., transient step) is performed. In the class **StdSolveStep** you can find the implementation of a standard step for 
each of the drivers. The interaction between the different classes for a static simulation is shown below.

![](/share/doc/developer/pages/pics/StaticDriverInteraction.png)

