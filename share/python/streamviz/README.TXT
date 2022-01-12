Summary
-------
streamviz visualizes transient/harmonic/optimization data.
It is a minimalistic webserver, receiving data from openCFS via
 <output>
   <streaming host="localhost" port="5000" protocol="http" sendMesh="true" /> 
 </output>
a public server is monster.mi.uni-erlangen.de

Installation
------------
The most important python packages are:
- flask
- matplotlib
- svgwrite
- objgraph

The only non python depency is graphziv

The python paraview package for paraview catalyst is not necessary, as it does not run that stable.

Running
-------
Start with python streamviz.py. See python streamviz.py --help for options.

External Access
----------------
To access streamviz from extern, you need to start streamwith with the real host ip.
Also open your firewall. Yast for openSUSE, sudo ufw allow 5000 for Ubuntu.

Streaming
---------
With the streaming output element in a background thread all data written to .cfs and .info.xml are
sent to the streaming target. This shall not have a performance penalty. If the streaming target is
not reachable, cfs still runs.

The file protocol is more for development reasons.

Authors
-------
The first and full implementation was by Simon Michalke. Extensions by Fabian Wein

