ANSYS
=====

Installation
------------

### ANSYS
Choose a directory installation directory `/path/to/installation/`.
The installer installs in the sub-direcory with the correct version number (e.g. v170 for v17.0).
You can not move the installation due to absolute paths in a number of files!

### mkmesh.mlib CFS++ extension
To activate the `mkmesh.mlib` extensions the modified start script `/path/to/mkmesh/start.ans` must be run when ANSYS starts.
For this purpose rename the example ANSYS start script `/path/to/installation/v<XXX>/ansys/apdl/start.ans` and replace it with a link to the one from `mkmesh`.
```
cd /path/to/installation/v170/ansys/apdl
mv start.ans start.ans.original
ln -s -T /path/to/mkmesh/start.ans start.ans
```
The `mkmesh.mlib` extensions should be now automatically available when ANSYS starts.
If this is not the case one could also execute both commands contained in the start script manually or append them to the top of the *.in file after ANSYS has already been started.
For more information see also `/path/to/mkmesh/start.ans`.

CFS++ mkmesh plugin
-------------------
In order to use the `mkmesh` awk script directly in ANSYS (*.in file) it has to be set as shell variable first.
One could set a link or directly make a copy of the script in one of the directories in $PATH.
'''
echo $PATH
cd /path/to/one/of/the/listed/directories/
ln -s /path/to/mkmesh mkmesh
'''

Command reference
-----------------
For more information see also http://cfs-doc.mdmt.tuwien.ac.at/mediawiki/index.php/Creating_Meshes
* init: Initializes the mesh-extension. 
* setelems,type,[order]: Select element type for meshing. Valid types are '2d-line', '3d-line', 'triangle', 'quadr', 'tetra', 'brick', 'pyramid' and 'wedge'. Valid orders are 'linear' (=default) and 'quad'.
* welems,name: Write selected elements with name attribute name.
* wregelem,name: Deprecated.
* wnodes: Write selected nodes.
* wnodbc,name: Write selected nodal boundary conditions with name attribute name.
* wsavnod,name: Write selected nodes for saving nodal results with attribute name.
* wsavelem,name: Write selected elements for saving element results with attribute name.
* mkmesh: Write *.mesh file.
* mkhdf5: Deprecated. Needs cplreader to be installed.






