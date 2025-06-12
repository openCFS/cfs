# Testing grid conversion with CGNS Reader

This test checks if the grid conversion using the CGNS Reader works as expected. Therefore, dummy simulation xml file is used and the grid conversion only is executed using `cfs -gG ShellElementConversion_CGNS`.

1. The grid is read by the from the cgns input
2. The grid is written to a hdf5 '.cfs' file
3. Additionally, the grid is written to the '.info.xml' file
4. All attributes in the 'regions' tag in the '.info.xml' are compared to the 'ref.info.xml', thus checking if all region names, ids, types, ... are identical.
5. All attributes in the 'grid' tag in the '.info.xml' are compared to the 'ref.info.xml', thus checking if the location of all nodes and connectivity and assignment of all elements inside the corresponding region are identical