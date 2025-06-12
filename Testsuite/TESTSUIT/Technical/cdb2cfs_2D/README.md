# Testing the grid conversion with cfs -g [--printGrid] option without material input

This test checks if the grid conversion using the option `cfs -g` without material input works as expected. Therefore, dummy simulation xml file and a geometry `*.cdb` file are used and the grid conversion is executed using `cfs -g cdb2cfs_2D`.

1. The grid is read by the `cfs -g` input
2. The grid is written to a hdf5 `*.cfs` file
3. Addionally, the grid is written to the `*.info.xml` file
4. All attributes in the `<regions>` tag in the `*.info.xml` are compared to the `ref.info.xml`, thus checking if all region names, ids, types, ... are identical.
5. All attributes in the `<grid>` tag in the `*.info.xml` are compared to the `ref.info.xml`, thus checking if the location of all nodes and connectivity and assignment of all elements inside the corresponding region are identical
