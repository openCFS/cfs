Developer's Manual                                                {#devmanual}
==================

\section introduction Introduction

### History of CFS++

The Finite-Element  (FE) software  has been named  CFS++ according  to Coupled
Field Simulation  and the ++ indicates  that as programming languange  C++ has
been chosen.   Originally, the first  struture of  CFS++ has been  designed by
M. Kaltenbacher and  S. Reitzinger during a  stay of M. Kaltenbacher  in Linz,
Austria in  the year  2000 at  the SFB013  (Numerical and  Symbolic Scientific
Computing). The  idea was to set  up a new simulation  software for scientific
computing capable  for multi-field applications.   Since then on the  group of
M.Kaltenbacher  at  the  Department   of  Sensor  Technology,  University  of
Erlangen-Nuremberg in Germany  was growing and the ideas were  realized in the
FE  software CFS++.   In 2008,  M. Kaltenbacher and  many co-workers  moved to
University Klagenfurt, Austria. Furthermore,  within the Cluster of Excellence
for   Engineering  of   Advanced  Materials   (University  Erlangen-Nuremberg,
Germany), the research group of M.Stingl contribute to CFS++ concerning shape
and topology optimization.

Among many scientists,  who have contributed to CFS++, we  want to acknowledge
the following:

-   M. Escobar: aeroacoustics

-   J. Grabinger: non-matching grids, aeroacoustics

-   A. Hauck: direct and iterative coupled PDEs, higher order FEM,
    revisions of class structure (function space, BDB-operator, ..)

-   M. Hofer: forms class

-   A. Hüppe: spectral FEM, revisions of class structure (function
    space, BDB-operator, ..)

-   M. Kaltenbacher: nonlinear PDEs, hysteresis, PML

-   M. Mohr: algebraic system and XML-structure

-   T. Lahmer: inverse schemes

-   G. Link: fluid flow and coupling to mechanics

-   S. Triebenbacher: non-matching grids, CFSDEPS, CMake environment,
    ctest, ..

-   F. Wein: optimization

-   E. Zhelezina: first implementations. hFEM

-   S. Zörner: fluid-structure-acoustic, XFEM

### Structure of CFS++ and Programming Environment

[Sec:CFSstructure]

From the beginning  on, the main focus  was and still is  on numerical schemes
for solving the arising system of PDEs within multi-field problems. Therewith,
we  have no  capabilities concerning  geometric modeling,  mesh generation  or
post-processing. Instead,  we provide  a variety of  interfaces to  obtain the
appropriate input data  from a number of supported pre-processors,  as well as
to write different file formats to feed a variety of post-processing tools.

#### Data I/O and Class Structure 

To obtain the necessary input information  for our numerical simulation tool a
number of input files need to be provided:

This data  provides the coordinates and  mesh topology of all  3D (volume), 2D
(surface/interface) and  1D (interface)  elements. Due to  the fact,  that all
elements are tagged with unique ids,  we can differ between different physical
regions where different material properties may be applied, interfaces for the
specification  of  couplings  as  well as  surfaces  within  the  computation.
Furthermore, it provides  special elements and nodes on  which boundary or/and
load conditions  can be specified or  which are used as  monitoring points for
saving special physical data (e.g.,  acoustic pressure over time). In addition
to  the  mesh information,  field  information  may be  read  as  input for  a
simulation (e.g.  acoustic  source terms). By default, this data  is stored in
the  hierarchical file  format HDF5  \cite hdf5.  However, we  also support  a
number of different pre-processors and data sources .

To have a standardized  and easy extendable file format, we  use for this task
the XML-format. By  specifying schemas \cite xmlschema for  our expected input
formats, it is possible  to validate input files written by  humans as well as
automatically generated ones.

##### Listing: XML for piezoelectric-acoustic coupling \anchor xml-piezo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.xml}
  <couplingList>
    <direct>
      <piezoDirect>
        <regionList>
          <region name="solid"/>
        </regionList>
      </piezoDirect>
      <acouMechDirect>
        <surfRegionList>
          <surfRegion name="interface"/>
        </surfRegionList>               
      </acouMechDirect>
    </direct>                        
  </couplingList>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to that when writing the input  files in a XML editor, tags can be
suggested by the  editor and the documentation of tags  inside the schema also
assists  users, so  that  they  seldomly need  to  consult  manuals. Also  the
availability of  a vast range  of tools for  transforming XML data  (e.g. XSLT
\cite xslt) simplifies conducting parameter studies by a great margin. We list
the main  commands for the  coupling section within  a [piezoelectric-acoustic
example](#xml-piezo).

The XML-file  format is also  used for  the material description  allowing for
automatic  processing  entries with  3rd  party  tools  as well  as  integrity
checking.  As an  example, we  show the  specification of  the electromagnetic
properties of aluminum (see Fig. [fig:xml-mat](#xml-mat)).

##### Listing: XML material file specification   of the electromagnetic  properties of aluminum \anchor xml-mat
\htmlinclude pygments_test.html
\latexonly
\input{pygments_test.tex}
\endlatexonly

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.xml}
<material name="alu">
  <magnetic>
    <electricConductivity> 1.0e7 </electricConductivity>
    <magneticPermeability>
      <linear>
        <isotropic> 1.25664e-6 </isotropic>
      </linear>
    </magneticPermeability>
  </magnetic>
</material>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

\section franz Franz

Concerning the result data, our standard file format is again HDF5. This
enables us to do chained simulations, since the output of a simulation
may be used as an input for a subsequent simulation. We may also use
third-party tools like Matlab in this chain, due to the fact that the
HDF5 format is supported by many tools. Furthermore, as displayed in
Fig. [fig:cfs~I~O], we provide many tool specific input- and output-file
formats to allow a use of different post-processing tools.


![Caption text](images/cfs-logo.png "Image title")

\tableofcontents