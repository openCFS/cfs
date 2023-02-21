XML Reference Manual{#xmlref}
====================

\section xmlbasics XML Basics

What is XML?
------------

XML is the *Extensible Markup Language*. It is designed to improve the
functionality of the Web by providing more flexible and adaptable
information identification.

It is called extensible because it is not a fixed format like HTML (a
single, predefined markup language). Instead, XML is actually a
*metalanguage* - a language for describing other languages - which lets
you design your own customized markup languages for limitless different
types of documents. XML can do this because it's written in *SGML*, the
international standard metalanguage for text markup systems (ISO 8879).

XML: Data Elements
------------------

In XML a single data item is enclosed by two tags:

\<node~n~ame\>ground\</node~n~ame\>

an *opening tag* and a *closing tag*. The whole assembly is called an
element. In XML, an element can have

ll a value &\
 attributes &\
 children &

10cm

\<nodelist\> \<node1\>....\</node1\> \<node2\>....\</node2\>
\</nodelist\>

\

or a combination of the above.

XML: Empty Elements
-------------------

The XML specification allows for empty elements, i.e. elements that
contain neither a value nor any childred. Empty, however, does not
necessarily mean that the element does not carry any information, but
its name, since an empty element may still possess attributes. Indeed
many of the elements used in the CFS++ parameter file are empty with
information encoded as attribute values. Note also that, if an element
is empty,

\<problem name="demo3d"\>\</problem\>

then the opening and the closing tag can, for brevity, be combined into
one:

\<problem name="demo3d"/\>

\tableofcontents

\section strucconfigfile Structure of Config-File

[sec:general] The basic structure of the config-file for CFS++ in our
XML syntax is as follows:

\<?xml version="1.0"?\> \<cfsSimulation xmlns="http://www.cfs++.org"\>

... Parameters ...

\</cfsSimulation\>

The first line gives version information for the XML parser, while the
second line is the starting tag for the main enclosing element called .
The last line then is the closing tag for this main element. All
parameters must appear between the opening and the closing tag.

Only a few tags are allowed in this main body. They will be documented
in the remaining sections. In the following, you can find a simple
example of a harmonic piezoelectric analysis.

\<?xml version="1.0"?\>

\<cfsSimulation xmlns="http://www.cfs++.org"\>

\<fileFormats\> \<output\> \<text/\> \<info/\> \<hdf5/\> \</output\>
\<materialData file="mat.xml" format="xml"/\> \</fileFormats\>

\<domain geometryType="axi"\> \<regionList\> \<region name="piezo~a~x"
material="PZT-4"/\> \</regionList\> \<surfRegionList\> \<surfRegion
name="surf"/\> \</surfRegionList\> \<nodeList\> \<nodes
name="ep-lower"/\> \<nodes name="ep-top"/\> \<nodes name="axisBoun"/\>
\<nodes name="save-u"/\> \</nodeList\> \</domain\>

\<sequenceStep\> \<analysis\> \<harmonic\> \<numFreq\> 5 \</numFreq\>
\<startFreq\> 2.0e6 \</startFreq\> \<stopFreq\> 6.0e6 \</stopFreq\>
\</harmonic\> \</analysis\>

\<pdeList\>

\<mechanic subType="axi"\>

\<regionList\> \<region name="piezo~a~x" complexMaterial="yes"/\>
\</regionList\>

\<bcsAndLoads\> \<fix name="axisBoun"\> \<comp dof="r"/\> \</fix\>
\</bcsAndLoads\>

\<storeResults\> \<nodeResult type="mechDisplacement"\> \<allRegions/\>
\</nodeResult\> \<elemResult type="mechStrain"\> \<allRegions/\>
\</elemResult\> \</storeResults\> \</mechanic\>

\<electrostatic\> \<regionList\> \<region name="piezo~a~x"
complexMaterial="yes"/\> \</regionList\> \<bcsAndLoads\> \<ground
name="ep-lower" /\> \<potential name="ep-top" value="10" phase="30"/\>
\</bcsAndLoads\>

\<storeResults\> \<elemResult type="elecFieldIntensity"\>
\<allRegions/\> \</elemResult\> \</storeResults\> \</electrostatic\>

\</pdeList\>

\<couplingList\> \<direct\> \<piezoDirect\>

\<regionList\> \<region name="piezo~a~x" complexMaterial="yes"/\>
\</regionList\>

\<storeResults\>

\<elemResult type="mechStress"\> \<allRegions/\> \</elemResult\>

\<surfRegionResult type="elecCharge"\> \<surfRegionList\> \<surfRegion
name="surf" neighborRegion="piezo~a~x" /\> \</surfRegionList\>
\</surfRegionResult\>

\</storeResults\> \</piezoDirect\> \</direct\> \</couplingList\>

\</sequenceStep\>

\</cfsSimulation\>