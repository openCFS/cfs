<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <!-- This stylesheet extracts the last design set from a density.xml file
    (optimization!) as single line such that it can be read in octave/matlab
    for delta-design in parameter studies 
  -->
  <xsl:output method="text" encoding="UTF-8" indent="no"/>
  <xsl:template match="/">
    <xsl:for-each select="/cfsErsatzMaterial/set[last()]/element">
      <xsl:value-of select="@design"/><xsl:text>  </xsl:text>
    </xsl:for-each>
    <xsl:text>
    </xsl:text>
  </xsl:template>
</xsl:stylesheet>