<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <!-- This stylesheet extracts the subiteratiosn from a meta optimization 
    for multiple excitations to a gnuplot file
  -->
  <xsl:output method="text" encoding="UTF-8" indent="no"/>
  <xsl:template match="/">
    <xsl:text>#iter </xsl:text>
    <xsl:for-each select="/cfsInfo/optimization/process/iteration[last()]">
      <xsl:for-each select="excitation">
        <xsl:value-of select="@index"/><xsl:text>_cost </xsl:text>
      </xsl:for-each>
      <xsl:for-each select="excitation">
        <xsl:value-of select="@index"/><xsl:text>_nor_weight </xsl:text>
      </xsl:for-each>
    </xsl:for-each>
<xsl:text>
</xsl:text>
    <xsl:for-each select="/cfsInfo/optimization/process/iteration">
      <xsl:value-of select="@number"/><xsl:text>  </xsl:text>
      <xsl:for-each select="excitation">
        <xsl:value-of select="@objective"/><xsl:text>  </xsl:text>
      </xsl:for-each>
      <xsl:for-each select="excitation">
        <xsl:value-of select="@normalized_weight"/><xsl:text>  </xsl:text>
      </xsl:for-each>
<xsl:text>
</xsl:text>
    </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>


