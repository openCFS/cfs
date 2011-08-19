<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <!-- This stylesheet extracts the results within info.xml from the first sequence step.
    In case of several steps (frequencies, timesteps, iterations) the last is extracted 
  -->
  <xsl:output method="text" encoding="UTF-8" indent="no"/>
  <xsl:template match="/">
    <xsl:text>#step </xsl:text>
    <xsl:for-each select="/cfsInfo/calculation/process/sequence">
      <!-- write the header with a combination of data and location, e.g. mechDeformedVolume_hot_se -->
      <xsl:for-each select="result">
        <xsl:value-of select="@data"/><xsl:text>_</xsl:text>
        <xsl:value-of select="@location"/><xsl:text> </xsl:text>
      </xsl:for-each>
    </xsl:for-each>
<xsl:text>
</xsl:text>
    <!-- for every result take the last item -->
    <!-- start with the last item step_nr from any result -->
    <xsl:value-of select="/cfsInfo/calculation/process/sequence/result/item[last()]/step/@analyis_id"/><xsl:text> </xsl:text> 
    <xsl:for-each select="/cfsInfo/calculation/process/sequence">
      <xsl:for-each select="result">
        <xsl:for-each select="item[last()]">
          <xsl:value-of select="@value"/><xsl:text> </xsl:text>
          <xsl:value-of select="@x"/><xsl:text> </xsl:text>
          <xsl:value-of select="@y"/><xsl:text> </xsl:text>
          <xsl:value-of select="@z"/><xsl:text> </xsl:text>
        </xsl:for-each>
      </xsl:for-each>
    </xsl:for-each>
<xsl:text>
</xsl:text>
  </xsl:template>
</xsl:stylesheet>  
