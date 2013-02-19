<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <!-- makes a gnuplot file from the optimization iterations where "homogenizedTensor/isotropy" is added if present -->
  <xsl:output method="text" encoding="UTF-8" indent="no"/>
  <xsl:template match="/">
    <!-- print header,-->
    <xsl:for-each select="/cfsInfo/optimization/process/iteration[1]/@*">
      <xsl:value-of select="name()"/><xsl:text>  </xsl:text>    
    </xsl:for-each>
    <xsl:for-each select="/cfsInfo/optimization/process/iteration[1]/homogenizedTensor/isotropy/@*">
      <xsl:text>  </xsl:text><xsl:value-of select="name()"/>    
    </xsl:for-each>
<xsl:text>
</xsl:text>
    <!-- print content -->
    <xsl:for-each select="/cfsInfo/optimization/process/iteration">
      <xsl:for-each select="@*">
        <xsl:value-of select="."/><xsl:text>  </xsl:text>
      </xsl:for-each>
      <xsl:for-each select="homogenizedTensor/isotropy/@*">
        <xsl:value-of select="."/><xsl:text>  </xsl:text>
      </xsl:for-each>
<xsl:text>
</xsl:text>
    </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>


