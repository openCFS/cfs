<!-- This stylesheet extracts the results form an eigenfrequency analysis 
     out of the info.xml file to a single line with a sequence of frequencies.
     This is meant collecting multiple eigenfrequencies calculations 
     when doing a parameter study.
     Example: xsltproc eigen.xsl problem.info.xml 
     Output: 322.346  612.249  612.249  993.428  1143.81
 -->
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="text" encoding="UTF-8" indent="no"/>
  <xsl:template match="/">
 <xsl:for-each select="/cfsInfo/analysis/eigenFrequency/result">
   <xsl:value-of select="@frequency"/><xsl:text>  </xsl:text>
  </xsl:for-each>
<xsl:text>
</xsl:text>
  </xsl:template>
  <xsl:template mode="copy-no-ns" match="*">
    <xsl:element name="{name(.)}">
      <xsl:copy-of select="@*"/>
      <xsl:apply-templates mode="copy-no-ns"/>
    </xsl:element>
  </xsl:template>
</xsl:stylesheet>


