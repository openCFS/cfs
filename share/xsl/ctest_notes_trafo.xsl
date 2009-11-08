<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
		xmlns:wh="whatever"
		version="1.0">
  
  <xsl:namespace-alias result-prefix="xsl"
                       stylesheet-prefix="wh"/>
  
  <xsl:output method="xml" encoding="UTF-8" indent="yes"/> 
  
  <!-- This template  will be applied to the  Testing/Configure.xml file which
       is generated  during a nightly CTest  run and generates  a new template
       which can  be applied to  the Testing/NotesTemp.xml which  is generated
       during the  nightly builds.  From  NotesTemp.xml it will produce  a new
       Notes.xml  where  all  machine  specific  entries  have  been  replaced
       corresponding to the local machine -->
  
  
  <xsl:template match="Site">
    <wh:stylesheet version="1.0">
      
      <wh:output method="xml" encoding="UTF-8" indent="yes"/> 
      
      <wh:variable name="generator">
	<xsl:value-of select="@Generator"/>
      </wh:variable>
      <wh:variable name="site">
	<xsl:value-of select="@Name"/>
      </wh:variable>
      <wh:variable name="buildname">
	<xsl:value-of select="@BuildName"/>
      </wh:variable>
      <wh:variable name="buildstamp">
	<xsl:value-of select="@BuildStamp"/>
      </wh:variable>

      <wh:template match="*|@*">
	<wh:copy>
	  <wh:apply-templates select="@*|node()"/>
	</wh:copy>
      </wh:template>
      
      <wh:template match="@Generator">
	<wh:attribute name="Generator">
	  <wh:value-of select="$generator"/>
	</wh:attribute>
      </wh:template>
      
      <wh:template match="Site/@Name">
        <wh:attribute name="Name">
          <wh:value-of select="$site"/>
        </wh:attribute>
      </wh:template>
      
      <wh:template match="@BuildName"> 
	<wh:attribute name="BuildName">
	  <wh:value-of select="$buildname"/>
	</wh:attribute>
      </wh:template>
      
      <wh:template match="@BuildStamp">
	<wh:attribute name="BuildStamp">
	  <wh:value-of select="$buildstamp"/>
	</wh:attribute> 
      </wh:template>

    </wh:stylesheet>
  </xsl:template>
</xsl:stylesheet>
