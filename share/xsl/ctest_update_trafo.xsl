<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
		xmlns:wh="whatever"
		version="1.0">

  <xsl:namespace-alias result-prefix="xsl"
                       stylesheet-prefix="wh"/>

  <xsl:output method="xml" encoding="UTF-8" indent="yes"/> 

  <!-- This template  will be applied to the  Testing/Configure.xml file which
       is generated  during a nightly CTest  run and generates  a new template
       which can be applied to the Testing/Update.xml which has been generated
       during updating of the working  copy on lse10.  From Update.xml it will
       produce a new  Update.xml where all machine specific  entries have been
       replaced corresponding to the local machine -->


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

      <wh:variable name="wcdir"/>

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
      
      <wh:template match="Site">
	<Site><wh:value-of select="$site"/></Site>
      </wh:template>
      
      <wh:template match="BuildName"> 
	<BuildName>
	  <wh:value-of select="$buildname"/>
	</BuildName>
      </wh:template>
      
      <wh:template match="BuildStamp">
	<BuildStamp>
	  <wh:value-of select="$buildstamp"/>
	</BuildStamp>
      </wh:template>

      <wh:template match="Name">
	<wh:choose>
	  <wh:when test="parent::Directory">
	    <Name>
              <wh:value-of select="$wcdir"/>/<wh:value-of select="node()"/>
	    </Name>
	  </wh:when>
	  <wh:otherwise>
	    <Name>
              <wh:value-of select="node()"/>
	    </Name>
	  </wh:otherwise>
	</wh:choose>
      </wh:template>
      
      <wh:template match="Updated/Directory">
	<Directory>
	  <wh:value-of select="$wcdir"/>/<wh:value-of select="node()"/>
	</Directory>
      </wh:template>

      <wh:template match="Updated/FullName">
	<FullName>
	  <wh:value-of select="$wcdir"/>/<wh:value-of select="node()"/>
	</FullName>
      </wh:template>
      
      <wh:template match="File/@Directory">
	<wh:attribute name="Directory">
	  <wh:value-of select="$wcdir"/>/<wh:value-of select="."/>
	</wh:attribute>
      </wh:template>      
      
    </wh:stylesheet>
  </xsl:template>

</xsl:stylesheet>
