<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:output method="text" encoding="UTF-8" indent="no"/>
	<xsl:template match="/">
		<xsl:text>#freq #step </xsl:text>
		<xsl:for-each select="/cfsInfo/calculation/process/sequence">
		<!-- write the header with a combination of data and location, e.g. mechDeformedVolume_hot_se -->
			<xsl:for-each select="result">
				<xsl:value-of select="@data"/><xsl:text>_</xsl:text>
				<xsl:value-of select="@location"/><xsl:text> </xsl:text>
			</xsl:for-each>
		</xsl:for-each>
		<xsl:text>
</xsl:text>
	<xsl:call-template name="loop">
		<xsl:with-param name="maxcount" select="/cfsInfo/calculation/process/sequence/result/item[last()]/@step_nr+1"/>
		<xsl:with-param name="incriment-factor" select="1"/>
		<xsl:with-param name="initial-value" select="/cfsInfo/calculation/process/sequence/result/item[1]/@step_nr"/>
	</xsl:call-template>
</xsl:template>
<xsl:template name="loop">
	<xsl:param name="maxcount"/>
	<xsl:param name="incriment-factor"/>
	<xsl:param name="initial-value"/>
	<xsl:if test="$initial-value &lt; $maxcount">
		<xsl:value-of select="/cfsInfo/analysis/harmonic/process/solveStep[$initial-value+1]/@currentFrequency"/><xsl:text> </xsl:text>
		<xsl:value-of select="/cfsInfo/calculation/process/sequence/result/item[$initial-value+1]/@step_nr"/><xsl:text> </xsl:text> 
		<xsl:value-of select="/cfsInfo/calculation/process/sequence/result[@data='acouPower']/item[$initial-value+1]/@value"/><xsl:text> </xsl:text>
		<xsl:value-of select="/cfsInfo/calculation/process/sequence/result[@data='mechDeformedVolume']/item[$initial-value+1]/@value"/>
		<xsl:text> 
</xsl:text>
		<xsl:call-template name="loop">
			<xsl:with-param name="maxcount" select="$maxcount"/>
			<xsl:with-param name="initial-value" select="$initial-value+$incriment-factor"/>
			<xsl:with-param name="incriment-factor" select="$incriment-factor"/>
		</xsl:call-template>
		</xsl:if>
</xsl:template>
</xsl:stylesheet>