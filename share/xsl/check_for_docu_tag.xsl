<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    version="1.0"
    xmlns:cfs="http://www.cfs++.org"
    exclude-result-prefixes="cfs">
    
    <xsl:output method="text" encoding="UTF-8" indent="no"/>
    
    <xsl:template match="/cfs:cfsSimulation">
        <xsl:choose>
            <xsl:when test="cfs:documentation">
                <xsl:text>SUCCESS: Documentation tag is present.&#10;</xsl:text>
                <xsl:apply-templates select="cfs:documentation/cfs:title"/>
                <xsl:apply-templates select="cfs:documentation/cfs:authors"/>
                <xsl:apply-templates select="cfs:documentation/cfs:date"/>
                <xsl:apply-templates select="cfs:documentation/cfs:keywords"/>
                <xsl:apply-templates select="cfs:documentation/cfs:references"/>
                <xsl:apply-templates select="cfs:documentation/cfs:description"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:text>WARNING: Documentation tag is missing!&#10;</xsl:text>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>
    
    <xsl:template match="cfs:title">
        <xsl:choose>
            <xsl:when test="normalize-space(.)">
                <xsl:text>SUCCESS: Title tag: </xsl:text> <xsl:value-of select="node()"/><xsl:text>&#10;</xsl:text>
            </xsl:when>
            <xsl:otherwise>
                <xsl:text>FAILURE: Title tag does not contain text!&#10;</xsl:text>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>
    
    <xsl:template match="cfs:authors">
        <xsl:for-each select="cfs:author">
            <xsl:choose>
                <xsl:when test="normalize-space(.)">
                    <xsl:text>SUCCESS: Author tag: </xsl:text> <xsl:value-of select="node()"/><xsl:text>&#10;</xsl:text>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:text>FAILURE: Author tag does not contain text!&#10;</xsl:text>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:for-each>
    </xsl:template>
    
    <xsl:template match="cfs:date">
        <xsl:choose>
            <xsl:when test="normalize-space(.)">
                <xsl:text>SUCCESS: Date tag: </xsl:text> <xsl:value-of select="node()"/><xsl:text>&#10;</xsl:text>
            </xsl:when>
            <xsl:otherwise>
                <xsl:text>FAILURE: Date tag does not contain text!&#10;</xsl:text>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>
    
    <xsl:template match="cfs:keywords">
        <xsl:for-each select="cfs:keyword">
            <xsl:choose>
                <xsl:when test="normalize-space(.)">
                    <xsl:text>SUCCESS: Keyword tag: </xsl:text> <xsl:value-of select="node()"/><xsl:text>&#10;</xsl:text>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:text>FAILURE: Keyword tag does not contain text!&#10;</xsl:text>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:for-each>
    </xsl:template>
    
    <xsl:template match="cfs:references">
        <xsl:choose>
            <xsl:when test="normalize-space(.)">
                <xsl:text>SUCCESS: References tag: </xsl:text> <xsl:value-of select="node()"/><xsl:text>&#10;</xsl:text>
            </xsl:when>
            <xsl:otherwise>
                <xsl:text>FAILURE: References tag does not contain text!&#10;</xsl:text>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>
    
    <xsl:template match="cfs:description">
        <xsl:choose>
            <xsl:when test="normalize-space(.)">
                <xsl:text>SUCCESS: Description tag: </xsl:text> <xsl:value-of select="node()"/><xsl:text>&#10;</xsl:text>
            </xsl:when>
            <xsl:otherwise>
                <xsl:text>FAILURE: Description tag does not contain text!&#10;</xsl:text>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>
    
</xsl:stylesheet>
