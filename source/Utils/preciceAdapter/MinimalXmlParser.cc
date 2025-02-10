#include "MinimalXmlParser.hh"
#include "General/Environment.hh"
#include <fstream>
#include <sstream>
#include <cctype>
#include <string>

namespace CoupledField {


    // Constructor: reads the entire file into a string.
    MinimalXmlParser::MinimalXmlParser(const std::string &filename)
    : pos(0)
    {
    std::ifstream ifs(filename);
    if (!ifs)
        EXCEPTION("Cannot open file: " + filename);
    std::stringstream ss;
    ss << ifs.rdbuf();
    xmlContent = ss.str();
    }

    // Skip whitespace characters.
    void MinimalXmlParser::skipWhitespace() {
    while (pos < xmlContent.size() && std::isspace(xmlContent[pos])) {
        pos++;
    }
    }

    // Skip over an XML comment.
    // Assumes that 'pos' is at the '<' character that starts the comment.
    void MinimalXmlParser::skipComment() {
    // We expect the comment to start with "<!--"
    pos += 4; // skip "<!--"
    while (pos + 2 < xmlContent.size() &&
            !(xmlContent[pos] == '-' && xmlContent[pos+1] == '-' && xmlContent[pos+2] == '>')) {
        pos++;
    }
    pos += 3; // skip "-->"
    }

    // Parse a tag name (until a whitespace, '/', '>' or '<' is encountered).
    std::string MinimalXmlParser::parseTagName() {
    skipWhitespace();
    std::string tagName;
    while (pos < xmlContent.size() && !std::isspace(xmlContent[pos]) &&
            xmlContent[pos] != '>' && xmlContent[pos] != '/' && xmlContent[pos] != '<') {
        tagName.push_back(xmlContent[pos]);
        pos++;
    }
    return tagName;
    }

    // Parse attributes until we hit '>' or '/'.
    std::map<std::string, std::string> MinimalXmlParser::parseAttributes() {
    std::map<std::string, std::string> attrs;
    skipWhitespace();
    while (pos < xmlContent.size() && xmlContent[pos] != '>' && xmlContent[pos] != '/') {
        // Parse attribute name.
        std::string attrName;
        while (pos < xmlContent.size() && !std::isspace(xmlContent[pos]) &&
            xmlContent[pos] != '=' && xmlContent[pos] != '>' && xmlContent[pos] != '/') {
        attrName.push_back(xmlContent[pos]);
        pos++;
        }
        skipWhitespace();
        if (pos >= xmlContent.size() || xmlContent[pos] != '=')
        break; // Malformed attribute; break out.
        pos++; // Skip '='.
        skipWhitespace();
        // Expect a quote character.
        if (pos >= xmlContent.size() || (xmlContent[pos] != '"' && xmlContent[pos] != '\'')) {
        EXCEPTION("Expected quote for attribute value at position " + std::to_string(pos));
        }
        char quote = xmlContent[pos];
        pos++; // Skip opening quote.
        std::string attrValue;
        while (pos < xmlContent.size() && xmlContent[pos] != quote) {
        attrValue.push_back(xmlContent[pos]);
        pos++;
        }
        if (pos < xmlContent.size())
        pos++; // Skip closing quote.
        attrs[attrName] = attrValue;
        skipWhitespace();
    }
    return attrs;
    }

    // Parse an element, including its tag, attributes, children, and content.
    XmlNode MinimalXmlParser::parseElement() {
    XmlNode node;
    skipWhitespace();
    if (pos >= xmlContent.size() || xmlContent[pos] != '<')
        EXCEPTION("Expected '<' at position " + std::to_string(pos));
    
    // Check for comment.
    if (pos + 3 < xmlContent.size() && xmlContent.substr(pos, 4) == "<!--") {
        skipComment();
        skipWhitespace();
        // After skipping the comment, we expect a new element.
        if (pos >= xmlContent.size() || xmlContent[pos] != '<')
        EXCEPTION("Expected '<' after comment at position " + std::to_string(pos));
    }
    
    pos++; // Skip '<'.
    node.name = parseTagName();
    node.attrs = parseAttributes();
    skipWhitespace();
    bool selfClosing = false;
    if (pos < xmlContent.size() && xmlContent[pos] == '/') {
        selfClosing = true;
        pos++;
    }
    if (pos >= xmlContent.size() || xmlContent[pos] != '>')
        EXCEPTION("Expected '>' at position " + std::to_string(pos));
    pos++; // Skip '>'.
    
    if (selfClosing)
        return node;
    
    // Parse content and/or child elements.
    std::string textContent;
    while (pos < xmlContent.size()) {
        skipWhitespace();
        if (pos < xmlContent.size() && xmlContent[pos] == '<') {
        // Check for comment.
        if (pos + 3 < xmlContent.size() && xmlContent.substr(pos, 4) == "<!--") {
            skipComment();
            continue;
        }
        // Check if this is a closing tag.
        if (pos + 1 < xmlContent.size() && xmlContent[pos + 1] == '/') {
            pos += 2; // Skip "</".
            std::string closeName = parseTagName();
            if (closeName != node.name)
            EXCEPTION("Mismatched closing tag for " + node.name);
            while (pos < xmlContent.size() && xmlContent[pos] != '>')
            pos++;
            if (pos < xmlContent.size())
            pos++; // Skip '>'.
            break;
        }
        // Otherwise, parse a child element.
        XmlNode child = parseElement();
        node.children.push_back(child);
        } else {
        // Accumulate text content.
        textContent.push_back(xmlContent[pos]);
        pos++;
        }
    }
    node.content = textContent;
    return node;
    }

    // Entry point: parse the whole document and return the root element.
    XmlNode MinimalXmlParser::parse() {
    skipWhitespace();
    return parseElement();
    }
}