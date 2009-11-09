#!/usr/bin/python
##########################################
##
##########################################

# this program does the following
# 1) copy tex-file to /html-subdirectory
# 2) load copied tex file
# 3) look for \begin{verbatim} .. \end{verbatim} environment
# 4) pygmentize content -> $pygoutput
# 5) Replace content by \HCode{ $pygoutput }
# 6) run htlatex

import os,  re,  sys
from optparse import OptionParser
from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import HtmlFormatter,  LatexFormatter

# global definitions
# -------------------

# environments, which should be parsed
envNames = ["xmllisting"]

# names of commands, which should get replaced
cmdNames = ["xml"]

# substitution for html (as the html code will be embedded 
# in TeX again
htmlSubst = { 
  "&#39;" : "&prime;",
  "}" : "\}", 
  "\n": "<br>\n", 
  "%": "\%",
  "#": "\#",
}


# parse command line
usage = "usage: %prog [-o newFile ] format file \n \n\
       Perform syntax highlightning of tex file and \n\
       write the output either in tex format or enclose it \n\
       in \HCode{} elements to be parsed by tex4ht."

parser = OptionParser(usage)
newFileName = ""
parser.set_defaults(newFileName="")
parser.add_option("-o", "",  dest="newFileName")
(options, args) = parser.parse_args()
if len(args) != 2:
        parser.error("please provide format and file")

parserType = args[0]
fileName = args[1]

# if no newFileName was given, just overwrtite the old one
if newFileName == "":
  newFileName = fileName


# callback function for html output
def codeReplaceHtml( txt ):
  result = highlight(txt.group(2), lexer, formatter)

  # substitute special characters which the tex-subsystem
  # would replace otherwise
  for k, s in htmlSubst.iteritems():
    result = result.replace(k,  s)
  ret =""
  
  # in addition we have to replace all literal spaces
  # by &nbsp;, as the \HCode{} statement removes all
  # trailing whitespaces, which woud ruin the indentation
  n = 0
  for i in result.split("\n"):
    n = n+1
    ret = ret + "\HCode{"
    
    # if we are in the first line, we omit everything, until
    # the <pre> element and just take 
    beginPos = 0
    if n == 1:
       beginPos = i.find("<pre>",0)
       if beginPos > 0:
         ret = ret + i[:beginPos+5] + "&nbsp;"
         i = i[beginPos+5:]
    
    temp = ""
    index = 0
    while index < len(i) and i[index] == " ":
      ret = ret + "&nbsp;"
      index = index+1
    ret=  ret + i + "}\n"
    
  return ret

# callback function for latex output
def codeReplaceTeX( txt ):
  result = highlight(txt.group(2), lexer, formatter)
  return result


# pygments stuff:
# ----------------
# define lexer ...
lexer = get_lexer_by_name("xml", stripall=False)
callback = lambda x: x

# and syntax highglighter depending on format html / tex
# and define callback function
formatter = LatexFormatter()
inlineFormatter = LatexFormatter()
if parserType == "tex":
  formatter = LatexFormatter()
  inlineFormatter = LatexFormatter()
  callback = codeReplaceTeX
  formatFile ="style.tex"
elif parserType == "html" :
  formatter = HtmlFormatter(style="default", linenos=False, 
                            cssclass="xmlcode", cssFilefull="False")
  inlineFormatter = HtmlFormatter(style="default", linenos=False, nowrap=True, nobackground=True,
                                  noclasses=True)
                            
  callback = codeReplaceHtml
  formatFile = "style.css"
else:
  print "Valid options for 'format' are 'tex' and 'html' "

# load file and omit all comments
f = open(fileName)
code = ""
commentExp = re.compile("^[ ]*\%")
for line in f:
  m = re.match(commentExp,  line)
  if m == None:
    code = code + line

f.close()

# iterate over all environments and do syntax highlightning
for actEnv in envNames:
  actRegEx = "(\\\\begin\{"+actEnv+"\})(.*)(?<!%.)(\\\end\{"+actEnv+"\})"
  myexp = re.compile(actRegEx,re.DOTALL)
  beginPos = 0
  endPos = 0
  while endPos < len(code) and endPos != -1:

    # return substrings with content between 
    # \begin{envName} ... and \end{envName}
    beginPos = code.find("\\begin{"+actEnv+"}", endPos)
    endPos = code.find("\\end{"+actEnv+"}", beginPos)
    #check if environment was found
    if not(endPos == -1 and beginPos == -1):
      endPos = endPos + 6 + len(actEnv)
      
      # perform syntax highlightning and replace original code
      temp = re.sub(myexp, callback, code[beginPos:endPos])
      code = code.replace(code[beginPos:endPos], temp)

# write html-css / tex-header to style file, before treating
# the inline case, which requires no external style definition
f = open(formatFile, 'w')
f.write(formatter.get_style_defs() )
f.close()


# perform inline syntax coloring hjus
if parserType == "html":
  # iterate over all commands and do syntax highlightning
  # but first we have to use the inlineFormatter
  formatter = inlineFormatter
  for actCmd in cmdNames:
    actRegEx = "(\\\\"+actCmd+"\{)(.*)(?<!%.)(})"
    myexp = re.compile(actRegEx,re.DOTALL)
    beginPos = 0
    endPos = 0
    while endPos < len(code) and endPos != -1:
      #print "looking for ","\\"+actCmd+"{"
      beginPos = code.find("\\"+actCmd+"{", endPos)
      endPos = code.find("}", beginPos)
   
      # return substrings with content between 
      # \cmdName{ ... and }
      if not(endPos == -1 and beginPos == -1):
        endPos = endPos + 1
        
        # perform syntax highlightning and replace original code
        temp = re.sub(myexp, callback, code[beginPos:endPos])
        code = code.replace(code[beginPos:endPos], temp)

# write highlighted code to file
f = open(newFileName, 'w')
f.write(code)
f.close()

