#!/usr/bin/env python

import libxml2
import os
import string
import sys
from datetime import datetime, date, time


# the directory where the testsuite is located
testsuitedir = "/home/fschury/cfs/cfs-test/TESTSUIT"
# the directory where the documentation is going to be created
htmldocdir = "/home/fschury/cfs/htmldoc"

# when converting the testsuite directory structure to files, we 
# need to substitute the directory delimiter '/' by another char
splitchar = '~'


# ---------------------------------------------------------------
# ---------------------------------------------------------------


# parse the CMakeLists.txt in the TESTSUITE root directory
# to obtain a list of testcases
# at the same time do some checks
def getTestcases():
  os.chdir(testsuitedir)
  count = 0
  ret = []
  with open("CMakeLists.txt") as cmakelist:
    for line in cmakelist:
      l = line.lstrip()
      if(not l.startswith('#') and (l.find('SUBDIRS') > -1) and (not l.find('cfstool') > -1)):
        candidate = l.split('(')[1].strip().rstrip(')')
        # test if the corresponding xml file exists
        if(os.path.isfile(testsuitedir + '/' + candidate + '/' + candidate.rsplit('/', 1)[1] + '.xml')):
          ret.append(candidate)
          count += 1
  print 'found {0} testcases'.format(count)
  ret.sort()
  return ret

# function that creates the html representation of the docu tags
def docTagToHtml(tags):
  html = []
  for v in tags:
    #if(v.name == 'keywords'):
      #print v.name + ': ' 
      #print v.content.strip()
      #continue
    if(v.name == 'title'):
      #html.append('<div id="' + v.name + '">' + v.content.strip().title() + '</div>')
      html.append('<h2>' + v.content.strip().title() + '</h2>')
      continue
    html.append('<div id="' + v.name + '"><strong>' + v.name.title() + '</strong>: ' + v.content.strip() + '</div>')
  return html

def writeHtmlHeader(file, css):
  file.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"\n\
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n\
  <html xmlns="http://www.w3.org/1999/xhtml" xml:lang="de">\n\n')
  file.write('<head>\n  <link rel="stylesheet" href="' + css + '" type="text/css"/>\n\
    <meta http-equiv="Content-Type" content= "application/xhtml+xml; charset=iso-8859-1"/>\n\
    <meta name="generator" content= "HTML, see www.w3.org" />\n\
    <meta name="author-personal" content="Fabian Schury" />\n\
    <title> CFS Testsuite Documentation </title>\n\
  </head>\n\n')

# write the css file
def createCSSfile():
  os.chdir(htmldocdir)
  file = open('main.css', "w")
  file.write('\
html, body, ul.nav, ul.nav li, ul.subnav {\n\
  margin:  0;\n\
  padding: 0;\n\
}\n\
  \n\
html, body {\n\
  height: 90%;\n\
  width: 90%;\n\
}\n\
  \n\
body {\n\
  background-color: #FFF;\n\
  background-repeat: repeat-y;\n\
  background-position: top center;\n\
  font-family: Arial, Helvetica, sans-serif;\n\
  font-size: 14px;\n\
  font-weight: normal;\n\
  color: #036;\n\
  margin: 10px 0px 10px 20px;\n\
  padding: 10px 0px 10px 20px;\n\
}\n\
  \n\
body div {\n\
  margin: 10px 0px 10px 20px;\n\
  padding: 10px 0px 10px 20px;\n\
}\n\
  \n\
#title {\n\
  font-size: 18px;\n\
}\n\
  \n\
h1, h2, h3 {\n\
  margin: 10px 0px 10px 10px;\n\
  padding: 10px 0px 10px 10px;\n\
  color: #FFF;\n\
  background-color: #036;\n\
}\n\
table {\n\
  border-spacing: 10px;\n\
}\n\
thead {\n\
  font-weight: bold;\n\
}\n\
')
  file.close()


# creates the index.html in the htmldocdir root directory
def createIndex(ind, info):
  os.chdir(htmldocdir)
  filename = 'index.html'
  file = open(filename, "w")
  writeHtmlHeader(file, 'main.css')
  file.write('<body>\n')
  file.write('<div> <a href="index_tags.html">See testcases by tags</a></div>\n')
  file.write('<h1>CFS Testsuite Documentation</h1>\n')
  file.write('<table>\n')
  file.write('<thead><tr><td>#</td><td>Category</td><td>Name</td><td>Dimension</td><td>Elements</td>\
  <td>Wall Clock</td></tr></thead>\n')
  file.write('<tbody>\n')

  count = 0
  strings = []
  # create the proper strings here, put them in the strings-array
  # so we can sort them later
  for t in ind:
    cand = t.split(splitchar, 1)
    string1 = '<td>' + cand[0].title() + '</td>\n <td> <a href="docu/' + t + '">' \
    + cand[1].replace(splitchar, ' - ').title().rsplit('.', 1)[0] + ' </a> </td>\n'

    string2 = ' ' + info[count] + '</tr>\n'
    count += 1
    strings.append([string1, string2])

  # write results into file
  for i in range(0, len(strings)):
    # insert different tr depending on i%2...
    file.write('<tr>\n')
    #if(i % 2 == 0):
      #file.write('<tr>\n')
    #else:
      #file.write('<tr style="background-color: #9cf;">\n')
    file.write('<td>' + str(i+1) + '</td>\n')
    file.write(strings[i][0])
    file.write(strings[i][1])

  # file body
  file.write('</tbody>\n')
  file.write('</table>\n')
  file.write('<h3>created: ' + str(datetime.now()).rsplit('.', 1)[0] + '</h3>\n')
  file.write('</body>\n</html>\n')
  file.close()


def getInfoXmlContent(infoxmlname, infoxml):
  ret = ''
  # check if info.xml is present
  if(os.path.exists(infoxmlname)):
    doc = libxml2.parseFile(infoxmlname)
    xml = doc.xpathNewContext()

    domain = xml.xpathEval('//domain/@*')
    timer = xml.xpathEval('//summary/timer/@*')

    if(not domain == []):
      ret = ret + '<td>' + str(domain[0]).split('"')[1] + '</td>\
      <td>' + str(domain[1]).split('"')[1] + '</td>'
    else:
      ret = ret + '<td>n. a.</td><td>n. a.</td>'

    if(not timer == []):
      ret = ret + '<td>' + str(timer[0]).split('"')[1] + '</td>'
    else:
      ret = ret + '<td>n. a.</td>'

  else:
    ret = ret + '<td>-</td><td>-</td><td>no info.xml found</td>'

  ret = ret + '\n'
  infoxml.append(ret)


####################
# main starts here
####################

print 'Testsuite in directory ' + testsuitedir
print 'Documentation will be created in directory ' + htmldocdir

# make sure the directories exist
if(not os.path.exists(testsuitedir)):
  print 'testsuite directory "' + testsuitedir + '" not found, exiting'
  exit(0)

if(not os.path.exists(htmldocdir)):
  print 'directory "' + htmldocdir + '" not found, exiting'
  exit(0)


list = getTestcases()
createCSSfile()

os.chdir(htmldocdir)
if(not os.path.exists('docu')):
  os.mkdir('docu')
os.chdir('docu')

# count the documented and undocumentated testcases
count = 0
count_undoc = 0
# contains file names for creating the index.html
index = []
# contains info from info.xml
infoxml = []

for l in list:
  # xmlname is the absolute file name of the corresponding xml file
  xmlname = testsuitedir + '/' + l + '/' + l.rsplit('/', 1)[1] + '.xml'
  # name of info.xml for further information of the testcase
  infoxmlname = xmlname.rsplit('.', 1)[0] + '.info.xml'
  # something like: /home/fschury/cfs/cfs-test/TESTSUIT/coupledfield/piezo-direct/preisach2d/preisach2d.info.xml

  doc = libxml2.parseFile(xmlname)
  xml = doc.xpathNewContext()
  xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')

  tags = xml.xpathEval('//cfs:documentation/*')

  # check if documentation tag was found
  if(not tags == []):
    htmltags = docTagToHtml(tags)

    # create files
    filename = l.replace('/', splitchar) + '.html'
    file = open(filename, "w")
    writeHtmlHeader(file, '../main.css')
    file.write('<body>\n')
    file.write('<div><a href="../index.html">Back to Index</a></div>\n')
    for t in htmltags:
      file.write('  ' + t + '\n')
    file.write('</body>\n</html>\n')
    file.close()

    # add filename to index
    index.append(filename)

    # check for info.xml and extract further data
    getInfoXmlContent(infoxmlname, infoxml)

    count += 1
  else:
    count_undoc += 1

# create index.html for browsing the docu
createIndex(index, infoxml)

# print some useful information
print 'found {0} documentation tags (and {1} undocumentated \
testcases, sum = {2})'.format(count, count_undoc, count + count_undoc)



# produce second index page where testcases are sorted by tags
# autotags
autotags = ['Homogenization', 'Scpip', 'Snopt', 'OptimalityCondition',
            'Coupled', 'Acoustic', 'Mechanic',
            'ExportLinSys', 'MultiObjective',
            'Cholmod', 'Pardiso']

os.chdir(htmldocdir)
filename = 'index_tags.html'
file = open(filename, "w")
writeHtmlHeader(file, 'main.css')
file.write('<body>\n')
file.write('<div> <a href="index.html">See testcases by alphabet</a></div>\n')
file.write('<h1>Testsuite by tags</h1>\n')


for tag in autotags:
  file.write('<h2>' + tag + '</h2>\n<ul>\n')
  foundone = False
  for l in list:
    filename = l.replace('/', splitchar) + '.html'
    xmlname = testsuitedir + '/' + l + '/' + l.rsplit('/', 1)[1] + '.xml'
    doc = libxml2.parseFile(xmlname)
    xml = doc.xpathNewContext()
    xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')

    # check if the tag applies
    r = []
    if(tag == 'Homogenization'):
      r = xml.xpathEval('//cfs:multipleExcitation[@type=\'homogenizationTestStrains\']')

    if(tag == 'Scpip'):
      r = xml.xpathEval('//cfs:optimizer[@type=\'scpip\']')

    if(tag == 'Snopt'):
      r = xml.xpathEval('//cfs:optimizer[@type=\'snopt\']')

    if(tag == 'OptimalityCondition'):
      r = xml.xpathEval('//cfs:optimizer[@type=\'optimalityCondition\']')

    if(tag == 'Coupled'):
      r = xml.xpathEval('//cfs:couplingList')

    if(tag == 'Mechanic'):
      r = xml.xpathEval('//cfs:pdeList/cfs:mechanic')

    if(tag == 'Acoustic'):
      r = xml.xpathEval('//cfs:pdeList/cfs:acoustic')

    if(tag == 'ExportLinSys'):
      r = xml.xpathEval('//cfs:exportLinSys')

    if(tag == 'MultiObjective'):
      r = xml.xpathEval('//cfs:costFunction[@type=\'multiObjective\']')

    if(tag == 'Cholmod'):
      r = xml.xpathEval('//cfs:solver[@type=\'cholmod\']')

    if(tag == 'Pardiso'):
      r = xml.xpathEval('//cfs:solver[@type=\'pardiso\']')

    if(not r == [] and os.path.isfile(htmldocdir + '/docu/' + filename)):
      nicename = filename.replace(splitchar, ' - ').title().rsplit('.', 1)[0]
      file.write('<li><a href="docu/' + filename + '">' + nicename + '</a></li>')
      foundone = True

  if(not foundone):
    file.write('no (documented) testcases found')


  file.write('</ul>\n')


file.write('<h3>created: ' + str(datetime.now()).rsplit('.', 1)[0] + '</h3>\n')
file.write('</body>\n</html>\n')
file.close()
