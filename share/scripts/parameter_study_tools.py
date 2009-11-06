import libxml2
import os

# replace a single xpath value -> must exsit once!
# xml is a xpathContext: doc = libxml2.parseFile("params.xml") -> xml = doc.xpathNewContext()
# save your doc later via doc.saveFile("param_tmp.xml")
def replace(xml, path, value):
  res = xml.xpathEval(path)
  if  len(res) == 0:
    raise RuntimeError(path + " not found")
  if len(res) > 1:
    raise RuntimeError(path + " has " + len(res) + " hits")
  data = res[0]
  data.setContent(value)
  return    

# returns an xpath value
def xpath(xml, path):
  res = xml.xpathEval(path)
  if  len(res) == 0:
    raise RuntimeError(path + " not found")
  if len(res) > 1:
    raise RuntimeError(path + " has " + len(res) + " hits")
  data = res[0]
  return data.getContent()    

  
# dump a xml node
def dump(xml, path):
  res = xml.xpathEval(path)
  if  len(res) == 0:
    raise RuntimeError(path + " empty")
  for i in res:
    print str(i) + ": " + i.getContent()
  
  
  
# mimic conditional operator
def cond(test, trueval, falseval):
  if test:
    return trueval
  else:
    return falseval

# execute cmd and rais error when not 0
def execute(cmd):
 ret = os.system(cmd) <> 0
 if(ret <> 0):
   raise RuntimeError("execution of '" + cmd + "' -> " + str(ret))
 return     

# return the first line of a file
def first_line(file_name):
   file = open(file_name, "r")
   line = file.readline()
   file.close() 
   if len(line) == 0:
     raise RuntimeError("file '" + file_name + "' is empty")
   return line

# return the second line of a file, this skips the gnuplot header
def second_line(file_name):
   file = open(file_name, "r")
   lines = file.readlines()
   file.close() 
   if len(lines) == 0:
     raise RuntimeError("file '" + file_name + "' is empty")
   return lines[1]


# return the last line of a file
def last_line(file_name):
   file = open(file_name, "r")
   lines = file.readlines()
   file.close() 
   if len(lines) == 0:
     raise RuntimeError("file '" + file_name + "' is empty")
   return lines[len(lines)-1]



