# -*- coding: utf-8 -*-
import locale
import string 
import os,sys
import subprocess

def do():
    print "sys.path[0]:   %s" % sys.path[0]
    print "os.getcwd():   %s" % os.getcwd()

class CFS:
    def __init__(self):
        self.executable = "cfsbin"
        self.root_dir = reduce (lambda l,r: l + os.path.sep + r, os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep )[:-2] )
        self.script_name = os.path.realpath( __file__ ).split( os.path.sep )[-1]
        self.script_bin_dir = os.path.join( self.root_dir, "bin" )
        self.scripts_dir = os.path.join( self.root_dir, "share/scripts" )
        self.distro_script = os.path.join( self.scripts_dir, "distro.sh" )
        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
        self.distro_out=subprocess.check_output([self.distro_script, "-p"])
        
        lib_options = {
            'I386'   : 'lib',
            'X86_64' : 'lib',
            'IA64'   : 'lib',
            }
        self.lib_suffix = lib_options[self.getArch()]
        self.exe_bin_dir = os.path.join( self.root_dir, "bin", self.getArchStr() )
        self.lib_dir = os.path.join( self.root_dir, self.lib_suffix, self.getArchStr() )

    def getArch(self):
        exec(self.distro_out)
        return ARCH

    def getArchStr(self):
        exec(self.distro_out)
        return str.join('_', [DIST, REV, ARCH])

    def getDist(self):
        exec(self.distro_out)
        return DIST

    def getRev(self):
        exec(self.distro_out)
        return REV


### EXAMPLE PYTHON MODULE
# Define some variables:
numberone = 1
ageofqueen = 78


# define some functions
def printhello():
    print "hello"
    
def timesfour(input):
    print input * 4
    
# define a class
class Piano:
    def __init__(self):
        self.type = raw_input("What type of piano? ")
        self.height = raw_input("What height (in feet)? ")
        self.price = raw_input("How much did it cost? ")
        self.age = raw_input("How old is it (in years)? ")
	
    def printdetails(self):
        print "This piano is a/an " + self.height + " foot",
        print self.type, "piano, " + self.age, "years old and costing\
 " + self.price + " dollars."
