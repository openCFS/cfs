#! /usr/bin/python 

""" converts selected testsuite examples to MediaWiki format and writes a tar-archive of the files """

testsuitepath = '../TESTSUIT'

# tests to use
wikitests = [\
  'Singlefield/Acoustics/Abc2d',\
  'Singlefield/Acoustics/Pml3dTrans',\
  'Singlefield/Acoustics/ConsInterp3d',\
  'Singlefield/Acoustics/Mortar3d',\
  'Singlefield/Acoustics/Nitsche2d',\
  'Singlefield/Acoustics/Periodic2d',\
  'Singlefield/AcousticMixed/channel3d_sfem_abc',\
  'Singlefield/WaterWaves/PML2d',
  'Singlefield/Electrostatics/CylCapAxi',\
  'Singlefield/Electrostatics/LShapedDomain',\
  'Singlefield/Electrostatics/CylCap3D',\
  'Singlefield/Electrostatics/Cube2dMortar',\
  'Singlefield/ElecConduction/OhmsLaw',\
  'Singlefield/Heat/SlabVolSource',\
  'Singlefield/Heat/HeatNLcondAndCap3d',\
  'Singlefield/Heat/HeatMortar',\
  'Singlefield/Flow/Channel2D',\
  'Singlefield/Flow/Channel2dMeanFlowPerturbed',\
  'Singlefield/Flow/LidDrivenCavity',\
  'Singlefield/Magnetics/Coil3DEdgeNL',\
  'Singlefield/Magnetics/CylindricCoilEdge',\
  'Singlefield/Magnetics/Coil3DVoltEdge',\
  'Singlefield/Magnetics/BlockConstFlux',\
  'Singlefield/Magnetics/ThinYoke',\
  'Singlefield/Mechanics/Plate3D',\
  'Singlefield/Mechanics/PreStress2dMS_EV',\
  'Singlefield/Mechanics/PreStress2dMS_Harm',\
  'Singlefield/Mechanics/MassCMEHarm',\
  'Singlefield/Mechanics/StiffCMEHarm',\
  'Singlefield/Mechanics/Torque3dNitsche',\
  'Singlefield/Mechanics/MechThermStatic',\
  'Coupledfield/Piezo/ShearActuatorStatic2D',\
  'Coupledfield/Piezo/ExtensionActuatorEV2D',\
  'Coupledfield/Piezo/PiezoTrans2D',\
  'Coupledfield/MechAcou/RingHarm2d',\
  'Coupledfield/MechAcou/RingHarm2dRestart',\
  'Coupledfield/FlowMech/OsciPipe2d',\
  'Coupledfield/MagMech/SolenoidAxi',\
  'Coupledfield/MagMech/PermaMS',\
  'Coupledfield/MagMechAcou/PlateAboveCoilFluid2D',\
  'Coupledfield/ElecTherm/ElTh_Static'\
]

# read the data for the defined tests
from os import path
from lxml import etree
data = dict()
from sys import stdout as f
for test in wikitests: 
    names = test.split('/')
    if names[0] not in data.keys() : # insert top level name 
        data[names[0]] = dict()
    if names[1] not in data[names[0]].keys() : # insert mid level name
        data[names[0]][names[1]] = dict ()

    xmlfile = path.join(*[testsuitepath,test,path.split(test)[-1]+'.xml'])
    tree = etree.parse(xmlfile)
    root = tree.getroot()
    
    txtdict = {'title':'NO TITLE','description':'NO DESCRIPTION'}
    for tag in txtdict.keys() :
        txt = root.xpath('w:documentation/w:%s/text()'%tag,\
        		namespaces={'w':'http://www.cfs++.org'})
        if len(txt) == 1 : # we found a text for the tag
            txt = txt[0]
            lines = txt.split('\n')
            txtdict[tag] = ('\n'.join([l.strip() for l in lines])).strip()
    data[names[0]][names[1]][names[2]] = txtdict

#  format the data
import sys
def writeHeadingWiki(txt,lvl=1,fid=sys.stdout) :
    """ wites a heading """
    fid.write("="*lvl+' %s '%txt+"="*lvl+'\n')

ni = 2 # initial heading level
for h1 in data.keys() :
    writeHeadingWiki(h1,ni,f)
    for h2 in data[h1].keys() :
        writeHeadingWiki(h2,ni+1,f)
        for h3 in data[h1][h2].keys() :
            writeHeadingWiki(h3,ni+2,f)
            info = data[h1][h2][h3]
            f.write("''%s''\n"%(info['title']))
            f.write("\n")
            f.write(info['description']+'\n')
            f.write("\n")

# make tar file
#print('writing archive of selected examples ...')
import tarfile
tar = tarfile.open('testsuite_examples.tar.gz','w|gz')
for name in wikitests :
    tar.add(path.join(testsuitepath,name),arcname=name,recursive=True,\
		    exclude=lambda x: 'CMakeLists.txt' in x)
tar.close()
