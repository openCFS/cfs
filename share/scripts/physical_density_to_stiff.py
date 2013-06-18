#!/usr/bin/env python
from optimization_tools import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", help="density.xml file where pyhsical is read")
parser.add_argument("output", help="density.xml file where pyhsical is written as stiff1 and stiff2")
args = parser.parse_args()


d = read_density_as_vector(args.input, "physical")
data = numpy.zeros((len(d), 2))
data[:,0] = d
data[:,1] = d
write_multi_design_file(args.output, data, "stiff1", "stiff2")