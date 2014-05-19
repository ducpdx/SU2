#!/usr/bin/env python

## \file shape_optimization.py
#  \brief Python script for writing a list of Hicks-Henne bumps for SU2.
#  \author T. Economon, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.1
#
# Stanford University Unstructured (SU2) Code
# Copyright (C) 2012 Aerospace Design Laboratory
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os, time
from optparse import OptionParser
from numpy import *

parser = OptionParser()
parser.add_option("-i", "--iDegree", dest="iDegree", default=4,
                  help="i degree of the FFD box", metavar="IDEGREE")
parser.add_option("-j", "--jDegree", dest="jDegree", default=4,
                  help="j degree of the FFD box", metavar="JDEGREE")
parser.add_option("-k", "--kDegree", dest="kDegree", default=1,
                  help="k degree of the FFD box", metavar="KDEGREE")
parser.add_option("-b", "--ffdid", dest="ffd_id", default=0,
                  help="ID of the FFD box", metavar="FFD_ID")
parser.add_option("-m", "--marker", dest="marker",
                  help="marker name of the design surface", metavar="MARKER")
parser.add_option("-s", "--scale", dest="scale", default=1.0,
                  help="scale factor for the bump functions", metavar="SCALE")

(options, args)=parser.parse_args()

# Process options
options.iOrder  = int(options.iDegree) + 1
options.jOrder  = int(options.jDegree) + 1
options.kOrder  = int(options.kDegree) + 1
options.ffd_id  = str(options.ffd_id)
options.marker = str(options.marker)
options.scale  = float(options.scale)

print " "
print "FFD_CONTROL_POINT"

iVariable = 0
dvList = "DEFINITION_DV= "
for kIndex in range(options.kOrder):
  for jIndex in range(options.jOrder-1):
    for iIndex in range(options.iOrder):
      iVariable = iVariable + 1
      dvList = dvList + "( 7, " + str(options.scale) + " | " + options.marker + " | "
      dvList = dvList + options.ffd_id + ", " + str(iIndex) + ", " + str(jIndex+1) + ", " + str(kIndex) + ", 0.0, 0.0, 1.0 )"
      if iVariable < (options.iOrder*(options.jOrder-1)*options.kOrder):
        dvList = dvList + "; "


print dvList

print " "
print "FFD_CAMBER & FFD_THICKNESS"

iVariable = 0
dvList = "DEFINITION_DV= "
for jIndex in range(options.jOrder-1):
  for iIndex in range(options.iOrder):
    iVariable = iVariable + 1
    dvList = dvList + "( 11, " + str(options.scale) + " | " + options.marker + " | "
    dvList = dvList + options.ffd_id + ", " + str(iIndex) + ", " + str(jIndex+1) + " )"
    dvList = dvList + "; "
iVariable = 0
for jIndex in range(options.jOrder-1):
  for iIndex in range(options.iOrder):
    iVariable = iVariable + 1
    dvList = dvList + "( 12, " + str(options.scale) + " | " + options.marker + " | "
    dvList = dvList + options.ffd_id + ", " + str(iIndex) + ", " + str(jIndex+1) + " )"
    if iVariable < (options.iOrder*(options.jOrder-1)):
      dvList = dvList + "; "

print dvList




