#!/usr/bin/env python 

## \file parse_config.py
#  \brief Builds a worksheet of all SU2.cpp options
#  \author Aniket Aranake, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.9
#
# Stanford University Unstructured (SU2).
# Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import os,sys,xlwt

# note: requires xlwt for spreadsheet output
# http://pypi.python.org/pypi/xlwt

class config_option:
  name             = ""
  type             = ""
  category         = ""
  available_values = []
  default          = ""
  description      = ""

  def __init__(self,name_in,type_in,category_in,available_values_in,default_in,description_in):
    self.name                  = name_in
    self.type                  = type_in
    self.category              = category_in
    self.available_values      = available_values_in
    self.default               = default_in
    self.description           = description_in

  def print_data(self):
    print 'Option Name: %s '%        self.name
    print 'Option Type: %s '%        self.type
    print 'Option Category: %s '%    self.category
    print 'Option values: ',         self.available_values
    print 'Option default: %s'%      self.default
    print 'Option description: %s '% self.description
    print ''

def parse_config(config_cpp, config_hpp):

  # List of option types
  option_types = ['AddEnumOption', 'AddMathProblem', 'AddSpecialOption', 'AddScalarOption', 'AddMarkerOption', 'AddMarkerPeriodic', 'AddMarkerDirichlet', 'AddMarkerInlet', 'AddMarkerOutlet', 'AddMarkerDisplacement', 'AddMarkerLoad', 'AddMarkerFlowLoad', 'AddArrayOption', 'AddListOption', 'AddConvectOption', 'AddEnumListOption', 'AddDVParamOption']

  # List of parser tags
  parser_tags = ['CONFIG_CATEGORY','DESCRIPTION','AVAILABLE_VALUES','DEFAULT']

  # Parse the hpp file
  # The following code builds a dictionary of enumerated options
  # By default, an AddEnumOption type will consider the options list in the CCreateMap as available
  # This behavior is overwritten for options which include an AVAILABLE_OPTIONS tag in option_structure.cpp
  enum_options = {}
  f = open(config_hpp,'r')
  while(1):  # Find beginning of enum definitions
    s = f.readline()
    if s.find('BEGIN_CONFIG_ENUMS') >-1:
      break
  while(1):
    s = f.readline()
    if s.find('END_CONFIG_ENUMS')>-1:
      break   # Reached end

    if s.find('CCreateMap')>-1:
      dict_key = (s.split('=')[0]).split('>')[1].strip()
      dict_val = []
      while(1):
        s2 = f.readline()
        thisval = s2.split('"')[1]
        dict_val.append(thisval)
        if s2.find(';')>-1:
          break;
      enum_options[dict_key] = dict_val
  f.close()

  # Temporary: For now, build a list of all of the schemes
  # 8/22/2013 - This is going to change in the main code!
  scheme_list   = enum_options['Upwind_Map']
  scheme_list.extend(enum_options['Centered_Map'][1:])

  # Read the Options section of config_structure.cpp into a list of strings
  # No parsing is done yet; this just isolates the code between the BEGIN_CONFIG_OPTIONS and END_CONFIG_OPTIONS tags
  lines = []
  f = open(config_cpp,'r')
  while(1): 
    s = f.readline()
    if s.find('BEGIN_CONFIG_OPTIONS')>-1:
      break
  while(1):
    s = f.readline()
    # Check if we've reached the end
    if s.find('END_CONFIG_OPTIONS')>-1:
      break
    lines.append(s)
  f.close()

  # The primary list which will be output by this routine
  # Members of this list are of type config_option, which is defined above in this file
  option_list = []

  #----- Main text parsing loop -----
  present_category = "None"
  default_override = ""
  descripion_override = ""
  available_values_override = []
  for j,line in enumerate(lines):

    # Check for a category description
    if line.find('CONFIG_CATEGORY')>-1:
      present_category = line.split(':')[1].strip().strip('*/').strip()

    # Check for a default value tag, which would override the default as we would search for it below
    if line.find('DEFAULT:')>-1:
      default_override = line.split(":")[1].strip()
      default_override = default_override.strip('/*')
    elif line.find('default_vec_3d[0]')>-1:
      pieces = line.strip().strip(';').split(';')
      default_override = "(" + ', '.join([piece.split('=')[1] for piece in pieces])+")"

    # Check for a description tag, which would override the default description of "None"
    # This description may span several lines
    # The description ends when the comment ends or another tag is found
    if line.find('DESCRIPTION')>-1:
      description_override = line.split(":")[1].strip()
      j2 = j
      if line.find('*/')==-1:
        while(True):
          j2+=1

          found_tag = False
          for tag in parser_tags:
            if (lines[j2].find(tag)>-1):   # Found another tag
              found_tag = True
              break
          if found_tag:
            break

          if (lines[j2].find('*/')>-1):    # Found end of comment
            description_override += "\n" + lines[j2].split('*/')[0].strip('*').strip()
            break
          else:                            # Found nothing; the next line is also part of the description
            description_override += "\n" + lines[j2].strip().strip('*')
      description_override = description_override.strip().strip('/*')
          
    # Check for an available_values override
    if line.find('AVAILABLE_VALUES')>-1:
      available_values_override = line.split(':')[1].strip('/*').strip().split(',')

    # Check for an option type
    for option_type in option_types:
      if line.find(option_type)>-1:  # Found an option

        # Get option name. This is the first parameter of the Add____ function, and is always surrounded by double quotes (")
        name = line.split('"')[1]

        # Decide what the available values are. As a default this is set to Yes or No, but depending on the option type
        # it may be necessary to draw from an enum map or parse the values from comments in preceeding lines
        if available_values_override:
          available_values = [entry.strip() for entry in available_values_override]
          available_values_override = []
        else:
          available_values = ['YES','NO']
          if option_type=='AddEnumOption':
            try:
              enum_mapname = line.split(',')[2].strip()
              available_values = enum_options[enum_mapname]
            except KeyError:
              print "KeyError, key=%s"%enum_mapname 
              print "enum_options: ",enum_options
              sys.exit(1)
            except TypeError:
              print "TypeError, key=%s"%enum_mapname 
              print "enum_options: ",enum_options
              sys.exit(1)
          elif option_type=='AddMathProblem':
            available_values = ['DIRECT','ADJOINT','LINEARIZED']
          elif option_type=='AddScalarOption':
            available_values = ['A scalar constant']
          elif option_type in ('AddMarkerOption', 'AddMarkerPeriodic', 'AddMarkerDirichlet', 'AddMarkerInlet', 'AddMarkerOutlet', 'AddMarkerDisplacement', 'AddMarkerLoad', 'AddMarkerFlowLoad'):
            available_values = ['Valid marker name from grid file']
          elif option_type == 'AddArrayOption':
            available_values = ['Array']
          elif option_type == 'AddListOption':
            available_values = ['List']
          elif option_type == 'AddConvectOption':
            available_values = scheme_list
          elif option_type == 'AddEnumListOption':
            available_values = ['Enum list'] 
          elif option_type == 'AddDVParamOption':
            available_values = ['DV Param']

        # If a default tag is specified, use it
        if not default_override=="":
          default = default_override
        
          # Clear default_override before continuing to parse
          default_override = ""

        # Otherwise, search for the default value
        else:

          # Check the last item in parenthesis, accounting for the fact that the function call
          # might span multiple lines. 
          # This will fail if the semicolon (;) is not on the same line as the last param!
          jdefault = j
          while(lines[jdefault].find(';')==-1):
            jdefault = jdefault+1
          default = lines[jdefault].strip().strip(');').split(',')[-1].strip().strip('"')

          # Correction for boolean options
          # true/false in c++ are replaced by YES and NO in config file
          if default=='false':
            default='NO'
          elif default=='true':
            default='YES'

          # Correction for filenames, which are sloppily parsed at this stage...
          if option_type=='AddScalarOption' and default.find('string')>-1:
            default = default.split('"')[1]

        # If a description was found, use it
        if not description_override=="":
          description = description_override
          description_override = ""
        else:
          description = "No description"

        # Add a new option to the list
        option_list.append(config_option(name,option_type[3:],present_category,available_values,default,description))

        break
  return option_list

def print_all(option_list):
  # Dumps the option list to screen
  for option in option_list:
    option.print_data()

def make_spreadsheet(filename, option_list):

  # note: package xlwt is required for spreadsheet output
  # http://pypi.python.org/pypi/xlwt
  # On Debian/Ubuntu: $ sudo apt-get install python-xlwt

  import xlwt
  wbk   = xlwt.Workbook()

  sheet_name = ""
  jp = 0
  for j,opt in enumerate(option_list):
    jp = jp+1
    if not sheet_name==opt.category:
      # Create new sheet for new category
      sheet_name = opt.category
      sheet = wbk.add_sheet(sheet_name[:31].replace('/','-'))

      # Write spreadsheet header
      sheet.write(0,0,'Option Name')
      sheet.write(0,1,'Option Type')
      sheet.write(0,2,'Option Category')
      sheet.write(0,3,'Option Values')
      sheet.write(0,4,'Option Default')
      sheet.write(0,5,'Option Description')
      jp = 1


    sheet.write(jp,0,opt.name)
    sheet.write(jp,1,opt.type)
    sheet.write(jp,2,opt.category)
    sheet.write(jp,3,(',').join(opt.available_values))
    sheet.write(jp,4,opt.default)
    sheet.write(jp,5,opt.description)

  wbk.save(filename)

def write_config_template(option_list, outfile):
  '''Routine to write out a template configuration file using descriptions gathered from C++ code. The parameter outfile is a filename. If the file already exists it will be overwritten.'''

  f = open(outfile,'w')

  # Start with a comment header
  f.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
  f.write('%                                                                              %\n')
  f.write('% Stanford University Unstructured (SU2) configuration file                    %\n')
  f.write('% Case description: _________________________________________________________  %\n')
  f.write('% Author: ___________________________________________________________________  %\n')
  f.write('% Institution: ______________________________________________________________  %\n')
  f.write('% Date: __________                                                             %\n')
  f.write('% File Version 2.0.6                                                           %\n')
  f.write('%                                                                              %\n')
  f.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')

  max_line_length = 80

  present_category = ""
  # Loop through options and write each to file
  for opt in option_list:

    if not opt.category==present_category:
      # Time for a new category
      present_category = opt.category

      # Figure out how many dashes to pad the line with to achieve max_line_length characters total
      total_dashes  = max_line_length - 5 - len(present_category)
      dashes_ahead  = total_dashes/2
      dashes_behind = total_dashes/2 + total_dashes%2
  
      # Assemble category string
      category_string = '\n% '
      for j in range(dashes_ahead):
        category_string += '-'
      category_string += ' '+present_category+' '
      for j in range(dashes_behind):
        category_string += '-'
      category_string += '%\n'
      f.write(category_string)

    # Write the lines for this option
    description_string = opt.description.replace('\n','\n% ')

    opt_string  = '%\n% ' + description_string + ' (' + ', '.join(opt.available_values) + ')\n'
    opt_string += opt.name + '= ' + opt.default + '\n'
    f.write(opt_string)

if __name__=="__main__":

  # These variables should point to the configuration files
  su2_home = os.environ['SU2_HOME']
  config_cpp = os.path.join(su2_home,'Common/src/config_structure.cpp')
  config_hpp = os.path.join(su2_home,'Common/include/option_structure.hpp')

  # Check that files exist
  if not (os.path.isfile(config_cpp) and os.path.isfile(config_hpp)):
    sys.exit('ERROR: Could not find source file, please check that $SU2_HOME environment variable is set correctly\n  (Currently set to: %s)'%su2_home)
 
  # Run the parser
  option_list = parse_config(config_cpp, config_hpp)

  # Dump parsed data to screen
  #print_all(option_list)

  # make_spreadsheet('out.xls',option_list)

  outfile = 'test.cfg'
  write_config_template(option_list, outfile)
