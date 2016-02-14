#!/usr/bin/env python
# ====================================================================== #
# fix_split_filenames.py
# ====================================================================== #
# There was a bug with the first version of split_data which 
# lead to files being named ending .SAC.SAC where the stream ID was 
# blank (and not putting the .. in the filename after the station 
# to signify this).
#
# This script renames files output by this buggy version of 
# split_data to the correct naming convention.
#
# Jamie Barron
# February 2011
# ====================================================================== #

import os

fixlist = [ f for f in os.listdir(os.getcwd()) if f.find("SAC.SAC") != -1 ]

for filename in fixlist:
    filesplit = filename.split('.')
    first_part = '.'.join(filesplit[0:8])
    second_part = '.'.join(filesplit[8:11])
    newname = "%s..%s" % (first_part, second_part)
    os.rename(filename,newname)
