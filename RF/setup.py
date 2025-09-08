#!/usr/bin/env python

###########################################################################
##    Copyright 2010 Rahul Suri and Tandy Warnow.
##    This file is part of spruce.
##
##    spruce is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    spruce is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with spruce.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################

from distutils.core import setup

setup(name = "spruce",
      version = "1.0",
      description = "A tidy package for supertree analysis.",
      packages = ["spruce"],

      url = "http://www.cs.utexas.edu/~phylo/software/spruce",
      author = "Rahul Suri",
      author_email = "rsuri@cs.utexas.edu",

      license = "General Public License (GPL)",
      requires = ["newick_modified (>= 1.3.1)"],
      provides = ["spruce"],

      scripts = ["spruce/scripts/checkIfPlenary.py",
                 "spruce/scripts/getFpFn.py",
                 "spruce/scripts/getGreedyConsensus.py",
                 "spruce/scripts/getNumSharedBipartitions.py",
                 "spruce/scripts/getResolution.py",
                 "spruce/scripts/getSumOfFpFnRf.py",
                 "spruce/scripts/makeRatchetFile.py",
                 "spruce/scripts/printPolytomyInfo.py"],

      classifiers = ["Environment :: Console",
                     "Intended Audience :: Developers",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: GNU General Public License (GPL)",
                     "Natural Language :: English",
                     "Operating System :: OS Independent",
                     "Programming Language :: Python",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"])
