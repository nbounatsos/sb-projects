#!/usr/bin/env python

from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='alignment.ali',
              knowns='2LRU', sequence='5G3Q',
              assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()