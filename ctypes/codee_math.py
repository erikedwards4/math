#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in libmath.so.

Do:
import codee_math.py as cm

And then:
x1 = cm.ones((2, 3))
x2 = cm.randn((2, 3))
y = cm.plus(x1, x2)
etc.
"""

# The project is organized by the following directory structure:

# All: Generate Matmanip Elementwise1 Elementwise2 Vec2scalar Vec2vec Vecs2scalar Complex Linalg
# Generate: Constants Other_Gen Rand
# Matmanip: Construct Matsel Rearrange Split_Join
# Elementwise1: Operators Trig Exp_Log Round Special Complex Nonlin
# Elementwise2: Arithmetic Trig2 Complex2
# Vec2scalar: Sums Prctiles Iprctiles Ranges Norms Moments Other_Means Other_Spreads Other_Stats
# Vecs2scalar: Similarity Dist
# Vec2vec: Center Scale Normalize Reorder Other_Vec2vec
# Linalg: Matmul Transform Sim_Mat Dist_Mat Other_Linalg

# Within each terminal directory, the functions are found in functions.py,
# and the testing code within unittest.py, doctest.py, etc.

from all.generate.constants import *
from all.generate.other_gen import *
from all.generate.rand import *

from all.matmanip.construct import *
from all.matmanip.matsel import *
from all.matmanip.rearrange import *
from all.matmanip.split_join import *

from all.elementwise1.operators import *
from all.elementwise1.trig import *
from all.elementwise1.exp_log import *
from all.elementwise1.round import *
from all.elementwise1.special import *
from all.elementwise1.nonlin import *

from all.elementwise2.arithmetic import *
from all.elementwise2.trig2 import *
from all.elementwise2.complex2 import *

from all.vec2scalar.sums import *
from all.vec2scalar.prctiles import *
from all.vec2scalar.iprctiles import *
from all.vec2scalar.ranges import *
from all.vec2scalar.norms import *
from all.vec2scalar.moments import *
from all.vec2scalar.other_means import *
from all.vec2scalar.other_spreads import *
from all.vec2scalar.other_stats import *

from all.vecs2scalar.similarity import *
from all.vecs2scalar.dist import *

from all.vec2vec.center import *
from all.vec2vec.scale import *
from all.vec2vec.normalize import *
from all.vec2vec.reorder import *
from all.vec2vec.other_vec2vec import *

from all.linalg.matmul import *
from all.linalg.transform import *
from all.linalg.sim_mat import *
from all.linalg.other_linalg import *
