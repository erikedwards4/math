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

from all.generate.constants.functions import *
from all.generate.other_gen.functions import *
from all.generate.rand.functions import *

from all.matmanip.construct.functions import *
from all.matmanip.matsel.functions import *
from all.matmanip.rearrange.functions import *
from all.matmanip.split_join.functions import *

from all.elementwise1.operators.functions import *
from all.elementwise1.trig.functions import *
from all.elementwise1.exp_log.functions import *
from all.elementwise1.round.functions import *
from all.elementwise1.special.functions import *
from all.elementwise1.nonlin.functions import *

from all.elementwise2.arithmetic.functions import *
from all.elementwise2.trig2.functions import *
from all.elementwise2.complex2.functions import *

from all.vec2scalar.sums.functions import *
from all.vec2scalar.prctiles.functions import *
from all.vec2scalar.iprctiles.functions import *
from all.vec2scalar.ranges.functions import *
from all.vec2scalar.norms.functions import *
from all.vec2scalar.moments.functions import *
from all.vec2scalar.other_means.functions import *
from all.vec2scalar.other_spreads.functions import *
from all.vec2scalar.other_stats.functions import *

from all.vecs2scalar.similarity.functions import *
from all.vecs2scalar.dist.functions import *

from all.vec2vec.center.functions import *
from all.vec2vec.scale.functions import *
from all.vec2vec.normalize.functions import *
from all.vec2vec.reorder.functions import *
from all.vec2vec.other_vec2vec.functions import *

from all.linalg.matmul.functions import *
from all.linalg.transform.functions import *
from all.linalg.sim_mat.functions import *
from all.linalg.other_linalg.functions import *
