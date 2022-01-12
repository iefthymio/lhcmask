#
# 

from cpymad.madx import Madx


md = Madx()

md.call(file='ForGuido_test_ibs.madx')

# -- verify the cpymad has correcly taken the beam
md.sequence['lhcb1'].beam

# -- now try to run a twiss inside cpymad

md.use('lhcb1')
md.twiss()