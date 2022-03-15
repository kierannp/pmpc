from numpy import pi

import mbuild as mb

from mbuild.lib.moieties.silane import Silane
from mbuild.lib.moieties.ch3 import CH3
from mpc import MPC
from initiator import Initiator


class Brush(mb.Compound):
    """ """
    def __init__(self, chain_length=4, alpha=pi/4):
        super(Brush, self).__init__()

        # Add parts
        self.add(Silane(), label='silane')
        self.add(Initiator(), label='initiator')
        chain = mb.recipes.Polymer()
        chain.add_monomer(MPC(alpha=alpha),indices=[1, 37], separation = .2, replace=False, orientation=[[0,1,0],[0,-1,0]])
        chain.build(n = chain_length)
        self.add(chain, label='pmpc')
        # self.add(CH3(), label='methyl')

        self.add(mb.Port(anchor=list(self['pmpc'].particles())[1]), label='pmpc_up')
        self.add(mb.Port(anchor=list(self['pmpc'].particles())[37]), label='pmpc_down')

        # print(self['methyl'])
        mb.force_overlap(self['initiator'], self['initiator']['down'], self['silane']['up'])
        mb.force_overlap(self['pmpc'], self['pmpc_up'], self['initiator']['up'])
        # mb.force_overlap(self['methyl'], self['methyl']['up'], self['pmpc_down'])

        # Make self.port point to silane.bottom_port
        self.add(self['silane']['down'], label='down', containment=False)

if __name__ == "__main__":
    pmpc = Brush()
    print(pmpc)
