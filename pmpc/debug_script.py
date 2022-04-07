from numpy import pi

import mbuild as mb

from mbuild.lib.atoms import H
from mbuild.lib.surfaces import Betacristobalite
from brush import Brush


class PMPCLayer(mb.lib.recipes.Monolayer):
    """Create a layer of grafted pMPC brushes on a beta-cristobalite surface."""
    def __init__(self, pattern, tile_x=1, tile_y=1, chain_length=4, alpha=pi / 4):
        surface = Betacristobalite()
        brush = Brush(chain_length=chain_length, alpha=alpha)
        hydrogen = H()
        super(PMPCLayer, self).__init__(surface, brush, backfill=hydrogen,
                                        pattern=pattern, tile_x=tile_x,
                                        tile_y=tile_y)

pattern = mb.Random2DPattern(10)
pmpc_layer = PMPCLayer(pattern=pattern, chain_length=3, alpha=pi / 4, tile_x=1, tile_y=1)
# pmpc_layer.visualize(show_ports=True)

pmpc_layer.save('pmpc_layer.top', forcefield_files='../utils/oplsaa_imodels.xml', overwrite=True)
pmpc_layer.save('pmpc_layer.lammps', forcefield_files='../utils/oplsaa_imodels.xml', overwrite=True)
# !cat pmpc_layer.top