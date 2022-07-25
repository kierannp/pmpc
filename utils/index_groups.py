from __future__ import division

import numpy as np


def generate_index_groups(system, freeze_thickness=0.5):
    bounding_box = system.boundingbox
    bot_of_box = bounding_box.mins[2]
    top_of_box = bounding_box.maxs[2]
    middle = (bot_of_box + top_of_box) / 2

    bottom_frozen = []
    top_frozen = []
    bottom_surface = []
    top_surface = []
    bottom_chains = []
    top_chains = []
    for i, particle in enumerate(system.particles()):
        z = particle.pos[2]
        ancestors = [ancestor.name for ancestor in particle.ancestors()]
        if z > middle:
            if 'Alkylsilane' in ancestors:
                top_chains.append(i + 1)
            else:
                top_surface.append(i + 1)
                if z > top_of_box - freeze_thickness:
                    top_frozen.append(i + 1)
        else:
            if 'Alkylsilane' in ancestors:
                bottom_chains.append(i + 1)
            else:
                bottom_surface.append(i + 1)
                if z < bot_of_box + freeze_thickness:
                    bottom_frozen.append(i + 1)

    bottom_frozen = np.asarray(bottom_frozen)
    print('bottom_frozen: {}'.format(len(bottom_frozen)))
    top_frozen = np.asarray(top_frozen)
    print('top_frozen: {}'.format(len(top_frozen)))

    bottom_surface = np.asarray(bottom_surface)
    print('bottom_surface: {}'.format(len(bottom_surface)))
    top_surface = np.asarray(top_surface)
    print('top_surface: {}'.format(len(top_surface)))
    surfaces = np.hstack((bottom_surface, top_surface))
    print('surfaces: {}'.format(len(surfaces)))

    bottom_chains = np.asarray(bottom_chains)
    print('bottom_chains: {}'.format(len(bottom_chains)))
    top_chains = np.asarray(top_chains)
    print('top_chains: {}'.format(len(top_chains)))
    chains = np.hstack((bottom_chains, top_chains))
    print('chains: {}'.format(len(chains)))

    bottom = np.hstack((bottom_surface, bottom_chains))
    print('bottom: {}'.format(len(bottom)))
    top = np.hstack((top_surface, top_chains))
    print('top: {}'.format(len(top)))
    system = np.hstack((bottom, top))
    print('System: {}'.format(len(system)))

    index_groups = {'System': system,
                    'bottom': bottom,
                    'top': top,
                    'bottom_frozen': bottom_frozen,
                    'top_frozen': top_frozen,
                    'surfaces': surfaces,
                    'bottom_surface': bottom_surface,
                    'top_surface': top_surface,
                    'chains': chains,
                    'bottom_chains': bottom_chains,
                    'top_chains': top_chains}

    return index_groups
