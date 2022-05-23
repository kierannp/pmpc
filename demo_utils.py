import mbuild as mb
import gmso
import pandas as pd
import unyt 
import networkx as nx
import numpy as np

import molecule_parts

class AceticAcid(mb.Compound):
    def __init__(self):
        super(AceticAcid, self).__init__()
        CH3 = mb.lib.moieties.CH3()
        carb = molecule_parts.Carboxylate()
        carb.visualize(show_ports=True)
        self.add(carb)
        mb.force_overlap(CH3,CH3['up'],carb['anterior'],CH3)
        self.add(CH3)
        
class Phenylcarbacid(mb.Compound):
    def __init__(self):
        super(Phenylcarbacid, self).__init__()
        
        C3 = mb.lib.atoms.C3()
        self.add(C3,'C0')
        for i in np.arange(5):
            C3 =  mb.lib.atoms.C3()
            mb.force_overlap(C3,C3['down'],self['C' + str(i)]['left'])
            H = mb.lib.atoms.H()
            mb.force_overlap(H,H['up'],C3['up'])
            self.add(C3,'C'+str(i+1))
            self.add(H,'H'+str(i+1))
            
        mb.force_overlap(self['C5'],self['C5']['left'],self['C0']['down'])
        
        CH2 = mb.lib.moieties.CH2()
        mb.force_overlap(CH2,CH2['down'],self['C0']['up'])
        self.add(CH2, "CH21")
        
        CH2 = mb.lib.moieties.CH2()
        mb.force_overlap(CH2,CH2['down'],self['CH21']['up'])
        #CH2.remove(CH2[1])
        self.add(CH2, "CH22")
        
        carb = molecule_parts.Carboxylate()
        mb.force_overlap(carb,carb['anterior'],self['CH22']['up'])
        self.add(carb,"carb")
        
        self.energy_minimize()
        
def atomtypes_to_datatables(graph, labels=None, atom_objects = False):
    if not labels:
        labels = []
    df = pd.DataFrame()
    df['index'] = np.arange(0,len(graph.nodes),1)
    df['atom_type'] = list(node.atom_type.name for node in graph.nodes)
    df['names'] = list(node.name for node in graph.nodes)
    df['charge'] = list((node.charge/unyt.electron_charge).round(4) for node in graph.nodes)
    for label in labels:
        if label == 'position':
            df['x'] = list(getattr(node,label)[0] for node in graph.nodes)
            df['y'] = list(getattr(node,label)[1] for node in graph.nodes)
            df['z'] = list(getattr(node,label)[2] for node in graph.nodes)
        else:
            df[label] = list(str(getattr(node,label)) for node in graph.nodes)
    if atom_objects:
        df['atom_id'] = list(node.__hash__() for node in graph.nodes)
    return df

def bondtypes_to_datatables(graph,topology,labels=None,atom_objects = False):
    df = pd.DataFrame()
    if not labels:
        labels = []
    list_of_nodes = list(node for node in graph.nodes)

    df['index'] = np.arange(0,len(graph.edges),1)
    df['Atom1'] = list(str(edge[0].name) + '(' + str(list_of_nodes.index(edge[0])) + ')' for edge in graph.edges)
    df['Atom2'] = list(str(edge[1].name) + '(' + str(list_of_nodes.index(edge[1])) + ')' for edge in graph.edges)
    df['Parameter 1 (k): ' + str(topology.bonds[0].bond_type.parameters['k'].units)] = (
        list(bond.bond_type.parameters['k'].round(3) for bond in topology.bonds))
    df['Parameter 2 (r_eq): ' + str(topology.bonds[0].bond_type.parameters['r_eq'].units)] = (
       list(bond.bond_type.parameters['r_eq'].round(3) for bond in topology.bonds))
    for label in labels:
        df[label] = list(getattr(node,label) for node in graph.nodes)
    if atom_objects:
        df['atom1_id'] = list(edge[0].__hash__() for edge in graph.edges)
        df['atom2_id'] = list(edge[1].__hash__() for edge in graph.edges)
    return df

def angletypes_to_datatables(graph,topology,labels=None,atom_objects=False):
    df = pd.DataFrame()
    if not labels:
        labels = []
    list_of_nodes = list(node for node in graph.nodes)

    df['index'] = np.arange(0,len(topology.angles),1)
    df['Atom1'] = list(str(angle.connection_members[0].name) +
                       '(' + 
                       str(list_of_nodes.index(angle.connection_members[0])) +
                       ')' for angle in topology.angles)
    df['Atom2'] = list(str(angle.connection_members[1].name) +
                       '(' + 
                       str(list_of_nodes.index(angle.connection_members[1])) +
                       ')' for angle in topology.angles)
    df['Atom3'] = list(str(angle.connection_members[2].name) +
                       '(' + 
                       str(list_of_nodes.index(angle.connection_members[2])) +
                       ')' for angle in topology.angles)
    df['Parameter 1 (k): ' + str(topology.angles[0].angle_type.parameters['k'].units)] = (
        list(angle.angle_type.parameters['k'].round(3) for angle in topology.angles))
    df['Parameter 2 (theta_eq): ' + str(topology.angles[0].angle_type.parameters['theta_eq'].units)] = (
       list(angle.angle_type.parameters['theta_eq'].round(3) for angle in topology.angles))
    for label in labels:
        df[label] = list(getattr(node,label) for node in graph.nodes)
    if atom_objects:
        df['atom1_id'] = list(angle.connection_members[0].__hash__() for angle in topology.angles)
        df['atom2_id'] = list(angle.connection_members[1].__hash__() for angle in topology.angles)
        df['atom3_id'] = list(angle.connection_members[2].__hash__() for angle in topology.angles)
    return df

def dihedraltypes_to_datatables(graph,topology,labels=None,atom_objects=False):
    df = pd.DataFrame()
    if not labels:
        labels = []
    list_of_nodes = list(node for node in graph.nodes)

    df['index'] = np.arange(0,len(topology.dihedrals),1)
    df['Atom1'] = list(str(dihedral.connection_members[0].name) +
                       '(' + 
                       str(list_of_nodes.index(dihedral.connection_members[0])) +
                       ')' for dihedral in topology.dihedrals)
    df['Atom2'] = list(str(dihedral.connection_members[1].name) +
                       '(' + 
                       str(list_of_nodes.index(dihedral.connection_members[1])) +
                       ')' for dihedral in topology.dihedrals)
    df['Atom3'] = list(str(dihedral.connection_members[2].name) +
                       '(' + 
                       str(list_of_nodes.index(dihedral.connection_members[2])) +
                       ')' for dihedral in topology.dihedrals)
    df['Atom4'] = list(str(dihedral.connection_members[3].name) +
                       '(' + 
                       str(list_of_nodes.index(dihedral.connection_members[3])) +
                       ')' for dihedral in topology.dihedrals)
    df['Parameter 1 (c0): ' + str(topology.dihedrals[0].dihedral_type.parameters['c0'].units)] = (
        list(dihedral.dihedral_type.parameters['c0'].round(3) for dihedral in topology.dihedrals))
    df['Parameter 2 (c1): ' + str(topology.dihedrals[0].dihedral_type.parameters['c1'].units)] = (
        list(dihedral.dihedral_type.parameters['c1'].round(3) for dihedral in topology.dihedrals))
    df['Parameter 3 (c2): ' + str(topology.dihedrals[0].dihedral_type.parameters['c2'].units)] = (
        list(dihedral.dihedral_type.parameters['c2'].round(3) for dihedral in topology.dihedrals))
    df['Parameter 4 (c3): ' + str(topology.dihedrals[0].dihedral_type.parameters['c3'].units)] = (
        list(dihedral.dihedral_type.parameters['c3'].round(3) for dihedral in topology.dihedrals))
    df['Parameter 5 (c4): ' + str(topology.dihedrals[0].dihedral_type.parameters['c4'].units)] = (
        list(dihedral.dihedral_type.parameters['c4'].round(3) for dihedral in topology.dihedrals))
    df['Parameter 6 (c5): ' + str(topology.dihedrals[0].dihedral_type.parameters['c5'].units)] = (
        list(dihedral.dihedral_type.parameters['c5'].round(3) for dihedral in topology.dihedrals))
    for label in labels:
        df[label] = list(getattr(node,label) for node in graph.nodes)
    if atom_objects:
        df['atom1_id'] = list(dihedral.connection_members[0].__hash__() for dihedral in topology.dihedrals)
        df['atom2_id'] = list(dihedral.connection_members[1].__hash__() for dihedral in topology.dihedrals)
        df['atom3_id'] = list(dihedral.connection_members[2].__hash__() for dihedral in topology.dihedrals)
        df['atom4_id'] = list(dihedral.connection_members[3].__hash__() for dihedral in topology.dihedrals)
    return df