#Contains all of the building blocks for basic atoms building to do forcefield testing on

import mbuild as mb
import numpy as np
import warnings
from mbuild.lib.recipes import Polymer
from mbuild.lib.moieties import Silane

class Methyl(mb.Compound):
    """A methyl group with one port labeled 'down' """
    def __init__(self):
        super(Methyl, self).__init__()
        
        #Add CH3
        CH3=mb.lib.moieties.CH3()
        self.add(CH3)
        self.add(self.all_ports()[0],'up',containment=False)
        
class Ethyl(mb.Compound):
    """A methyl group with one port labeled 'down' """
    def __init__(self):
        super(Ethyl, self).__init__()
        
        #Add CH3
        CH3=mb.lib.moieties.CH3()
        self.add(CH3)
        CH3_2=mb.lib.moieties.CH3()
        CH3_2.translate([0.2,0,0])
        mb.force_overlap(CH3_2, CH3_2.all_ports()[0], self.all_ports()[0])
        self.add(CH3_2)
        #remove a hydrogen to add a port for connecting to the main carbon
        self.remove(self[1])
        self.add(self.all_ports()[0],'up',containment=False)
        
class Carboxylate(mb.Compound):
    """A carboxylate group with one extra port labeled ['C']['down'] """
    def __init__(self):
        super(Carboxylate, self).__init__()
        
        #Add central carbonyl with two extra ports
        self.add(mb.Particle(name='C'))
        
        self.add(mb.Port(anchor=self[0]), 'double')
        self['double'].translate(np.array([0, 0.062, 0]))

        self.add(mb.Port(anchor=self[0]), 'posterior')
        self['posterior'].translate(np.array([0, 0.07, 0]))
        self['posterior'].rotate(np.pi * 120/180, [0, 0, 1])
        
        self.add(mb.Port(anchor=self[0]), 'anterior')
        self['anterior'].translate(np.array([0, 0.07, 0]))
        self['anterior'].rotate(np.pi * 120/180, [0, 0, -1])
        
        self.add(mb.Particle(name='O'))
        self.add(mb.Port(anchor=self[1]), label='Odouble')
        self['Odouble'].translate(np.array([0, 0.062, 0]))
        

        mb.force_overlap(self[1], self['Odouble'], self['double'], add_bond=True)
        
        
        #Add posterior OH
        self.add(mb.Particle(name='O'))

        self.add(mb.Port(anchor=self[2]), label='Oanterior')
        self['Oanterior'].translate(np.array([0, .07, 0]))
        mb.force_overlap(self[2], self['Oanterior'], self['posterior'])
        
        self.add(mb.Port(anchor=self[2]), label='Oposterior')
        self['Oposterior'].translate(np.array([0, -.0485, 0]))
        
        self.add(mb.Particle(name='H'))
        self.add(mb.Port(anchor=self[3]), label='Hbond')
        self['Hbond'].translate(np.array([0, .0485, 0]))
        mb.force_overlap(self[3], self['Hbond'], self['Oposterior'])
        #self.add(self.all_ports()[0],'anterior',containment=False) 
        
        #Add CH2
        #CH2=mb.lib.moieties.CH2()
        #mb.force_overlap(CH2, CH2['up'], self['anterior'])
        #self.add(CH2)
        
class AlkylSilane(mb.Compound):
    """A silane functionalized alkane chain with two Ports. One is labeled ['down'] 
    and the other is ['alkane']['up']"""
    def __init__(self, chain_length):
        super(AlkylSilane, self).__init__()
        
        self.add(mb.Particle(name='H'))
        self.add(mb.Port(anchor=self[0]), 'oH')
        self['oH'].spin(np.pi, [0, 0, 1])
        self['oH'].translate(np.array([0, 0.097/2, 0]))
        
        self.add(mb.Particle(name='O'))
        self.add(mb.Port(anchor=self[1]), label='Oh')
        self['Oh'].spin(np.pi, [0, 0, 1])
        self['Oh'].translate(np.array([0, 0.097/2, 0]))
        mb.force_overlap(self[0], self['oH'], self['Oh'])
        self.add(mb.Port(anchor=self[1]), label='alcohol')
        self['alcohol'].spin(np.pi, [0, 0, 1])
        self['alcohol'].translate(np.array([0, -0.097/2, 0]))
        
        CH2=mb.lib.moieties.CH2()
        alkane = Polymer(CH2, chain_length, port_labels=['up','down'])
        silane = Silane()
        mb.force_overlap(silane, silane['down'],self['alcohol'])
        self.add(silane, 'silane')
        mb.force_overlap(alkane, alkane['down'], silane['up'])
        self.add(alkane, 'alkane')
        
        
        # Hoist silane port to AlkylSilane level.
        #self.add(silane['down'], 'down', containment=False)
        
        #remove a hydrogen to add a port for the second functionalgroup
        self.remove(self['alkane'][3*chain_length-1])
        self.add(self.all_ports()[0],'secondary',containment=False)
        self.add(self.all_ports()[1],'primary',containment=False)
        
class H(mb.Compound):
    """A hydrogen atom with two overlayed ports."""
    def __init__(self):
        super(H, self).__init__()
        self.add(mb.Particle(name='H'))

        self.add(mb.Port(anchor=self[0]), 'down')
        self['down'].spin(np.pi, [0, 0, 1])
        self['down'].translate(np.array([0, 0.07, 0]))
        
class H_Cap(mb.Compound):
    """A hydrogen atom with two overlayed ports."""
    def __init__(self,molecule):
        super(H_Cap, self).__init__()
        hydrogen=H()
        for port in molecule.all_ports():
            mb.force_overlap(hydrogen,hydrogen['down'], port)
            self.add(hydrogen)
            hydrogen=H()
        self.add(molecule)
        
class FunctionalizedAlkylSilane(mb.Compound):
    """A silane functionalized alkane chain with one Port labeled ['down']. """
    def __init__(self, chain_length, group1, group2,carbonyl_bool=False,chain_group='carbonyl'):
        super(FunctionalizedAlkylSilane, self).__init__()
        
        #Create functional group and main chain
        if group1 == 1:
            group_1=Carboxylate()
        if group2 == 1:
            group_2 = Methyl()
        if group2 == 2: 
            group_2=Ethyl()
        if group2 == 3: 
            group_2=Propyl()
        if group2 == 3.5: 
            group_2=IPA()
        if group2 == 4: 
            group_2=Butyl3() 
            
        alkylsilane=AlkylSilane(chain_length)
        
        #Check carbonyl_bool value to see if there should be an O double bond before the linker
        if carbonyl_bool:
            alkylsilane=ModifiedChain(chain_length,chain_group)
        
        #Overlap these at the top port, leaving the silane port available
        mb.force_overlap(group_1, group_1['anterior'], 
                         alkylsilane['primary'])
        self.add(alkylsilane)
        self.add(group_1,'Carboxylate')
        
        
        #Overlap the methyl group onto the secondary port of the chain
        mb.force_overlap(group_2, group_2['up'], 
                         alkylsilane['secondary'])
        self.add(group_2,'Methyl')
        
        #hoist alkylsilane port to the molecule level.
        #self.add(alkylsilane['down'], 'down', containment=False)
        
class ModifiedChain(mb.Compound):
    def __init__(self,chain_length,group):
        super(ModifiedChain,self).__init__()
        
        if group == 'ether':
            alkylsilane=AlkylSilane(chain_length)
            list(alkylsilane.particles())[-5].name = list(alkylsilane.particles())[-5].name.replace('C','O') 
            print('You have added an ether!')
            alkylsilane.remove(list(alkylsilane.particles())[-3])
            alkylsilane.remove(alkylsilane.all_ports()[0:1])
            alkylsilane.remove(list(alkylsilane.particles())[-3])
            alkylsilane.remove(alkylsilane.all_ports()[0:1])
        else:    
            alkylsilane=AlkylSilane(chain_length)
            alkylsilane.remove(list(alkylsilane.particles())[-6])
            alkylsilane.remove(alkylsilane.all_ports()[0:1])
            list(alkylsilane.particles())[-6].name = list(alkylsilane.particles())[-6].name.replace('H','O')
            if group == 'ester':
                list(alkylsilane.particles())[-5].name = list(alkylsilane.particles())[-5].name.replace('C','O') 
                print('You have added an ester!')
                alkylsilane.remove(list(alkylsilane.particles())[-3])
                alkylsilane.remove(alkylsilane.all_ports()[0:1])
                alkylsilane.remove(list(alkylsilane.particles())[-3])
                alkylsilane.remove(alkylsilane.all_ports()[0:1])
            elif group == 'amide':
                list(alkylsilane.particles())[-5].name = list(alkylsilane.particles())[-5].name.replace('C','N') 
                print('You have added an amide!')
                alkylsilane.remove(list(alkylsilane.particles())[-3])
                alkylsilane.remove(alkylsilane.all_ports()[0:1])
            else: 
                print('You have added a carbonyl!')
        
        self.add(alkylsilane)
        #hoist alkylsilane port to the molecule level.
        #self.add(alkylsilane['down'], 'down', containment=False)
        self.add(alkylsilane['primary'], 'primary', containment=False)
        self.add(alkylsilane['secondary'], 'secondary', containment=False)
        

        
   