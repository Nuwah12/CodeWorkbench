#### Object for handling bonds
import numpy as np
class bondUpdater(object):

    def __init__(self, LEFpositions):
        """
        Initialize a bondUpdater object

        :param LEFpositions: numpy array of extruder positions wrt polymer position
        """
        self.LEFpositions = LEFpositions
        self.curtime  = 0
        self.allBonds = []

    def setParams(self, activeParamDict, inactiveParamDict):
        """
        A method to set parameters for bonds.
        It is a separate method because you may want to have a Simulation object already existing

        :param activeParamDict: a dict (argument:value) of addBond arguments for active bonds
        :param inactiveParamDict:  a dict (argument:value) of addBond arguments for inactive bonds

        """
        self.activeParamDict = activeParamDict
        self.inactiveParamDict = inactiveParamDict


    def setup(self, bondForce,  blocks=100, smcStepsPerBlock=1):
        """
        A method that milks smcTranslocator object
        and creates a set of unique bonds, etc.

        :param bondForce: a bondforce object (new after simulation restart!)
        :param blocks: number of blocks to precalculate
        :param smcStepsPerBlock: number of smcTranslocator steps per block
        :return:
        """


        if len(self.allBonds) != 0:
            raise ValueError("Not all bonds were used; {0} sets left".format(len(self.allBonds)))

        self.bondForce = bondForce # force_dict from simulation object (bondForce obj)

        #precalculating all bonds
        allBonds = []
        
        loaded_positions  = self.LEFpositions[self.curtime : self.curtime+blocks] # Get all extruder positions from curtime to curtime+blocks
        allBonds = [[(int(loaded_positions[i, j, 0]), int(loaded_positions[i, j, 1])) 
                        for j in range(loaded_positions.shape[1])] for i in range(blocks)] # Get all positions for both legs of extruder, i.e. location of the 'bonds'
        self.allBonds = allBonds
        #print("bonds set to {}".format(self.allBonds))
        self.uniqueBonds = list(set(sum(allBonds, []))) # unlist the bonds and get unique
        #print("Unique bonds: {}".format(self.uniqueBonds))
        #adding forces and getting bond indices
        self.bondInds = []
        self.curBonds = allBonds.pop(0) 
        #print("Cur Bonds: {}".format(self.curBonds))
        for bond in self.uniqueBonds: # Loop thru all bonds, if bonds are in curBonds (first positions for legs, since this is the setup func)
            paramset = self.activeParamDict if (bond in self.curBonds) else self.inactiveParamDict # Determine if bond is active and get parameters
            ind = bondForce.addBond(bond[0], bond[1], **paramset) # This is where we add the bond to the actual bondForce object 
            self.bondInds.append(ind)
        self.bondToInd = {i:j for i,j in zip(self.uniqueBonds, self.bondInds)} # Dict of {bond : bond index}
        
        self.curtime += blocks # Advance blocks
        
        return self.curBonds,[]


    def step(self, context, verbose=True):
        """
        Update the bonds to the next step.
        It sets bonds for you automatically!
        :param context:  context
        :return: (current bonds, previous step bonds); just for reference
        """
        if len(self.allBonds) == 0:
            raise ValueError("No bonds left to run; you should restart simulation and run setup  again")

        pastBonds = self.curBonds
        self.curBonds = self.allBonds.pop(0)  # getting current bonds
        bondsRemove = [i for i in pastBonds if i not in self.curBonds] # ID bonds to remove
        bondsAdd = [i for i in self.curBonds if i not in pastBonds] # Bonds to add
        bondsStay = [i for i in pastBonds if i in self.curBonds] # Bonds to stay
        if verbose:
            print("{0} bonds stay, {1} new bonds, {2} bonds removed".format(len(bondsStay),
                                                                            len(bondsAdd), len(bondsRemove)))
        bondsToChange = bondsAdd + bondsRemove # Total num of bonds to change
        bondsIsAdd = [True] * len(bondsAdd) + [False] * len(bondsRemove) # Flag for if a bond should be added
        for bond, isAdd in zip(bondsToChange, bondsIsAdd): # loop thru all bonds
            ind = self.bondToInd[bond] # Get index of this bond
            paramset = self.activeParamDict if isAdd else self.inactiveParamDict # Fetch parameters
            self.bondForce.setBondParameters(ind, bond[0], bond[1], **paramset)  # actually updating bonds
        self.bondForce.updateParametersInContext(context)  # now run this to update things in the context
        return self.curBonds, pastBonds