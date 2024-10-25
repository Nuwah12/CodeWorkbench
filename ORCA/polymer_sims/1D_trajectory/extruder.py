import numpy as np
import sys

class Extruder():
    def __init__(self, extruder_index, leg1, leg2, left_blockers_capture, right_blockers_capture, left_blockers_release, right_blockers_release, extrusion_occupancy, lifetime=100, lifetime_stalled=10, targeted_loading=False, loading_spots=None, loading_dist=None):
        """
        Class defining a generic Loop Extrusion Factor (LEF)
        Parameters:
            leftpos - int, (initial) position of left leg of extruder
            rightpos - int, (initial) position of right leg of extruder
            ebf1 - boolean, is ebf1 present?
            left/right_ blockers - dicts: each dict of form {pos:prob}, where pos is position of blocker and prob is capture probability [0,1]. Left-facing blockers and right-facing in respective dicts
            extrusion_occupancy - list of ints, denotes positions on polymer that are occupied by an extruder
            targeted_loading - boolean, are extruders to be loaded at predetermined positions?
            loading_spots - list, all loading index positions for all extruders
            loading_dist - list, probabilities of loading at each spot 
            random_with_targeted - boolean, do random (uniform) loading in addition to targeted loading
        """
        self.ex_index = extruder_index 
        self.targeted_loading = targeted_loading
        # Checks
        if (not targeted_loading and loading_spots != None) or (targeted_loading and loading_spots == None):
            print('ERROR - Please specify both targeted loading as true AND a list of ALL loading spots for ALL extruders.')
            sys.exit(1)
        if targeted_loading: # If we have specified loading, we will first load the extruders at their defined spots 
            self.loading_spots = loading_spots
            leg1 = loading_spots[self.ex_index] # int of leg position
            leg2 = leg1+1 # int of other leg position
            self.nLEF = len(loading_spots)
            if loading_dist == None:
                print('{} extruders'.format(self.nLEF))
                self.loading_dist = np.linspace(1/self.nLEF, 1, num=self.nLEF)
                print(self.loading_dist)
            else:
                # make probabilities cumulative
                loading_distribution = []
                s = 0
                for i,x in enumerate(loading_dist):
                    s += x
                    loading_distribution.append(s)
                    
                self.loading_dist = loading_distribution
            
        self.leg1 = self.ExtruderLeg(leg1, -1)
        self.leg2 = self.ExtruderLeg(leg2, 1)
        self.blockers_capture = {-1:left_blockers_capture, 1:right_blockers_capture}
        self.blockers_release = {-1:left_blockers_release, 1:right_blockers_release}
        self.occupied = extrusion_occupancy
        #print(self.occupied)
        self.occupied[self.leg1.pos] = 1
        self.occupied[self.leg2.pos] = 1
        #print(self.occupied)

        self.lifetime = lifetime
        self.lifetime_stalled = lifetime_stalled

    def _any(self, attribute):
        """
        Returns true if either leg is true for attribute
        """
        return self.leg1.attrs[attribute] or self.leg2.attrs[attribute]
    def _all(self, attribute):
        """
        Returns true if both legs are true for attribute
        """
        return self.leg1.attrs[attribute] and self.leg2.attrs[attribute]
    def capture(self):
        """
        Attempt to 'capture' a LEF with a blocker
        """
        #print('Starting capture, leg info:\n')
        #self.printLegInfo()
        legs = [self.leg1, self.leg2]
        for i in range(len(legs)):
            p = np.random.random()
            #print('Capture prob. for leg {} at pos {}: {}'.format(leg.side, leg.pos, self.blockers_capture[leg.side].get(leg.pos, 0)))
            if p < self.blockers_capture[legs[i].side].get(legs[i].pos, 0):
                #print('Leg {} captured at pos {} with prob. {}'.format(legs[i].side, legs[i].pos, p))
                legs[i].setAttribute('captured',True)
    def release(self):
        """
        And attempt to release an LEF captured by a blocker
        """
        if not self._any('captured'):
            #print('No captured legs, not trying to release...')
            return
        for leg in [self.leg1, self.leg2]:
            p = np.random.random()
            if (p < self.blockers_release[leg.side].get(leg.pos, 0)):
                #print('Leg {} released at pos {} with prob {}'.format(leg.side, leg.pos, p))
                leg.attrs['captured'] = False
    def translocate(self, occupied):
        """
        The main function. Performs 3 main functions:
            1. Attempts to unload LEFs with prob. 1/lifetime (1/stalled_lifetime if stalled)
            2. Attempts to capture extrusion blockers and release them
            3. Translocate (move/extrude) the LEF and determine if it is stalled
        """
        # 1 - attempt to unload LEF
        unload_prob = self.getUnloadProb()
        if np.random.random() < unload_prob:
            #print('Unloading cohesin at pos {}, {}'.format(self.leg1.pos, self.leg2.pos))
            self.occupied[self.leg1.pos] = 0
            self.occupied[self.leg2.pos] = 0
            self.loadNew()
        # 2 - attemp to capture and release blockers
        self.capture()
        self.release()
        # 3 - translocate cohesin
        #print(self.occupied)
        for leg in [self.leg1, self.leg2]:
            #print('Leg {}: captured = {}'.format(leg.pos, leg.attrs['captured']))
            if not leg.attrs['captured']:
                if (leg.pos + leg.side) >= len(self.occupied) or self.occupied[leg.pos + leg.side] != 0:
                    #print('Leg {} stalled at pos {} because it is marked as {}'.format(leg.side, leg.pos, self.occupied[leg.pos+leg.side]))
                    leg.attrs['stalled'] = True
                else:
                    #print('Leg at {} translocating'.format(leg.pos))
                    leg.attrs['stalled'] = False
                    self.occupied[leg.pos] = 0
                    self.occupied[leg.pos + leg.side] = 1
                    #print(self.occupied)
                    leg.pos += leg.side
        #print('Extruder now at ({},{})'.format(self.leg1.pos,self.leg2.pos))
        return self.occupied
    def getUnloadProb(self):
        if self._any('stalled'):
            return 1/self.lifetime_stalled
        return 1/self.lifetime
    def loadNew(self):
        """
        TODO: function to initialize new leg positions after LEF is unloaded
        """
        while True:
            if self.targeted_loading: # This is where the ""magic"" happens for targeted loading!
                p = np.random.random()
                #print(p)
                pos = self.loading_spots[len(self.loading_dist)-1]
                # find first index greater than p
                for i,x in enumerate(self.loading_dist):
                    if p > x:
                        continue
                    elif p < x:
                        pos = self.loading_spots[i]
                        break
                pos = self.loading_spots[i]
            else:    
                pos = np.random.randint(len(self.occupied)-1)
            if self.occupied[pos] != 1 and self.occupied[pos+1] != 1:
                #print('loading extruder at spot {}'.format(pos))
                break
        self.leg1.pos = pos
        self.leg2.pos = pos+1
        # Reset leg attributes
        for leg in [self.leg1, self.leg2]:
            leg.attrs['captured'] = False
            leg.attrs['stalled'] = False
        self.occupied[pos] = 1
        self.occupied[pos+1] = 1
    def printLegInfo(self):
        for leg in [self.leg1, self.leg2]:
            s = "Leg {}:\nPosition: {}\n Captured: {}\n Stalled: {}\n".format(leg.side, leg.pos, leg.attrs['captured'], leg.attrs['stalled'])
            #print(s)
    class ExtruderLeg():
        """
        Class defining one side / 'leg' of a loop extruder
        Conceptually, the extruder will be contacting two points on the polymer as it pulls two distal locations closer (like a cohesin ring)
        """
        def __init__(self, pos, side, attrs = None):
            self.pos = pos
            self.side = side
            if attrs is None:
                self.attrs = {'stalled': False, 'captured': False}
        def setAttribute(self, attribute, value):
            self.attrs[attribute] = value

if __name__ == "__main__":
    pass

