import numpy as np

class Extruder():
    def __init__(self, leg1, leg2, ebf1, left_blockers, right_blockers, extrusion_occupancy, lifetime=100, lifetime_stalled=10):
        """
        Class defining a generic Loop Extrusion Factor (LEF)
        Parameters:
            leftpos - int, (initial) position of left leg of extruder
            rightpos - int, (initial) position of right leg of extruder
            ebf1 - boolean, is ebf1 present?
            left/right_ blockers - dicts: each dict of form {pos:prob}, where pos is position of blocker and prob is capture probability [0,1]. Left-facing blockers and right-facing in respective dicts
            extrusion_occupancy - list of ints, denotes positions on polymer that are occupied by an extruder
        """
        self.leg1 = self.ExtruderLeg(leg1, -1)
        self.leg2 = self.ExtruderLeg(leg2, 1)
        self.ebf1 = ebf1
        self.blockers = {-1:left_blockers, 1:right_blockers}
        self.stalled = stalled
        self.occupied = extrusion_occupancy

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
        for leg in [self.leg1, self.leg2]:
            if np.random.random() < self.blockers[leg.side].get(leg.pos, 0):
                leg.attrs['captured'] = True
    def release(self):
        """
        And attempt to release an LEF captured by a blocker
        """
        if not self._any('captured'):
            return
        for leg in [self.leg1, self.leg2]:
            if (np.random.random() < self.blockers[leg.side].get(leg.pos, 0)) and (leg.attrs['captured']):
                leg.attrs['captured'] = False
    def translocate(self):
        """
        The main function. Performs 3 main functions:
            1. Attempts to unload LEFs with prob. 1/lifetime (1/stalled_lifetime if stalled)
            2. Attempts to capture extrusion blockers and release them
            3. Translocate (move/extrude) the LEF and determine if it is stalled
        """
        # 1 - attempt to unload LEF
        unload_prob = self.getUnloadProb()
        # 2 - attemp to capture and release blockers
        self.capture()
        self.release()
        # 3 - translocate cohesin
        for leg in [self.leg1, self.leg2]:
            if not leg.attrs['captured']:
                if self.occupied[leg.pos + leg.side] != 0:
                    leg.attrs['stalled'] = True
                else:
                    leg.attrs['stalled'] = False
                    self.occupied[leg.pos] = 0
                    self.occupied[led.pos + leg.side] = 1
                    leg.pos += leg.side
    def getUnloadProb(self):
        if self._any('stalled'):
            return 1/self.lifetime_stalled
        return 1/self.lifetime
    class ExtruderLeg():
        """
        Class defining one side / 'leg' of a loop extruder
        Conceptually, the extruder will be contacting two points on the polymer as it pulls two distal locations closer (like a cohesin ring)
        """
        def __init__(self, pos, side, attrs = {'stalled':False,'captured':False}):
            self.pos = pos
            self.side = side
            self.attrs = attrs

if __name__ == "__main__":
    pass

