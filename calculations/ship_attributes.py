from calculations.coefficients import *

class ShipAttributes:
    """
    Base class for common ship attributes and methods.
    """
    def __init__(self, L: float, B: float, TLC: float, TSC: float, Cb: float,
                 kr: float, GM: float, bilge_keel: bool):
        self.L = L                                  # Rule Length (m)
        self.L0 = max(110.0, self.L)                # Rule Length but not less than 110 (m)
        self.B = B                                  # Breadth Moulded (m)
        self.TLC = TLC                              # Loading Condition Draught (m)
        self.TSC = TSC                              # Scantling Draught (m)
        self.Cb = Cb                                # Block Coefficient (-)
        self.kr = kr                                # Roll Radius of Gyration (m)
        self.GM = GM                                # Metacentric height (m)
        self.fBK = calculate_fBK(bilge_keel)        # Coefficient fBK

        # Constants
        self.rho = 1.025                            # Seawater density (t/m³)
        self.g = 9.81                               # Gravity Acceleration (m/s²)