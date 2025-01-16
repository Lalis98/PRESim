import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from calculations.data_processing import *
from calculations.coefficients import *
from calculations.pressure_calculations import *
from calculations.plotting import *

class ExternalSeaPressureCalc:
    """
    This class calculates the pressures for dynamic load cases according to Common Structural Rules,
    based on various parameters like ship dimensions, loading conditions, and draughts. It is used for
    bulk carriers and tankers.

    :param coordinates: 2D numpy array with shape (n, 3), where n is the number of coordinates
    :type coordinates: np.ndarray
    :param L: Rule Length of the ship (m).
    :type L: float
    :param B: Breadth Moulded of the ship (m).
    :type B: float
    :param TLC: Loading Condition Draught (m).
    :type TLC: float
    :param TSC: Scantling Draught (m).
    :type TSC: float
    :param Cb: Block Coefficient (-).
    :type Cb: float
    :param Bx_distribution: Moulded Breadth at the Waterline positions and values
    :type Bx_distribution: np.ndarray
    :param kr: Roll Radius of Gyration (m).
    :type kr: float
    :parameter GM: Metacentric height (m).
    :type GM: float

    """

    def __init__(self, coordinates: np.ndarray, L: float, B: float, TLC: float, TSC: float, Cb: float,
                 Bx_distribution: np.ndarray, kr: float, GM:float, bilge_keel: bool):


        if not isinstance(coordinates, np.ndarray):   # Ensure coordinates is a 2D numpy array
            raise ValueError("coordinates must be a numpy array")


        if coordinates.ndim != 2 or coordinates.shape[1] != 3:  # Ensure dimensions are correct
            raise ValueError("coordinates must be a 2D numpy array with 3 columns")

        # Set the attributes
        self.coordinates = coordinates              # Coordinates (x,y,z) (m)
        self.L = L                                  # Rule Length (m)
        self.L0 = max(110.0, self.L)                # Rule Length but not taken less than 110 (m)
        self.B = B                                  # Breadth Moulded (m)
        self.TLC = TLC                              # Loading Condition Draught (m)
        self.TSC = TSC                              # Scantling Draught (m)
        self.Cb = Cb                                # Block Coefficient (-)
        self.Bx_positions = Bx_distribution[0, :]   # Moulded Breath at the Waterline (m)
        self.Bx_values = Bx_distribution[1, :]      # Moulded Breath at the Waterline (m)
        self.kr = kr                                # Roll Radius of Gyration (m)
        self.GM = GM                                # Metacentric height (m)
        self.fBK = calculate_fBK(bilge_keel)        # Coefficient fBK

        # Constants
        self.rho = 1.025                            # Seawater density (t/m^3)
        self.g = 9.81                               # Gravity Acceleration (m/s^2)

        # Internal Calculations
        self.fT = max(0.5, TLC / TSC)                   # Ratio of Loading to Scantling Draught
        self.Cw = calculate_wave_coefficient(self.L)    # Wave coefficient (-)

        # Save the results with condition
        self.last_pressure_data = {
            "data": np.array([]),
            "case": ""
        }


    def calculate_load_case(self, base_case: str, sub_case: str, fps: float, load_scenario: str = "Extreme Sea Loads",
                            fb: float = None, design_load: str = "S+D"):
        """
        Calculates the load case based on the specified base case, subcase, and parameters.

        :param base_case: The base case identifier. Options include: "HSM", "FSM", "HSA",
         "BSR", "BSP", "OST", "OSA".
        :type base_case: str
        :param sub_case: The subcase identifier used to refine the load calculation. Options include: "1", "2", "1P",
         "1S", "2P", "2S".
        :type sub_case: str
        :param fps: Coefficient for strength assessments.
        :type fps: float
        :param load_scenario: Description of the load scenario (if needed), "Ballast Water Exchange" or "Extreme Sea Loads".
        :type load_scenario: str
        :param fb: Heading correction factor (if needed).
        :type fb: float
        :param design_load: Load scenario ("S", "D", or "S+D"). Defaults to "S+D".
        :type design_load: str

        :return: The result of the load case calculation, specific to the base case and subcase.
        :rtype: Any
        """

        if base_case == "HSM":
            return self._calculate_hsm_case(sub_case, load_scenario, fb, fps, design_load)
        elif base_case == "FSM":
            return self._calculate_fsm_case(sub_case, load_scenario, fb, fps, design_load)
        elif base_case == "HSA":
            return self._calculate_hsa_case(sub_case, load_scenario, fps, design_load)
        elif base_case == "BSR":
            return self._calculate_bsr_case(sub_case, fb, fps, design_load)
        elif base_case == "BSP":
            return self._calculate_bsp_case(sub_case, load_scenario, fb, fps, design_load)
        elif base_case == "OST":
            return self._calculate_ost_case(sub_case, load_scenario, fps, design_load)
        elif base_case == "OSA":
            return self._calculate_osa_case(sub_case, load_scenario, fps, design_load)  # fb not used
        else:
            raise ValueError(f"Unknown base case: {base_case}")

    def print_data_summary(self):
        """
        Print the dimensions and important attributes of the load analysis in a table format with units in parentheses.
        """
        print("=" * 70)  #f
        print(f"{'EXTERNAL SEA LOADS: DATA SUMMARY':^70}")
        print("=" * 70)
        print(f"{'Parameter':<40}{'Unit':<10}{'Value':>20}")
        print("-" * 70)
        print(f"{'Rule Length (L)':<40}{'(m)':<10}{f'{self.L:.3f}':>20}")
        print(f"{'Breadth (B)':<40}{'(m)':<10}{f'{self.B:.2f}':>20}")
        print(f"{'Loading Condition Draught (TLC)':<40}{'(m)':<10}{f'{self.TLC:.3f}':>20}")
        print(f"{'Scantling Draught (TSC)':<40}{'(m)':<10}{f'{self.TSC:.3f}':>20}")
        print(f"{'Block Coefficient (Cb)':<40}{'(-)':<10}{f'{self.Cb:.3f}':>20}")
        print(f"{'Roll Radius of Gyration (kr)':<40}{'(m)':<10}{f'{self.kr:.3f}':>20}")
        print(f"{'Metacentric Height (GM)':<40}{'(m)':<10}{f'{self.GM:.3f}':>20}")
        print(f"{'Wave Coefficient (Cw)':<40}{'(-)':<10}{f'{self.Cw:.3f}':>20}")
        print(f"{'Loading/Scantling Draught Ratio (fT)':<40}{'(-)':<10}{f'{self.fT:.3f}':>20}")
        print(f"{'Seawater Density (ρ)':<40}{'(t/m³)':<10}{f'{self.rho:.3f}':>20}")
        print(f"{'Gravity Acceleration (g)':<40}{'(m/s²)':<10}{f'{self.g:.2f}':>20}")
        print("=" * 70)

    def _calculate_hsm_case(self, sub_case: str, load_scenario: str, fb: float, fps: float, design_load: str)-> np.ndarray:
        """
        Calculates Head Sea Equivalent Design Waves that minimize and maximize the vertical wave bending moment amidships.

        :param sub_case: Identifier for the Head Sea case, e.g., "1" or "2".
        :type sub_case: str
        :param load_scenario: Description of the load scenario, e.g., "Ballast Water Exchange" or "Extreme Sea Loads".
        :type load_scenario: str
        :param fb: Heading correction factor.
        :type fb: float
        :param fps: Coefficient for strength assessments.
        :type fps: float
        :param design_load: Load scenario ("S", "D", or "S+D"). Defaults to "S+D".
        :type design_load: str

        :return: Calculated pressure based on the design load.
        :rtype: np.ndarray
        """

        lamda = 0.6 * (1 + self.fT) * self.L  # Wave length for the HSM [m]

        def _get_fnl_arrays(load_scenario_):

            if load_scenario_ == "Ballast Water Exchange":

                fnl_fxL_array = np.array([
                    0.0,
                    0.3,
                    0.7,
                    1.0
                ])

                fnl_array = np.array([
                    0.85,
                    0.95,
                    0.95,
                    0.80
                ])

            elif load_scenario_ == "Extreme Sea Loads":

                fnl_fxL_array = np.array([
                    0.0,
                    0.3,
                    0.7,
                    1.0
                ])

                fnl_array = np.array([
                    0.7,
                    0.9,
                    0.9,
                    0.6
                ])
            else:
                raise ValueError(f"Unrecognized load_scenario: {load_scenario_}")

            return fnl_fxL_array, fnl_array

        def _calculate_coefficients(x, y, z):

            Bx = interpolate(x, self.Bx_positions, self.Bx_values)
            fyB = np.abs(2 * y) / Bx

            fxL = x / self.L

            fyz = z / self.TLC + fyB + 1

            fh = 3.0 * (1.21 - 0.66 * self.fT)

            fnl_fxL_array, fnl_array = _get_fnl_arrays(load_scenario_=load_scenario)
            fnl = interpolate(fxL, fnl_fxL_array, fnl_array)

            if fxL < 0.15:
                ka = (0.5 + self.fT) * ((3 - 2 * np.sqrt(fyB)) - 20 / 9 * fxL * (7 - 6 * np.sqrt(fyB))) + 2 / 3 * (1 - self.fT)
            elif 0.15 <= fxL <= 0.7:
                ka = 1.0
            else:
                ka = 1 + (fxL - 0.7) * ((40 / 3 * self.fT - 5) + 2 * (1 - fyB)
                                        * (18 / self.Cb * self.fT * (fxL - 0.7) - 0.25 * (2 - self.fT)))

            kp_fxL_array = np.array([0, 0.3 - 0.1 * self.fT, 0.35 - 0.1 * self.fT, 0.8 - 0.2 * self.fT, 0.9 - 0.2 * self.fT, 1.0])
            kp_array = np.array([-0.25 * self.fT * (1 + fyB), -1.0, 1.0, 1.0, -1.0, -1.0])
            kp = interpolate(fxL, kp_fxL_array, kp_array)

            return fnl, fh, ka, kp, fyz

        def _pressure_hs(fb_, fps_, fnl_, fh_, ka_, kp_, fyz_, Cw_, L0_, lamda_, L_):
            return fb_ * fps_ * fnl_ * fh_ * ka_ * kp_ * fyz_ * Cw_ * np.sqrt((L0_ + lamda_ - 125) / L_)

        def _calculate_total_pressure(data_points):

            pressure_data_ = []

            for i, point in enumerate(data_points):

                x, y, z = point
                fnl, fh, ka, kp, fyz = _calculate_coefficients(x, y, z)
                Bx = interpolate(x, self.Bx_positions, self.Bx_values)
                fnl_wl, fh_wl, ka_wl, kp_wl, fyz_wl = _calculate_coefficients(x, Bx / 2,
                                                                              self.TLC)  # Pw,wl=Pw for y=B/2 and z=TLC

                Phs = _pressure_hs(fb, fps, fnl, fh, ka, kp, fyz, self.Cw, self.L0, lamda, self.L)
                if sub_case == "1":  # HSM-1
                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -_pressure_hs(fb, fps, fnl_wl, fh_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )
                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Phs)

                elif sub_case == "2":

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        _pressure_hs(fb, fps, fnl_wl, fh_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )

                    hw = Pw_wl / (self.rho * self.g)
                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Phs)

                else:
                    raise ValueError("Wrong Dynamic Load Condition. Please check sub_case parameter!")

                Ps = calculate_hydrostatic_pressure(z, self.TLC, self.rho, self.g)
                Pex = calculate_pressure(Ps, Pw, design_load)  # Pex = Pw + Ps
                pressure_data_.append(np.array([x, y, z, Pex]))

            return pressure_data_

        pressure_data = _calculate_total_pressure(data_points=self.coordinates)
        pressure_data = np.array(pressure_data)  # Convert to numpy array

        self.update_pressure_data(pressure_data, f"Dynamic Load Case: HSM-{sub_case} ({design_load})")  # Update last pressure data

        return pressure_data

    def _calculate_fsm_case(self, sub_case: str, load_scenario: str, fb: float, fps: float, design_load: str)-> np.ndarray:
        """
        Calculates the Following Sea Equivalent Design Waves that minimize and maximize the vertical wave bending moment amidships.

        :param sub_case: Identifier for the Head Sea case, e.g., "1" or "2".
        :type sub_case: str
        :param load_scenario: Description of the load scenario, e.g., "Ballast Water Exchange" or "Extreme Sea Loads".
        :type load_scenario: str
        :param fb: Heading correction factor.
        :type fb: float
        :param fps: Coefficient for strength assessments.
        :type fps: float
        :param design_load: Load scenario ("S", "D", or "S+D"). Defaults to "S+D".
        :type design_load: str

        :return: Calculated pressure based on the design load.
        :rtype: np.ndarray
        """

        lamda = 0.6 * (1 + 2 / 3 * self.fT) * self.L  # Wave length for FSM [m]
        if load_scenario == "Extreme Sea Loads":
            fnl = 0.9
        elif load_scenario == "Ballast Water Exchange":
            fnl = 0.95
        else:
            raise ValueError(f"Unrecognized load_scenario: {load_scenario}")
        fh = 2.6

        def _calculate_coefficients(x, y, z):

            fxL = x / self.L

            Bx = interpolate(x, self.Bx_positions, self.Bx_values)
            fyB = np.abs(2 * y) / Bx

            fyz = z / self.TLC + fyB + 1

            if fxL < 0.20:
                ka = 1 + (3.75 - 2 * self.fT) * (1 - 5 * fxL) * (1 - fyB)
            elif 0.2 <= fxL <= 0.9:
                ka = 1.0
            else:
                ka = 1 + 20 * (1 - fyB) * (fxL - 0.9)

            kp_fxL_array = np.array([
                0,
                0.35 - 0.1 * self.fT,
                0.50 - 0.2 * self.fT,
                0.75,
                0.80,
                1.0
            ])

            kp_array = np.array([
                -0.75 - 0.25 * fyB,
                -1.0,
                1.0,
                1.0,
                -1.0,
                -0.75 - 0.25 * fyB
            ])

            kp = interpolate(fxL, kp_fxL_array, kp_array)

            return ka, kp, fyz

        def _pressure_fs(fb_, fps_, fnl_, fh_, ka_, kp_, fyz_, Cw_, L0_, lamda_, L_):
            return fb_ * fps_ * fnl_ * fh_ * ka_ * kp_ * fyz_ * Cw_ * np.sqrt((L0_ + lamda_ - 125) / L_)

        def _calculate_total_pressure(data_points):

            pressure_data_ = []

            for i, point in enumerate(data_points):

                x, y, z = point
                ka, kp, fyz = _calculate_coefficients(x, y, z)
                Bx = interpolate(x, self.Bx_positions, self.Bx_values)
                ka_wl, kp_wl, fyz_wl = _calculate_coefficients(x, Bx / 2, self.TLC)
                fnl_wl, fh_wl = fnl, fh

                Pfs = _pressure_fs(fb, fps, fnl, fh, ka, kp, fyz, self.Cw, self.L0, lamda, self.L)

                if sub_case == "1":  # FSM-1

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -_pressure_fs(fb, fps, fnl_wl, fh_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )
                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Pfs)

                elif sub_case == "2":

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        _pressure_fs(fb, fps, fnl_wl, fh_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )
                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Pfs)

                else:
                    raise ValueError("Wrong Dynamic Load Condition. Please check sub_case parameter!")

                Ps = calculate_hydrostatic_pressure(z, self.TLC, self.rho, self.g)
                Pex = calculate_pressure(Ps, Pw, design_load)  # Pex = Pw + Ps
                pressure_data_.append(np.array([x, y, z, Pex]))

            return pressure_data_

        pressure_data = _calculate_total_pressure(data_points=self.coordinates)
        pressure_data = np.array(pressure_data)  # Convert to numpy array

        self.update_pressure_data(pressure_data, f"Dynamic Load Case: FSM-{sub_case} ({design_load})") # Update last pressure data

        return pressure_data

    def _calculate_hsa_case(self, sub_case: str, load_scenario: str, fps: float, design_load: str) -> np.ndarray:
        """
        Calculates the Head Sea Equivalent Design Waves that maximize and minimize the
        head sea vertical acceleration at FP.

        :param sub_case: Identifier for the Head Sea case, e.g., "1" or "2".
        :type sub_case: str
        :param load_scenario: Description of the load scenario, e.g., "Ballast Water Exchange" or "Extreme Sea Loads".
        :type load_scenario: str
        :param fps: Coefficient for strength assessments.
        :type fps: float
        :param design_load: Load scenario ("S", "D", or "S+D"). Defaults to "S+D".
        :type design_load: str

        :return: Calculated pressure based on the design load.
        :rtype: np.ndarray
        """
        lamda = 0.6 * (1 + self.fT) * self.L  # Wave length for FSM [m]
        fh = 2.4 * (1.21 - 0.66 * self.fT)  # Coefficient

        def _pressure_hs(fps_, fnl_, fh_, ka_, kp_, fyz_, Cw_, L0_, lamda_, L_):
            return fps_ * fnl_ * fh_ * ka_ * kp_ * fyz_ * Cw_ * np.sqrt((L0_ + lamda_ - 125) / L_)

        def _get_load_scenario_arrays(load_scenario_):

            if load_scenario_ == "Ballast Water Exchange":

                fnl_fxL_array = np.array([
                    0.0,
                    0.3,
                    0.7,
                    1.0
                ])

                fnl_array = np.array([
                    0.85,
                    0.95,
                    0.95,
                    0.80
                ])

            elif load_scenario_ == "Extreme Sea Loads":

                fnl_fxL_array = np.array([
                    0.0,
                    0.3,
                    0.7,
                    1.0
                ])

                fnl_array = np.array([
                    0.7,
                    0.9,
                    0.9,
                    0.6
                ])
            else:
                raise ValueError(f"Unrecognized load_scenario: {load_scenario_}")

            return fnl_fxL_array, fnl_array

        def _calculate_coefficients(x, y, z):

            fxL = x / self.L

            Bx = interpolate(x, self.Bx_positions, self.Bx_values)
            fyB = np.abs(2 * y) / Bx

            fyz = z / self.TLC + fyB + 1

            if fxL < 0.15:
                ka = (0.5 + self.fT) * ((3 - 2 * np.sqrt(fyB)) - 20 / 9 * fxL * (7 - 6 * np.sqrt(fyB))) + 2 / 3 * (1 - self.fT)
            elif 0.15 <= fxL <= 0.7:
                ka = 1.0
            else:
                ka = 1 + (fxL - 0.7) * ((40 / 3 * self.fT - 5) + 2 * (1 - fyB)
                                        * (18 / self.Cb * self.fT * (fxL - 0.7) - 0.25 * (2 - self.fT)))

            kp_fxL_array = np.array([
                0,
                0.30 - 0.10 * self.fT,
                0.50 - 0.20 * self.fT,
                0.80 - 0.20 * self.fT,
                0.90 - 0.20 * self.fT,
                1.0
            ])

            kp_array = np.array([
                1.5 - self.fT - 0.5 * fyB,
                -1.0,
                1.0,
                1.0,
                -1.0,
                -1.0
            ])

            kp = interpolate(fxL, kp_fxL_array, kp_array)

            fnl_fxL_array, fnl_array = _get_load_scenario_arrays(load_scenario_=load_scenario)
            fnl = interpolate(fxL, fnl_fxL_array, fnl_array)

            return ka, kp, fyz, fnl

        def _calculate_total_pressure(data_points):

            pressure_data_ = []

            for i, point in enumerate(data_points):

                x, y, z = point

                ka, kp, fyz, fnl = _calculate_coefficients(x, y, z)

                Bx = interpolate(x, self.Bx_positions, self.Bx_values)

                ka_wl, kp_wl, fyz_wl, fnl_wl = _calculate_coefficients(x, Bx / 2,
                                                                       self.TLC)  # Pw,wl=Pw for y=B/2 and z=TLC

                fh_wl = fh

                Pfs = _pressure_hs(fps, fnl, fh, ka, kp, fyz, self.Cw, self.L0, lamda, self.L)

                if sub_case == "1":  # Dynamic Load Case HSA-1

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -_pressure_hs(fps, fnl_wl, fh_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )
                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Pfs)

                elif sub_case == "2":  # Dynamic Load Case HSA-2

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        _pressure_hs(fps, fnl_wl, fh_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )
                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Pfs)

                else:
                    raise ValueError("Wrong Dynamic Load Condition. Please check sub_case parameter!")

                Ps = calculate_hydrostatic_pressure(z, self.TLC, self.rho, self.g)
                Pex = calculate_pressure(Ps, Pw, design_load)  # Pex = Pw + Ps
                pressure_data_.append(np.array([x, y, z, Pex]))

            return pressure_data_

        pressure_data = _calculate_total_pressure(data_points=self.coordinates)
        pressure_data = np.array(pressure_data)  # Convert to numpy array

        self.update_pressure_data(pressure_data, f"Dynamic Load Case: HSA-{sub_case} ({design_load})")  # Update last pressure data

        return pressure_data

    def _calculate_bsr_case(self, sub_case: str, fb: float, fps: float, design_load: str) -> np.ndarray:
        """
        Calculates Beam Sea Equivalent Design Waves that minimize and maximize roll motion on the port or starboard side.

        :param sub_case: Beam Sea case identifier, e.g., "1P", "1S", "2P", or "2S".
        :type sub_case: str
        :param fb: Description of the load scenario, e.g., "Ballast Water Exchange" or "Extreme Sea Loads".
        :type fb: float
        :param fps: Coefficient for strength assessments.
        :type fps: float
        :param design_load: Load scenario ("S", "D", or "S+D"). Defaults to "S+D".
        :type design_load: str

        :return: Calculated pressure based on the design load.
        :rtype: np.ndarray
        """

        fnl = 1  # For both "Extreme Sea Loads" or "Ballast Water Exchange" fnl = 1.0

        def _get_roll_angle(kr_, GM_):
            T_thita_ = 2.3 * np.pi * kr_ / np.sqrt(self.g * GM_)
            thita_ = 9000 * (1.25 - 0.025 * T_thita_) * fps * self.fBK / ((self.B + 75) * np.pi)  # in degrees
            thita_rad_ = thita_ * np.pi / 180.0

            return T_thita_, thita_rad_

        def _pressure_bsr(fb_, fps_, fnl_, fyB1_, Cw_, L0_, lamda_, L_, y_, thita_, condition):

            if condition in ["1P", "2P"]:
                return fb_ * fnl_ * (10 * y_ * np.sin(thita_) + 0.88 * fps_ * Cw_ *
                                     np.sqrt((L0_ + lamda_ - 125) / L_) * (fyB1_ + 1))

            elif condition in ["1S", "2S"]:
                return fb_ * fnl_ * (-10 * y_ * np.sin(thita_) + 0.88 * fps_ * Cw_ *
                                     np.sqrt((L0_ + lamda_ - 125) / L_) * (fyB1_ + 1))

            else:
                raise ValueError("Wrong Dynamic Load Condition. Please check sub_case parameter!")

        def _calculate_coefficients(x, y, z):

            T_thita_, thita_ = _get_roll_angle(self.kr, self.GM)
            fyB1_ = np.abs(2 * y) / self.B
            lamda_ = self.g / (2 * np.pi ** 2) * T_thita_ ** 2  # Wave length for BSR dynamic Load Case [m]

            return fyB1_, lamda_, thita_

        def _calculate_total_pressure(data_points):
            # Pw = max(-Pfs, rho * g * (z - TLC))

            pressure_data_ = []

            for i, point in enumerate(data_points):

                x, y, z = point
                fyB1, lamda, thita = _calculate_coefficients(x, y, z)
                Bx = interpolate(x, self.Bx_positions, self.Bx_values)
                fyB1_wl, lamda_wl, thita_wl = _calculate_coefficients(x, Bx / 2, self.TLC)

                Pbsr = _pressure_bsr(fb, fps, fnl, fyB1, self.Cw, self.L0, lamda, self.L, y, thita, sub_case)

                if sub_case == "1P":  # Dynamic Load Case BSR-1P

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        _pressure_bsr(fb, fps, fnl, fyB1_wl, self.Cw, self.L0, lamda_wl, self.L, y, thita_wl, sub_case)
                    )
                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Pbsr)

                elif sub_case == "2P":  # Dynamic Load Case BSR-2P

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -_pressure_bsr(fb, fps, fnl, fyB1_wl, self.Cw, self.L0, lamda_wl, self.L, y, thita_wl, sub_case)
                    )
                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Pbsr)

                elif sub_case == "1S":  # Dynamic Load Case BSR-1S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        _pressure_bsr(fb, fps, fnl, fyB1_wl, self.Cw, self.L0, lamda_wl, self.L, y, thita_wl, sub_case)
                    )
                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Pbsr)

                elif sub_case == "2S":  # Dynamic Load Case BSR-2S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -_pressure_bsr(fb, fps, fnl, fyB1_wl, self.Cw, self.L0, lamda_wl, self.L, y, thita_wl, sub_case)
                    )
                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Pbsr)

                else:
                    raise ValueError(
                        "Wrong Dynamic Load Condition. Please check condition for BSR-1P, BSR-2P, BSR-1S, BSR-2S!")

                Ps = calculate_hydrostatic_pressure(z, self.TLC, self.rho, self.g)
                Pex = calculate_pressure(Ps, Pw, design_load)  # Pex = Pw + Ps
                pressure_data_.append(np.array([x, y, z, Pex]))

            return pressure_data_

        pressure_data = _calculate_total_pressure(data_points=self.coordinates)
        pressure_data = np.array(pressure_data)  # Convert to numpy array

        self.update_pressure_data(pressure_data, f"Dynamic Load Case: BSR-{sub_case} ({design_load})")  # Update last pressure data

        return pressure_data

    def _calculate_bsp_case(self, sub_case: str, load_scenario: str, fb: float, fps: float, design_load: str) -> np.ndarray:
        """
        Calculates Beam Sea Equivalent Design Waves that maximize and minimize hydrodynamic pressure at the waterline amidships.

        :param sub_case: Beam Sea case identifier, e.g., "1P", "1S", "2P", or "2S".
        :type sub_case: str
        :param load_scenario: Description of the load scenario, e.g., "Ballast Water Exchange" or "Extreme Sea Loads".
        :type load_scenario: str
        :param fb: Heading correction factor, used to adjust calculations for wave heading.
        :type fb: float
        :param fps: Coefficient for strength assessments.
        :type fps: float
        :param design_load: Load scenario ("S", "D", or "S+D"). Defaults to "S+D".
        :type design_load: str

        :return: Calculated pressure based on the design load.
        :rtype: np.ndarray
        """

        lamda = 0.2 * (1 + 2 * self.fT) * self.L  # Wave length for BSP[m]

        def _get_load_scenario_arrays(load_scenario_):

            if load_scenario_ == "Ballast Water Exchange":

                fnl_fxL_array = np.array([
                    0.0,
                    0.3,
                    0.7,
                    1.0
                ])

                fnl_array = np.array([
                    0.6,
                    0.8,
                    0.8,
                    0.6
                ])

            elif load_scenario_ == "Extreme Sea Loads":

                fnl_fxL_array = np.array([
                    0.0,
                    0.3,
                    0.7,
                    1.0
                ])

                fnl_array = np.array([
                    0.6,
                    0.8,
                    0.8,
                    0.6
                ])
            else:
                raise ValueError(f"Unrecognized load_scenario: {load_scenario_}")

            return fnl_fxL_array, fnl_array

        def _pressure_bsp(fb_, fps_, fnl_, fyz_, Cw_, L0_, lamda_, L_):
            return 4.5 * fb_ * fps_ * fnl_ * fyz_ * Cw_ * np.sqrt((L0_ + lamda_ - 125) / L_)

        def _calculate_fyz(y_, z_, TLC_, fyB1_, condition_):

            if condition_ in ["1P", "2P"]:
                if y_ >= 0.0:
                    fyz = 2 * z_ / TLC_ + 2.5 * fyB1_ + 0.5
                else:
                    fyz = 2 / 3 * z_ / TLC_ + 1 / 2 * fyB1_ + 0.5
            elif condition_ in ["1S", "2S"]:
                if y_ >= 0.0:
                    fyz = 2 / 3 * z_ / TLC_ + 1 / 2 * fyB1_ + 0.5
                else:
                    fyz = 2 * z_ / TLC_ + 2.5 * fyB1_ + 0.5
            else:
                raise ValueError(f"Invalid condition: {condition_}. Must be one of BSP-1P, BSP-2P, BSP-1S, BSP-2S.")

            return fyz

        def _calculate_coefficients(x, y, z):

            fxL = x / self.L

            fnl_fxL_array, fnl_array = _get_load_scenario_arrays(load_scenario_=load_scenario)
            fnl = interpolate(fxL, fnl_fxL_array, fnl_array)

            fyB1 = np.abs(2 * y) / self.B

            fyz = _calculate_fyz(y, z, self.TLC, fyB1, sub_case)

            return fnl, fyz

        def _calculate_total_pressure(data_points):

            pressure_data_ = []

            for i, point in enumerate(data_points):

                x, y, z = point
                Bx = interpolate(x, self.Bx_positions, self.Bx_values)
                fnl, fyz = _calculate_coefficients(x, y, z)
                fnl_wl, fyz_wl = _calculate_coefficients(x, Bx / 2, self.TLC)  # Pw,wl=Pw for y=B/2 and z=TLC
                Pbsp = _pressure_bsp(fb, fps, fnl, fyz, self.Cw, self.L0, lamda, self.L)

                if sub_case == "1P":  # Dynamic Load Case BSP-1S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        _pressure_bsp(fb, fps, fnl_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Pbsp)

                elif sub_case == "2P":  # Dynamic Load Case BSP-2P

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -_pressure_bsp(fb, fps, fnl_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Pbsp)

                elif sub_case == "1S":  # Dynamic Load Case BSP-1S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        _pressure_bsp(fb, fps, fnl_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Pbsp)

                elif sub_case == "2S":  # Dynamic Load Case BSP-2S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -_pressure_bsp(fb, fps, fnl_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Pbsp)

                else:
                    raise ValueError(
                        "Wrong Dynamic Load Condition. Please check condition for BSP-1P, BSP-2P, BSP-1S, BSP-2S!")

                Ps = calculate_hydrostatic_pressure(z, self.TLC, self.rho, self.g)
                Pex = calculate_pressure(Ps, Pw, design_load)  # Pex = Pw + Ps
                pressure_data_.append(np.array([x, y, z, Pex]))

            return pressure_data_

        pressure_data = _calculate_total_pressure(data_points=self.coordinates)
        pressure_data = np.array(pressure_data)  # Convert to numpy array

        self.update_pressure_data(pressure_data, f"Dynamic Load Case: BSP-{sub_case} ({design_load})")  # Update last pressure data

        return pressure_data

    def _calculate_ost_case(self, sub_case: str, load_scenario: str, fps: float, design_load: str) -> np.ndarray:
        """
        Calculates Oblique Sea Equivalent Design Waves that minimize and maximize torsional moment
        at 0.25L from the aft end.

        :param sub_case: Oblique Sea case identifier, e.g., "1P", "1S", "2P", or "2S".
        :type sub_case: str
        :param load_scenario: Description of the load scenario, e.g., "Ballast Water Exchange" or "Extreme Sea Loads".
        :type sub_case: str
        :param fps: Coefficient for strength assessments.
        :type sub_case: float
        :param design_load: Load scenario ("S", "D", or "S+D"). Defaults to "S+D".
        :type sub_case: str

        :return: Calculated pressure based on the design load.
        :rtype: np.ndarray
        """

        lamda = 0.45 * self.L  # Wave length for OST [m]

        def pressure_ost(fps_, fnl_, ka_, kp_, fyz_, Cw_, L0_, lamda_, L_):
            return 1.38 * fps_ * fnl_ * ka_ * kp_ * fyz_ * Cw_ * np.sqrt((L0_ + lamda_ - 125) / L_)

        def _calculate_ka(y_, fxL_, fyB_, fT_, condition_):

            if condition_ in ["1P", "2P"]:

                if y_ >= 0.0:

                    if fxL_ <= 0.2:
                        ka_ = 1.0 + 3.5 * (1 - fyB_) * (1 - 5 * fxL_)
                    elif 0.2 < fxL_ <= 0.8:
                        ka_ = 1.0
                    else:  # fxL > 0.8
                        ka_ = 1.0

                else:

                    if fxL_ <= 0.2:
                        ka_ = 1.0 + (3.5 - (4 * fT_ - 0.5) * fyB_) * (1 - 5 * fxL_)
                    elif 0.2 < fxL_ <= 0.8:
                        ka_ = 1.0
                    else:  # fxL > 0.8
                        ka_ = 1.0 + 4 * (1 - fT_) * (5 * fxL_ - 4) * fyB_

            elif condition_ in ["1S", "2S"]:

                if y_ >= 0.0:

                    if fxL_ <= 0.2:
                        ka_ = 1.0 + (3.5 - (4 * fT_ - 0.5) * fyB_) * (1 - 5 * fxL_)
                    elif 0.2 < fxL_ <= 0.8:
                        ka_ = 1.0
                    else:  # fxL > 0.8
                        ka_ = 1.0 + 4 * (1 - fT_) * (5 * fxL_ - 4) * fyB_

                else:

                    if fxL_ <= 0.2:
                        ka_ = 1.0 + 3.5 * (1 - fyB_) * (1 - 5 * fxL_)
                    elif 0.2 < fxL_ <= 0.8:
                        ka_ = 1.0
                    else:  # fxL > 0.8
                        ka_ = 1.0

            else:
                raise ValueError(f"Invalid condition: {condition_}. Must be one of OST-1P, OST-2P, OST-1S, OST-2S.")

            return ka_

        def _get_kp_arrays(y_, fyB_, fT_, condition_):

            if condition_ in ["1P", "2P"]:

                if y_ >= 0:

                    kp_fxL_array_ = np.array([
                        0.0,
                        0.2,
                        0.4,
                        0.5,
                        0.7,
                        0.9,
                        1.0
                    ])

                    kp_array_ = np.array([
                        1.0,
                        1.0,
                        -1.0,
                        -1.0,
                        -0.1 + (1.6 * fT_ - 1.5) * fyB_,
                        0.8 + 0.2 * fyB_,
                        -1.0 + fyB_
                    ])

                else:  # y < 0

                    kp_fxL_array_ = np.array([
                        0.0,
                        0.2,
                        0.4,
                        0.5,
                        0.7,
                        0.9,
                        1.0
                    ])

                    kp_array_ = np.array([
                        1.0,
                        1.0 + (0.75 - 1.5 * fT_) * fyB_,
                        -1.0 + (1.75 - 0.5 * fT_) * fyB_,
                        -1.0 + (1.75 - 0.5 * fT_) * fyB_,
                        -0.1 + (0.25 - 0.3 * fT_) * fyB_,
                        0.8 - (0.9 * fT_ + 0.85) * fyB_,
                        -1.0 + (0.5 - 0.5 * fT_) * fyB_
                    ])

            elif condition_ in ["1S", "2S"]:

                if y_ >= 0:

                    kp_fxL_array_ = np.array([
                        0.0,
                        0.2,
                        0.4,
                        0.5,
                        0.7,
                        0.9,
                        1.0
                    ])

                    kp_array_ = np.array([
                        1.0,
                        1.0 + (0.75 - 1.5 * fT_) * fyB_,
                        -1.0 + (1.75 - 0.5 * fT_) * fyB_,
                        -1.0 + (1.75 - 0.5 * fT_) * fyB_,
                        -0.1 + (0.25 - 0.3 * fT_) * fyB_,
                        0.8 - (0.9 * fT_ + 0.85) * fyB_,
                        -1.0 + (0.5 - 0.5 * fT_) * fyB_
                    ])

                else:  # y < 0

                    kp_fxL_array_ = np.array([
                        0.0,
                        0.2,
                        0.4,
                        0.5,
                        0.7,
                        0.9,
                        1.0
                    ])

                    kp_array_ = np.array([
                        1.0,
                        1.0,
                        -1.0,
                        -1.0,
                        -0.1 + (1.6 * fT_ - 1.5) * fyB_,
                        0.8 + 0.2 * fyB_,
                        -1.0 + fyB_
                    ])

            else:
                raise ValueError(f"Invalid condition: {condition_}. Must be one of OST-1P, OST-2P, OST-1S, OST-2S.")

            return kp_fxL_array_, kp_array_

        def _calculate_fyz(y_, z_, TLC_, fyB_, condition_):

            if condition_ in ["1P", "2P"]:

                if y_ >= 0.0:
                    fyz = 5 * z_ / TLC_ + 3.5 * fyB_ + 0.5
                else:
                    fyz = 1.5 * z_ / TLC_ + 1.5
            elif condition_ in ["1S", "2S"]:
                if y_ >= 0.0:
                    fyz = 1.5 * z_ / TLC_ + 1.5
                else:
                    fyz = 5 * z_ / TLC_ + 3.5 * fyB_ + 0.5
            else:

                raise ValueError(f"Invalid condition: {condition_}. Must be one of OST-1P, OST-2P, OST-1S, OST-2S.")

            return fyz

        def _calculate_coefficients(x, y, z):

            fxL = x / self.L

            Bx_ = interpolate(x, self.Bx_positions, self.Bx_values)
            fyB = np.abs(2 * y) / Bx_

            if load_scenario == "Extreme Sea Loads":
                fnl = 0.8
            elif load_scenario == "Ballast Water Exchange":
                fnl = 0.9
            else:
                raise ValueError(f"Invalid load scenario: {load_scenario}")

            ka = _calculate_ka(y, fxL, fyB, self.fT, sub_case)

            kp_fxL_array_, kp_array_ = _get_kp_arrays(y, fyB, self.fT, sub_case)
            kp = interpolate(fxL, kp_fxL_array_, kp_array_)

            fyz = _calculate_fyz(y, z, self.TLC, fyB, sub_case)

            return fnl, ka, kp, fyz

        def _calculate_total_pressure(data_points):

            pressure_data_ = []

            for i, point in enumerate(data_points):

                x, y, z = point
                Bx_ = interpolate(x, self.Bx_positions, self.Bx_values)
                fnl, ka, kp, fyz = _calculate_coefficients(x, y, z)
                fnl_wl, ka_wl, kp_wl, fyz_wl = _calculate_coefficients(x, Bx_ / 2,
                                                                       self.TLC)  # Pw,wl=Pw for y=B/2 and z=TLC

                Posp = pressure_ost(fps, fnl, ka, kp, fyz, self.Cw, self.L0, lamda, self.L)

                if sub_case == "1P":  # Dynamic Load Case BSP-1S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        pressure_ost(fps, fnl_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Posp)

                elif sub_case == "2P":  # Dynamic Load Case BSP-2P

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -pressure_ost(fps, fnl_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Posp)

                elif sub_case == "1S":  # Dynamic Load Case BSP-1S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        pressure_ost(fps, fnl_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Posp)

                elif sub_case == "2S":  # Dynamic Load Case BSP-2S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -pressure_ost(fps, fnl_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Posp)

                else:
                    raise ValueError(
                        "Wrong Dynamic Load Condition. Please check condition for OST-1P, OST-2P, OST-1S, OST-2S!")

                Ps = calculate_hydrostatic_pressure(z, self.TLC, self.rho, self.g)
                Pex = calculate_pressure(Ps, Pw, design_load)  # Pex = Pw + Ps
                pressure_data_.append(np.array([x, y, z, Pex]))

            return pressure_data_

        pressure_data = _calculate_total_pressure(data_points=self.coordinates)
        pressure_data = np.array(pressure_data)  # Convert to numpy array

        self.update_pressure_data(pressure_data, f"Dynamic Load Case: OST-{sub_case} ({design_load})")  # Update last pressure data

        return pressure_data

    def _calculate_osa_case(self, sub_case: str, load_scenario: str, fps: float, design_load: str) -> np.ndarray:
        """
        Calculates Oblique Sea Equivalent Design Waves that maximize and minimize pitch acceleration with waves from the port or starboard side.

        :param sub_case: Oblique Sea case identifier, e.g., "1P", "1S", "2P", or "2S".
        :type sub_case: str
        :param load_scenario: Description of the load scenario, e.g., "Ballast Water Exchange" or "Extreme Sea Loads".
        :type load_scenario: str
        :param fps: Coefficient for strength assessments.
        :type fps: float
        :param design_load: Load scenario ("S", "D", or "S+D"). Defaults to "S+D".
        :type design_load: str

        :return: Calculated pressure based on the design load.
        :rtype: np.ndarray
        """

        lamda = 0.70 * self.L  # Wave length for OSA [m]

        def _pressure_osa(fps_, fnl_, ka_, kp_, fyz_, Cw_, L0_, lamda_, L_, fT_):
            return 0.81 * fps_ * fnl_ * ka_ * kp_ * fyz_ * Cw_ * np.sqrt((L0_ + lamda_ - 125) / L_) * (1 + 0.5 * fT_)

        def _calculate_ka(y_, fxL_, fyB_, fT_, condition_):

            A_ = 22 - 15 * fT_ + 3 * (22 * (fxL_ - 0.8) - 0.25 * (2 - fT_))

            if condition_ in ["1P", "2P"]:

                if y_ >= 0.0:

                    if fxL_ <= 0.2:
                        ka_ = 1 + 3 * (2 - fT_) * (1 - 5 * fxL_) * (1 - fyB_) * A_
                    elif 0.2 < fxL_ <= 0.5:
                        ka_ = 1.0
                    elif 0.5 < fxL_ <= 0.8:
                        ka_ = 1.0
                    else:  # fxL > 0.8
                        ka_ = 1.0 + (fxL_ - 0.8) * (1 - fyB_) * A_

                else:

                    if fxL_ <= 0.2:
                        ka_ = 1.0 + 3 * (2 - fT_) * (1 - 5 * fxL_) + ((28 * fxL_ - 5) + 3 * fT_ * (1 - 5 * fxL_)) * fyB_
                    elif 0.2 < fxL_ <= 0.5:
                        ka_ = 1.0 + (1 - 2 * fxL_) * fyB_
                    elif 0.5 < fxL_ <= 0.8:
                        ka_ = 1.0 + 1.5 * (2 * fxL_ - 1) * fyB_
                    else:  # fxL > 0.8
                        ka_ = 1.0 + (1.5 * (2 * fxL_ - 1) - (fxL_ - 0.8) * A_) * fyB_ + (fxL_ - 0.8) * A_

            elif condition_ in ["1S", "2S"]:

                if y_ >= 0.0:

                    if fxL_ <= 0.2:
                        ka_ = 1.0 + 3 * (2 - fT_) * (1 - 5 * fxL_) + ((28 * fxL_ - 5) + 3 * fT_ * (1 - 5 * fxL_)) * fyB_
                    elif 0.2 < fxL_ <= 0.5:
                        ka_ = 1.0 + (1 - 2 * fxL_) * fyB_
                    elif 0.5 < fxL_ <= 0.8:
                        ka_ = 1.0 + 1.5 * (2 * fxL_ - 1) * fyB_
                    else:  # fxL > 0.8
                        ka_ = 1.0 + (1.5 * (2 * fxL_ - 1) - (fxL_ - 0.8) * A_) * fyB_ + (fxL_ - 0.8) * A_

                else:

                    if fxL_ <= 0.2:
                        ka_ = 1.0 + 3 * (2 - fT_) * (1 - 5 * fxL_) + ((28 * fxL_ - 5) + 3 * fT_ * (1 - 5 * fxL_)) * fyB_
                    elif 0.2 < fxL_ <= 0.5:
                        ka_ = 1.0 + (1 - 2 * fxL_) * fyB_
                    elif 0.5 < fxL_ <= 0.8:
                        ka_ = 1.0 + 1.5 * (2 * fxL_ - 1) * fyB_
                    else:  # fxL > 0.8
                        ka_ = 1.0 + (1.5 * (2 * fxL_ - 1) - (fxL_ - 0.8) * A_) * fyB_ + (fxL_ - 0.8) * A_

            else:
                raise ValueError(f"Invalid condition: {condition_}. Must be one of OSA-1P, OSA-2P, OSA-1S, OSA-2S.")

            return ka_

        def _get_kp_arrays(y_, fyB_, fT_, condition_):

            if condition_ in ["1P", "2P"]:

                if y_ >= 0:

                    kp_fxL_array_ = np.array([
                        0.0,
                        0.2,
                        0.4,
                        0.5,
                        0.6,
                        0.85,
                        1.0
                    ])

                    kp_array_ = np.array([
                        0.75 - 0.5 * fyB_,
                        fT_ - 0.25 + (1.25 - fT_) * fyB_,
                        1.0,
                        1.25 - 0.5 * fT_ + (0.5 * fT_ - 0.25) * fyB_,
                        1.5 - fT_ + (fT_ - 1.07) * fyB_,
                        0.5 * fT_ - 1.25 + (0.25 - 0.5 * fT_) * fyB_,
                        0.5 * fT_ - 1.25 + (0.25 - 0.5 * fT_) * fyB_
                    ])

                else:  # y < 0

                    kp_fxL_array_ = np.array([
                        0.0,
                        0.2,
                        0.4,
                        0.5,
                        0.6,
                        0.85,
                        1.0
                    ])

                    kp_array_ = np.array([
                        0.75,
                        fT_ - 0.25 + (0.35 * fT_ - 0.47) * fyB_,
                        1.0 + (2.7 * fT_ - 3.2) * fyB_,
                        1.25 - 0.5 * fT_ + (2.7 * fT_ - 3.2) * fyB_,
                        1.5 - fT_ + (2.68 * fT_ - 3.19) * fyB_,
                        0.5 * fT_ - 1.25 + (0.2 - 0.1 * fT_) * fyB_,
                        0.5 * fT_ - 1.25 + (0.2 - 0.1 * fT_) * fyB_
                    ])

            elif condition_ in ["1S", "2S"]:

                if y_ >= 0:

                    kp_fxL_array_ = np.array([
                        0.0,
                        0.2,
                        0.4,
                        0.5,
                        0.6,
                        0.85,
                        1.0
                    ])

                    kp_array_ = np.array([
                        0.75,
                        fT_ - 0.25 + (0.35 * fT_ - 0.47) * fyB_,
                        1.0 + (2.7 * fT_ - 3.2) * fyB_,
                        1.25 - 0.5 * fT_ + (2.7 * fT_ - 3.2) * fyB_,
                        1.5 - fT_ + (2.68 * fT_ - 3.19) * fyB_,
                        0.5 * fT_ - 1.25 + (0.2 - 0.1 * fT_) * fyB_,
                        0.5 * fT_ - 1.25 + (0.2 - 0.1 * fT_) * fyB_
                    ])

                else:  # y < 0

                    kp_fxL_array_ = np.array([
                        0.0,
                        0.2,
                        0.4,
                        0.5,
                        0.6,
                        0.85,
                        1.0
                    ])

                    kp_array_ = np.array([
                        0.75 - 0.5 * fyB_,
                        fT_ - 0.25 + (1.25 - fT_) * fyB_,
                        1.0,
                        1.25 - 0.5 * fT_ + (0.5 * fT_ - 0.25) * fyB_,
                        1.5 - fT_ + (fT_ - 1.07) * fyB_,
                        0.5 * fT_ - 1.25 + (0.25 - 0.5 * fT_) * fyB_,
                        0.5 * fT_ - 1.25 + (0.25 - 0.5 * fT_) * fyB_
                    ])

            else:
                raise ValueError(f"Invalid condition: {condition_}. Must be one of OSA-1P, OSA-2P, OSA-1S, OSA-2S.")

            return kp_fxL_array_, kp_array_

        def _get_fnl_arrays(load_scenario_):

            if load_scenario_ == "Extreme Sea Loads":  # "Ballast Water Exchange"

                fnl_fxL_array_ = np.array([
                    0.0,
                    0.3,
                    0.7,
                    1.0
                ])

                fnl_array_ = np.array([
                    0.5,
                    0.8,
                    0.8,
                    0.6
                ])

            elif load_scenario_ == "Ballast Water Exchange":

                fnl_fxL_array_ = np.array([
                    0.0,
                    0.3,
                    0.7,
                    1.0
                ])

                fnl_array_ = np.array([
                    0.75,
                    0.9,
                    0.9,
                    0.8
                ])
            else:
                raise ValueError(f"Unrecognized load_scenario: {load_scenario_}")

            return fnl_fxL_array_, fnl_array_

        def _calculate_fyz(y_, z_, TLC_, fyB_, condition_):

            if condition_ in ["1P", "2P"]:

                if y_ >= 0.0:
                    fyz = 5.5 * z_ / TLC_ + 5.3 * fyB_ + 2.2
                else:
                    fyz = 0.9 * z_ / TLC_ + 0.4 * fyB_ + 2.2
            elif condition_ in ["1S", "2S"]:
                if y_ >= 0.0:
                    fyz = 0.9 * z_ / TLC_ + 0.4 * fyB_ + 2.2
                else:
                    fyz = 5.5 * z_ / TLC_ + 5.3 * fyB_ + 2.2
            else:
                raise ValueError(f"Invalid condition: {condition_}. Must be one of OSA-1P, OSA-2P, OSA-1S, OSA-2S.")

            return fyz

        def _calculate_coefficients(x, y, z):

            fxL = x / self.L

            Bx_ = interpolate(x, self.Bx_positions, self.Bx_values)
            fyB = np.abs(2 * y) / Bx_

            fnl_fxL_array_, fnl_array_ = _get_fnl_arrays(load_scenario)
            fnl = interpolate(fxL, fnl_fxL_array_, fnl_array_)

            ka = _calculate_ka(y, fxL, fyB, self.fT, sub_case)

            kp_fxL_array_, kp_array_ = _get_kp_arrays(y, fyB, self.fT, sub_case)
            kp = interpolate(fxL, kp_fxL_array_, kp_array_)

            fyz = _calculate_fyz(y, z, self.TLC, fyB, sub_case)

            return fnl, ka, kp, fyz

        def _calculate_total_pressure(data_points):
            # Pw = max(-Pfs, rho * g * (z - TLC))

            pressure_data_ = []

            for i, point in enumerate(data_points):

                x, y, z = point
                Bx_ = interpolate(x, self.Bx_positions, self.Bx_values)
                fnl, ka, kp, fyz = _calculate_coefficients(x, y, z)
                fnl_wl, ka_wl, kp_wl, fyz_wl = _calculate_coefficients(x, Bx_ / 2, self.TLC)

                Posa = _pressure_osa(fps, fnl, ka, kp, fyz, self.Cw, self.L0, lamda, self.L, self.fT)

                if sub_case == "1P":  # Dynamic Load Case BSP-1S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        _pressure_osa(fps, fnl_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L, self.fT)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Posa)

                elif sub_case == "2P":  # Dynamic Load Case BSP-2P

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -_pressure_osa(fps, fnl_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L, self.fT)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Posa)

                elif sub_case == "1S":  # Dynamic Load Case BSP-1S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        _pressure_osa(fps, fnl_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L, self.fT)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, Posa)

                elif sub_case == "2S":  # Dynamic Load Case BSP-2S

                    Pw_wl = max(
                        0.0,  # rho * g * (TLC - TLC),
                        -_pressure_osa(fps, fnl_wl, ka_wl, kp_wl, fyz_wl, self.Cw, self.L0, lamda, self.L, self.fT)
                    )

                    hw = Pw_wl / (self.rho * self.g)

                    Pw = calculate_wave_pressure(z, self.TLC, Pw_wl, self.rho, self.g, hw, -Posa)

                else:
                    raise ValueError(
                        "Wrong Dynamic Load Condition. Please check condition for OSA-1P, OSA-2P, OSA-1S, OSA-2S!")

                Ps = calculate_hydrostatic_pressure(z, self.TLC, self.rho, self.g)
                Pex = calculate_pressure(Ps, Pw, design_load)  # Pex = Pw + Ps
                pressure_data_.append(np.array([x, y, z, Pex]))

            return pressure_data_

        pressure_data = _calculate_total_pressure(data_points=self.coordinates)
        pressure_data = np.array(pressure_data)  # Convert to numpy array

        self.update_pressure_data(pressure_data, f"Dynamic Load Case: OSA-{sub_case} ({design_load})") # Update last pressure data

        return pressure_data

    def update_pressure_data(self, pressure_data, case):
        """
        Update the last pressure data with new data.

        :param pressure_data: The pressure data to update the last data with.
        :type pressure_data: np.ndarray
        :param case: Identifier or description of the case for the pressure data.
        :type case: str

        :return: None
        """

        self.last_pressure_data["data"] = pressure_data  # Save the provided pressure data last called
        self.last_pressure_data["case"] = case  # Save the provided case method last called

    def plot_last_pressure_data(self, color: str = 'coolwarm', size: int = 5, file_path: str = None,
                                file_name: str = None, show_plot: bool = True):
        """
        Plots the last calculated pressure data.

        :param color: Colormap for the plot (default is 'coolwarm').
        :type color: str
        :param size: Size of the markers in the plot (default is 5).
        :type size: int
        :param show_plot: Whether to display the plot or not (default is False).
        :rtype: False
        :param file_name: Name of the image i.e. "image.png"
        :type file_name: str
        :param file_path: Path of the plot to be saved.
        :type file_path: str

        :return: None
        """

        plot_pressure_data(
            pressure_data=self.last_pressure_data["data"],
            title=self.last_pressure_data["case"],
            color=color,
            size=size,
            file_path=file_path,
            file_name=file_name,
            show_plot=show_plot
        )


    def save_to_csv(self, path: str, file_name: str):
        """
        Save the pressure data to a CSV file.

        :param path: (str) Directory where the CSV file will be saved.
        :param file_name: (str) Name of the CSV file.

        :return: None
        """

        # Extract values and convert to a DataFrame
        x_vals = self.last_pressure_data["data"][:, 0]
        y_vals = self.last_pressure_data["data"][:, 1]
        z_vals = self.last_pressure_data["data"][:, 2]
        Pressure_vals = self.last_pressure_data["data"][:, 3]

        # Create the DataFrame and set it as an attribute
        pressure_data = pd.DataFrame({
            'x (m)': x_vals,
            'y (m)': y_vals,
            'z (m)': z_vals,
            'Pressure (N/m^2)': Pressure_vals * 1E3
        })

        full_file_path = f"{path.rstrip('/')}/{file_name}"  # Construct the full file path

        pressure_data.to_csv(full_file_path, index=False)  # Save the DataFrame to CSV

        print(f"Data exported successfully to {full_file_path}")

