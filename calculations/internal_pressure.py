import numpy as np
from calculations.coefficients import *
from calculations.pressure_calculations import *
from calculations.load_combination_factors import *
from calculations.accelerations import *
from calculations.plotting import *

class IntCargoPressureCalc:

    def __init__(self, L: float, B: float, D: float, TLC: float, TSC: float, Cb: float,
                  kr: float, GM: float, hdb: float, cargo_type: str,
                 cargo_density: float, cargo_density_effective: float, load_scenario: str, bilge_keel: bool):
        """
        A class to calculate internal cargo pressures and associated parameters for ship designs.

        :param L: Rule length of the ship (m).
        :type L: float
        :param B: Breadth moulded of the ship (m).
        :type B: float
        :param D: Depth moulded of the ship (m).
        :type D: float
        :param TLC: Loading condition draught (m).
        :type TLC: float
        :param TSC: Scantling draught (m).
        :type TSC: float
        :param Cb: Block coefficient of the ship (-).
        :type Cb: float
        :param kr: Roll radius of gyration (m).
        :type kr: float
        :param GM: Metacentric height (m).
        :type GM: float
        :param hdb: Height of the double bottom (m).
        :type hdb: float
        :param cargo_type: Type of cargo, used to determine specific parameters (e.g., "oil", "grain").
        :type cargo_type: str
        :param cargo_density: Density of the cargo (t/m^3).
        :type cargo_density: float
        :param load_scenario: Load scenario determining specific coefficients.
        :type load_scenario: str
        :param bilge_keel: Indicates if the ship has a bilge keel.
        :type bilge_keel: bool
        """

        # Set the attributes
        # self.coordinates = coordinates
        self.L = L                                  # Rule Length (m)
        self.L0 = max(110.0, self.L)                # Rule Length but not taken less than 110 (m)
        self.B = B                                  # Breadth Moulded (m)
        self.D = D                                  # Depth Moulded (m)
        self.TLC = TLC                              # Loading Condition Draught (m)
        self.TSC = TSC                              # Scantling Draught (m)
        self.Cb = Cb                                # Block Coefficient (-)
        self.kr = kr                                # Roll Radius of Gyration (m)
        self.GM = GM                                # Metacentric height (m)
        self.hdb = hdb                              # Height of Double Hull (m)
        self.rho_c = cargo_density                  # Cargo Density (t/m^3)
        self.rho_eff = cargo_density_effective      # Effective Cargo Density (t/m^3)
        self.load_scenario = load_scenario          # Load scenario determining the fp (fps)
        self.fBK = calculate_fBK(bilge_keel)        # Coefficient for Bilge Keel

        # Constants
        self.rho = 1.025                            # Seawater density (t/m^3)
        self.g = 9.81                               # Gravity Acceleration (m/s^2)

        # Internal Calculations
        self.fT = max(0.5, self.TLC / self.TSC)         # Ratio of Loading to Scantling Draught
        self.Cw = calculate_wave_coefficient(self.L)    # Wave coefficient (-)
        self.R = min(                                   # Vertical coordinate of ship rotation centre (m)
            D / 4 + TLC / 2,
            D / 2
        )

        # Angle of Repose of Cargo
        self.psi = calculate_angle_psi(cargo_type)  # Angle ψ (rad)

        # Coefficient fp (fps)
        self.fps = calculate_fps_coefficient(load_scenario)

        # Pitch (Τφ) and Roll (Τθ) periods
        self.T_thita = calculate_roll_period(self.kr, self.g, self.GM)  # Pitch Τθ period (s)
        self.T_phi = calculate_pitch_period(self.fT, self.L, self.g)    # Roll Τφ period (s)

        # Pitch (φ) and Roll (θ) angles
        self.phi = calculate_pitch_angle(self.fps, self.L, self.g)                      # Pitch φ angle (rad)
        self.thita = calculate_roll_angle(self.T_thita, self.fps, self.fBK, self.B)     # Pitch θ angle (rad)

        # Pitch (φ) and Roll (θ) angles with fp = 1
        fp = 1.0  # Will be used for calculating a_pitch and a_roll
        self.phi_fp = calculate_pitch_angle(fp, self.L, self.g)                     # Pitch φ angle (rad) for fp=1
        self.thita_fp = calculate_roll_angle(self.T_thita, fp, self.fBK, self.B)    # Pitch θ angle (rad) for fp=1

        # Accelerations
        self.a0 = calculate_a0(self.Cb, self.L)  # Acceleration a0 (m/s^2)
        self.a_surge = calculate_acceleration_surge(self.fps, self.a0, self.g)  # Acceleration surge (m/s^2)
        self.a_sway = calculate_acceleration_sway(self.fps, self.a0, self.g)    # Acceleration sway (m/s^2)
        self.a_heave = calculate_acceleration_heave(self.fps, self.a0, self.g)  # Acceleration heave (m/s^2)
        self.a_pitch = calculate_acceleration_pitch(self.fps, self.g, self.L, self.phi, self.T_phi) # Acceleration pitch (rad/s^2)
        self.a_roll = calculate_acceleration_roll(self.fps, self.thita, self.T_thita) # Acceleration roll (rad/s^2)

        # Save the results with condition
        self.last_pressure_data = {
            "data": np.array([]),
            "case": ""
        }

        self.all_pressure_data = np.array([])

    def print_data_summary(self):
        """
        Print the dimensions and important attributes of the load analysis in a table format with units in parentheses.
        """
        print("=" * 70)
        print(f"{'INTERNAL CARGO PRESSURE: DATA SUMMARY':^70}")
        print("=" * 70)
        print(f"{'Parameter':<40}{'Unit':<10}{'Value':>20}")
        print("-" * 70)
        print(f"{'Rule Length (L)':<40}{'(m)':<10}{f'{self.L:.3f}':>20}")
        print(f"{'Breadth (B)':<40}{'(m)':<10}{f'{self.B:.3f}':>20}")
        print(f"{'Depth (D)':<40}{'(m)':<10}{f'{self.D:.3f}':>20}")
        print(f"{'Loading Condition Draught (TLC)':<40}{'(m)':<10}{f'{self.TLC:.3f}':>20}")
        print(f"{'Scantling Draught (TSC)':<40}{'(m)':<10}{f'{self.TSC:.3f}':>20}")
        print(f"{'Block Coefficient (Cb)':<40}{'(-)':<10}{f'{self.Cb:.3f}':>20}")
        print(f"{'Roll Radius of Gyration (kr)':<40}{'(m)':<10}{f'{self.kr:.3f}':>20}")
        print(f"{'Metacentric Height (GM)':<40}{'(m)':<10}{f'{self.GM:.3f}':>20}")
        print(f"{'Height of Double Hull (hdb)':<40}{'(m)':<10}{f'{self.hdb:.3f}':>20}")
        print(f"{'Seawater Density (ρ)':<40}{'(t/m³)':<10}{f'{self.rho:.3f}':>20}")
        print(f"{'Cargo Density (ρc)':<40}{'(t/m³)':<10}{f'{self.rho_c:.3f}':>20}")
        print(f"{'Effective Cargo Density (ρeff)':<40}{'(t/m³)':<10}{f'{self.rho_eff:.3f}':>20}")
        print(f"{'Gravity Acceleration (g)':<40}{'(m/s²)':<10}{f'{self.g:.2f}':>20}")
        print(f"{'Wave Coefficient (Cw)':<40}{'(-)':<10}{f'{self.Cw:.3f}':>20}")
        print(f"{'Loading/Scantling Draught Ratio (fT)':<40}{'(-)':<10}{f'{self.fT:.3f}':>20}")
        print(f"{'Angle ψ (Repose of Cargo)':<40}{'(°)':<10}{f'{np.degrees(self.psi):.2f}':>20}")
        print(f"{'Pitch Angle (φ)':<40}{'(°)':<10}{f'{np.degrees(self.phi):.2f}':>20}")
        print(f"{'Roll Angle (θ)':<40}{'(°)':<10}{f'{np.degrees(self.thita):.2f}':>20}")
        print(f"{'Surge Acceleration (a_surge)':<40}{'(m/s²)':<10}{f'{self.a_surge:.3f}':>20}")
        print(f"{'Sway Acceleration (a_sway)':<40}{'(m/s²)':<10}{f'{self.a_sway:.3f}':>20}")
        print(f"{'Heave Acceleration (a_heave)':<40}{'(m/s²)':<10}{f'{self.a_heave:.3f}':>20}")
        print(f"{'Pitch Acceleration (a_pitch)':<40}{'(rad/s²)':<10}{f'{self.a_pitch:.3f}':>20}")
        print(f"{'Roll Acceleration (a_roll)':<40}{'(rad/s²)':<10}{f'{self.a_roll:.3f}':>20}")
        print("=" * 70)


    def calculate_load_case_full(self, base_case: str, sub_case: str, coordinates: np.ndarray, fb: float, center_of_gravity: np.ndarray,
                                 hHPU: float, S0: float, Bh: float, Vhc: float, lh: float, angle_alpha: np.ndarray | float,
                                 shear_load_hopper: bool = False, shear_load_inner_bottom: bool = False,
                                  fdc: float = 1.0, design_load: str = "S+D") -> np.ndarray:

        """
        Calculates the load case based on the specified base case and subcase for full loading of cargo hold.

        :param base_case: The base case identifier. Options include: "HSM", "FSM", "HSA",
                          "BSR", "BSP", "OST", "OSA".
        :type base_case: str
        :param sub_case: The subcase identifier used to refine the load calculation. Options include: "1", "2", "1P",
                         "1S", "2P", "2S".
        :type sub_case: str
        :param coordinates: Array of coordinates representing the geometry (m).
        :type coordinates: np.ndarray
        :param fb: Heading correction factor.
        :type fb: float
        :param center_of_gravity: The center of gravity coordinates (xG, yG, zG) of the cargo hold (m).
        :type center_of_gravity: np.ndarray
        :param hHPU: Height hHPU (m).
        :type hHPU: float
        :param S0: Area S0 (m²).
        :type S0: float
        :param Bh: Beam height (m).
        :type Bh: float
        :param Vhc: Hatch Coaming volume (m³).
        :type Vhc: float
        :param lh: Cargo length (m).
        :type lh: float
        :param angle_alpha: Angle a (deg).
        :type angle_alpha: np.ndarray | float
        :param shear_load_hopper: Whether to include shear load (for the hopper tank or lower stool only). Default is False.
        :type shear_load_hopper: bool
        :param shear_load_inner_bottom: Whether to include shear load (for the inner bottom only). Default is False.
        :type shear_load_inner_bottom: bool
        :param fdc: Dry Cargo factor. Default is 1.0.
        :type fdc: float
        :param design_load: Load scenario ("S", "D", or "S+D"). Defaults to "S+D".
        :type design_load: str

        :return: The result of the load case calculation, specific to the base case and subcase.
        :rtype: np.ndarray
        """

        def calculate_total_pressure(data_points):
            pressure_data_ = []

            for i, point in enumerate(data_points):
                x, y, z = point

                alpha = get_angle_for_z(z, angle_alpha)
                Kc = calculate_Kc_coefficient(alpha, self.psi)

                xG = center_of_gravity[0]  # x-coordinate for center of gravity
                flp = calculate_flp_coefficient(xG, self.L)
                flp_ost = calculate_flp_ost(xG, self.L)
                flp_osa = calculate_flp_osa(xG, self.L, self.fT)
                lcf = calculate_lcf(base_case, sub_case, flp, self.fT, flp_ost, flp_osa)

                hc = calculate_hc_full(hHPU, S0, Bh, Vhc, lh)  # hc for full load

                zc = self.hdb + hc

                ax = calculate_acceleration_x(lcf["c_xg"], self.g, self.phi, lcf["c_xs"],  # ax
                                              self.a_surge, lcf["c_xp"], self.a_pitch, z, self.R)
                ay = calculate_acceleration_y(lcf["c_yg"], self.g, self.thita, lcf["c_ys"],  # ay
                                              self.a_sway, lcf["c_yr"], self.a_roll, z, self.R)
                az = calculate_acceleration_z(lcf["c_zh"], self.a_heave, lcf["c_zr"],  # az
                                              self.a_roll, y, lcf["c_zp"], self.a_pitch, x, self.L)

                a = np.array([ax, ay, az])

                Pbs = calculate_static_cargo_pressure(self.rho_c, self.g, Kc, zc, z)  # Static Pressure

                Pbd = calculate_dynamic_cargo_pressure(fb, fdc, Kc, self.rho_c, a,  # Dynamic Pressure
                                                       center_of_gravity, point, zc)

                if shear_load_hopper:
                    Pbs_s = calculate_static_shear_load_on_hopper(self.rho_eff, self.g, Kc, zc,
                                          z, alpha)
                    Pbs_d = calculate_dynamic_shear_load_on_hopper(fb, az, self.rho_eff, Kc, zc,
                    z, alpha)
                    Pbs = Pbs + Pbs_s
                    Pbd = Pbd + Pbs_d

                if shear_load_inner_bottom:
                    Pbs_dx, Pbs_dy = calculate_dynamic_shear_load_on_inner_bottom(fb, ax, ay, self.rho_eff, hc)
                    Pbd = Pbd + Pbs_dx + Pbs_dy

                Pin = calculate_pressure(Pbs, Pbd, design_load)  # Pin = Pbs + Pbd
                pressure_data_.append(np.array([x, y, z, Pin]))

            return pressure_data_

        pressure_data = calculate_total_pressure(data_points=coordinates)
        pressure_data = np.array(pressure_data)  # Convert to numpy array

        self.update_pressure_data(pressure_data,
                                  f"Dynamic Load Case: {base_case}-{sub_case} ({design_load}) for Full Load")

        return pressure_data

    def calculate_load_case_partial(self, base_case: str, sub_case: str, coordinates: np.ndarray, fb: float, center_of_gravity: np.array,
                                    hHPL: float, Bh: float, Bib: float, M:float,
                                    Vts: float, lh: float, angle_alpha: np.array,
                                    shear_load_hopper: bool = False, shear_load_inner_bottom: bool = False,
                                    fdc=1.0, design_load="S+D"):

        """
        Calculates the load case based on the specified base case and subcase for partial loading of cargo hold.

        :param base_case: The base case identifier. Options include: "HSM", "FSM", "HSA",
                          "BSR", "BSP", "OST", "OSA".
        :type base_case: str
        :param sub_case: The subcase identifier used to refine the load calculation. Options include: "1", "2", "1P",
                         "1S", "2P", "2S".
        :type sub_case: str
        :param coordinates: Array of coordinates representing the geometry (m).
        :type coordinates: np.ndarray
        :param fb: Heading correction factor.
        :type fb: float
        :param center_of_gravity: The center of gravity coordinates (xG, yG, zG) of the cargo hold (m).
        :type center_of_gravity: np.ndarray
        :param hHPL: The vertical distance hHPL(m).
        :type hHPL: float
        :param Bh: The breadth of the cargo hold (m).
        :type Bh: float
        :param Bib: The breadth of the inner bottom (m).
        :type Bib: float
        :param M: The mass of the bulk cargo being considered (t).
        :type M: float
        :param Vts: The total volume of the lower bulkhead stools within the cargo hold length
                    and inboard of the hopper tanks (m³).
        :type Vts: float
        :param lh: The length of the cargo hold (m).
        :type lh: float
        :param angle_alpha: An array or single value representing the inclination angles (deg)
                            for the hopper side slope or similar components.
        :type angle_alpha: np.ndarray | float
        :param shear_load_hopper: Whether to include shear load (for the hopper tank or lower stool only). Default is False.
        :type shear_load_hopper: bool
        :param shear_load_inner_bottom: Whether to include shear load (for the inner bottom only). Default is False.
        :type shear_load_inner_bottom: bool
        :param design_load: Load scenario ("S", "D", or "S+D"). Defaults to "S+D".
        :type design_load: str
        :param fdc: Dry Cargo factor. Default is 1.0.
        :type fdc: float

        :return: The result of the load case calculation, specific to the base case and subcase.
        :rtype: np.ndarray
        """

        def calculate_total_pressure(data_points):

            pressure_data_ = []

            for i, point in enumerate(data_points):
                x, y, z = point

                alpha = get_angle_for_z(z, angle_alpha)
                Kc = calculate_Kc_coefficient(alpha, self.psi)

                xG = center_of_gravity[0]  # x-coordinate for center of gravity
                flp = calculate_flp_coefficient(xG, self.L)
                flp_ost = calculate_flp_ost(xG, self.L)
                flp_osa = calculate_flp_osa(xG, self.L, self.fT)
                lcf = calculate_lcf(base_case, sub_case, flp, self.fT, flp_ost, flp_osa)

                hc, hc_cl = calculate_hc_partial(Bh, Bib, self.psi, y, M, self.rho_c, lh, hHPL, Vts)  # hc partial load

                center_of_gravity[2] = self.hdb + hc_cl / 2  # zG

                zc = self.hdb + hc

                ax = calculate_acceleration_x(lcf["c_xg"], self.g, self.phi, lcf["c_xs"],  # ax
                                              self.a_surge, lcf["c_xp"], self.a_pitch, z, self.R)
                ay = calculate_acceleration_y(lcf["c_yg"], self.g, self.thita, lcf["c_ys"],  # ay
                                              self.a_sway, lcf["c_yr"], self.a_roll, z, self.R)
                az = calculate_acceleration_z(lcf["c_zh"], self.a_heave, lcf["c_zr"],  # az
                                              self.a_roll, y, lcf["c_zp"], self.a_pitch, x, self.L)

                a = np.array([ax, ay, az])

                Pbs = calculate_static_cargo_pressure(self.rho_c, self.g, Kc, zc, z) # Static Pressure

                Pbd = calculate_dynamic_cargo_pressure(fb, fdc, Kc, self.rho_c, a,   # Dynamic Pressure
                                                       center_of_gravity, point, zc)

                if shear_load_hopper:
                    Pbs_s = calculate_static_shear_load_on_hopper(self.rho_eff, self.g, Kc, zc, z, alpha)
                    Pbs_d = calculate_dynamic_shear_load_on_hopper(fb, az, self.rho_eff, Kc, zc, z, alpha)
                    Pbs = Pbs + Pbs_s
                    Pbd = Pbd + Pbs_d

                if shear_load_inner_bottom:
                    Pbs_dx, Pbs_dy = calculate_dynamic_shear_load_on_inner_bottom(fb, ax, ay, self.rho_eff, hc)
                    Pbd = Pbd + Pbs_dx + Pbs_dy

                Pin = calculate_pressure(Pbs, Pbd, design_load)  # Pin = Pbs + Pbd
                pressure_data_.append(np.array([x, y, z, Pin]))

            return pressure_data_

        pressure_data = calculate_total_pressure(data_points=coordinates)
        pressure_data = np.array(pressure_data)  # Convert to numpy array

        self.update_pressure_data(pressure_data,
                                  f"Dynamic Load Case: {base_case}-{sub_case} ({design_load}) for Partial Load")

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

    def add_pressure_data(self):
        if self.all_pressure_data.size == 0:  # Check if array is empty
            self.all_pressure_data = self.last_pressure_data["data"]
        else:
            self.all_pressure_data = np.vstack((self.all_pressure_data, self.last_pressure_data["data"]))

    def clear_all_pressure_data(self):

        self.all_pressure_data = np.array([])


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

    def plot_all_pressure_data(self, color: str = 'coolwarm', size: int = 5, file_path: str = None,
                                file_name: str = None, show_plot: bool = True, title: str = ""):
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
            pressure_data=self.all_pressure_data,
            title=title,
            color=color,
            size=size,
            file_path=file_path,
            file_name=file_name,
            show_plot=show_plot
        )




