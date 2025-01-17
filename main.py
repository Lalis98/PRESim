from calculations.external_pressure import *
from calculations.internal_pressure import *
from utils.utils import *

def main():

    # ==================================================================================================================
    #                                          CONFIGURATION
    # ==================================================================================================================

    # Configurations for LaTeX rendering
    configure_latex_font(
        use_latex=True
    )

    # ==================================================================================================================
    #                                         IMPORT GEOMETRIES
    # ==================================================================================================================

    # Hull Geometry
    file_path = "geometry_data/outter_hull/outter_hull.csv"
    coordinates = np.loadtxt(file_path, delimiter=',', skiprows=1)    # Load the CSV file into a 2D numpy array

    # Cargo Hold Geometry
    file_path_2 = "geometry_data/cargo_hold/"
    coordinates_bulk = np.loadtxt(file_path_2 + "bulkheads.csv", delimiter=',', skiprows=1)
    coordinates_hopper = np.loadtxt(file_path_2 + "hopper_tank.csv", delimiter=',', skiprows=1)
    coordinates_inner = np.loadtxt(file_path_2 + "inner_bottom.csv", delimiter=',', skiprows=1)
    coordinates_stool = np.loadtxt(file_path_2 + "lower_stool.csv", delimiter=',', skiprows=1)
    coordinates_side = np.loadtxt(file_path_2 + "side_shell.csv", delimiter=',', skiprows=1)

    # ==================================================================================================================
    #                                                DATA
    # ==================================================================================================================

    rule_length = 218.372                   # Rule Length (according to CSR) (m)
    breadth = 32.240                        # Breadth Moulded (m)
    depth = 20.200                          # Depth Moulded (m)
    double_hull = 1.74                      # Double Hull height from baseline (m)
    loading_draught = 14.555                # Loading Condition Draught (m)
    scantling_draught = 14.555              # Scantling Draught (m)
    block_coefficient = 0.797               # Block Coefficient (-)
    breadth_distribution = np.array([       # Distribution of Breadth at waterline relative to x-coordinate
            [82.27, 158.71],                # x-coordinates along the ship's length (m)
            [32.24, 32.24]                  # Corresponding breadth values at those positions at the waterline (m)
        ])                                  # (2, N): 2 rows (x-coordinates and breadth), N columns (positions NUMBER)
    roll_radius = 0.4 * breadth             # Roll Radius of Gyration (m)
    metacentric_height = 0.2 * breadth      # Metacentric height (m)
    cog = np.array([95.62, 0.0, 11.19]) # 95.62 / 145.36      # Center of Gravity of Cargo Hold (x, y, z) (m)

    angle_distribution_side = np.array([
        [5.68, 15.20, 90.0],                # 5.68 m < z < 15.20 m, angle is 90.0 deg
        [15.20, 23.90, 147.688]             # 15.20 m < z < 23.90 m, angle is 147.0 deg
    ])

    hHPU = 13.46                            # Height hHPU (m)
    S0 = 125.840                            # Area S0 (m²)
    Bh = 32.240                             # Breadth of Cargo Hold (m)
    Vhc = 360.78                            # Volume of Hatch Coaming (m³)
    lh = 25.480                             # Length of Cargo (m)

    hHPL = 3.95                             # Height hHPL (m)
    Bib = 23.80                             # Breadth of Inner Bottom (m)
    M = 21_401.5                            # Considered Mass of the Cargo (t)
    Vts = 700.0                             # Volume of lower stools and hopper tanks within lh (m³)

    # ==================================================================================================================
    #                                    INTERNAL CARGO PRESSURES - FULL LOADING
    # ==================================================================================================================

    load_cases = [
        ["HSM", "1", 1.05],
        # ["HSM", "2", 1.05],
        # ["FSM", "1", 1.05],
        # ["FSM", "2", 1.05],
        # ["HSA", "1", 1.0],
        # ["HSA", "2", 1.0],
        # ["BSR", "1P", 0.8],
        # ["BSR", "1S", 0.8],
        # ["BSR", "2P", 0.8],
        # ["BSR", "2S", 0.8],
        # ["BSP", "1P", 0.8],
        # ["BSP", "1S", 0.8],
        # ["BSP", "2P", 0.8],
        # ["BSP", "2S", 0.8],
        # ["OST", "1P", 1.0],
        # ["OST", "1S", 1.0],
        # ["OST", "2P", 1.0],
        # ["OST", "2S", 1.0],
        # ["OSA", "1P", 1.0],
        # ["OSA", "1S", 1.0],
        # ["OSA", "2P", 1.0],
        # ["OSA", "2S", 1.0],
    ]

    comp_full_load = [
        [coordinates_bulk, cog, hHPU, S0, Bh, Vhc, lh, 90.0, False, False, 1.0, "S+D"],
        [coordinates_hopper, cog, hHPU, S0, Bh, Vhc, lh, 51.249, True, False, 1.0, "S+D"],
        [coordinates_inner, cog, hHPU, S0, Bh, Vhc, lh, 0.0, False, True, 1.0, "S+D"],
        [coordinates_stool, cog, hHPU, S0, Bh, Vhc, lh, 48.7, True, False, 1.0, "S+D"],
        [coordinates_side, cog, hHPU, S0, Bh, Vhc, lh, angle_distribution_side, False, False, 1.0, "S+D"]
    ]

    internal_analysis = IntCargoPressureCalc(
        L=rule_length,
        B=breadth,
        D=depth,
        TLC=loading_draught,
        TSC=scantling_draught,
        Cb=block_coefficient,
        kr=roll_radius,
        GM=metacentric_height,
        hdb=double_hull,
        cargo_type="General",  # Allowed values are: 'General', 'Iron Ore', or 'Cement'
        cargo_density=1.452,
        cargo_density_effective=1.535,
        load_scenario="Extreme Sea Loads",  # Allowed values are: "Extreme Sea Loads", "Water Exchange", "Accidental Flooded", "Harbour / Sheltered Water"
        bilge_keel=True
    )

    process_load_cases_full(internal_analysis, load_cases, comp_full_load)

    # ==================================================================================================================
    #                                 INTERNAL CARGO PRESSURES - PARTIAL LOADING
    # ==================================================================================================================

    internal_analysis = IntCargoPressureCalc(
        L=rule_length,
        B=breadth,
        D=depth,
        TLC=loading_draught,
        TSC=scantling_draught,
        Cb=block_coefficient,
        kr=roll_radius,
        GM=metacentric_height,
        hdb=double_hull,
        cargo_type="General",
        cargo_density=3.0,
        cargo_density_effective=3.17,
        load_scenario="Extreme Sea Loads",
        bilge_keel=True
    )

    internal_analysis.print_data_summary()

    comp_partial_load = [
        [coordinates_bulk, cog, hHPL, Bh, Bib, M, Vts, lh, 90.0, False, False, 1.0, "S+D"],
        [coordinates_hopper, cog, hHPL, Bh, Bib, M, Vts, lh, 51.249, True, False, 1.0, "S+D"],
        [coordinates_inner, cog, hHPL, Bh, Bib, M, Vts, lh, 0.0, False, True, 1.0, "S+D"],
        [coordinates_stool, cog, hHPL, Bh, Bib, M, Vts, lh, 48.7, True, False, 1.0, "S+D"],
        [coordinates_side, cog, hHPL, Bh, Bib, M, Vts, lh, angle_distribution_side, False, False, 1.0, "S+D"]
    ]

    # ==================================================================================================================
    #                                         EXTERNAL SEA PRESSURES
    # ==================================================================================================================

    load_cases = [
        ["HSM", "1", 1.0, "Extreme Sea Loads", 1.05, "S+D"],
        ["HSM", "2", 1.0, "Extreme Sea Loads", 1.05, "S+D"],
        ["FSM", "1", 1.0, "Extreme Sea Loads", 1.05, "S+D"],
        ["FSM", "2", 1.0, "Extreme Sea Loads", 1.05, "S+D"],
        ["HSA", "1", 1.0, "Extreme Sea Loads", None, "S+D"],
        ["HSA", "2", 1.0, "Extreme Sea Loads", None, "S+D"],
        ["BSR", "1P", 1.0, None, 0.8, "S+D"],
        ["BSR", "1S", 1.0, None, 0.8, "S+D"],
        ["BSR", "2P", 1.0, None, 0.8, "S+D"],
        ["BSR", "2S", 1.0, None, 0.8, "S+D"],
        ["BSP", "1P", 1.0, "Extreme Sea Loads", 0.8, "S+D"],
        ["BSP", "1S", 1.0, "Extreme Sea Loads", 0.8, "S+D"],
        ["BSP", "2P", 1.0, "Extreme Sea Loads", 0.8, "S+D"],
        ["BSP", "2S", 1.0, "Extreme Sea Loads", 0.8, "S+D"],
        ["OST", "1P", 1.0, "Extreme Sea Loads", None, "S+D"],
        ["OST", "1S", 1.0, "Extreme Sea Loads", None, "S+D"],
        ["OST", "2P", 1.0, "Extreme Sea Loads", None, "S+D"],
        ["OST", "2S", 1.0, "Extreme Sea Loads", None, "S+D"],
        ["OSA", "1P", 1.0, "Extreme Sea Loads", None, "S+D"],
        ["OSA", "1S", 1.0, "Extreme Sea Loads", None, "S+D"],
        ["OSA", "2P", 1.0, "Extreme Sea Loads", None, "S+D"],
        ["OSA", "2S", 1.0, "Extreme Sea Loads", None, "S+D"],
    ]

    external_analysis = ExternalSeaPressureCalc(        # Call class "ExternalSeaPressureCalc" and initialize the attributes
        coordinates=coordinates,
        L=rule_length,
        B=breadth,
        TLC=loading_draught,
        TSC=scantling_draught,
        Cb=block_coefficient,
        Bx_distribution= breadth_distribution,
        kr=roll_radius,
        GM=metacentric_height,
        bilge_keel=True
    )

    external_analysis.print_data_summary()

    process_load_cases_external(
        external_analysis,
        load_cases,
        color='coolwarm',
        # file_path="C:/Users/user/Downloads/figures",
        size=5,
        show_plot=True)


if __name__ == "__main__":
    main()

