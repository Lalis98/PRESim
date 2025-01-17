import matplotlib.pyplot as plt

def process_load_cases_full(internal_analysis, load_cases, load_cases_full_load,
                       show_plot: bool = True, file_path: str = None, file_name: str = None,
                       color: str = 'coolwarm', size: int = 5):
    """
    Processes all load cases and calculates full loads, then generates and plots the pressure data.

    :param internal_analysis: The internal analysis object to perform calculations.
    :param load_cases: List of load case names and their parameters.
    :param load_cases_full_load: List of full load cases with necessary parameters.
    :param show_plot: Whether to display the plot after calculation (default is True).
    :param file_path: Path where the plot should be saved (default is None).
    :param file_name: Name of the file to save the plot (default is None).
    :param color: Colormap for the plot (default is 'coolwarm').
    :param size: Size of the markers in the plot (default is 5).

    :return: None
    """
    design_load = None

    for name_case in load_cases:
        for load_case in load_cases_full_load:
            internal_analysis.calculate_load_case_full(
                base_case=name_case[0],
                sub_case=name_case[1],
                coordinates=load_case[0],
                fb=name_case[2],
                center_of_gravity=load_case[1],
                hHPU=load_case[2],
                S0=load_case[3],
                Bh=load_case[4],
                Vhc=load_case[5],
                lh=load_case[6],
                angle_alpha=load_case[7],
                shear_load_hopper=load_case[8],
                shear_load_inner_bottom=load_case[9],
                fdc=load_case[10],
                design_load=load_case[11]
            )
            design_load = load_case[11]
            internal_analysis.add_pressure_data()

        internal_analysis.plot_all_pressure_data(
            title=f"Dynamic Load Case: {name_case[0]}-{name_case[1]} ({design_load})",
            show_plot=show_plot,
            file_path=file_path,
            file_name=file_name,
            color=color,
            size=size
        )
        internal_analysis.clear_all_pressure_data()


def process_load_cases_partial(internal_analysis, load_cases, load_cases_partial_load,
                       show_plot: bool = True, file_path: str = None, file_name: str = None,
                       color: str = 'coolwarm', size: int = 5):
    """
    Processes all load cases and calculates partial loads, then generates and plots the pressure data.

    :param internal_analysis: The internal analysis object to perform calculations.
    :param load_cases: List of load case names and their parameters.
    :param load_cases_partial_load: List of partial load cases with necessary parameters.
    :param show_plot: Whether to display the plot after calculation (default is True).
    :param file_path: Path where the plot should be saved (default is None).
    :param file_name: Name of the file to save the plot (default is None).
    :param color: Colormap for the plot (default is 'coolwarm').
    :param size: Size of the markers in the plot (default is 5).

    :return: None
    """
    design_load = None

    for name_case in load_cases:
        for load_case in load_cases_partial_load:
            internal_analysis.calculate_load_case_partial(
                base_case=name_case[0],
                sub_case=name_case[1],
                coordinates=load_case[0],
                fb=name_case[2],
                center_of_gravity=load_case[1],
                hHPL=load_case[2],
                Bh=load_case[3],
                Bib=load_case[4],
                M=load_case[5],
                Vts=load_case[6],
                lh=load_case[7],
                angle_alpha=load_case[8],
                shear_load_hopper=load_case[9],
                shear_load_inner_bottom=load_case[10],
                fdc=load_case[11],
                design_load=load_case[12]
            )

            design_load = load_case[11]
            internal_analysis.add_pressure_data()

        # After processing all partial loads for the case
        internal_analysis.plot_all_pressure_data(
            title=f"Dynamic Load Case: {name_case[0]}-{name_case[1]} ({design_load})",
            show_plot=show_plot,
            file_path=file_path,
            file_name=file_name,
            color=color,
            size=size,
        )
        internal_analysis.clear_all_pressure_data()


def process_load_cases_external(
        external_analysis,
        load_cases,
        file_path=None,
        color='coolwarm',
        size=5,
        show_plot=True
):
    """
    Calculates load cases and plots the last pressure data for each case.

    :param external_analysis: Instance of the ExternalSeaPressureCalc class.
    :param load_cases: List of load cases to process.
    :param file_path: Directory where the plot images will be saved (default is "C:/Users/micha/Downloads/figures").
    :param color: Colormap for the plot (default is 'coolwarm').
    :param size: Size of the markers in the plot (default is 5).
    :param show_plot: Whether to display the plot or not (default is True).
    """
    for case in load_cases:
        base_case, sub_case, fps, load_scenario, fb, design_load = case

        try:
            # Calculate load case
            external_analysis.calculate_load_case(
                base_case=base_case,
                sub_case=sub_case,
                fps=fps,
                load_scenario=load_scenario,
                fb=fb,
                design_load=design_load
            )
            print(f"Load case {base_case}-{sub_case} calculated successfully.")

            external_analysis.plot_last_pressure_data(
                color=color,
                size=size,
                file_path=file_path,
                file_name=f"{base_case}-{sub_case}.png",
                show_plot=show_plot,
            )
        except Exception as e:
            print(f"Error calculating load case {base_case}-{sub_case}: {e}")

def configure_latex_font(use_latex=True):
    """
    Configures the global settings for Matplotlib font.

    Parameters:
        use_latex (bool): Enables or disables LaTeX rendering for Matplotlib.
    """
    if use_latex:
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"],
        })
    else:
        plt.rcParams.update({
            "text.usetex": False,
            "font.family": "sans-serif",
        })
