import numpy as np

def calculate_wave_coefficient(L_: float):
    """
    Calculate the wave coefficient based on the given length.

    :param L_: Rule length of the ship (m).
    :type L_: float

    :return: The calculated wave coefficient.
    :rtype: float
    """
    if 90.0 <= L_ <= 300.0:
        wave_coefficient = 10.75 - ((300 - L_) / 100) ** 1.5
    elif 300.0 < L_ <= 350.0:
        wave_coefficient = 10.75
    elif 350.0 < L_ <= 500.0:
        wave_coefficient = 10.75 - ((L_ - 350) / 150) ** 1.5
    else:
        raise ValueError("Length L is out of bounds. Valid range is [90, 500].")

    return wave_coefficient

def calculate_fBK(bilge_keel_: bool) -> float:
    """
    Calculates the fBK coefficient based on the presence of bilge keels.

    :param bilge_keel_: Indicates whether bilge keels are present (True) or not (False).
    :type bilge_keel_: bool

    :return: The calculated fBK coefficient.
    :rtype: float
    """
    if bilge_keel_:
        fBK_ = 1.0
    elif not bilge_keel_:
        fBK_ = 1.2
    else:
        raise ValueError("Invalid value for bilge_keel. Expected True or False.")
    return fBK_


def calculate_flp_coefficient(x_: float , L_: float):
    """
    Calculate the flp coefficient.

    :param x_: The x-coordinate of the point.
    :type x_: float
    :param L_: The rule length.
    :type L_: float

    :return: The flp coefficient.
    :rtype: float
    """

    if x_ / L_ <= 0.5:
        flp_ = 1.0
    else:
        flp_ = -1.0

    return flp_

def calculate_flp_ost(x_: float, L_: float):
    """
    Calculate the flp_ost coefficient.

    :param x_: The x-coordinate of the point.
    :type x_: float
    :param L_: The rule length.
    :type L_: float

    :return: The flp_ost coefficient.
    :rtype: float
    """

    fxL_ = x_ / L_

    if x_ / L_ < 0.2:
        flp_ost_ = 5 * fxL_
    elif 0.2 <= x_ / L_ < 0.4:
        flp_ost_ = 1.0
    elif 0.4 <= x_ / L_ < 0.65:
        flp_ost_ = -7.6 * fxL_ + 4.04
    elif 0.65 <= x_ / L_ < 0.85:
        flp_ost_ = -0.9
    else:
        flp_ost_ = 6 * fxL_ - 6

    return flp_ost_


def calculate_flp_osa(x_: float, L_: float, fT_: float):
    """
    Calculate the flp_osa coefficient.

    :param x_: The x-coordinate of the point.
    :type x_: float
    :param L_: The rule length.
    :type L_: float
    :param fT_: The draught loading to scantling ratio.
    :type fT_: float

    :return: The flp_osa coefficient.
    :rtype: float
    """
    fxL_ = x_ / L_

    if x_ / L_ < 0.4:
        flp_osa_ = -(0.2 + 0.3 * fT_)
    elif 0.4 <= x_ / L_ < 0.6:
        flp_osa_ = -(0.2 + 0.3 * fT_) * (5.6 - 11.5 * fxL_)
    else:
        flp_osa_ = 1.3 * (0.2 + 0.3 * fT_)

    return flp_osa_

def calculate_fps_coefficient(load_scenario_: str):
    """
    Calculates the Strength assessment coefficient fps for different load scenarios.

    :param load_scenario_: Load scenario (Allowed values are "Extreme Sea Loads", "Water Exchange", "Accidental Flooded", "Harbour / Sheltered Water").
    :type load_scenario_: str

    :return: The fps coefficient corresponding to the specified load scenario.
    :rtype: float
    """

    if load_scenario_ == "Extreme Sea Loads":
        fps = 1.0
    elif load_scenario_ == "Water Exchange":
        fps = 0.8
    elif load_scenario_ == "Accidental Flooded":
        fps = 0.8
    elif load_scenario_ == "Harbour / Sheltered Water":
        fps = 0.4
    else:
        raise ValueError(f"Invalid load scenario: {load_scenario_}")

    return fps

def calculate_Kc_coefficient(alpha_: float, psi_: float):
    """
    Calculates the Kc coefficient based on the angles α (alpha) and ψ (psi).

    :param alpha_: The angle α (rad).
    :type alpha_: float
    :param psi_: The angle ψ (rad).
    :type psi_: float

    :return: The calculated Kc coefficient.
    :rtype: float
    """

    if alpha_ <= (90 * np.pi / 180.0): # a <= 90 deg

        Kc_ = np.cos(alpha_) ** 2 + (1 - np.sin(psi_)) * np.sin(alpha_) ** 2

    elif (90 * np.pi / 180.0) < alpha_ <= (120 * np.pi / 180.0): # 90 deg < a <= 120 deg

        Kc_ = (1 - np.sin(psi_)) * np.sin(alpha_) ** 2

    elif alpha_ > (120 * np.pi / 180.0) and (alpha_ + psi_) < (180 * np.pi / 180.0):  # a > 120 deg + a + ψ < 180 deg

        Kc_ = 0.75 * (1 - np.sin(psi_)) * (1 - (alpha_ - 120.0 * np.pi / 180.0) / (60 * np.pi / 180.0 - psi_))

    else:   # a + ψ >= 180 deg

        Kc_ = 0.0

    return Kc_