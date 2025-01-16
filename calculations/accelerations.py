import numpy as np

def calculate_a0(Cb_: float, L_: float):
    """
    Calculates the acceleration parameter a0.

    :param Cb_: The block coefficient.
    :type Cb_: float
    :param L_: The rule length of the ship (m).
    :type L_: float

    :return: The calculated acceleration parameter a0.
    :rtype: float
    """
    a0_= (1.58 - 0.47 * Cb_) * (2.4 / np.sqrt(L_) + 34 / L_ + 600 / L_ ** 2)

    return a0_

def calculate_acceleration_surge(fp_: float, a0_: float, g_: float):
    """
    Calculates the surge acceleration.

    :param fp_: The strength assessment coefficient.
    :type fp_: float
    :param a0_: The acceleration parameter a0.
    :type a0_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float

    :return: The calculated surge acceleration (m/s²).
    :rtype: float
    """
    return 0.2 * fp_ * a0_ * g_

def calculate_acceleration_sway(fp_: float, a0_: float, g_: float):
    """
    Calculates the sway acceleration.

    :param fp_: The strength assessment coefficient.
    :type fp_: float
    :param a0_: The acceleration parameter a0.
    :type a0_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float

    :return: The calculated sway acceleration (m/s²).
    :rtype: float
    """
    return 0.3 * fp_ * a0_ * g_

def calculate_acceleration_heave(fp_: float, a0_: float, g_: float):
    """
    Calculates the heave acceleration.

    :param fp_: The strength assessment coefficient.
    :type fp_: float
    :param a0_: The acceleration parameter a0.
    :type a0_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float

    :return: The calculated heave acceleration (m/s²).
    :rtype: float
    """
    return fp_ * a0_ * g_

def calculate_acceleration_pitch(fp_: float, g_: float, L_: float, phi_: float, T_phi_: float):
    """
    Calculates the pitch acceleration.

    :param fp_: The strength assessment coefficient.
    :type fp_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float
    :param L_: The rule length of the ship (m).
    :type L_: float
    :param phi_: The pitch angle φ (rad).
    :type phi_: float
    :param T_phi_: The period of pitch motion Τφ (s).
    :type T_phi_: float

    :return: The calculated pitch acceleration (rad/s²).
    :rtype: float
    """

    a_pitch_ = fp_ * (3.1 / np.sqrt(g_ * L_) + 1) * phi_ * ((2 * np.pi) / T_phi_) ** 2
    return a_pitch_

def calculate_acceleration_roll(fp_: float, thita_: float, T_thita_: float):
    """
    Calculates the roll acceleration.

    :param fp_: The strength assessment coefficient.
    :type fp_: float
    :param thita_: The roll angle θ (rad).
    :type thita_: float
    :param T_thita_: The period of roll motion Τθ (seconds).
    :type T_thita_: float

    :return: The calculated roll acceleration (rad/s²).
    :rtype: float
    """

    a_roll_ = fp_ * thita_ * (2 * np.pi / T_thita_) ** 2

    return a_roll_

def calculate_roll_angle(T_thita_: float, fp_: float, fBK_: float, B_: float):
    """
    Calculates the roll angle.

    :param T_thita_: The period of roll motion (seconds).
    :type T_thita_: float
    :param fp_: The strength assessment coefficient.
    :type fp_: float
    :param fBK_: The bilge keel coefficient (1.2 if the ship has a bilge keel, 1.0 otherwise).
    :type fBK_: float
    :param B_: The moulded breadth of the ship (m).
    :type B_: float

    :return: The calculated roll angle θ (rad).
    :rtype: float
    """

    thita_ = 9000 * (1.25 - 0.025 * T_thita_) * fp_ * fBK_ / ((B_ + 75) * np.pi)  # in degrees
    thita_rad_ = thita_ * np.pi / 180.0  # in rad

    return thita_rad_

def calculate_roll_period(kr_: float, g_: float, GM_: float):
    """
    Calculates the roll period of the ship.

    :param kr_: The roll radius of gyration (m).
    :type kr_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float
    :param GM_: The metacentric height (m).
    :type GM_: float

    :return: The calculated roll period Tθ (s).
    :rtype: float
    """
    return 2.3 * np.pi * kr_ / np.sqrt(g_ * GM_)

def calculate_pitch_angle(fp_: float, L_: float, g_: float):
    """
    Calculates the pitch angle of the ship.

    :param fp_: The strength assessment coefficient.
    :type fp_: float
    :param L_: The rule length of the ship (m).
    :type L_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float

    :return: The calculated pitch angle φ (rad).
    :rtype: float
    """

    phi_ = 1350 * fp_ * L_ ** (-0.94) * (     # in deg
        1 + (2.57 / np.sqrt(g_ * L_) ** 1.2)
    )

    phi_rad_ = phi_ * np.pi / 180.0   # in rad

    return phi_rad_

def calculate_pitch_period(fT_: float, L_: float, g_: float):
    """
    Calculates the pitch period of the ship.

    :param fT_: The draught loading to scantling ratio.
    :type fT_: float
    :param L_: The rule length of the ship (m).
    :type L_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float

    :return: The calculated pitch period Τφ (s).
    :rtype: float
    """
    lamda_phi_ = 0.6 * (1 + fT_) * L_

    T_phi_ = np.sqrt(
        (2 * np.pi * lamda_phi_) / g_
    )

    return T_phi_

def calculate_acceleration_x(c_xg_: float, g_,phi_: float, c_xs_: float, a_surge_: float,
                             c_xp_: float, a_pitch_: float, z_: float, R_: float):
    """
    Calculates the acceleration in the x-direction (ax) for a given dynamic system.

    :param c_xg_: The coefficient C_XG.
    :type c_xg_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float
    :param phi_: The roll angle φ (rad).
    :type phi_: float
    :param c_xs_: The coefficient C_XS.
    :type c_xs_: float
    :param a_surge_: The surge acceleration (m/s²).
    :type a_surge_: float
    :param c_xp_: The coefficient C_XP.
    :type c_xp_: float
    :param a_pitch_: The pitch acceleration (rad/s²).
    :type a_pitch_: float
    :param z_: The z-coordinate of the point (m).
    :type z_: float
    :param R_: The vertical coordinate of ship rotation centre (m).
    :type R_: float

    :return: The calculated acceleration in the x-direction (ax) (m/s²).
    :rtype: float
    """
    ax_ = -c_xg_ * g_ * np.sin(phi_) + c_xs_ * a_surge_ + c_xp_ * a_pitch_ * (z_ - R_)

    return ax_

def calculate_acceleration_y(c_yg_: float, g_: float, thita_: float, c_ys_: float, a_sway_: float,
                             c_yr_: float, a_roll_: float, z_: float, R_: float):
    """
    Calculates the acceleration in the y-direction (ay) for a given dynamic system.

    :param c_yg_: The coefficient C_YG.
    :type c_yg_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float
    :param thita_: The pitch angle θ (rad).
    :type thita_: float
    :param c_ys_: The coefficient C_YS.
    :type c_ys_: float
    :param a_sway_: The sway acceleration (m/s²).
    :type a_sway_: float
    :param c_yr_: The coefficient C_YR.
    :type c_yr_: float
    :param a_roll_: The roll acceleration (rad/s²).
    :type a_roll_: float
    :param z_: The z-coordinate of the point (m).
    :type z_: float
    :param R_: The vertical coordinate of ship rotation centre (m).
    :type R_: float

    :return: The calculated acceleration in the y-direction (ay) (m/s²).
    :rtype: float
    """
    ay_ = c_yg_ * g_ * np.sin(thita_) + c_ys_ * a_sway_ - c_yr_ * a_roll_ * (z_ - R_)

    return ay_

def calculate_acceleration_z(c_zh_: float, a_heave_: float, c_zr_: float, a_roll_: float,
                             y_: float, c_zp_: float, a_pitch_: float, x_: float, L_: float):
    """
    Calculates the acceleration in the z-direction (az) for a given dynamic system.

    :param c_zh_: The coefficient C_ZH.
    :type c_zh_: float
    :param a_heave_: The heave acceleration (m/s²).
    :type a_heave_: float
    :param c_zr_: The coefficient C_ZR.
    :type c_zr_: float
    :param a_roll_: The roll acceleration (rad/s²).
    :type a_roll_: float
    :param y_: The y-coordinate of the point (m).
    :type y_: float
    :param c_zp_: The coefficient C_ZP.
    :type c_zp_: float
    :param a_pitch_: The pitch acceleration (rad/s²).
    :type a_pitch_: float
    :param x_: The x-coordinate of the point (m).
    :type x_: float
    :param L_: The ship's rule length (m).
    :type L_: float

    :return: The calculated acceleration in the z-direction (az) (m/s²).
    :rtype: float
    """
    az_ = c_zh_ * a_heave_ + c_zr_ * a_roll_ * y_ - c_zp_ * a_pitch_ * (x_ - 0.45 * L_)

    return az_