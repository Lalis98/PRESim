from typing import Tuple

import numpy as np

def calculate_pressure(Ps: float, Pd: float, design_load: str):
    """
    Calculates the pressure based on the design load type.

    :param Ps: Value of Static Pressure.
    :type Ps: float
    :param Pd: Value of Dynamic Pressure.
    :type Pd: float
    :param design_load: Load scenario ("S", "D", or "S+D").
    :type design_load: str

    :return: Calculated pressure based on the design load.
    :rtype: float
    """
    design_load = design_load.strip().upper()

    if design_load == "S":
        return max(0.0, Ps)
    elif design_load == "D":
        return Pd
    elif design_load == "S+D":
        return max(0.0, Ps + Pd)
    else:
        raise ValueError("Invalid design_load. Expected 'S', 'D', or 'S+D'.")

def calculate_wave_pressure(z_, TLC_, Pw_wl_, rho_, g_, hw_, P_case_) -> float:
    """
    Calculates the wave pressure Pw based on the given conditions.

    :param z_: The vertical position of the point.
    :type z_: float
    :param TLC_: The total vertical clearance.
    :type TLC_: float
    :param Pw_wl_: The wave pressure at the waterline.
    :type Pw_wl_: float
    :param rho_: The density.
    :type rho_: float
    :param g_: The gravitational acceleration.
    :type g_: float
    :param hw_: The wave height at the waterline.
    :type hw_: float
    :param P_case_: Pressure related to wave pressure case.
    :type P_case_: float

    :return: The calculated wave pressure Pw.
    :rtype: float
    """
    if z_ <= TLC_:
        return max(P_case_, rho_ * g_ * (z_ - TLC_))  # For FSM-1 or below TLC
    elif TLC_ < z_ <= hw_ + TLC_:
        return Pw_wl_ - rho_ * g_ * (z_ - TLC_)  # For wave height above TLC
    else:
        return 0.0  # Outside of wave height range

def calculate_hydrostatic_pressure(z_, TLC_, rho_, g_) -> float:
    """
    Calculates the hydrostatic pressure Ps based on the given conditions.

    :param z_: The vertical position of the point.
    :type z_: float
    :param TLC_: The total vertical clearance.
    :type TLC_: float
    :param rho_: The density of the fluid.
    :type rho_: float
    :param g_: The gravitational acceleration.
    :type g_: float

    :return: The calculated hydrostatic pressure Ps.
    :rtype: float
    """
    if z_ <= TLC_:
        return rho_ * g_ * (TLC_ - z_)  # Hydrostatic pressure below the waterline (TLC)
    else:
        return 0.0  # No hydrostatic pressure above the waterline


def get_angle_for_z(z: float, sections: np.ndarray | float) -> float:
    """
    * Retrieves the angle corresponding to a given z-coordinate or directly returns the angle if a float is provided.

    :param z: The z-coordinate of the point.
    :type z: float
    :param sections: Either: (1) A numpy array where each row is [z_start, z_end, angle], representing sections
     with start and end z-coordinates and their respective angles (degrees).
      (2) A float value representing a direct angle in degrees.
    :type sections: np.ndarray | float

    :return: The angle corresponding to the z-coordinate (rad).
    :rtype: float
    """
    if isinstance(sections, float):  # If sections is a float, treat it as the angle
        return sections * np.pi / 180.0

    # If sections is an array, find the row where z is within [z_start, z_end]
    row = sections[(sections[:, 0] <= z) & (z <= sections[:, 1])]

    if row.size > 0:
        return row[0, 2] * np.pi / 180.0  # Return the angle
    else:
        raise ValueError(f"z-coordinate {z} is out of range.")


def calculate_dynamic_cargo_pressure(fb_: float, fdc_: float, Kc_: float, rho_c_: float,
                                     a_: np.ndarray, center_of_gravity_: np.ndarray,
                                     coordinates_: np.ndarray, zc_: float):
    """
    Calculates the dynamic cargo pressure on a given point in the ship's hold.

    :param fb_: Heading correction factor.
    :type fb_: float
    :param fdc_: Dry cargo factor.
    :type fdc_: float
    :param Kc_: Coefficient related to cargo.
    :type Kc_: float
    :param rho_c_: The density of the cargo ρc (t/m³).
    :type rho_c_: float
    :param a_: The accelerations in the x, y, z directions [ax, ay, az] (m/s²).
    :type a_: np.ndarray
    :param center_of_gravity_: The coordinates of the ship's center of gravity [xg, yg, zg] (m).
    :type center_of_gravity_: np.ndarray
    :param coordinates_: The coordinates of the specific point in the cargo hold [x, y, z] (m).
    :type coordinates_: np.ndarray
    :param zc_: The height of the upper surface of the cargo above the baseline in way of the
     load point (m).
    :type zc_: float

    :return: The calculated dynamic cargo pressure Pbd (kN/m²).
    :rtype: float
    """

    ax_, ay_, az_ = a_
    xg_, yg_, zg_ = center_of_gravity_
    x_, y_, z_ = coordinates_

    if z_ > zc_:
        Pbd_ = 0
    else:
        Pbd_ = fb_ * rho_c_ * (
                0.25 * ax_ * (xg_ - x_) + 0.25 * ay_ * (yg_ - y_) + fdc_ * Kc_ * az_ * (zc_ - z_)
        )

    return Pbd_

def calculate_static_cargo_pressure(rho_c_: float, g_: float, Kc_: float, zc_: float, z_: float):
    """
    Calculates the static cargo pressure on a given point in the ship's hold.

    :param rho_c_: The density of the cargo ρc (t/m³).
    :type rho_c_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float
    :param Kc_: Coefficient related to cargo.
    :type Kc_: float
    :param zc_: The height of the upper surface of the cargo above the baseline in way of the
                 load point (m).
    :type zc_: float
    :param z_: The z-coordinate of the point (m).
    :type z_: float

    :return: The calculated static cargo pressure Pbs (kN/m²).
    :rtype: float
    """
    Pbs_ = max(
        0.0,
        rho_c_ * g_ * Kc_ * (zc_ - z_)
    )

    return Pbs_

def calculate_hc_full(hHPU_: float, S0_: float, Bh_: float, Vhc_: float, lh_: float):
    """
    Calculates the vertical distance hc from the inner bottom to the top of the hatch coaming.

    :param hHPU_: The vertical distance from the inner bottom at the centerline to
                   the lower intersection of the topside tank and side shell or inner side for double side
                   bulk carriers, determined at mid-length of the cargo hold at midship (m).
    :type hHPU_: float
    :param S0_: The shaded area above the lower intersection of the topside tank and
                side shell or inner side (m²).
    :type S0_: float
    :param Bh_: The breadth of the cargo hold (m).
    :type Bh_: float
    :param Vhc_: The volume of the hatch coaming (m³).
    :type Vhc_: float
    :param lh_: The length of the cargo hold (m).
    :type lh_: float

    :return: The calculated vertical distance hc (m).
    :rtype: float
    """

    SA_ = S0_ + Vhc_ / lh_

    h0_ = SA_ / Bh_
    hc_ = hHPU_ + h0_

    return hc_


def calculate_angle_psi(cargo_type_: str):
    """
    Calculates the angle ψ based on the cargo type.

    :param cargo_type_: The type of cargo being handled.
                        Accepted values are "General", "Iron Ore", or "Cement".
    :type cargo_type_: str

    :return: The calculated angle ψ in radians.
    :rtype: float
    """
    if cargo_type_ == "General":
        psi_ = 30 * np.pi / 180.0  # rad
    elif cargo_type_ == "Iron Ore":
        psi_ = 35 * np.pi / 180.0  # rad
    elif cargo_type_ == "Cement":
        psi_ = 25 * np.pi / 180.0  # rad
    else:
        raise ValueError(f"Invalid cargo type: {cargo_type_}. Allowed values are 'General', 'Iron Ore', or 'Cement'.")

    return psi_

def calculate_hc_partial(Bh_: float, Bib_: float, psi_: float, y_: float, M_: float, rho_c_: float,
                         lh_: float, hHPL_: float, Vts_: float):
    """
    Calculates the vertical distance (hc) and the corrected vertical distance (hc_cl) for a given bulk cargo hold.

    :param Bh_: The breadth of the cargo hold (m).
    :type Bh_: float
    :param Bib_: The breadth of the inner bottom (m).
    :type Bib_: float
    :param psi_: The angle ψ related to the cargo type (rad).
    :type psi_: float
    :param y_: The transverse position of the specific point in the cargo hold (m).
    :type y_: float
    :param M_: The mass of the bulk cargo being considered (t).
    :type M_: float
    :param rho_c_: The density of the cargo ρc (t/m³).
    :type rho_c_: float
    :param lh_: The length of the cargo hold (m).
    :type lh_: float
    :param hHPL_: The vertical distance from the inner bottom at centerline to the upper intersection
     of the hopper tank and side shell (m).
    :type hHPL_: float
    :param Vts_: The total volume of the portion of the lower bulkhead stools within the cargo hold
     length and inboard of the hopper tanks (m³).
    :type Vts_: float

    :return: A tuple containing vertical distances hc and hc_cl (m).
    :rtype: (float, float)
    """

    h1_ = (M_ / (rho_c_ * Bh_ * lh_) - (Bh_ + Bib_) / (2 * Bh_) * hHPL_ - 3 / 16 *
           Bh_ * np.tan(psi_ / 2) + Vts_ / (Bh_ * lh_))

    # Calculate hc_cl
    if h1_ >= 0:

        B2_ = Bh_
        h2_ = Bh_ / 4 * np.tan(psi_ / 2)

        hc_cl_ = hHPL_ + h1_ + h2_

    else:  # h1 < 0

        B2_ = np.sqrt(
            (1 / lh_ * (M_ / rho_c_ + Vts_) + 1 / 2 * (hHPL_ * Bib_ ** 2) / (Bh_ - Bib_) + Bh_ ** 2 / 16 * np.tan(psi_ / 2))
            /
            (1 / 2 * (hHPL_ / (Bh_ - Bib_)) + 1 / 2 * np.tan(psi_ / 2))
        )

        h11_ = hHPL_ * (B2_ - Bib_) / (Bh_ - Bib_)
        h22_ = (B2_ / 2 - Bh_ / 4) * np.tan(psi_ / 2)

        hc_cl_ = h11_ + h22_

    # Calculate hc
    if np.abs(y_) <= Bh_ / 4:

        hc_ = hc_cl_

    elif Bh_ / 4 <= np.abs(y_) < B2_ / 2:

        hc_ = hc_cl_ - (np.abs(y_) - Bh_ / 4) * np.tan(psi_ / 2)

    else:

        hc_ = 0

    return hc_, hc_cl_


def calculate_static_shear_load_on_hopper(rho_c_: float, g_: float, Kc_: float, zc_: float,
                                          z_: float, alpha_: float) -> float:

    """
    Calculates the static shear load pressure on a given point in lower stool or hopper tank.

    :param rho_c_: The density of the cargo ρc (t/m³).
    :type rho_c_: float
    :param g_: The acceleration due to gravity (m/s²).
    :type g_: float
    :param Kc_: Coefficient related to cargo.
    :type Kc_: float
    :param zc_: The height of the upper surface of the cargo above the baseline in way of the load point (m).
    :type zc_: float
    :param z_: The x-coordinate of the point.
    :type z_: float
    :param alpha_: The angle α (rad).
    :type alpha_: float

    :return: The calculated static shear load pressure Pbs_s (kN/m²).
    :rtype: float
    """

    Pbs_s_ = rho_c_ * g_ * (1 - Kc_) * (zc_ - z_) / np.tan(alpha_)

    return Pbs_s_


def calculate_dynamic_shear_load_on_hopper(fb_:float, az_:float, rho_c_: float, Kc_: float, zc_: float,
                                          z_: float, alpha_: float) -> float:
    """
    Calculates the dynamic shear load pressure on a given point in lower stool or hopper tank.

    :param fb_: Head Correction factor.
    :type fb_: float
    :param az_: The acceleration on z-direction (m/s²).
    :type az_: float
    :param rho_c_: The density of the cargo ρc (t/m³).
    :type rho_c_: float
    :param Kc_: Coefficient related to cargo.
    :type Kc_: float
    :param zc_: The height of the upper surface of the cargo above the baseline in way of the load point (m).
    :type zc_: float
    :param z_: The x-coordinate of the point.
    :type z_: float
    :param alpha_: The angle α (rad).
    :type alpha_: float

    :return: The calculated dynamic shear load pressure Pbs_d (kN/m²).
    :rtype: float
    """

    Pbs_d_ = fb_ * rho_c_ * az_ * (1 - Kc_) * (zc_ - z_) / np.tan(alpha_)

    return Pbs_d_


def calculate_dynamic_shear_load_on_inner_bottom(fb_: float, ax_: float, ay_: float,
                                                 rho_c_: float, hc_: float) -> tuple[float, float]:
    """
    Calculates the vertical distance (hc) for a given bulk cargo hold.

    :param fb_: Head Correction factor.
    :type fb_: float
    :param ax_: The acceleration on x-direction (m/s²).
    :type ax_: float
    :param ay_: The acceleration on y-direction (m/s²).
    :type ay_: float
    :param rho_c_: The density of the cargo ρc (t/m³).
    :type rho_c_: float
    :param hc_: Coefficient related to cargo.
    :type hc_: float

    :return: The calculated dynamic shear load pressure (Pbs_dx, Pbs_dy)  (kN/m²).
    :rtype: (float, float)
    """

    Pbs_dx_ = -0.75 * fb_ * rho_c_ * ax_ * hc_
    Pbs_dy_ = -0.75 * fb_ * rho_c_ * ay_ * hc_

    return Pbs_dx_, Pbs_dy_


