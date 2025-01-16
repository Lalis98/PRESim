
def calculate_lcf(case_: str, sub_case_: str, flp_: float, fT_: float, flp_ost_: float, flp_osa_: float):
    """
    Calculates the Load Combination Factors (LCF) for HSM Equivalent Design Waves.

    :param sub_case_: The identifier for the HSM case ("1" or "2").
    :param case_: The identifier for the HSM case ("1" or "2").
    :type case_: str
    :param flp_: The flp coefficient.
    :type flp_: float
    :param fT_: The draught loading to scantling ratio.
    :type fT_: float

    :return: A dictionary containing the calculated Load Combination Factors (LCF) for the selected case.
    :rtype: dict
    """

    load_cases_ = {
        "HSM":
            {
                "1":
                    {
                        "c_wv": -1.0,
                        "c_qv": -1.0 * flp_,
                        "c_wh": 0.0,
                        "c_wt": 0.0,
                        "c_xs": 0.3 - 0.2 * fT_,
                        "c_xp": -0.7,
                        "c_xg": 0.6,
                        "c_ys": 0.0,
                        "c_yr": 0.0,
                        "c_yg": 0.0,
                        "c_zh": 0.5 * fT_ - 0.15,
                        "c_zr": 0.0,
                        "c_zp": -0.7
                    },
                "2":
                    {
                        "c_wv": 1.0,
                        "c_qv": 1.0 * flp_,
                        "c_wh": 0.0,
                        "c_wt": 0.0,
                        "c_xs": 0.2 - 0.3 * fT_,
                        "c_xp": 0.7,
                        "c_xg": -0.6,
                        "c_ys": 0.0,
                        "c_yr": 0.0,
                        "c_yg": 0.0,
                        "c_zh": 0.15 - 0.5 * fT_,
                        "c_zr": 0.0,
                        "c_zp": 0.7
                    }
            },
        "HSA":
            {
                "1":
                    {
                        "c_wv": -0.7,
                        "c_qv": -0.6 * flp_,
                        "c_wh": 0.0,
                        "c_wt": 0.0,
                        "c_xs": 0.2,
                        "c_xp": -0.4 * fT_ - 0.4,
                        "c_xg": 0.4 * fT_ + 0.4,
                        "c_ys": 0.0,
                        "c_yr": 0.0,
                        "c_yg": 0.0,
                        "c_zh": 0.4 * fT_ - 0.1,
                        "c_zr": 0.0,
                        "c_zp": -0.4 * fT_ - 0.4
                    },
                "2":
                    {
                        "c_wv": 0.7,
                        "c_qv": 0.6 * flp_,
                        "c_wh": 0.0,
                        "c_wt": 0.0,
                        "c_xs": -0.2,
                        "c_xp": 0.4 * fT_ + 0.4,
                        "c_xg": -0.4 * fT_ - 0.4,
                        "c_ys": 0.0,
                        "c_yr": 0.0,
                        "c_yg": 0.0,
                        "c_zh": 0.1 - 0.4 * fT_,
                        "c_zr": 0.0,
                        "c_zp": 0.4 * fT_ + 0.4
                    }
            },
        "FSM":
            {
                "1":
                    {
                        "c_wv": -0.4 * fT_ - 0.6,
                        "c_qv": -1.0 * flp_,
                        "c_wh": 0.0,
                        "c_wt": 0.0,
                        "c_xs": 0.2 - 0.4 * fT_,
                        "c_xp": 0.15,
                        "c_xg": -0.2,
                        "c_ys": 0.0,
                        "c_yr": 0.0,
                        "c_yg": 0.0,
                        "c_zh": 0.0,
                        "c_zr": 0.0,
                        "c_zp": 0.15
                },
                "2":
                    {
                        "c_wv": 0.4 * fT_ + 0.6,
                        "c_qv": 1.0 * flp_,
                        "c_wh": 0.0,
                        "c_wt": 0.0,
                        "c_xs": 0.4 * fT_ - 0.2,
                        "c_xp": -0.15,
                        "c_xg": 0.2,
                        "c_ys": 0.0,
                        "c_yr": 0.0,
                        "c_yg": 0.0,
                        "c_zh": 0.0,
                        "c_zr": 0.0,
                        "c_zp": -0.15
                    }
            },
        "BSR":
            {
                "1P":
                    {
                        "c_wv": 0.1 - 0.2 * fT_,
                        "c_qv": (0.1 - 0.2 * fT_) * flp_,
                        "c_wh": 1.2 - 1.1 * fT_,
                        "c_wt": 0.0,
                        "c_xs": 0.0,
                        "c_xp": 0.0,
                        "c_xg": 0.0,
                        "c_ys": 0.2 - 0.2 * fT_,
                        "c_yr": 1.0,
                        "c_yg": -1.0,
                        "c_zh": 0.7 - 0.4 * fT_,
                        "c_zr": 1.0,
                        "c_zp": 0.0
                    },
                "2P":
                    {
                        "c_wv": 0.2 * fT_ - 0.1,
                        "c_qv": (0.2 * fT_ - 0.1) * flp_,
                        "c_wh": 1.1 * fT_ - 1.2,
                        "c_wt": 0.0,
                        "c_xs": 0.0,
                        "c_xp": 0.0,
                        "c_xg": 0.0,
                        "c_ys": 0.2 * fT_ - 0.2,
                        "c_yr": -1.0,
                        "c_yg": 1.0,
                        "c_zh": 0.4 * fT_ - 0.7,
                        "c_zr": -1.0,
                        "c_zp": 0.0
                    },
                "1S":
                    {
                        "c_wv": 0.1 - 0.2 * fT_,
                        "c_qv": (0.1 - 0.2 * fT_) * flp_,
                        "c_wh": 1.1 * fT_ - 1.2,
                        "c_wt": 0.0,
                        "c_xs": 0.0,
                        "c_xp": 0.0,
                        "c_xg": 0.0,
                        "c_ys": 0.2 * fT_ - 0.2,
                        "c_yr": -1.0,
                        "c_yg": 1.0,
                        "c_zh": 0.7 - 0.4 * fT_,
                        "c_zr": -1.0,
                        "c_zp": 0.0
                    },
                "2S":
                    {
                        "c_wv": 0.2 * fT_ - 0.1,
                        "c_qv": (0.2 * fT_ - 0.1) * flp_,
                        "c_wh": 1.2 - 1.1 * fT_,
                        "c_wt": 0.0,
                        "c_xs": 0.0,
                        "c_xp": 0.0,
                        "c_xg": 0.0,
                        "c_ys": 0.2 - 0.2 * fT_,
                        "c_yr": 1.0,
                        "c_yg": -1.0,
                        "c_zh": 0.4 * fT_ - 0.7,
                        "c_zr": 1.0,
                        "c_zp": 0.0
                    }

            },
        "BSP":
            {
                "1P":
                    {
                        "c_wv": 0.3 - 0.8 * fT_,
                        "c_qv": (0.3 - 0.8 * fT_) * flp_,
                        "c_wh": 0.7 - 0.7 * fT_,
                        "c_wt": 0.0,
                        "c_xs": 0.0,
                        "c_xp": 0.1 - 0.3 * fT_,
                        "c_xg": 0.3 * fT_ - 0.1,
                        "c_ys": -0.9,
                        "c_yr": 0.3,
                        "c_yg": -0.2,
                        "c_zh": 1.0,
                        "c_zr": 0.3,
                        "c_zp": 0.1 - 0.3 * fT_
                    },
                "2P":
                    {
                        "c_wv": 0.8 * fT_ - 0.3,
                        "c_qv": (0.8 * fT_ - 0.3) * flp_,
                        "c_wh": 0.7 * fT_ - 0.7,
                        "c_wt": 0.0,
                        "c_xs": 0.0,
                        "c_xp": 0.3 * fT_ - 0.1,
                        "c_xg": 0.1 - 0.3 * fT_,
                        "c_ys": 0.9,
                        "c_yr": -0.3,
                        "c_yg": 0.2,
                        "c_zh": -1.0,
                        "c_zr": -0.3,
                        "c_zp": 0.3 * fT_ - 0.1
                    },
                "1S":
                    {
                        "c_wv": 0.3 - 0.8 * fT_,
                        "c_qv": (0.3 - 0.8 * fT_) * flp_,
                        "c_wh": 0.7 * fT_ - 0.7,
                        "c_wt": 0.0,
                        "c_xs": 0.0,
                        "c_xp": 0.1 - 0.3 * fT_,
                        "c_xg": 0.3 * fT_ - 0.1,
                        "c_ys": 0.9,
                        "c_yr": -0.3,
                        "c_yg": 0.2,
                        "c_zh": 1.0,
                        "c_zr": -0.3,
                        "c_zp": 0.1 - 0.3 * fT_
                    },
                "2S":
                    {
                        "c_wv": 0.8 * fT_ - 0.3,
                        "c_qv": (0.8 * fT_ - 0.3) * flp_,
                        "c_wh": 0.7 - 0.7 * fT_,
                        "c_wt": 0.0,
                        "c_xs": 0.0,
                        "c_xp": 0.3 * fT_ - 0.1,
                        "c_xg": 0.1 - 0.3 * fT_,
                        "c_ys": -0.9,
                        "c_yr": 0.3,
                        "c_yg": -0.2,
                        "c_zh": -1.0,
                        "c_zr": 0.3,
                        "c_zp": 0.3 * fT_ - 0.1
                    }
            },
        "OST":
            {
                "1P":
                    {
                        "c_wv": -0.3 - 0.2 * fT_,
                        "c_qv": (-0.35 - 0.2 * fT_) * flp_,
                        "c_wh": -0.9,
                        "c_wt": -flp_ost_,
                        "c_xs": 0.1 * fT_ - 0.15,
                        "c_xp": 0.7 - 0.3 * fT_,
                        "c_xg": 0.2 * fT_ - 0.45,
                        "c_ys": 0.0,
                        "c_yr": 0.4 * fT_ - 0.25,
                        "c_yg": 0.1 - 0.2 * fT_,
                        "c_zh": 0.2 * fT_ - 0.05,
                        "c_zr": 0.4 * fT_ - 0.25,
                        "c_zp": 0.7 - 0.3 * fT_
                    },
                "2P":
                    {
                        "c_wv": 0.3 + 0.2 * fT_,
                        "c_qv": (0.35 + 0.2 * fT_) * flp_,
                        "c_wh": 0.9,
                        "c_wt": flp_ost_,
                        "c_xs": 0.15 - 0.1 * fT_,
                        "c_xp": 0.3 * fT_ - 0.7,
                        "c_xg": 0.45 - 0.2 * fT_,
                        "c_ys": 0.0,
                        "c_yr": 0.25 - 0.4 * fT_,
                        "c_yg": 0.2 * fT_ - 0.1,
                        "c_zh": 0.05 - 0.2 * fT_,
                        "c_zr": 0.25 - 0.4 * fT_,
                        "c_zp": 0.3 * fT_ - 0.7
                    },
                "1S":
                    {
                        "c_wv": -0.3 - 0.2 * fT_,
                        "c_qv": (-0.35 - 0.2 * fT_) * flp_,
                        "c_wh": 0.9,
                        "c_wt": flp_ost_,
                        "c_xs": 0.1 * fT_ - 0.15,
                        "c_xp": 0.7 - 0.3 * fT_,
                        "c_xg": 0.2 * fT_ - 0.45,
                        "c_ys": 0.0,
                        "c_yr": 0.25 - 0.4 * fT_,
                        "c_yg": 0.2 * fT_ - 0.1,
                        "c_zh": 0.2 * fT_ - 0.05,
                        "c_zr": 0.25 - 0.4 * fT_,
                        "c_zp": 0.7 - 0.3 * fT_
                    },
                "2S":
                    {
                        "c_wv": 0.3 + 0.2 * fT_,
                        "c_qv": (0.35 + 0.2 * fT_) * flp_,
                        "c_wh": -0.9,
                        "c_wt": -flp_ost_,
                        "c_xs": 0.15 - 0.1 * fT_,
                        "c_xp": 0.3 * fT_ - 0.7,
                        "c_xg": 0.45 - 0.2 * fT_,
                        "c_ys": 0.0,
                        "c_yr": 0.4 * fT_ - 0.25,
                        "c_yg": 0.1 - 0.2 * fT_,
                        "c_zh": 0.05 - 0.2 * fT_,
                        "c_zr": 0.4 * fT_ - 0.25,
                        "c_zp": 0.3 * fT_ - 0.7
                    }
            },
        "OSA":
            {
                "1P":
                    {
                        "c_wv": 0.75 - 0.5 * fT_,
                        "c_qv": (0.6 - 0.4 * fT_) * flp_,
                        "c_wh": 0.55 + 0.2 * fT_,
                        "c_wt": -flp_osa_,
                        "c_xs": 0.1 * fT_ - 0.45,
                        "c_xp": 1.0,
                        "c_xg": -1.0,
                        "c_ys": -0.2 - 0.1 * fT_,
                        "c_yr": 0.3 - 0.2 * fT_,
                        "c_yg": 0.1 * fT_ - 0.2,
                        "c_zh": -0.2 * fT_,
                        "c_zr": 0.3 - 0.2 * fT_,
                        "c_zp": 1.0
                    },
                "2P":
                    {
                        "c_wv": -0.75 + 0.5 * fT_,
                        "c_qv": (-0.6 + 0.4 * fT_) * flp_,
                        "c_wh": -0.55 - 0.2 * fT_,
                        "c_wt": flp_osa_,
                        "c_xs": 0.45 - 0.1 * fT_,
                        "c_xp": -1.0,
                        "c_xg": 1.0,
                        "c_ys": 0.2 + 0.1 * fT_,
                        "c_yr": 0.2 * fT_ - 0.3,
                        "c_yg": 0.2 - 0.1 * fT_,
                        "c_zh": 0.2 * fT_,
                        "c_zr": 0.2 * fT_ - 0.3,
                        "c_zp": -1.0
                    },
                "1S":
                    {
                        "c_wv": 0.75 - 0.5 * fT_,
                        "c_qv": (0.6 - 0.4 * fT_) * flp_,
                        "c_wh": -0.55 - 0.2 * fT_,
                        "c_wt": flp_osa_,
                        "c_xs": -0.45 + 0.1 * fT_,
                        "c_xp": 1.0,
                        "c_xg": -1.0,
                        "c_ys": 0.2 + 0.1 * fT_,
                        "c_yr": 0.2 * fT_ - 0.3,
                        "c_yg": 0.2 - 0.1 * fT_,
                        "c_zh": -0.2 * fT_,
                        "c_zr": 0.2 * fT_ - 0.3,
                        "c_zp": 1.0
                    },
                "2S":
                    {
                        "c_wv": -0.75 + 0.5 * fT_,
                        "c_qv": (-0.6 + 0.4 * fT_) * flp_,
                        "c_wh": 0.55 + 0.2 * fT_,
                        "c_wt": -flp_osa_,
                        "c_xs": 0.45 - 0.1 * fT_,
                        "c_xp": -1.0,
                        "c_xg": 1.0,
                        "c_ys": -0.2 - 0.1 * fT_,
                        "c_yr": 0.3 - 0.2 * fT_,
                        "c_yg": 0.1 * fT_ - 0.2,
                        "c_zh": 0.2 * fT_,
                        "c_zr": 0.3 - 0.2 * fT_,
                        "c_zp": -1.0
                    }
            }
    }

    return load_cases_[case_][sub_case_]