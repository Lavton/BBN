reactions = [
    [
        "p to n, elements.H_1.p_to_n, elements.H_1.n_to_p",
    ],
    [
        "p + n to D, elements.H_2.H_2_forw_rate, elements.H_2.H_2_backward_rate",
    ],
    [
        "p + D to He_3, elements.He_3.d_pg_he3, elements.He_3.he3_gp_d",
        "D + D to n + He_3, elements.He_3.dd_nhe3, elements.He_3.nhe3_dd"
    ],
    [
        "n + D to T, elements.H_3.nd_tg, elements.H_3.tg_nd",
        "n + He_3 to p + T, elements.H_3.nhe3_pt, elements.H_3.pt_nhe3",
        "D + D to p + T, elements.H_3.dd_pt, elements.H_3.pt_dd"
    ],
    [
        "p + T to He_4, elements.He_4.pt_he4g, elements.He_4.he4g_pt",
        "n + He_3 to He_4, elements.He_4.nhe3_he4g, elements.He_4.he4g_nhe3",
        "D + D to He_4, elements.He_4.dd_he4g, elements.He_4.he4g_dd",
        "D + He_3 to He_4 + p, elements.He_4.dhe3_he4p, elements.He_4.he4p_dhe3",
        "D + T to He_4 + n, elements.He_4.dt_he4n, elements.He_4.he4n_dt",
        "He_3 + He_3 to He_4 + p + p, elements.He_4.he3he3_he4pp, elements.He_4.he4pp_he3he3",
        "T + T to He_4 + n + n, elements.He_4.tt_he4nn, elements.He_4.he4nn_tt",
        "He_3 + T to He_4 + p + n, elements.He_4.he3t_he4pn, elements.He_4.he4pn_he3t",
        "He_3 + T to He_4 + D, elements.He_4.he3t_he4d, elements.He_4.he4d_he3t"
    ],
    [
        "He_3 + He_4 to Be_7, elements.Be_7.he3he4_be7g, elements.Be_7.be7g_he3he4",
        "n + Be_7 to He_4 + He_4, elements.Be_7.nBe7_He4He4, elements.Be_7.He4He4_nBe7",
    ],
    [
        "T + He_4 to Li_7, elements.Li_7.the4_li7g, elements.Li_7.li7g_the4",
        "n + Be_7 to p + Li_7, elements.Li_7.nbe7_pli7, elements.Li_7.pli7_nbe7",
        "p + Li_7 to He_4 + He_4, elements.Li_7.pli7_he4he4, elements.Li_7.he4he4_pli7"
    ],
    [
        # "He_4 + D to Li_6, elements.Li_6.he4d_li6g, elements.Li_6.li6g_he4d",
        # "He_4 + n + p to Li_6, elements.Li_6.he4np_li6g, elements.Li_6.li6g_he4np",
        # "He_4 + T + n to Li_6, elements.Li_6.he4t_li6n, elements.Li_6.be7g_li6p",
        # "Li_6 + p to Be_7, elements.Li_6.li6p_be7g, elements.Li_6.li6n_he4t",
        # "Li_6 + p to He_4 + He_3, elements.Li_6.li6p_he3he4, elements.Li_6.he3he4_li6p"
    ]
]
