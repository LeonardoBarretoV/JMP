# ============================================================
# Time series functions: To feed in the transition path
# ============================================================
function make_time_series(c, e_init, e_final, b_path_in, Λ_bstar_f)
    Tsim = length(b_path_in)

    # -----------------------
    # Unemployment Markov paths
    # -----------------------
    Pee_path = fill(c.Pee, Tsim)
    Pee_path[1] = 0.89
    Pee_path[2] = 0.96
    Pee_path[3] = 0.962
    Pee_path[3] = 0.964
    Pee_path[4] = 0.962
    Pee_path[5] = 0.964
    Pee_path[6] = 0.971

    Puu_path = fill(c.Puu, Tsim)
    Puu_path[1] = 0.55

    # -----------------------
    # Inflation path
    # -----------------------
    π_path = fill(g.π_bar, Tsim) 
    π_path[1] = (1.0065)^0.25-1; π_path[2] = (1.0137)^0.25-1; π_path[3] = (1.0136)^0.25-1; #2020Q2-2020Q4
    π_path[4] = (1.0262)^0.25-1; π_path[5] = (1.0539)^0.25-1; π_path[6] = (1.0539)^0.25-1; π_path[7] = (1.0704)^0.25-1; #2021
    π_path[8] = (1.0854)^0.25-1; π_path[9] = (1.0906)^0.25-1; π_path[10] = (1.0820)^0.25-1; π_path[11] = (1.0645)^0.25-1; #2022
    π_path[12] = (1.0498)^0.25-1; π_path[13] = (1.0297)^0.25-1; π_path[14] = (1.0370)^0.25-1; π_path[15] = (1.0335)^0.25-1; #2023
    π_path[16] = (1.0348)^0.25-1; π_path[17] = (1.0297)^0.25-1; π_path[18] = (1.0244)^0.25-1; π_path[19] = (1.0289)^0.25-1; #2023
 
    # -----------------------
    # Debt path (qb_path in your call)
    # Keep your exact logic: start from steady final and override early quarters
    # -----------------------
    b_path = fill(e_final.q * e_final.b, Tsim) ./ y(e_final)
    b_path[1] = e_init.b ; b_path[2] = 1.24 * 4.0 + 0.0085 * 4.0; b_path[3] = 1.26 * 4.0 # 2020Q2-2020Q4
    b_path[4] = 1.24 * 4.0; b_path[5] = 1.22 * 4.0; b_path[6] = 1.19 * 4.0 ; b_path[7] = 1.20 * 4.0 # 2021
    b_path[8] = 1.21 * 4.0; b_path[9] = 1.18 * 4.0; b_path[10] = 1.18 * 4.0; b_path[11] = 1.18 * 4.0 # 2022
    b_path[12] = 1.16 * 4.0; b_path[13] = 1.177 * 4.0; b_path[14] = 1.186 * 4.0; b_path[15] = 1.20 * 4.0 # 2023
    b_path[16] = 1.208 * 4.0; b_path[17] = 1.20 * 4.0; b_path[18] = 1.207 * 4.0; b_path[19] = 1.218 * 4.0 # 2024
    b_path[20] = 1.208 * 4.0 # 2025

    # -----------------------
    # Covid / lockdowns
    # -----------------------
    covid_path = fill(c.covid, Tsim)
    covid_path[1] = 1.0

    d_path = fill(c.d_, Tsim)

    Γ_n = 0.9023011517193031
    Γ_c = 0.9229604357071445

    # -----------------------
    # Fiscal policy paths
    # -----------------------
    T_path = fill(e_init.T, Tsim)
    T_path[1] = 1.4 * 4 * 0.051

    Tu_path = fill(e_init.lump_sum_Tu, Tsim) ./ y(e_init)
    Tu_path[1] = 1.4 * 4 * 0.046 + 4 * 0.005

    τ_k_path = fill(e_init.τ_k, Tsim)
    τ_k_path[1] = -0.47

    G_path = fill(e_init.G, Tsim) ./ y(e_init)
    G_path[1]  = e_init.G * 1.1; G_path[2]  = 0.0533;   G_path[3]  = 0.0526;
    G_path[4]  = 0.0548;    G_path[5]  = 0.0523;   G_path[6]  = 0.0503;
    G_path[7]  = 0.0494;    G_path[8]  = 0.0481;   G_path[9]  = 0.0471;
    G_path[10] = 0.0470;    G_path[11] = 0.0474;   G_path[12] = 0.0478;
    G_path[13] = 0.0473;    G_path[14] = 0.0478;   G_path[15] = 0.0476;
    G_path[16] = 0.0478;    G_path[17] = 0.0478;   G_path[18] = 0.0484;
    G_path[19] = 0.0487;

    # -----------------------
    # bstar paths
    # -----------------------
    bstar_path = fill(Λ_bstar_f, Tsim)
    bstar_path[1] = e_init.bstar ; bstar_path[2] = 0.326 * 4.0; bstar_path[3] = 0.320 * 4.0 # 2020Q2-2020Q4
    bstar_path[4] = 0.310 * 4.0; bstar_path[5] = 0.322 * 4.0; bstar_path[6] = 0.316 * 4.0 ; bstar_path[7] = 0.312 * 4.0 # 2021
    bstar_path[8] = 0.30 * 4.0; bstar_path[9] = 0.29 * 4.0; bstar_path[10] = 0.276 * 4.0; bstar_path[11] = 0.274 * 4.0 # 2022
    bstar_path[12] = 0.275 * 4.0; bstar_path[13] = 0.275 * 4.0; bstar_path[14] = 0.268 * 4.0; bstar_path[15] = 0.28 * 4.0 # 2023
    bstar_path[16] = 0.283 * 4.0; bstar_path[17] = 0.284 * 4.0; bstar_path[18] = 0.296 * 4.0; bstar_path[19] = 0.287 * 4.0 # 2024

    bstar_path_r = fill(e_final.bstar, Tsim) ./ e_final.b
    bstar_path_r[1]  = e_init.bstar;     # Note b_path[2] = q[]
    bstar_path_r[2]  = 0.266;    bstar_path_r[3]  = 0.262;    bstar_path_r[4]  = 0.255;    bstar_path_r[5]  = 0.250;
    bstar_path_r[6]  = 0.264;    bstar_path_r[7]  = 0.266;    bstar_path_r[8]  = 0.261;    bstar_path_r[9]  = 0.250;
    bstar_path_r[10] = 0.243;    bstar_path_r[11] = 0.234;    bstar_path_r[12] = 0.233;    bstar_path_r[13] = 0.237;
    bstar_path_r[14] = 0.234;    bstar_path_r[15] = 0.226;    bstar_path_r[16] = 0.234;    bstar_path_r[17] = 0.234;
    bstar_path_r[18] = 0.237;    bstar_path_r[19] = 0.246;    bstar_path_r[20] = 0.236;

    return (
        Pee_path = Pee_path,
        Puu_path = Puu_path,
        π_path = π_path,
        qb_path = b_path,
        covid_path = covid_path,
        d_path = d_path,
        Γ_n = Γ_n,
        Γ_c = Γ_c,
        T_path = T_path,
        Tu_path = Tu_path,
        τ_k_path = τ_k_path,
        G_path = G_path,
        bstar_path = bstar_path,
        bstar_path_r = bstar_path_r,
    )
end