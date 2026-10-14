using LinearAlgebra, ProgressMeter, Roots, NLsolve, StructArrays, StatsBase, Tullio,
      Polyester, PrettyTables, CairoMakie, LaTeXStrings, DelimitedFiles, XLSX, CSV,
      DataFrames, ColorTypes, Printf, JLD2, MAT
const ROOT = normpath(joinpath(@__DIR__, ".."))
for f in ("types","preferences","technology_and_fiscal","calibration","auxiliary_functions",
          "household_problem_and_aggregation","stationary_equilibria","transition",
          "union_block","transition_sw","counterfactuals_sw","statistics","time_series")
    include(joinpath(ROOT, "src", f * ".jl"))
end
Jf(f) = joinpath(ROOT,"output","jld2_sw",f*".jld2")
_i = load_object(Jf("steady_states")); e = _i[1]; ef = _i[2]
h = e.h; t = e.t; g = e.g
@printf("--- technology ---\n")
@printf("alpha            = %.10f\n", t.α)
@printf("delta            = %.12f   (ann %.4f%%)\n", t.δ, 100*(1-(1-t.δ)^4))
@printf("A (Z)            = %.10f\n", t.A)
@printf("mu (markup)      = %.10f\n", t.μ)
@printf("--- preferences ---\n")
@printf("beta             = %.16f\n", get_β(h.u))
@printf("1/ies (xi_inv)   = %.10f  -> ies xi = %.4f\n", get_inverse_ies(h.u), 1/get_inverse_ies(h.u))
@printf("theta            = %.10f\n", h.v.θ)
@printf("nu (Frisch)      = %.10f\n", h.v.ν)
@printf("--- fiscal ---\n")
@printf("tau (progressiv) = %.10f\n", g.τ)
@printf("tau_n ss         = %.10f\n", e.τ_n)
@printf("tau_k ss         = %.10f\n", e.τ_k)
@printf("--- assets / prices, initial SS ---\n")
@printf("r ss             = %.12f  (ann %.4f%%)\n", e.r, 100*((1+e.r)^4-1))
@printf("w ss             = %.10f\n", e.w)
@printf("pi_bar           = %.12f  (ann %.4f%%)\n", e.π_, 100*((1+e.π_)^4-1))
@printf("q ss             = %.10f\n", e.q)
@printf("k ss             = %.10f    b ss = %.10f   bstar ss = %.10f\n", e.k, e.b, e.bstar)
@printf("y ss             = %.10f\n", y(e))
@printf("n ss (N)         = %.10f\n", e.n)
@printf("a_min            = %.6f    a_max = %.2f    n_a = %d    n_z = %d\n",
        first(h.a_grid), last(h.a_grid), length(h.a_grid), length(h.z_grid))
@printf("z_grid[1]        = %.4f  (unemployment state)   z_grid[end] = %.4f\n", h.z_grid[1], h.z_grid[end])
@printf("--- terminal SS (baseline) ---\n")
@printf("r_final          = %.12f  (ann %.4f%%)\n", ef.r, 100*((1+ef.r)^4-1))
i_ss = (1+ef.r)*(1+e.π_)-1
@printf("i* = (1+r_f)(1+pi_bar)-1 = %.12f  (ann %.4f%%)\n", i_ss, 100*((1+i_ss)^4-1))
@printf("--- union block normalisation ---\n")
DY = union_aggregates_(h, g, e.ws; s=1.0, w=e.w, τ_n=e.τ_n, n_ss=e.ws.n_pol[1,:])
@printf("D_ss = %.10f   Y_ss = %.10f   w_ss*Y_ss = %.10f\n", DY.D, DY.Y, e.w*DY.Y)
@printf("--- lockdown, from the rule-2 time series ---\n")
Base.@kwdef struct Cal
    θ=1.0; ν=1.0; ies=1.0; crra=4.0; β=0.9880431080390226
    a_min=-0.02; a_max=120.0; n_a=1000; z_bar=12; Pee=0.975; Puu=0.4
    Λ_t=0.0; Λ_tu=0.005; Λ_a=5.06*4.0; Λ_g=0.05; Λ_f=0.04; Λ_b=1.05*4.0
    Λ_bstar=0.31*4.0; Λ_rk=0.33; r_target=(1.01)^(1/4)-1; π_bar=(1.02)^(1/4)-1
    mc=20.0; τ=0.181; nu_target=0.7; A=0.035; μ=1.0
    d_=1.0; Γ_n=1.0; Γ_c=1.0; covid=0.0; χ=0.44
end
c = Cal()
bt = e.b/y(e); T=700; H=100; ρB=0.08
bl = Array{Float64,1}(undef,T+H); bl[1]=bt; bl[2]=bt; bl[T:end].=bt
for i in 3:T-1; bl[i]=bl[1]*ρB^(i-1)+(1-ρB^(i-1))*bt; end
ts = make_time_series_sw(c, e, ef, bl, 0.287*4.0, g)
@printf("Gamma_n          = %.14f\n", ts.Γ_n)
@printf("Gamma_c          = %s\n", ts.Γ_c isa Number ? @sprintf("%.14f", ts.Γ_c) : @sprintf("path, [1]=%.14f [2]=%.14f", ts.Γ_c[1], ts.Γ_c[2]))
@printf("covid_path[1:3]  = %s\n", string(ts.covid_path[1:3]))
@printf("d_path[1:3]      = %s\n", string(round.(ts.d_path[1:3], digits=6)))
@printf("chi (foreign)    = %.4f\n", c.χ)
@printf("Pee = %.4f  Puu = %.4f\n", c.Pee, c.Puu)
