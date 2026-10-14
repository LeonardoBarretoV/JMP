# ============================================================
# Preferences: Consumption and Labor
# ============================================================
# # Disutility of labor methods
# Fixed Labor
labor(v::FixedLabor) = v.n  # labor supply
labor_income(v::FixedLabor; w) = w * v.n  # labor income
disutility(::FixedLabor{R}) where {R} = zero(R)  # disutility
disutility_given_w(::FixedLabor{R}; w) where {R} = zero(R)  # disutility

# GHH with HSV
labor(v::GHH, g::TaxSys; w , τ_n) = ((1 - g.τ) * τ_n * w^(1 - g.τ) / v.θ)^(v.ν / (1 + g.τ * v.ν)) # labor supply
labor_income(v::GHH, g::TaxSys; w , τ_n) = w * ((1 - g.τ) * τ_n * w^(1 - g.τ) / v.θ)^(v.ν / (1 + g.τ * v.ν))  # labor income
disutility(v::GHH; n) = v.θ * n^(1 + 1/v.ν) / (1 + 1/v.ν)
disutility_given_w(v::GHH, g::TaxSys; w , τ_n) = v.θ * ((1 - g.τ) * τ_n * w^(1 - g.τ) / v.θ)^((1 + v.ν) / (1 + g.τ * v.ν)) / (1 + 1/v.ν)

# GHH with Pandemic Constraints
labor_(v::GHH, g::TaxSys; w , τ_n, n_, Γ_n, covid) = covid == 1 ? min(Γ_n * n_ , ((1 - g.τ) * τ_n * w^(1 - g.τ) / v.θ)^(v.ν / (1 + g.τ * v.ν))) : ((1 - g.τ) * τ_n * w^(1 - g.τ) / v.θ)^(v.ν / (1 + g.τ * v.ν)) # labor supply
labor_income_(v::GHH, g::TaxSys; w , τ_n, n_, Γ_n, covid) = covid == 1 ? w * min(Γ_n * n_ , ((1 - g.τ) * τ_n * w^(1 - g.τ) / v.θ)^(v.ν / (1 + g.τ * v.ν))) : w * ((1 - g.τ) * τ_n * w^(1 - g.τ) / v.θ)^(v.ν / (1 + g.τ * v.ν))  # labor income
disutility_(v::GHH; n) = v.θ * n^(1 + 1/v.ν) / (1 + 1/v.ν)
disutility_given_w_(v::GHH, g::TaxSys; w , τ_n, n_, Γ_n, covid) = covid == 1 ? v.θ * (min(Γ_n * n_ , ((1 - g.τ) * τ_n * w^(1 - g.τ) / v.θ)^(v.ν / (1 + g.τ * v.ν)))) ^(1 + 1/v.ν) / (1 + 1/v.ν) : v.θ * ((1 - g.τ) * τ_n * w^(1 - g.τ) / v.θ)^((1 + v.ν) / (1 + g.τ * v.ν)) / (1 + 1/v.ν)

# GHH without HSV (lump_sum transfers)
labor_lump_sum(v::GHH; w) = (w / v.θ)^v.ν  # labor supply
labor_income_lump_sum(v::GHH; w) = w^(1 + v.ν) / (v.θ^v.ν)  # labor income
disutility_lump_sum(v::GHH; n) = v.θ * n^(1 + 1/v.ν) / (1 + 1/v.ν)
disutility_given_w_lump_sum(v::GHH; w) = v.θ * (w / v.θ)^(v.ν + 1) / (1 + 1/v.ν)
