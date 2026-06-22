# ---------------------------------------------------------------------------
# tangential_force.jl
#
# Crankpin tangential force and excitation moment for a single cylinder of
# an internal-combustion engine, as a function of crank angle alpha and
# angular velocity omega.
#
# Reference:
#   Mendes A.S., Meirelles P.S., Zampieri D.E. (2008),
#   "Analysis of torsional vibration in internal combustion engines:
#    modelling and experimental validation",
#   Proc. IMechE Vol. 222 Part K (J. Multi-body Dynamics), pp. 155-178.
#
# Implements equations (8)-(13). Target reproduction: Fig. 14
# (crankpin tangential forces at 2000 and 2550 r/min).
#
# Pressure model: parametric Gaussian approximation of Figs. 18 & 19.
# The paper uses MEASURED p(alpha) curves; this is a stand-in only and
# must be replaced with measured data for quantitative agreement.
# ---------------------------------------------------------------------------

using Plots

# === Engine geometry & oscillating mass (paper p.156) ======================
const D_PISTON = 0.105        # piston diameter dp [m]
const STROKE   = 0.137        # piston stroke s   [m]
const R_CRANK  = STROKE / 2   # crank radius r    [m]    (s = 2r)
const L_ROD    = 0.207        # con-rod length L  [m]
const M_OSC    = 2.521        # oscillating mass ma (piston + mab) [kg]
const LAMBDA   = R_CRANK / L_ROD          # lambda = r/L  ~ 0.331

rpm_to_omega(n_e_rpm) = 2pi * n_e_rpm / 60.0       # [rad/s]

# === Kinematics ============================================================
"Connecting-rod inclination beta from sin(beta) = lambda * sin(alpha)."
beta_angle(alpha) = asin(LAMBDA * sin(alpha))

"Geometric factor sin(alpha+beta)/cos(beta), shared by eqs (9) and (11)."
function geo_factor(alpha)
    b = beta_angle(alpha)
    return sin(alpha + b) / cos(b)
end

# === Cylinder pressure (parametric Gaussian approximation) =================
# Peak combustion pressure read off Fig. 19 of the paper [Pa]
function p_peak_fig19(n_e_rpm)
    rpm = [1000.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2550.0]
    pp  = [14.2,   14.5,   14.5,   14.8,   15.0,   15.0,   13.0,    6.0] .* 1e6
    n = clamp(n_e_rpm, rpm[1], rpm[end])
    i = something(findlast(<=(n), rpm), 1)
    i = min(i, length(rpm) - 1)
    t = (n - rpm[i]) / (rpm[i+1] - rpm[i])
    return pp[i] + t * (pp[i+1] - pp[i])
end

"""
    cylinder_pressure(alpha_rad, n_e_rpm; p_atm, alpha_peak_deg, sigma_deg) -> Pa

Approximate cylinder pressure for one 4-stroke cycle (0..720 deg).
Centre of combustion ~12 deg ATDC (compression TDC = 360 deg by convention
used in Fig. 18). Replace by measured p(alpha) for accurate analyses.
"""
function cylinder_pressure(alpha_rad, n_e_rpm;
                           p_atm=1.0e5,
                           alpha_peak_deg=372.0,
                           sigma_deg=25.0)
    alpha_deg = mod(rad2deg(alpha_rad), 720.0)
    p_peak = p_peak_fig19(n_e_rpm)
    return p_atm + (p_peak - p_atm) *
           exp(-0.5 * ((alpha_deg - alpha_peak_deg) / sigma_deg)^2)
end

# === Forces (eqs 8 - 12) ===================================================
"Gas load on piston, eq. (8): F_g = pi * dp^2 / 4 * p."
gas_load(p) = pi * D_PISTON^2 / 4 * p

"Tangential gas load, eq. (9): F_tg = F_g * sin(alpha+beta)/cos(beta)."
tangential_gas_load(alpha, p) = gas_load(p) * geo_factor(alpha)

"""
Oscillating inertial force, eq. (10):
F_ia = ma r w^2 ( cos a + lambda cos 2a - (lambda^3/4) cos 4a
                                        + (9 lambda^5 / 128) cos 6a )
"""
function oscillating_inertia(alpha, omega)
    return M_OSC * R_CRANK * omega^2 * (
        cos(alpha)
        + LAMBDA       * cos(2alpha)
        - (LAMBDA^3 / 4)            * cos(4alpha)
        + (9 * LAMBDA^5 / 128)      * cos(6alpha)
    )
end

"Tangential inertial force, eq. (11): F_ta = F_ia * sin(alpha+beta)/cos(beta)."
tangential_inertia(alpha, omega) =
    oscillating_inertia(alpha, omega) * geo_factor(alpha)

"""
    tangential_force(alpha, omega, p) -> (F_tg, F_ta, F_t)

Crankpin tangential forces at crank angle `alpha` [rad] and angular speed
`omega` [rad/s], with cylinder pressure `p` [Pa] at this angle.
Returns gas, inertia and resultant components (eq. 12).
"""
function tangential_force(alpha, omega, p)
    F_tg = tangential_gas_load(alpha, p)
    F_ta = tangential_inertia(alpha, omega)
    return F_tg, F_ta, F_tg + F_ta
end

"Excitation moment for one crank throw, eq. (13): M_t = F_t * r."
moment_from_Ft(F_t) = F_t * R_CRANK

# === Sweep over one 4-stroke cycle =========================================
"""
    sweep_cycle(n_e_rpm; alpha_deg=0:1:720, p_func=cylinder_pressure)

Compute the three tangential force components and the excitation moment
over one 4-stroke cycle at engine speed `n_e_rpm`.
Returns a NamedTuple with vectors `alpha_deg`, `F_tg`, `F_ta`, `F_t`, `M_t`.
"""
function sweep_cycle(n_e_rpm; alpha_deg=0:1:720, p_func=cylinder_pressure)
    omega = rpm_to_omega(n_e_rpm)
    a_rad = deg2rad.(collect(alpha_deg))
    F_tg = similar(a_rad); F_ta = similar(a_rad); F_t = similar(a_rad)
    for i in eachindex(a_rad)
        p_i = p_func(a_rad[i], n_e_rpm)
        F_tg[i], F_ta[i], F_t[i] = tangential_force(a_rad[i], omega, p_i)
    end
    # World-frame decomposition of the tangential force vector.
    # alpha from TDC (cylinder axis), rotation increases alpha.
    # crankpin position    : (r sin a, r cos a)
    # tangent unit vector  : (cos a, -sin a)   (direction of rotation)
    F_t_x =  F_t .* cos.(a_rad)        # horizontal component (sideways)
    F_t_y = -F_t .* sin.(a_rad)        # vertical component   (along cylinder axis, + up)
    return (alpha_deg = collect(alpha_deg),
            F_tg = F_tg, F_ta = F_ta, F_t = F_t,
            F_t_x = F_t_x, F_t_y = F_t_y,
            M_t  = moment_from_Ft.(F_t))
end

# === Plot reproducing Fig. 14 ==============================================
function plot_fig14(; savepath::Union{Nothing,String}=nothing)
    r1 = sweep_cycle(2000)
    r2 = sweep_cycle(2550)

    common = (xlabel="Crankshaft angle [deg]", ylabel="[N]",
              xticks=0:60:720, legend=:topleft, lw=1.5,
              framestyle=:box)

    p1 = plot(r1.alpha_deg, r1.F_tg, label="Tangential gas load",
              color=:gray, title="Tangential Forces at 2000 rpm"; common...)
    plot!(p1, r1.alpha_deg, r1.F_ta, label="Tangential inertia force",
          color=:black)
    plot!(p1, r1.alpha_deg, r1.F_t, label="Resultant tangential force",
          color=:black, ls=:dash)

    p2 = plot(r2.alpha_deg, r2.F_tg, label="Tangential gas load",
              color=:gray, title="Tangential Forces at 2550 rpm"; common...)
    plot!(p2, r2.alpha_deg, r2.F_ta, label="Tangential inertia force",
          color=:black)
    plot!(p2, r2.alpha_deg, r2.F_t, label="Resultant tangential force",
          color=:black, ls=:dash)

    fig = plot(p1, p2, layout=(2,1), size=(900, 800))
    savepath !== nothing && savefig(fig, savepath)
    return fig
end
# === plot for interactive query of tangential forces at a single speed =================
function plot_rpm_fig(speed_1; savepath::Union{Nothing,String}=nothing)
    r1 = sweep_cycle(speed_1)

    common = (xlabel="Crankshaft angle [deg]", ylabel="[N]",
              xticks=0:60:720, legend=:topleft, lw=1.5,
              framestyle=:box)

    p1 = plot(r1.alpha_deg, r1.F_tg, label="Tangential gas load",
              color=:gray, title="Tangential Forces at 2000 rpm"; common...)
    plot!(p1, r1.alpha_deg, r1.F_ta, label="Tangential inertia force",
          color=:black)
    plot!(p1, r1.alpha_deg, r1.F_t, label="Resultant tangential force",
          color=:black, ls=:dash)

    fig = plot(p1, size=(900, 800))
    #savepath !== nothing && savefig(fig, savepath)
    return fig
end

# Run with:
#   julia> include("tangential_force.jl")
#   julia> plot_fig14(savepath="imgs/tangential_force_fig14.png")
#
# Single-point query:
#   julia> tangential_force(deg2rad(380.0), rpm_to_omega(2000),
#                           cylinder_pressure(deg2rad(380.0), 2000))
#=If my sign convention for F_t,y (positive up) is the opposite of what you wanted, 
    just flip the sign in line 132 of tangential_force.jl (F_t_y = -F_t .* sin.(a_rad) → +F_t .* sin.(a_rad)).=#