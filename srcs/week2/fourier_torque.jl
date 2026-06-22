# ---------------------------------------------------------------------------
# fourier_torque.jl
#
# Fourier-series decomposition of the resultant crankpin tangential force
# F_t (and equivalently the per-cylinder excitation torque M_t = F_t * r)
# obtained from `sweep_cycle` in `tangential_force.jl`.
#
# Signal period:  one full 4-stroke cycle = 720 deg crank angle = 2 engine
#                 revolutions, so the fundamental Fourier frequency is
#                 omega_fund = omega_engine / 2.
# This is the same convention as eq. (17) of Mendes et al. (2008):
#
#     M_t(t) = A0/2 + sum_{n=1..N} [ C_n exp(+i n omega_fund t)
#                                  + conj(C_n) exp(-i n omega_fund t) ]
#
# Two-sided complex form (returned by this module):
#
#     F_t(t) = sum_{n = -N..N}  H_n * exp( i * omega_n * t )
#
# with omega_n = n * omega_fund  and  H_{-n} = conj(H_n).
# ---------------------------------------------------------------------------

#include("tangential_force.jl")

using FFTW
using Plots

# --- internal: STFT spectrogram of a periodic signal ------------------------
# Tiles the one-cycle signal `n_cycles` times so the sliding window has room.
# Returns (magnitude, crank_angle_deg_axis, engine_order_axis).
function _stft_spectrogram(sig, n_per_cycle;
                           n_cycles::Integer = 3,
                           win_frac::Real    = 0.5,
                           hop_frac::Real    = 0.025)
    tiled   = repeat(sig, n_cycles)
    N       = length(tiled)
    win_len = max(8, round(Int, n_per_cycle * win_frac))
    hop     = max(1, round(Int, n_per_cycle * hop_frac))
    w       = @. 0.5 * (1 - cos(2pi * (0:win_len-1) / max(win_len - 1, 1)))
    n_win   = (N - win_len) ÷ hop + 1
    nbins   = win_len ÷ 2 + 1
    spec    = zeros(nbins, n_win)
    for k in 0:n_win-1
        s = k * hop + 1
        seg = view(tiled, s:s+win_len-1) .* w
        spec[:, k+1] = abs.(fft(seg))[1:nbins]
    end
    crank_deg = ((0:n_win-1) .* hop .+ (win_len - 1)/2) .* (720.0 / n_per_cycle)
    orders    = (0:win_len÷2) .* (n_per_cycle / (2 * win_len))   # engine orders
    return spec, crank_deg, orders
end

"""
    fourier_decompose(rpm, order; signal=:F_t_y, n_samples=1440,
                      max_order_display=15) -> NamedTuple

Fourier-decompose one of the per-cycle force/moment traces from
`sweep_cycle(rpm)`.  Valid values for `signal`:
  :F_t_y  -- vertical component of the tangential force (default)
  :F_t_x  -- horizontal component of the tangential force
  :F_t    -- scalar tangential force (magnitude along tangent)
  :M_t    -- excitation moment (= F_t * r)

The cycle is sampled at `n_samples` uniform crank-angle points over
0..720 deg (excluding the duplicate endpoint), the DFT is computed, and
the first `order` harmonics (plus the DC term) are kept.

The returned `plot` field is a 3-panel figure:
  1. Top    -- original signal vs Fourier reconstruction (time domain).
  2. Middle -- one-sided amplitude spectrum |H_n| vs engine order.
  3. Bottom -- STFT spectrogram in dB (signal tiled 3 cycles).

Returned fields
  omegas        :: Vector{Float64}     two-sided frequencies omega_n  [rad/s]
                                       length = 2*order + 1, ordered -N..N
  coefficients  :: Vector{ComplexF64}  Fourier coefficients H_n matching omegas
  omega_fund    :: Float64             fundamental = omega_engine / 2  [rad/s]
  alpha_deg     :: Vector{Float64}     sample crank angles               [deg]
  signal        :: Vector{Float64}     original sampled signal           [N or N*m]
  reconstructed :: Vector{Float64}     Fourier reconstruction at same samples
  spectrogram   :: NamedTuple          (magnitude, crank_deg, orders) from STFT
  plot          :: Plots.Plot          3-panel composite figure

Reconstruction at arbitrary time `t [s]`:
  f_t = real(sum(coefficients .* exp.(1im .* omegas .* t)))
or use the helper `evaluate_fourier(t, omegas, coefficients)`.
"""
function fourier_decompose(rpm::Real, order::Integer;
                           signal::Symbol  = :F_t_y,
                           n_samples::Integer = 1440,
                           max_order_display::Real = 15.0)
    # --- sample one full cycle ---------------------------------------------
    alpha_deg = collect(range(0.0, 720.0, length = n_samples + 1))[1:end-1]
    res = sweep_cycle(rpm; alpha_deg = alpha_deg)
    sig = signal === :M_t   ? res.M_t   :
          signal === :F_t_x ? res.F_t_x :
          signal === :F_t_y ? res.F_t_y :
          signal === :F_t   ? res.F_t   :
          error("fourier_decompose: unknown signal $(signal); use :F_t_y, :F_t_x, :F_t, or :M_t")
    N = length(sig)

    omega_engine = rpm_to_omega(rpm)          # [rad/s]
    omega_fund   = omega_engine / 2           # 720 deg period -> half engine speed

    # --- DFT --> complex Fourier coefficients ------------------------------
    Y = fft(sig)
    c = Y ./ N                                # c_n,  n = 0,1,...,N-1
                                              # c_{N-n} = conj(c_n)  for real sig

    nmax = min(order, N ÷ 2)
    c_pos = c[1:nmax+1]                       # n = 0, 1, ..., nmax

    # two-sided arrays, ordered -nmax .. nmax
    omegas = collect(-nmax:nmax) .* omega_fund
    coeffs = ComplexF64[ n >= 0 ? c_pos[n + 1] : conj(c_pos[-n + 1])
                         for n in -nmax:nmax ]

    # --- reconstruct on the same crank-angle grid --------------------------
    t = deg2rad.(alpha_deg) ./ omega_engine   # time [s] at each sample
    recon = [ real(sum(coeffs .* exp.(1im .* omegas .* tt))) for tt in t ]

    # --- STFT spectrogram (tile 3 cycles) ----------------------------------
    spec, sg_crank_deg, sg_orders = _stft_spectrogram(sig, N; n_cycles = 3)
    spec_db = 20 .* log10.(spec ./ maximum(spec) .+ 1e-6)
    keep = findall(sg_orders .<= max_order_display)

    # --- plot panels -------------------------------------------------------
    ylab = signal === :M_t   ? "Excitation moment M_t [N m]" :
           signal === :F_t_x ? "F_t horizontal component  F_t,x [N]" :
           signal === :F_t_y ? "F_t vertical component  F_t,y [N]" :
                                "Resultant tangential force F_t [N]"

    lmarg = 12Plots.mm

    t_ms = t .* 1e3                                   # time axis in milliseconds
    ms_per_deg = 60_000.0 / (rpm * 360)               # 1° crank in ms

    fmt(x) = isinteger(x) ? string(Int(x)) :
                            string(round(x, digits = 1))

    # tick positions chosen on the crank-angle grid (multiples of 60°) and
    # labelled with both ms and deg so the user reads the same tick two ways.
    tick_deg_time  = 0:60:720
    tick_ms_time   = collect(tick_deg_time) .* ms_per_deg
    tick_labels_t  = ["$(fmt(ms)) ms\n$(d)°"
                      for (ms, d) in zip(tick_ms_time, tick_deg_time)]

    p_time = plot(t_ms, sig,
                  label = "Original", color = :black, lw = 1.5,
                  xlabel = "Time [ms]  /  Crank angle [deg]", ylabel = ylab,
                  xlims  = (0, t_ms[end]),
                  xticks = (tick_ms_time, tick_labels_t),
                  framestyle = :box,
                  left_margin = lmarg,
                  bottom_margin = 4Plots.mm,
                  title  = "Fourier reconstruction @ $(rpm) rpm,  order = $(nmax)")
    plot!(p_time, t_ms, recon,
          label = "Fourier (n = 0..$(nmax))",
          color = :red, lw = 1.5, ls = :dash)

    # one-sided amplitude spectrum (engine orders 0, 0.5, 1, ...)
    orders_x = collect(0:nmax) ./ 2
    mags = Float64[ n == 0 ? abs(c_pos[1]) : 2 * abs(c_pos[n + 1])
                    for n in 0:nmax ]
    p_spec = bar(orders_x, mags,
                 label = "|H_n|",
                 xlabel = "Engine order  (n / 2)", ylabel = "Amplitude",
                 color = :steelblue, lw = 0, bar_width = 0.35,
                 xticks = 0:1:max(1, ceil(orders_x[end])),
                 framestyle = :box,
                 left_margin = lmarg,
                 title = "One-sided amplitude spectrum")

    # dual tick labels for the spectrogram (x is crank angle, also show ms)
    sg_lo = first(sg_crank_deg);  sg_hi = last(sg_crank_deg)
    tick_deg_sg = filter(d -> sg_lo <= d <= sg_hi, 0:360:2880)
    tick_labels_sg = ["$(Int(d))°\n$(fmt(d * ms_per_deg)) ms" for d in tick_deg_sg]

    p_sg = heatmap(sg_crank_deg, sg_orders[keep], spec_db[keep, :],
                   xlabel = "Crank angle [deg]  /  Time [ms]   (3 cycles tiled)",
                   ylabel = "Engine order",
                   color  = :viridis, clim = (-60.0, 0.0),
                   xticks = (tick_deg_sg, tick_labels_sg),
                   colorbar_title = "dB (rel. max)",
                   framestyle = :box,
                   left_margin = lmarg,
                   bottom_margin = 4Plots.mm,
                   title = "STFT spectrogram")

    plt = plot(p_time, p_spec, p_sg, layout = (3, 1), size = (960, 1100))

    return (omegas        = omegas,
            coefficients  = coeffs,
            omega_fund    = omega_fund,
            alpha_deg     = alpha_deg,
            signal        = sig,
            reconstructed = recon,
            spectrogram   = (magnitude  = spec,
                             crank_deg  = sg_crank_deg,
                             orders     = sg_orders),
            plot          = plt)
end

"""
    evaluate_fourier(t, omegas, coefficients) -> Float64

Evaluate f(t) = sum_n H_n * exp(i * omega_n * t) at a single time `t [s]`,
returning the real part (the imaginary part is numerical noise for the
Hermitian-symmetric two-sided arrays returned by `fourier_decompose`).
"""
evaluate_fourier(t, omegas, coefficients) =
    real(sum(coefficients .* exp.(1im .* omegas .* t)))

"""
    spectrum(result) -> (orders, magnitudes)

Convenience: convert a `fourier_decompose` result into a single-sided
amplitude spectrum (engine orders vs |2 H_n| for n >= 1, |H_0| for n = 0).
`orders[n+1] = n / 2` because omega_fund = omega_engine / 2.
"""
function spectrum(result)
    nmax = (length(result.omegas) - 1) ÷ 2
    orders = collect(0:nmax) ./ 2                       # engine orders
    mags   = Float64[ n == 0 ? abs(result.coefficients[nmax + 1]) :
                                2 * abs(result.coefficients[nmax + 1 + n])
                      for n in 0:nmax ]
    return orders, mags
end

# --- example -----------------------------------------------------------------
# include("fourier_torque.jl")
# r = fourier_decompose(2000, 24)                  # 24 harmonics, like the paper
# display(r.plot)
# savefig(r.plot, "imgs/fourier_torque_2000rpm.png")
# evaluate_fourier(0.005, r.omegas, r.coefficients)   # F_t at t = 5 ms
# orders, mags = spectrum(r)
#fourier_decompose(2000, 24)                 # vertical component (default)
#fourier_decompose(2000, 24; signal=:F_t_x)  # horizontal component
#fourier_decompose(2000, 24; signal=:F_t)    # scalar tangential 
#fourier_decompose(2000, 24; signal=:M_t)    # excitation moment