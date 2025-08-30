# -*- coding: utf-8 -*-
# =============================================================================
# Emergent α (Path-I) — Paper-aligned naming with the closed root:
#     𝓕(α; Λ, {L_{2m}}) := 𝓓_C(α)  −  𝓖_ind^ref(Λ)  =  0
# -----------------------------------------------------------------------------
# What this file does (CLOSED; no fitting, no physical constants):
#   • Build once-only geometric blocks: UV/IR/OUT, TT–χ kernel, and spectrum {K, L_{2m}}
#   • Iterate:
#       [α-step]  Solve 𝓓_C(α) = 𝓖_ind^ref(Λ)            (bisection, no feedback locking)
#       [Λ-step]  Update Λ ← Λ_base + ∆Λ_out^(sync)(α)   (TT–χ, no-shear; optional ladders)
#   • Repeat until α converges.
#
# All names match the paper:
#   𝓓_C(α), 𝓖_ind(Λ), 𝓖_ind^ref(Λ)=𝓖_ind(Λ)+ε_Λ,  C_log = (π²/𝓓_C)·ζ·(1+ζ),  ζ=(K/2π²)Λ
# No “S3 stage” or “δ_*” names appear.
# =============================================================================

from __future__ import annotations
import mpmath as mp

pi = mp.pi

# =========================== USER KNOBS (safe defaults) =======================
mp.mp.dps = 90  # arbitrary precision (60–150 typical)

# Geometry (Path-I baseline: ε=1/√π, η=1/π, ℓ0=εη)
EPSILON = 1 / mp.sqrt(pi)
ETA0    = 1 / pi
ELL0    = EPSILON * ETA0

# OUT evaluator
OUT_MODE  = 'exact'   # 'exact' (preferred) | 'series' | 'dipole'
OUT_LMAX  = 19        # odd-only enforced; try 19–41 (odd)
GL_NODES  = 512       # Gauss–Legendre nodes for OUT 'exact' (256–1024 typical)

# Spectral depth / tolerances for {K, L_{2m}}
SPEC_M_MAX = 20                  # include L_{2m} up to m=SPEC_M_MAX (≥2)
SPEC_TOL   = mp.mpf('1e-40')     # spectral tail cutoff for K, L_{2m}

# Curvature series for γ_geom
CURV_SERIES_ON    = True
CURV_SERIES_ORDER = 12           # η²/6 + η⁴/120 + ... up to this even power (4–12 typical)

# Iteration / solve knobs
ALPHA_TOL_REL = mp.mpf('5e-25')  # relative α tolerance for convergence (1e-18–1e-22 typical)
MAX_ITERS     = 20
BISECT_TOL    = None             # None → auto from dps

# Optional α reference (ONLY for ppb logging; not used in solve)
ALPHA_REF = mp.mpf('7.297352564311e-3')  # CODATA2022; set None to disable

# Optional ladders in sync
CHI_LADDER_ON    = True
CHI_LADDER_MODE  = 'closed'      # 'closed' | 'series'
CHI_LADDER_TERMS = 10

SELF_LADDER_ON    = True
SELF_LADDER_MODE  = 'closed'     # 'closed' | 'series'
SELF_LADDER_TERMS = 10

# Paper’s remainder knobs (keep zero for CLOSED runs)
DELTA_LAMBDA_DYN = mp.mpf('0')   # ∆Λ_dyn (Eq. 56) — user-injected dynamic remainder
EPS_GIND         = mp.mpf('0')   # ε_Λ     (Eq. 86) — reference offset on 𝓖_ind^ref

# =========================== GEOMETRIC CONSTANTS ==============================
C0_UNI     = (1 / pi) * (mp.mpf('4') / 3 + 1 / (4 * pi**2))      # C0^uni
C0_GAUSS   = mp.mpf('0.5') * (mp.log(2) + mp.euler)              # ½(ln2+γ)
LAMBDA_IND = mp.log(8 * mp.sqrt(pi)) - 2                         # Λ_ind (a=R/√π)
ITOT       = mp.mpf('1') / 6 - mp.mpf('1') / (4 * pi**2)         # ∫_0^1 x^2 sin^2(πx) dx

# =========================== IR (TT–χ) KERNEL =================================
def f_swirl(x: mp.mpf, ell: mp.mpf) -> mp.mpf:
    """Collar swirl profile: f(x) = (1-x)^2 / [(1-x)^2 + ell^2]."""
    return ((1 - x)**2) / (((1 - x)**2) + ell**2)

def P_IR_chi(ell: mp.mpf) -> mp.mpf:
    """
    TT–χ overlap on the collar (Path-I), normalized by ITOT:
      P^(IR)_χ(ell) = [ ∫ w(x)·(1 - f/3)·e^{-((1-x)/ell)^2} dx ] / ITOT,  w=x^2 sin^2(πx).
    """
    w = lambda x: (x**2) * (mp.sin(pi * x)**2)
    num = mp.quad(lambda x: w(x) * (1 - mp.mpf('1')/3 * f_swirl(x, ell)) * mp.e**(-((1 - x)/ell)**2), [0, 1])
    return num / ITOT

# =========================== OUTER FIELD (Λ_OUT) ==============================
_GL_CACHE = {}

def _gauss_legendre(n: int):
    if n in _GL_CACHE: return _GL_CACHE[n]
    xs, ws = [], []
    tol = mp.mpf(10)**(-(mp.mp.dps - 8))
    for k in range(1, n + 1):
        x = mp.cos(pi * (k - mp.mpf('0.25')) / (n + mp.mpf('0.5')))
        for _ in range(25):
            Pn  = mp.legendre(n, x)
            dPn = n / (1 - x**2) * (mp.legendre(n - 1, x) - x * Pn)
            dx  = -Pn / dPn
            x  += dx
            if abs(dx) < tol: break
        w = 2 / ((1 - x**2) * (dPn**2))
        xs.append(x); ws.append(w)
    _GL_CACHE[n] = (xs, ws); return xs, ws

def _B_rho(rho: mp.mpf, z: mp.mpf) -> mp.mpf:
    rc = mp.sqrt((1 + rho)**2 + z**2)
    k2 = 4 * rho / ((1 + rho)**2 + z**2)
    if k2 <= 0: return mp.mpf('0')
    if k2 >= 1: k2 = mp.mpf('1') - mp.mpf('1e-30')
    K_ = mp.ellipk(k2); E_ = mp.ellipe(k2)
    denom = (1 - rho)**2 + z**2
    if rho == 0: return mp.mpf('0')
    return (z / (2 * pi * rho * rc)) * (-K_ + ((1 + rho**2 + z**2) / denom) * E_)

def _B_z(rho: mp.mpf, z: mp.mpf) -> mp.mpf:
    rc = mp.sqrt((1 + rho)**2 + z**2)
    k2 = 4 * rho / ((1 + rho)**2 + z**2)
    if k2 <= 0: return mp.mpf('1') / (2 * rc**3)
    if k2 >= 1: k2 = mp.mpf('1') - mp.mpf('1e-30')
    K_ = mp.ellipk(k2); E_ = mp.ellipe(k2)
    denom = (1 - rho)**2 + z**2
    return (1 / (2 * pi * rc)) * (K_ + ((1 - rho**2 - z**2) / denom) * E_)

def _Btheta_on_sphere_x(x: mp.mpf, rstar: mp.mpf) -> mp.mpf:
    s   = mp.sqrt(1 - x**2); rho = rstar * s; z = rstar * x
    return _B_rho(rho, z) * x - _B_z(rho, z) * s

def _dPdx_leg(l: int, x: mp.mpf) -> mp.mpf:
    if l == 0: return mp.mpf('0')
    Pl  = mp.legendre(l, x); Pl1 = mp.legendre(l - 1, x)
    return (l / (1 - x**2)) * (Pl1 - x * Pl)

def Lambda_OUT_exact(eta: mp.mpf, lmax: int = OUT_LMAX, gl_nodes: int = GL_NODES) -> mp.mpf:
    if OUT_MODE == 'dipole':  # coarse approximation
        return - (pi / 6) * (eta**3)
    rstar = 1 / eta
    xs, ws = _gauss_legendre(gl_nodes)
    Uout = mp.mpf('0')
    Lm = lmax if (lmax % 2 == 1) else (lmax - 1)
    for l in range(1, Lm + 1, 2):
        Il = l * (l + 1) * 2 / (2 * l + 1)
        s  = mp.mpf('0')
        for x, w in zip(xs, ws):
            s += _Btheta_on_sphere_x(x, rstar) * (-(1 - x**2) * _dPdx_leg(l, x)) * w
        a_l  = - (rstar**(l + 2)) * s / Il
        Uout += ((l + 1) / (2 * l + 1)) * (a_l**2) * (rstar**(-(2 * l + 1)))
    Uout *= (2 * pi)
    return -2 * Uout

def Lambda_OUT_series(eta: mp.mpf, lmax: int = OUT_LMAX) -> mp.mpf:
    s  = mp.mpf('0')
    Lm = lmax if (lmax % 2 == 1) else (lmax - 1)
    for n in range(1, Lm + 1, 2):
        s += pi / ((n + 1) * (2*n + 1)) * (eta**(2*n + 1))
    return -s

def Lambda_OUT(eta: mp.mpf, lmax: int = OUT_LMAX, mode: str = 'exact') -> mp.mpf:
    if mode == 'series': return Lambda_OUT_series(eta, lmax=lmax)
    return Lambda_OUT_exact(eta, lmax=lmax, gl_nodes=GL_NODES)

def Lambda_OUT_extrapolated(eta: mp.mpf, lbase: int = 19, mode: str = 'exact') -> mp.mpf:
    S1 = Lambda_OUT(eta, lmax=lbase,     mode=mode)
    S2 = Lambda_OUT(eta, lmax=lbase + 2, mode=mode)
    S3 = Lambda_OUT(eta, lmax=lbase + 4, mode=mode)
    denom = (S3 - 2 * S2 + S1)
    if denom == 0 or abs(denom) < mp.mpf('1e-30') * max(1, abs(S3)): return S3
    Sout = S1 - (S2 - S1)**2 / denom
    if abs(Sout) > 10 * max(abs(S1), abs(S2), abs(S3)): return S3
    return Sout

# =========================== SPECTRUM (K, L_{2m}) =============================
def In_m1(n: int) -> mp.mpf:
    n = mp.mpf(n)
    return (-1)**(n - 1) / (((n - 1) * pi)**2) + (-1)**n / (((n + 1) * pi)**2)

def series_K(tol: mp.mpf = SPEC_TOL) -> mp.mpf:
    S, n = mp.mpf('0'), 2
    while True:
        t = (2 * In_m1(n))**2 / (n**2 - 1)
        S += t
        if n > 120 and abs(t) < tol: break
        n += 1
    return (2 / pi**2) * S

def _CS_pair(m: int, a: mp.mpf):
    if a == 0: return mp.mpf('1') / (m + 1), mp.mpf('0')
    C = mp.sin(a) / a; S = (1 - mp.cos(a)) / a
    if m == 0: return C, S
    for k in range(1, m + 1):
        C, S = mp.sin(a) / a - (k / a) * S, (1 - mp.cos(a)) / a + (k / a) * C
    return C, S

def I_nm(n: int, m2: int) -> mp.mpf:
    C1, _ = _CS_pair(m2, (n - 1) * pi); C2, _ = _CS_pair(m2, (n + 1) * pi)
    return mp.mpf('0.5') * (C1 - C2)

def L_2m(m: int, tol: mp.mpf = SPEC_TOL, nmax: int = 1200) -> mp.mpf:
    S = mp.mpf('0')
    for n in range(2, nmax + 1):
        t = (2 * I_nm(n, 2 * m))**2 / (n**2 - 1)
        S += t
        if n > 80 and abs(t) < tol: break
    return (2 / pi**2) * S

def precompute_spectrum(M: int = SPEC_M_MAX, tol: mp.mpf = SPEC_TOL):
    """Return K and list of (m, L_{2m}) for m=2..M (fixed across the solve)."""
    K = series_K(tol=tol)
    L_list = []
    if M >= 2:
        for m in range(2, M + 1):
            L_list.append((m, L_2m(m, tol=tol)))
    return K, L_list

# =========================== 𝓓_C(α)  &  C_log(α; ζ) ===========================
def DC_of_alpha_fixedK(alpha: mp.mpf, K: mp.mpf, L_list, M: int = SPEC_M_MAX) -> mp.mpf:
    r"""
    𝓓_C(α) = (α/π) √(1 − ξ) − (α/π) (ξ/2) K − (α/π) Σ_{m≥2} (ξ/2)^m L_{2m},  ξ = 2 C0_UNI α.
    """
    a  = mp.mpf(alpha)
    xi = 2 * C0_UNI * a
    D  = (a / pi) * mp.sqrt(1 - xi) - (a / pi) * (xi / 2) * K
    if M >= 2:
        for (m, Lm) in L_list:
            D -= (a / pi) * ((xi / 2)**m) * Lm
    return D

def C_log_from_DC_zeta(DC: mp.mpf, zeta: mp.mpf) -> mp.mpf:
    """C_log = (π²/𝓓_C)·ζ·(1+ζ)."""
    return (pi**2 / DC) * zeta * (1 + zeta)

# =========================== 𝓖_ind(Λ)  &  𝓖_ind^ref(Λ) ========================
def G_ind_from_Lambda(Lambda: mp.mpf, K: mp.mpf) -> mp.mpf:
    """𝓖_ind(Λ) = (3/2) K Λ (1 + KΛ / 2π²)."""
    return (mp.mpf('3') / 2) * K * Lambda * (1 + (K * Lambda) / (2 * pi**2))

def G_ind_ref_from_Lambda(Lambda: mp.mpf, K: mp.mpf, eps_g: mp.mpf = EPS_GIND) -> mp.mpf:
    """𝓖_ind^ref(Λ) = 𝓖_ind(Λ) + ε_Λ  (ε_Λ default 0 for CLOSED runs)."""
    return G_ind_from_Lambda(Lambda, K) + eps_g

# =========================== γ_eff (curvature + map) ==========================
def curvature_series_eta(eta: mp.mpf, order: int = CURV_SERIES_ORDER) -> mp.mpf:
    if not CURV_SERIES_ON:
        return eta**2 / 6
    if order < 2:
        return mp.mpf('0')
    s = mp.mpf('0'); max_k = int(order // 2)
    for k in range(1, max_k + 1):
        s += (eta**(2 * k)) / mp.factorial(2 * k + 1)
    return s

def gamma_geom(eta: mp.mpf) -> mp.mpf:
    return mp.mpf('0.5') * (1 + curvature_series_eta(eta, order=CURV_SERIES_ORDER))

def gamma_eff(eta: mp.mpf, K: mp.mpf, DC_lock: mp.mpf, P_ir: mp.mpf) -> mp.mpf:
    """
    γ_eff = γ_geom + γ_map, with the map evaluated at the provided 𝓓_C.
    γ_map = (K / (2 𝓓_C)) · c0^Gauss · P^(IR)_χ.
    """
    return gamma_geom(eta) + (K / (2 * DC_lock)) * C0_GAUSS * P_ir

# ----- Optional ladder corrections (χ-map and self) -----
def deltaLambda_chi_ladder_extra(eta: mp.mpf, K: mp.mpf, DC_lock: mp.mpf, P_ir: mp.mpf,
                                 dLambda_OUT: mp.mpf, mode: str = 'closed', terms: int = 10) -> mp.mpf:
    k = mp.sinh(eta) / eta - 1
    x = (K / (2 * DC_lock)) * C0_GAUSS * P_ir
    if mode == 'closed':
        return (-k * x**2 / (1 + k * x)) * P_ir * dLambda_OUT
    s = mp.mpf('0')
    for n in range(2, max(2, int(terms)) + 1):
        s += ((-k)**(n - 1)) * (x**n)
    return s * P_ir * dLambda_OUT

def deltaLambda_self_ladder(eta: mp.mpf, K: mp.mpf, P_ir: mp.mpf, Lambda_eff: mp.mpf,
                            alpha: mp.mpf, mode: str = 'closed', terms: int = 10) -> mp.mpf:
    ggeom = gamma_geom(eta)
    k = mp.sinh(eta) / eta - 1
    ep = (alpha / mp.pi) * (K / (2 * mp.pi**2)) * P_ir * Lambda_eff
    if mode == 'closed':
        return - ggeom * P_ir * Lambda_eff * (ep / (1 + k * ep))
    s = mp.mpf('0')
    for n in range(2, max(2, int(terms)) + 1):
        s += ((-k)**(n - 2)) * (ep**(n - 1))
    return - ggeom * P_ir * Lambda_eff * s

# =========================== Λ_base & Λ_eff ===================================
def build_Lambda_base(ell: mp.mpf = ELL0, eta: mp.mpf = ETA0):
    """
    Λ_base = Λ_ind + ∆Λ_UV→IR + ∆Λ_OUT  (no sync terms).
    Returns: (Λ_base, P^(IR)_χ, ∆Λ_UV→IR, ∆Λ_OUT)
    """
    P_ir   = P_IR_chi(ell)
    dUVIR  = C0_GAUSS * P_ir
    dOUT   = Lambda_OUT_extrapolated(eta, lbase=(OUT_LMAX if OUT_LMAX % 2 == 1 else OUT_LMAX - 1), mode=OUT_MODE)
    Lambda_base = LAMBDA_IND + dUVIR + dOUT
    return Lambda_base, P_ir, dUVIR, dOUT

# =========================== Root solve helpers ===============================
def auto_bracket(F, a0, a_min: mp.mpf = mp.mpf('2e-3'), a_max: mp.mpf = mp.mpf('2e-2'),
                 w0: mp.mpf = mp.mpf('0.15'), grow: mp.mpf = mp.mpf('1.7'), max_expand: int = 48):
    lo = max(a_min, a0 * (1 - w0)); hi = min(a_max, a0 * (1 + w0))
    f_lo, f_hi = F(lo), F(hi)
    if f_lo == 0: return lo, lo
    if f_hi == 0: return hi, hi
    for _ in range(max_expand):
        if f_lo * f_hi < 0: return lo, hi
        lo = max(a_min, lo / grow); hi = min(a_max, hi * grow)
        f_lo, f_hi = F(lo), F(hi)
    raise RuntimeError("Auto-bracket failed.")

# Globals for logging convenience
_K_GLOBAL, _L_LIST = None, None

def DC_at_alpha(alpha, K_fixed=None):
    K = _K_GLOBAL if K_fixed is None else K_fixed
    D = DC_of_alpha_fixedK(alpha, K, _L_LIST, M=SPEC_M_MAX)
    return D, K

def bisection(F, lo, hi, tol=None, alpha_ref=None, Lambda=None, K=None):
    if tol is None:
        tol = mp.mpf(10)**(-min(60, mp.mp.dps - 12))
    print("\n[it]           alpha_mid         Δα(ppb)         F(mid)            𝓓_C(mid)          C_log(mid)")
    it = 0
    f_lo, f_hi = F(lo), F(hi)
    while (hi - lo) > tol * max(1, abs(lo), abs(hi)):
        it += 1
        mid   = (lo + hi) / 2
        f_mid = F(mid)
        # diagnostics:
        DC_mid = DC_of_alpha_fixedK(mid, K, _L_LIST, M=SPEC_M_MAX)
        zeta   = (K / (2 * pi**2)) * Lambda
        Cmid   = (pi**2 / DC_mid) * zeta * (1 + zeta)
        err_ppb = (mid - alpha_ref) / alpha_ref * mp.mpf('1e9') if alpha_ref is not None else mp.mpf('nan')
        print(f"[{it:02d}] {mp.nstr(mid, 18)} {mp.nstr(err_ppb, 12)} {mp.nstr(f_mid, 14)} {mp.nstr(DC_mid, 14)} {mp.nstr(Cmid, 14)}")
        if f_mid == 0:
            return mid
        if f_lo * f_mid < 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    return (lo + hi) / 2

# =========================== MAIN CLOSED ITERATION ============================
def main():
    global _K_GLOBAL, _L_LIST

    print("=== Emergent α — Paper-aligned naming (𝓕 = 𝓓_C − 𝓖_ind^ref = 0) ===")
    print(f"mp.mp.dps       = {mp.mp.dps}")
    print(f"Geometry: ε={mp.nstr(EPSILON, 12)}, η0={mp.nstr(ETA0, 12)}, ℓ0={mp.nstr(ELL0, 12)}")
    print(f"OUT: mode={OUT_MODE}, Lmax={OUT_LMAX} (odd-only), GL_NODES={GL_NODES}")
    print(f"Spectrum: M_max={SPEC_M_MAX}, tol={mp.nstr(SPEC_TOL, 2)}")
    print(f"Curvature: series_on={CURV_SERIES_ON}, order={CURV_SERIES_ORDER}")

    # Spectrum (fixed)
    K, L_list = precompute_spectrum(M=SPEC_M_MAX, tol=SPEC_TOL)
    _K_GLOBAL, _L_LIST = K, L_list
    print("\n-- Spectrum --")
    print(f"K (spectral)      = {mp.nstr(K, 22)}")
    if L_list:
        print(f"L_2m (m=2..{SPEC_M_MAX}) computed ({len(L_list)} terms)")

    # Λ_base (fixed)
    Lambda_base, P_ir, dUVIR, dOUT = build_Lambda_base(ell=ELL0, eta=ETA0)
    print("\n-- Geometry blocks --")
    print(f"P^(IR)_χ(ℓ0)      = {mp.nstr(P_ir, 22)}")
    print(f"∆Λ^(UV→IR)        = {mp.nstr(dUVIR, 22)}")
    print(f"∆Λ_OUT(η0,Lmax)   = {mp.nstr(dOUT, 22)}")
    print(f"Λ_ind             = {mp.nstr(LAMBDA_IND, 22)}")
    print(f"Λ_base            = {mp.nstr(Lambda_base, 22)}")

    # Initialize Λ_eff and α guess
    Lambda_eff = Lambda_base
    G_ind_ref  = G_ind_ref_from_Lambda(Lambda_eff, K, eps_g=EPS_GIND)
    alpha_guess = pi * G_ind_ref  # linearized guess α ≈ π·𝓖_ind^ref
    print(f"\n𝓖_ind^ref(Λ_eff)  = {mp.nstr(G_ind_ref, 22)}")
    print(f"α guess (≈π·𝓖_ref) = {mp.nstr(alpha_guess, 18)}")

    print("\n[it]    alpha_k              Δα(ppb)           𝓖_ind^ref(Λ_k)    𝓓_C(α_k)         ∆Λ_sync           Λ_eff(k)          C_log(α_k)")
    alpha_prev = None

    for it in range(1, MAX_ITERS + 1):
        # --- α-step: solve 𝓕(α; Λ_eff) = 0 with fixed Λ_eff
        G_ind_ref = G_ind_ref_from_Lambda(Lambda_eff, K, eps_g=EPS_GIND)
        F = lambda a: DC_of_alpha_fixedK(a, K, L_list, M=SPEC_M_MAX) - G_ind_ref
        lo, hi = auto_bracket(F, alpha_guess)
        alpha_k = bisection(F, lo, hi, tol=BISECT_TOL, alpha_ref=ALPHA_REF, Lambda=Lambda_eff, K=K)

        # Diagnostics at α_k
        DC_k   = DC_of_alpha_fixedK(alpha_k, K, L_list, M=SPEC_M_MAX)
        zeta_k = (K / (2 * pi**2)) * Lambda_eff
        Clog_k = C_log_from_DC_zeta(DC_k, zeta_k)
        err_ppb = (alpha_k - ALPHA_REF) / ALPHA_REF * mp.mpf('1e9') if ALPHA_REF is not None else mp.mpf('nan')

        # --- Λ-step: update Λ via sync (TT–χ, no-shear base, plus ladders if enabled)
        gamma_k = gamma_eff(ETA0, K, DC_k, P_ir)
        DeltaLambda_sync = gamma_k * P_ir * dOUT

        if CHI_LADDER_ON:
            dL_extra = deltaLambda_chi_ladder_extra(ETA0, K, DC_k, P_ir, dOUT, mode=CHI_LADDER_MODE, terms=CHI_LADDER_TERMS)
            print("∆Λ_sync  → χ-ladder extra", mp.nstr(dL_extra, 18))
            DeltaLambda_sync += dL_extra

        if SELF_LADDER_ON:
            dL_self = deltaLambda_self_ladder(ETA0, K, P_ir, Lambda_eff, alpha_k, mode=SELF_LADDER_MODE, terms=SELF_LADDER_TERMS)
            print("∆Λ_sync  → self-ladder   ", mp.nstr(dL_self, 18))
            DeltaLambda_sync += dL_self

        # New Λ_eff (add dynamic remainder only if user sets it explicitly)
        Lambda_eff = Lambda_base + DeltaLambda_sync
        if abs(DELTA_LAMBDA_DYN) > 0:
            Lambda_eff += DELTA_LAMBDA_DYN

        print(f"[{it:02d}] {mp.nstr(alpha_k, 18)} {mp.nstr(err_ppb, 12)} {mp.nstr(G_ind_ref, 14)} {mp.nstr(DC_k, 14)} "
              f"{mp.nstr(DeltaLambda_sync, 14)} {mp.nstr(Lambda_eff, 14)} {mp.nstr(Clog_k, 14)}")

        # Convergence check
        if alpha_prev is not None:
            if abs(alpha_k - alpha_prev) <= ALPHA_TOL_REL * max(1, abs(alpha_prev)):
                print("\nConverged: |∆α|/α < tol."); break
        alpha_prev  = alpha_k
        alpha_guess = alpha_k

    # Final report
    print("\n-- Final solution --")
    print(f"alpha_emergent     = {mp.nstr(alpha_k, 22)}")
    print(f"alpha_em^-1        = {mp.nstr(1 / alpha_k, 16)}")
    print(f"Λ_eff (final)      = {mp.nstr(Lambda_eff, 22)}")
    print(f"K (spectral)       = {mp.nstr(K, 22)}")
    print(f"P^(IR)_χ(ℓ0)       = {mp.nstr(P_ir, 22)}")
    print(f"∆Λ_OUT (η0)        = {mp.nstr(dOUT, 22)}")
    print(f"∆Λ^(UV→IR)         = {mp.nstr(dUVIR, 22)}")
    zeta_fin = (K / (2 * pi**2)) * Lambda_eff
    Clog_fin = C_log_from_DC_zeta(DC_of_alpha_fixedK(alpha_k, K, L_list, M=SPEC_M_MAX), zeta_fin)
    print(f"C_log(α_em)        = {mp.nstr(Clog_fin, 18)}   (∆ vs 1/3 = {mp.nstr(Clog_fin - mp.mpf('1')/3, 12)})")

    if ALPHA_REF is not None:
        err_ppb = (alpha_k - ALPHA_REF) / ALPHA_REF * mp.mpf('1e9')
        print(f"[Context] α_ref    = {mp.nstr(ALPHA_REF, 16)}  →  ∆α(ppb) = {mp.nstr(err_ppb, 12)}")

    if DELTA_LAMBDA_DYN:
        print(f"\n⚠️  DELTA_LAMBDA_DYN is nonzero ({DELTA_LAMBDA_DYN}). Results are not strictly CLOSED.")

if __name__ == "__main__":
    main()
