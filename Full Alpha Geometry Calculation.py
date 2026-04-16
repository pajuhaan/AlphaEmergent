
from __future__ import annotations

"""
Complete reviewer-facing geometric audit for the minimal electron branch:

    Λ_geom  ->  D_C^lock(Λ_geom)  ->  α_pred
             ->  pure-photonic g-2 coefficient sequence
             ->  reference-only QED and benchmark comparisons

Design principles
-----------------
1) Primary layer is geometry/scalar only.
   No benchmark α and no QED coefficient enters:
       - Λ_geom
       - D_C^lock(Λ_geom)
       - article/current scalar inversion α_pred
       - derived Relator g-2 coefficients

2) Reference-only layer is quarantined.
   Benchmark α and pure-photonic QED coefficients are used only in:
       - external comparison tables
       - QED-induced D_C cross-check
       - coefficient-error ledgers

3) One-file reproducibility.
   Every numerical table is computed inside this file and the stdout report is
   written byte-for-byte to a UTF-8 text file.

Manuscript convention note
--------------------------
The uploaded TeX source contains two slightly different refined-kernel
conventions across its scalar and g-2 sections. To preserve reproducibility,
this script keeps both conventions available:

    refined (scalar-table convention): kernel scale = λ_u
    refined (realized g-2 convention): kernel scale = λ_u²

with λ_u = (1 - ρ_dyn²)^(-1), matching the numerical table convention in the
uploaded manuscript. The current article lock α_pred remains on the current
rank-5 scalar evaluator and is unaffected by this choice.
"""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence
import argparse
import io
import os

import mpmath as mp
from tabulate import tabulate

# =============================================================================
# Precision, controls, and output settings
# =============================================================================

DEFAULT_DPS = 100
DEFAULT_MAX_ODD_MODE = 101
DEFAULT_GL_NODES = 512
DEFAULT_PRINT_DIGITS = 40
LOW_MODE_PRINT_MAX = 11

ARTICLE_RANK = 5
RANK_AUDIT_LIST = (1, 3, 5, 100)
QED_SERIES_ORDER = 5
PREDICTION_ORDER = 10

ROOT_TOL = mp.mpf("1e-95")
ROOT_MAX_ITERS = 1200
ALPHA_BRACKET_LEFT = mp.mpf("0.001")
ALPHA_BRACKET_RIGHT = mp.mpf("0.02")

# Refined-kernel powers used to mirror the manuscript numerics.
REFINED_SCALAR_TABLE_KERNEL_POWER = 1   # reproduces scalar c_n table
REFINED_REALIZED_KERNEL_POWER = 2       # reproduces refined g-2 tables

SENSITIVITY_EPS = mp.mpf("1e-8")

mp.mp.dps = DEFAULT_DPS

# =============================================================================
# Exact constants used by the primary derivation
# =============================================================================

ZERO = mp.mpf("0")
ONE = mp.mpf("1")
TWO = mp.mpf("2")
THREE = mp.mpf("3")
FOUR = mp.mpf("4")
FIVE = mp.mpf("5")
HALF = mp.mpf("0.5")
QUARTER = mp.mpf("0.25")

PI = mp.pi
SQRT_PI = mp.sqrt(PI)
GAMMA_E = mp.euler

EPSILON_0 = ONE / SQRT_PI
ETA_0 = ONE / PI
ELL_0 = EPSILON_0 * ETA_0

LAMBDA_IND = mp.log(8 * SQRT_PI) - TWO
C0_GAUSS = HALF * (mp.log(TWO) + GAMMA_E)
GAMMA_GEOM = HALF * mp.sinh(ETA_0) / ETA_0

A_UV = mp.log(TWO) * ETA_0**4
B_IR = ETA_0**4 / (8 * PI**2)
N_OMEGA = mp.sqrt(TWO / PI)

C_CHI = TWO * PI**2
S_UV = mp.log(2)
S_IR = ONE / (8 * PI**2)
C_UNI = (ONE / PI) * (mp.mpf("4") / THREE + ONE / (4 * PI**2))

# Orientation/benchmark point used only in comparison tables and rank audits.
ALPHA_INV_REF = mp.mpf("137.035999177")
ALPHA_INV_REF_SIGMA = mp.mpf("0.000000021")
ALPHA_REF = ONE / ALPHA_INV_REF
ALPHA_REF_SIGMA = ALPHA_INV_REF_SIGMA / (ALPHA_INV_REF**2)

# Pure-photonic universal QED coefficients A_1^(2n), used only in comparison.
QED_A1_UNIVERSAL: Dict[int, mp.mpf] = {
    1: mp.mpf("0.5"),
    2: mp.mpf("-0.328478965579193"),
    3: mp.mpf("1.181241456587"),
    4: mp.mpf("-1.9122457649264455741526471531265"),
    5: mp.mpf("5.891"),
}

# Optional full-electron coefficients, kept reference-only.
QED_AE_FULL: Dict[int, mp.mpf] = {
    1: mp.mpf("0.5"),
    2: mp.mpf("-0.328478444002617"),
    3: mp.mpf("1.181234016818"),
    4: mp.mpf("-1.91132213891344557"),
    5: mp.mpf("6.08"),
}

def geometric_shell_constant_k() -> mp.mpf:
    return (150 * PI**2 - 8 * PI**4 - 315) / (180 * PI**6)

K_GEOMETRIC = geometric_shell_constant_k()

# =============================================================================
# Formatting and table helpers
# =============================================================================

PRINT_DIGITS = DEFAULT_PRINT_DIGITS
SCI_DIGITS = 12

def fmt(x: mp.mpf | int | float, digits: int | None = None) -> str:
    return mp.nstr(mp.mpf(x), n=digits or PRINT_DIGITS)

def fmt_sci(x: mp.mpf | int | float, digits: int = SCI_DIGITS) -> str:
    return mp.nstr(mp.mpf(x), n=digits, min_fixed=0, max_fixed=0)

def fmt_pct(x: mp.mpf | int | float) -> str:
    return f"{float(x):+.6f}%"

def fmt_ppt(delta_alpha: mp.mpf, denominator: mp.mpf) -> str:
    value = mp.mpf("1e12") * delta_alpha / denominator
    return f"{float(value):+.6f}"

def render_table(title: str, headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    rendered = tabulate(
        rows,
        headers=headers,
        tablefmt="rounded_grid",
        stralign="left",
        numalign="right",
        disable_numparse=True,
    )
    return f"\n{title}\n{rendered}\n"

# =============================================================================
# Vector channel; TT–χ shell admission quotient
# =============================================================================

def f_swirl(x: mp.mpf, ell: mp.mpf = ELL_0) -> mp.mpf:
    return (ONE - x) ** 2 / ((ONE - x) ** 2 + ell**2)

def shell_weight(x: mp.mpf) -> mp.mpf:
    return x**2 * mp.sin(PI * x) ** 2

def p_ir_denominator_exact() -> mp.mpf:
    return ONE / 6 - ONE / (4 * PI**2)

def p_ir_numerator(ell: mp.mpf = ELL_0) -> mp.mpf:
    def integrand(x: mp.mpf) -> mp.mpf:
        return shell_weight(x) * (ONE - f_swirl(x, ell) / THREE) * mp.exp(-((ONE - x) / ell) ** 2)
    return mp.quad(integrand, [ZERO, ONE])

def p_ir_closed_rayleigh(ell: mp.mpf = ELL_0) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    numerator = p_ir_numerator(ell)
    denominator = p_ir_denominator_exact()
    return numerator, denominator, numerator / denominator

def q_admission_eigenvalue(j: int, ell: mp.mpf = ELL_0) -> mp.mpf:
    jj = mp.mpf(j * (j + 1))
    return N_OMEGA * (ONE - (ell**2 * jj) / (THREE * (ONE + ell**2 * jj)))

# =============================================================================
# Vector channel; exact Maxwell shell source and Λ_OUT
# =============================================================================

def b_rho_mp(rho: mp.mpf, z: mp.mpf) -> mp.mpf:
    if abs(rho) < mp.mpf("1e-60"):
        return ZERO
    rc = mp.sqrt((ONE + rho) ** 2 + z**2)
    k2 = (4 * rho) / ((ONE + rho) ** 2 + z**2)
    if k2 >= ONE:
        k2 = ONE - mp.mpf("1e-60")
    K = mp.ellipk(k2)
    E = mp.ellipe(k2)
    denom = (ONE - rho) ** 2 + z**2
    return (z / (TWO * PI * rho * rc)) * (-K + ((ONE + rho**2 + z**2) / denom) * E)

def b_z_mp(rho: mp.mpf, z: mp.mpf) -> mp.mpf:
    rc = mp.sqrt((ONE + rho) ** 2 + z**2)
    k2 = (4 * rho) / ((ONE + rho) ** 2 + z**2)
    if k2 <= ZERO:
        return ONE / (TWO * rc**3)
    if k2 >= ONE:
        k2 = ONE - mp.mpf("1e-60")
    K = mp.ellipk(k2)
    E = mp.ellipe(k2)
    denom = (ONE - rho) ** 2 + z**2
    return (ONE / (TWO * PI * rc)) * (K + ((ONE - rho**2 - z**2) / denom) * E)

def b_theta_on_shell_mp(x: mp.mpf, eta: mp.mpf = ETA_0) -> mp.mpf:
    r_star = ONE / eta
    s = mp.sqrt(max(ZERO, ONE - x * x))
    rho = r_star * s
    z = r_star * x
    return x * b_rho_mp(rho, z) - s * b_z_mp(rho, z)

def toroidal_basis_T(j: int, x: mp.mpf) -> mp.mpf:
    return mp.mpf(j) * (mp.legendre(j - 1, x) - x * mp.legendre(j, x))

_GL_CACHE: dict[int, tuple[list[mp.mpf], list[mp.mpf]]] = {}

def gauss_legendre_nodes(n: int) -> tuple[list[mp.mpf], list[mp.mpf]]:
    if n in _GL_CACHE:
        return _GL_CACHE[n]
    xs: list[mp.mpf] = []
    ws: list[mp.mpf] = []
    tol = mp.mpf(10) ** (-(mp.mp.dps - 12))
    for k in range(1, n + 1):
        x = mp.cos(PI * (mp.mpf(k) - mp.mpf("0.25")) / (mp.mpf(n) + mp.mpf("0.5")))
        for _ in range(80):
            pn = mp.legendre(n, x)
            dpn = n / (ONE - x**2) * (mp.legendre(n - 1, x) - x * pn)
            dx = -pn / dpn
            x += dx
            if abs(dx) < tol:
                break
        w = TWO / ((ONE - x**2) * dpn**2)
        xs.append(x)
        ws.append(w)
    _GL_CACHE[n] = (xs, ws)
    return xs, ws

def mode_integral(j: int, eta: mp.mpf = ETA_0, gl_nodes: int = DEFAULT_GL_NODES) -> mp.mpf:
    xs, ws = gauss_legendre_nodes(gl_nodes)
    return mp.fsum(b_theta_on_shell_mp(x, eta) * toroidal_basis_T(j, x) * w for x, w in zip(xs, ws))

def split_chunks(items: list[int], n_chunks: int) -> list[list[int]]:
    if n_chunks <= 1:
        return [items]
    chunks = [[] for _ in range(min(n_chunks, len(items)))]
    for i, item in enumerate(items):
        chunks[i % len(chunks)].append(item)
    return [chunk for chunk in chunks if chunk]

def _mode_worker_chunk(args: tuple[list[int], str, str, int, int]) -> list[tuple[int, str, str, str, str, str]]:
    js, eta_s, beta_s, dps, gl_nodes = args
    mp.mp.dps = dps
    eta = mp.mpf(eta_s)
    _ = mp.mpf(beta_s)  # beta is not needed in the projection itself; included for interface symmetry.
    r_star = ONE / eta
    xs, ws = gauss_legendre_nodes(gl_nodes)
    out: list[tuple[int, str, str, str, str, str]] = []
    for j in js:
        I_j = TWO * j * (j + 1) / (TWO * j + ONE)
        projection = mp.fsum(
            b_theta_on_shell_mp(x, eta) * toroidal_basis_T(j, x) * w
            for x, w in zip(xs, ws)
        )
        a_sh = (r_star ** (j + 2)) * projection / I_j
        a_hat = TWO * (j + 1) * (r_star ** (-(j + HALF))) * a_sh
        J = mp.sqrt(TWO * PI / ((j + 1) * (2 * j + 1))) * a_hat
        out.append(
            (
                j,
                mp.nstr(I_j, 120),
                mp.nstr(projection, 120),
                mp.nstr(a_sh, 120),
                mp.nstr(a_hat, 120),
                mp.nstr(J, 120),
            )
        )
    return out

@dataclass(frozen=True)
class ModeData:
    j: int
    I_j: mp.mpf
    projection: mp.mpf
    a_sh: mp.mpf
    a_hat: mp.mpf
    J: mp.mpf
    weight: mp.mpf
    q_j: mp.mpf
    phi_j: mp.mpf
    R_j: mp.mpf
    g_j: mp.mpf

@dataclass(frozen=True)
class SourceData:
    modes: list[ModeData]
    norm_sq: mp.mpf
    lambda_out: mp.mpf
    lambda_out_from_norm: mp.mpf
    norm_identity_error: mp.mpf
    weight_sum: mp.mpf

@dataclass(frozen=True)
class JacobiData:
    m1: mp.mpf
    m2: mp.mpf
    m3: mp.mpf
    a0: mp.mpf
    b1_sq: mp.mpf
    b1: mp.mpf
    a1: mp.mpf
    R_rank2: mp.mpf
    delta_R: mp.mpf
    delta_lambda: mp.mpf

@dataclass(frozen=True)
class LambdaGeomResult:
    p_num: mp.mpf
    p_den: mp.mpf
    p_ir: mp.mpf
    lambda_uv_to_ir: mp.mpf
    beta: mp.mpf
    source: SourceData
    F_chi: mp.mpf
    R_chi: mp.mpf
    exterior_memory: mp.mpf
    lambda_geom: mp.mpf
    jacobi: JacobiData

def beta_cur(p_ir: mp.mpf, lambda_uv_to_ir: mp.mpf) -> mp.mpf:
    return p_ir * GAMMA_GEOM + (p_ir * lambda_uv_to_ir) / (TWO * (ONE + A_UV))

def compute_source_data(
    eta: mp.mpf,
    max_odd_mode: int,
    beta: mp.mpf,
    workers: int,
    dps: int,
    gl_nodes: int,
) -> SourceData:
    odd_modes = list(range(1, max_odd_mode + 1, 2))
    chunks = split_chunks(odd_modes, workers)
    jobs = [(chunk, mp.nstr(eta, 120), mp.nstr(beta, 120), dps, gl_nodes) for chunk in chunks]

    if workers <= 1:
        nested = [_mode_worker_chunk(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            nested = list(executor.map(_mode_worker_chunk, jobs))

    raw_strings = [row for chunk in nested for row in chunk]
    raw = []
    norm_sq = ZERO
    lambda_out = ZERO
    for j, I_s, proj_s, a_sh_s, a_hat_s, J_s in sorted(raw_strings, key=lambda row: row[0]):
        I_j = mp.mpf(I_s)
        projection = mp.mpf(proj_s)
        a_sh = mp.mpf(a_sh_s)
        a_hat = mp.mpf(a_hat_s)
        J = mp.mpf(J_s)
        raw.append((j, I_j, projection, a_sh, a_hat, J))
        norm_sq += J**2
        lambda_out += -PI * a_hat**2 / ((j + 1) * (2 * j + 1))

    modes: list[ModeData] = []
    weight_sum = ZERO
    for j, I_j, projection, a_sh, a_hat, J in raw:
        weight = J**2 / norm_sq if norm_sq != ZERO else ZERO
        m = mp.mpf((j - 1) // 2)
        q_j = q_admission_eigenvalue(j, ELL_0)
        phi_j = (ONE + HALF * B_IR * m) / (ONE + A_UV * m)
        R_j = ONE + beta * phi_j
        g_j = ONE - ONE / R_j
        weight_sum += weight
        modes.append(ModeData(j, I_j, projection, a_sh, a_hat, J, weight, q_j, phi_j, R_j, g_j))

    lambda_out_from_norm = -HALF * norm_sq
    return SourceData(
        modes=modes,
        norm_sq=norm_sq,
        lambda_out=lambda_out,
        lambda_out_from_norm=lambda_out_from_norm,
        norm_identity_error=lambda_out - lambda_out_from_norm,
        weight_sum=weight_sum,
    )

def build_jacobi_audit(source: SourceData, R_chi: mp.mpf) -> JacobiData:
    m1 = mp.fsum(mode.weight * mode.g_j for mode in source.modes)
    m2 = mp.fsum(mode.weight * mode.g_j**2 for mode in source.modes)
    m3 = mp.fsum(mode.weight * mode.g_j**3 for mode in source.modes)
    a0 = m1
    b1_sq = m2 - m1**2
    b1 = mp.sqrt(max(ZERO, b1_sq))
    a1 = (m3 - TWO * m1 * m2 + m1**3) / b1_sq if b1_sq != ZERO else mp.nan
    R_rank2 = ONE / (ONE - a0 - b1_sq / (ONE - a1))
    delta_R = R_chi - R_rank2
    delta_lambda = source.lambda_out * delta_R
    return JacobiData(m1, m2, m3, a0, b1_sq, b1, a1, R_rank2, delta_R, delta_lambda)

def evaluate_lambda_geom(max_odd_mode: int, workers: int, dps: int, gl_nodes: int) -> LambdaGeomResult:
    p_num, p_den, p_ir = p_ir_closed_rayleigh(ELL_0)
    lambda_uv_to_ir = C0_GAUSS * p_ir
    beta = beta_cur(p_ir, lambda_uv_to_ir)
    source = compute_source_data(ETA_0, max_odd_mode, beta, workers, dps, gl_nodes)
    F_chi = mp.fsum(mode.weight * mode.phi_j for mode in source.modes)
    R_chi = ONE + beta * F_chi
    exterior_memory = source.lambda_out * R_chi
    lambda_geom = LAMBDA_IND + lambda_uv_to_ir + exterior_memory
    jacobi = build_jacobi_audit(source, R_chi)
    return LambdaGeomResult(
        p_num=p_num,
        p_den=p_den,
        p_ir=p_ir,
        lambda_uv_to_ir=lambda_uv_to_ir,
        beta=beta,
        source=source,
        F_chi=F_chi,
        R_chi=R_chi,
        exterior_memory=exterior_memory,
        lambda_geom=lambda_geom,
        jacobi=jacobi,
    )

# =============================================================================
# Scalar channel, static completion, and scenario builders
# =============================================================================

@dataclass(frozen=True)
class BaseShellGeometry:
    c11: mp.mpf
    theta1_core: mp.mpf
    theta1_coll: mp.mpf
    theta1_log: mp.mpf
    a_shell: mp.mpf
    s_uv: mp.mpf
    s_ir: mp.mpf
    eta0: mp.mpf
    ell0: mp.mpf
    gamma_clog_0: mp.mpf
    b11_glue: mp.mpf
    b11_glue_reconstructed: mp.mpf
    b11_glue_delta: mp.mpf

@dataclass(frozen=True)
class StaticCompletionData:
    n_c: mp.mpf
    n_l: mp.mpf
    rho_stat: mp.mpf
    lambda_h: mp.mpf
    c0: mp.mpf
    a_perp: mp.mpf
    b_perp: mp.mpf
    b11_3rho: mp.mpf
    theta1_3rho: mp.mpf

@dataclass(frozen=True)
class ScalarScenario:
    label: str
    family: str
    rank: int
    theta_eff: mp.mpf
    a_uv: mp.mpf
    a_ir: mp.mpf
    chi: mp.mpf
    rho_dyn: mp.mpf
    lambda_u_diag: mp.mpf
    kernel_scale: mp.mpf
    refined_kernel_power: int
    delta_a_mix: mp.mpf
    delta_a0: mp.mpf
    delta_a_diag: mp.mpf
    b11_effective: mp.mpf
    s_uv: mp.mpf
    s_ir: mp.mpf

@dataclass(frozen=True)
class PhysicalPointData:
    dc: mp.mpf
    zeta: mp.mpf
    lambda_pi: mp.mpf
    d_total: mp.mpf
    g_e_over_2: mp.mpf
    a_e: mp.mpf

@dataclass(frozen=True)
class ExactSensitivityData:
    d_lock_d_lambda: mp.mpf
    d_lock_d_k: mp.mpf
    d_alpha_d_d: mp.mpf
    d_alpha_d_lambda: mp.mpf
    d_alpha_d_k: mp.mpf
    log_alpha_log_lambda: mp.mpf
    log_alpha_log_k: mp.mpf
    d_alpha_d_theta: mp.mpf
    log_alpha_log_theta: mp.mpf
    log_alpha_log_b11: mp.mpf
    r_star: mp.mpf
    phi_dyn_star: mp.mpf
    phi_dyn_prime_star: mp.mpf
    radicand_prime_star: mp.mpf

@lru_cache(maxsize=None)
def f_p_cached(p: int, ell: mp.mpf) -> mp.mpf:
    p_mpf = mp.mpf(p)
    phase = 1j * p_mpf * PI * ell / TWO
    return (SQRT_PI * ell / TWO) * mp.e ** (-(p_mpf * PI * ell) ** 2 / FOUR) * mp.re(
        mp.e ** (1j * p_mpf * PI) * (mp.erf(ONE / ell + phase) - mp.erf(phase))
    )

def c11_closed_form(ell: mp.mpf) -> mp.mpf:
    term1 = (SQRT_PI * ell / TWO) * mp.erf(ONE / ell)
    term2 = (SQRT_PI * ell / TWO) * mp.e ** (-(PI**2) * (ell**2))
    term2 *= mp.re(mp.erf(ONE / ell - 1j * PI * ell) - mp.erf(-1j * PI * ell))
    return term1 + term2

def c_mn(m: int, n: int, ell: mp.mpf) -> mp.mpf:
    return m * n * (f_p_cached(abs(m - n), ell) + f_p_cached(m + n, ell))

def gamma_clog(d_value: mp.mpf, ell: mp.mpf) -> mp.mpf:
    root = mp.sqrt(d_value + QUARTER)
    return (SQRT_PI * ell / TWO) * mp.e ** ((ell**2 / FOUR) * (d_value + QUARTER)) * mp.erfc((ell / TWO) * root)

def reconstruct_b11_glue(theta1_coll: mp.mpf, theta1_log: mp.mpf, gamma0: mp.mpf) -> mp.mpf:
    k_vis = mp.matrix(
        [
            [ONE / theta1_coll, gamma0 / TWO],
            [gamma0 / TWO, ONE / theta1_log],
        ]
    )
    t_vis = mp.matrix([[ONE], [ONE]])
    return (t_vis.T * (k_vis ** -1) * t_vis)[0]

def build_base_shell_geometry() -> BaseShellGeometry:
    c11 = c11_closed_form(ELL_0)
    theta1_core = TWO * PI * C_UNI
    theta1_coll = (GAMMA_E / PI**2) * c11
    theta1_log = ETA_0**2 / (8 * (ONE - ETA_0**2 / FOUR))
    a_shell = PI**2 * c11 + GAMMA_E / FOUR
    gamma0 = gamma_clog(ZERO, ELL_0)
    b11_glue = (theta1_coll + theta1_log - gamma0 * theta1_coll * theta1_log) / (
        ONE - (gamma0**2 / FOUR) * theta1_coll * theta1_log
    )
    b11_glue_reconstructed = reconstruct_b11_glue(theta1_coll, theta1_log, gamma0)
    return BaseShellGeometry(
        c11=c11,
        theta1_core=theta1_core,
        theta1_coll=theta1_coll,
        theta1_log=theta1_log,
        a_shell=a_shell,
        s_uv=S_UV,
        s_ir=S_IR,
        eta0=ETA_0,
        ell0=ELL_0,
        gamma_clog_0=gamma0,
        b11_glue=b11_glue,
        b11_glue_reconstructed=b11_glue_reconstructed,
        b11_glue_delta=b11_glue_reconstructed - b11_glue,
    )

def build_static_completion(base: BaseShellGeometry) -> StaticCompletionData:
    n_c = mp.sqrt(PI / 8) * base.ell0
    n_l = ONE
    rho_stat = base.gamma_clog_0 / mp.sqrt(n_c * n_l)
    lambda_h = ONE / (ONE - rho_stat**2)

    a = base.theta1_coll
    b = base.theta1_log
    c0 = base.gamma_clog_0 * a * b / TWO
    a_perp = a - (c0**2) / b
    b_perp = b - (c0**2) / a

    t_perp = mp.matrix(
        [
            [mp.sqrt(a_perp / a)],
            [mp.sqrt(b_perp / b)],
            [ZERO],
        ]
    )
    k_stat = mp.matrix(
        [
            [ONE / a_perp, base.gamma_clog_0 / TWO, mp.sqrt(c0)],
            [base.gamma_clog_0 / TWO, ONE / b_perp, mp.sqrt(c0)],
            [mp.sqrt(c0), mp.sqrt(c0), lambda_h],
        ]
    )
    b11_3rho = (t_perp.T * (k_stat ** -1) * t_perp)[0]
    theta1_3rho = base.theta1_core + b11_3rho

    return StaticCompletionData(
        n_c=n_c,
        n_l=n_l,
        rho_stat=rho_stat,
        lambda_h=lambda_h,
        c0=c0,
        a_perp=a_perp,
        b_perp=b_perp,
        b11_3rho=b11_3rho,
        theta1_3rho=theta1_3rho,
    )

def higher_mode_sums(rank: int, ell: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    delta_a_mix = ZERO
    delta_a0 = ZERO
    if rank < 2:
        return delta_a_mix, delta_a0
    for n in range(2, rank + 1):
        c1n = c_mn(1, n, ell)
        cnn = c_mn(n, n, ell)
        gap = mp.mpf(n**2 - 1)
        delta_a_mix += (GAMMA_E / PI**2) * (c1n**2 / gap)
        delta_a0 += -(GAMMA_E / PI**2) * (c1n**2 * cnn / gap**3)
    return delta_a_mix, delta_a0

def update_dynamic_diagnostics(scenario: ScalarScenario) -> ScalarScenario:
    rho_dyn = scenario.chi / mp.sqrt(scenario.s_uv * scenario.s_ir) if scenario.chi != ZERO else ZERO
    lambda_u_diag = ONE / (ONE - rho_dyn**2) if rho_dyn != ZERO else ONE
    kernel_scale = ONE if scenario.family == "current" else lambda_u_diag ** scenario.refined_kernel_power
    return replace(
        scenario,
        rho_dyn=rho_dyn,
        lambda_u_diag=lambda_u_diag,
        kernel_scale=kernel_scale,
    )

def build_current_scalar_scenario(rank: int, base: BaseShellGeometry) -> ScalarScenario:
    delta_a_mix, delta_a0 = higher_mode_sums(rank, base.ell0)
    delta_a_diag = -(base.eta0**2) * delta_a0
    chi = delta_a_mix**2
    a_uv = base.a_shell + delta_a_mix + delta_a0 + delta_a_diag
    a_ir = base.a_shell - delta_a_mix + delta_a0 + delta_a_diag
    raw = ScalarScenario(
        label=f"current_rank_{rank}",
        family="current",
        rank=rank,
        theta_eff=base.theta1_core + base.b11_glue,
        a_uv=a_uv,
        a_ir=a_ir,
        chi=chi,
        rho_dyn=ZERO,
        lambda_u_diag=ONE,
        kernel_scale=ONE,
        refined_kernel_power=0,
        delta_a_mix=delta_a_mix,
        delta_a0=delta_a0,
        delta_a_diag=delta_a_diag,
        b11_effective=base.b11_glue,
        s_uv=base.s_uv,
        s_ir=base.s_ir,
    )
    return update_dynamic_diagnostics(raw)

def build_refined_scalar_scenario(
    rank: int,
    base: BaseShellGeometry,
    static_data: StaticCompletionData,
    kernel_power: int,
) -> ScalarScenario:
    delta_a_mix, delta_a0 = higher_mode_sums(rank, base.ell0)
    delta_a_diag = -(base.eta0**2) * delta_a0
    chi = delta_a_mix**2
    a_uv = base.a_shell + delta_a_mix + delta_a0 + delta_a_diag
    a_ir = base.a_shell - delta_a_mix + delta_a0 + delta_a_diag
    raw = ScalarScenario(
        label=f"refined_rank_{rank}_kernelpow_{kernel_power}",
        family="refined",
        rank=rank,
        theta_eff=static_data.theta1_3rho,
        a_uv=a_uv,
        a_ir=a_ir,
        chi=chi,
        rho_dyn=ZERO,
        lambda_u_diag=ONE,
        kernel_scale=ONE,
        refined_kernel_power=kernel_power,
        delta_a_mix=delta_a_mix,
        delta_a0=delta_a0,
        delta_a_diag=delta_a_diag,
        b11_effective=static_data.b11_3rho,
        s_uv=base.s_uv,
        s_ir=base.s_ir,
    )
    return update_dynamic_diagnostics(raw)

# =============================================================================
# Scalar evaluator and primary α inversion
# =============================================================================

def scalar_dynamic_numerator(d_value: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    return (
        scenario.a_uv * (ONE + scenario.s_ir * d_value)
        + scenario.a_ir * (ONE + scenario.s_uv * d_value)
        - TWO * scenario.chi * d_value * mp.sqrt(scenario.a_uv * scenario.a_ir)
    )

def scalar_dynamic_denominator(d_value: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    return (ONE + scenario.s_uv * d_value) * (ONE + scenario.s_ir * d_value) - (scenario.chi**2) * (d_value**2)

def scalar_dynamic_phi(d_value: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    return scenario.kernel_scale * scalar_dynamic_numerator(d_value, scenario) / scalar_dynamic_denominator(d_value, scenario)

def scalar_radicand(d_value: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    return ONE - scenario.theta_eff * d_value + d_value**2 * scalar_dynamic_phi(d_value, scenario)

def solve_scalar_dc_for_alpha(alpha: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    x_value = alpha / PI

    def residual(d_value: mp.mpf) -> mp.mpf:
        rad = scalar_radicand(d_value, scenario)
        if rad < 0:
            raise ValueError("Negative scalar radicand during root solve.")
        return d_value - x_value * mp.sqrt(rad)

    left = ZERO
    right = max(mp.mpf("0.003"), TWO * x_value)
    f_left = residual(left)
    f_right = residual(right)
    if not (f_left < 0 < f_right):
        raise RuntimeError("Failed to bracket scalar D_C root.")
    for _ in range(ROOT_MAX_ITERS):
        mid = (left + right) / TWO
        f_mid = residual(mid)
        if abs(f_mid) < ROOT_TOL or abs(right - left) < ROOT_TOL:
            return mid
        if f_mid > 0:
            right = mid
        else:
            left = mid
    return (left + right) / TWO

def dc_of_alpha(alpha: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    return solve_scalar_dc_for_alpha(alpha, scenario)

def zeta_from_lambda(lambda_value: mp.mpf) -> mp.mpf:
    return K_GEOMETRIC * lambda_value / (TWO * PI**2)

def dc_lock_from_lambda(lambda_value: mp.mpf) -> mp.mpf:
    return (THREE / TWO) * K_GEOMETRIC * lambda_value * (ONE + (K_GEOMETRIC * lambda_value) / (TWO * PI**2))

def lambda_from_dc_lock(dc_value: mp.mpf) -> mp.mpf:
    zeta = (mp.sqrt(ONE + (FOUR * dc_value) / (THREE * PI**2)) - ONE) / TWO
    return (TWO * PI**2 * zeta) / K_GEOMETRIC

def zeta_from_dc(dc_value: mp.mpf) -> mp.mpf:
    return HALF * (mp.sqrt(ONE + (FOUR * dc_value) / (THREE * PI**2)) - ONE)

def lambda_pi_from_dc(dc_value: mp.mpf) -> mp.mpf:
    return (TWO * PI**2 / K_GEOMETRIC) * zeta_from_dc(dc_value)

def d_total_from_dc_and_alpha(dc_value: mp.mpf, alpha: mp.mpf) -> mp.mpf:
    x = alpha / PI
    zeta = zeta_from_dc(dc_value)
    return dc_value - x * zeta + (x**2 * zeta**2) / (FOUR * dc_value)

def g_e_over_2_from_d_total(d_total: mp.mpf) -> mp.mpf:
    return ONE / mp.sqrt(ONE - d_total)

def c_log_from_dc_and_lambda(dc_value: mp.mpf, lambda_value: mp.mpf) -> mp.mpf:
    zeta = zeta_from_lambda(lambda_value)
    return (C_CHI / (TWO * dc_value)) * zeta * (ONE + zeta)

def solve_alpha_from_dc_model(dc_model: Callable[[mp.mpf], mp.mpf], dc_target: mp.mpf) -> mp.mpf:
    def residual(alpha_value: mp.mpf) -> mp.mpf:
        return dc_model(alpha_value) - dc_target

    left = ALPHA_BRACKET_LEFT
    right = ALPHA_BRACKET_RIGHT
    f_left = residual(left)
    f_right = residual(right)
    if not (f_left < 0 < f_right):
        raise RuntimeError("Failed to bracket alpha root.")
    for _ in range(ROOT_MAX_ITERS):
        mid = (left + right) / TWO
        f_mid = residual(mid)
        if abs(f_mid) < ROOT_TOL or abs(right - left) < ROOT_TOL:
            return mid
        if f_mid > 0:
            right = mid
        else:
            left = mid
    return (left + right) / TWO

# =============================================================================
# Power-series tools
# =============================================================================

def zero_series(order: int) -> List[mp.mpf]:
    return [ZERO for _ in range(order + 1)]

def one_series(order: int) -> List[mp.mpf]:
    out = zero_series(order)
    out[0] = ONE
    return out

def x_series(order: int) -> List[mp.mpf]:
    out = zero_series(order)
    if order >= 1:
        out[1] = ONE
    return out

def series_add(a: List[mp.mpf], b: List[mp.mpf]) -> List[mp.mpf]:
    return [ai + bi for ai, bi in zip(a, b)]

def series_sub(a: List[mp.mpf], b: List[mp.mpf]) -> List[mp.mpf]:
    return [ai - bi for ai, bi in zip(a, b)]

def series_scale(a: List[mp.mpf], c: mp.mpf) -> List[mp.mpf]:
    return [c * ai for ai in a]

def series_mul(a: List[mp.mpf], b: List[mp.mpf], order: int) -> List[mp.mpf]:
    out = zero_series(order)
    for n in range(order + 1):
        out[n] = mp.fsum(a[k] * b[n - k] for k in range(n + 1))
    return out

def series_inv(a: List[mp.mpf], order: int) -> List[mp.mpf]:
    if a[0] == 0:
        raise ZeroDivisionError("Series inversion requires nonzero constant term.")
    out = zero_series(order)
    out[0] = ONE / a[0]
    for n in range(1, order + 1):
        s = mp.fsum(a[k] * out[n - k] for k in range(1, n + 1))
        out[n] = -s / a[0]
    return out

def series_sqrt(a: List[mp.mpf], order: int) -> List[mp.mpf]:
    if a[0] <= 0:
        raise ValueError("Square-root series requires positive constant term.")
    out = zero_series(order)
    out[0] = mp.sqrt(a[0])
    for n in range(1, order + 1):
        s = mp.fsum(out[k] * out[n - k] for k in range(1, n))
        out[n] = (a[n] - s) / (TWO * out[0])
    return out

def series_eval(coeffs: Dict[int, mp.mpf], x: mp.mpf) -> mp.mpf:
    max_order = max(coeffs) if coeffs else 0
    return mp.fsum(coeffs.get(n, ZERO) * x**n for n in range(max_order + 1))

def coeff_dict_from_series(series: List[mp.mpf], start: int = 1, stop: int | None = None) -> Dict[int, mp.mpf]:
    if stop is None:
        stop = len(series) - 1
    return {n: series[n] for n in range(start, stop + 1)}

# =============================================================================
# Derived scalar/g-2 series from the current/refined evaluators
# =============================================================================

def build_dc_series_from_scenario(scenario: ScalarScenario, order: int) -> Dict[int, mp.mpf]:
    xs = x_series(order)
    dc = x_series(order)

    sqrt_aauv = mp.sqrt(scenario.a_uv * scenario.a_ir)
    n0 = scenario.a_uv + scenario.a_ir
    n1 = scenario.a_uv * scenario.s_ir + scenario.a_ir * scenario.s_uv - TWO * scenario.chi * sqrt_aauv
    q1 = scenario.s_uv + scenario.s_ir
    q2 = scenario.s_uv * scenario.s_ir - scenario.chi**2

    for _ in range(80):
        d2 = series_mul(dc, dc, order)
        n_series = series_add(series_scale(one_series(order), n0), series_scale(dc, n1))
        q_series = series_add(one_series(order), series_add(series_scale(dc, q1), series_scale(d2, q2)))
        phi_series = series_scale(series_mul(n_series, series_inv(q_series, order), order), scenario.kernel_scale)
        rad_series = series_add(
            one_series(order),
            series_sub(series_scale(dc, -scenario.theta_eff), series_scale(series_mul(d2, phi_series, order), -ONE)),
        )
        new_dc = series_mul(xs, series_sqrt(rad_series, order), order)
        if all(abs(new_dc[n] - dc[n]) < mp.mpf("1e-95") for n in range(order + 1)):
            dc = new_dc
            break
        dc = new_dc
    return coeff_dict_from_series(dc, start=1, stop=order)

def zeta_series_from_dc(dc_series: List[mp.mpf], order: int) -> List[mp.mpf]:
    inside = series_add(one_series(order), series_scale(dc_series, FOUR / (THREE * PI**2)))
    return series_scale(series_sub(series_sqrt(inside, order), one_series(order)), HALF)

def d_total_series_from_dc(dc_coeffs: Dict[int, mp.mpf], order: int) -> Dict[int, mp.mpf]:
    dc_series = zero_series(order)
    for n, value in dc_coeffs.items():
        if n <= order:
            dc_series[n] = value
    zeta = zeta_series_from_dc(dc_series, order)
    xs = x_series(order)
    u = [dc_series[n + 1] if n + 1 <= order else ZERO for n in range(order + 1)]  # D/x
    u_inv = series_inv(u, order)
    zeta_sq = series_mul(zeta, zeta, order)
    term3 = series_scale(series_mul(xs, series_mul(zeta_sq, u_inv, order), order), QUARTER)
    d_total = series_add(dc_series, series_sub(series_scale(series_mul(xs, zeta, order), -ONE), series_scale(term3, -ONE)))
    return coeff_dict_from_series(d_total, start=1, stop=order)

def g2_coefficients_from_dc(dc_coeffs: Dict[int, mp.mpf], order: int) -> Dict[int, mp.mpf]:
    d_total_coeffs = d_total_series_from_dc(dc_coeffs, order)
    d_total_series = zero_series(order)
    for n, value in d_total_coeffs.items():
        if n <= order:
            d_total_series[n] = value
    g_series = series_inv(series_sqrt(series_sub(one_series(order), d_total_series), order), order)
    return coeff_dict_from_series(g_series, start=1, stop=order)

def evaluate_physical_point(alpha: mp.mpf, scenario: ScalarScenario) -> PhysicalPointData:
    dc = dc_of_alpha(alpha, scenario)
    zeta = zeta_from_dc(dc)
    lambda_pi = lambda_pi_from_dc(dc)
    d_total = d_total_from_dc_and_alpha(dc, alpha)
    g2 = g_e_over_2_from_d_total(d_total)
    return PhysicalPointData(dc=dc, zeta=zeta, lambda_pi=lambda_pi, d_total=d_total, g_e_over_2=g2, a_e=g2 - ONE)

# =============================================================================
# QED-induced scalar series and QED anomaly helpers (reference-only)
# =============================================================================

def qed_terms(alpha: mp.mpf, coeffs: Dict[int, mp.mpf]) -> list[tuple[int, mp.mpf, mp.mpf]]:
    x = alpha / PI
    return [(n, coeffs[n], coeffs[n] * x**n) for n in sorted(coeffs)]

def qed_anomaly(alpha: mp.mpf, coeffs: Dict[int, mp.mpf]) -> mp.mpf:
    return mp.fsum(term for _, _, term in qed_terms(alpha, coeffs))

def build_d_total_series_from_qed_ae(ae_coeffs: Dict[int, mp.mpf], order: int) -> List[mp.mpf]:
    ae = zero_series(order)
    for n, value in ae_coeffs.items():
        if n <= order:
            ae[n] = value
    inv_one_plus_ae = series_inv(series_add(one_series(order), ae), order)
    inv_sq = series_mul(inv_one_plus_ae, inv_one_plus_ae, order)
    return series_sub(one_series(order), inv_sq)

def build_dc_series_from_qed_ae(ae_coeffs: Dict[int, mp.mpf], order: int) -> Dict[int, mp.mpf]:
    d_total_target = build_d_total_series_from_qed_ae(ae_coeffs, order)
    dc = x_series(order)
    xs = x_series(order)
    for _ in range(80):
        zeta = zeta_series_from_dc(dc, order)
        term2 = series_mul(xs, zeta, order)
        zeta_sq = series_mul(zeta, zeta, order)
        u = [dc[n + 1] if n + 1 <= order else ZERO for n in range(order + 1)]
        u_inv = series_inv(u, order)
        term3 = series_scale(series_mul(xs, series_mul(zeta_sq, u_inv, order), order), QUARTER)
        dc_new = series_add(d_total_target, series_sub(term2, term3))
        if all(abs(dc_new[n] - dc[n]) < mp.mpf("1e-95") for n in range(order + 1)):
            dc = dc_new
            break
        dc = dc_new
    return coeff_dict_from_series(dc, start=1, stop=order)

def dc_from_series(alpha: mp.mpf, series: Dict[int, mp.mpf]) -> mp.mpf:
    return series_eval(series, alpha / PI)

# =============================================================================
# Sensitivities and stability diagnostics
# =============================================================================

def scalar_dynamic_primes(scenario: ScalarScenario, d_value: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    n = scalar_dynamic_numerator(d_value, scenario)
    q = scalar_dynamic_denominator(d_value, scenario)
    n_prime = scenario.a_uv * scenario.s_ir + scenario.a_ir * scenario.s_uv - TWO * scenario.chi * mp.sqrt(scenario.a_uv * scenario.a_ir)
    q_prime = scenario.s_uv + scenario.s_ir + TWO * d_value * (scenario.s_uv * scenario.s_ir - scenario.chi**2)
    phi_prime = scenario.kernel_scale * (n_prime * q - n * q_prime) / (q**2)
    return n, q, n_prime, phi_prime

def compute_exact_sensitivities(lambda_geom: mp.mpf, dc_lock: mp.mpf, alpha_pred: mp.mpf, scenario: ScalarScenario) -> ExactSensitivityData:
    r_star = scalar_radicand(dc_lock, scenario)
    phi_dyn_star = scalar_dynamic_phi(dc_lock, scenario)
    _, _, _, phi_dyn_prime_star = scalar_dynamic_primes(scenario, dc_lock)
    radicand_prime_star = -scenario.theta_eff + TWO * dc_lock * phi_dyn_star + dc_lock**2 * phi_dyn_prime_star

    d_lock_d_lambda = (THREE / TWO) * K_GEOMETRIC * (ONE + (K_GEOMETRIC * lambda_geom) / (PI**2))
    d_lock_d_k = (THREE / TWO) * lambda_geom * (ONE + (K_GEOMETRIC * lambda_geom) / (PI**2))
    d_alpha_d_d = PI * (ONE - (dc_lock / TWO) * (radicand_prime_star / r_star)) / mp.sqrt(r_star)
    d_alpha_d_lambda = d_alpha_d_d * d_lock_d_lambda
    d_alpha_d_k = d_alpha_d_d * d_lock_d_k

    d_alpha_d_theta = alpha_pred * dc_lock / (TWO * r_star)
    log_alpha_log_theta = scenario.theta_eff * dc_lock / (TWO * r_star)
    log_alpha_log_b11 = scenario.b11_effective * dc_lock / (TWO * r_star)

    return ExactSensitivityData(
        d_lock_d_lambda=d_lock_d_lambda,
        d_lock_d_k=d_lock_d_k,
        d_alpha_d_d=d_alpha_d_d,
        d_alpha_d_lambda=d_alpha_d_lambda,
        d_alpha_d_k=d_alpha_d_k,
        log_alpha_log_lambda=lambda_geom * d_alpha_d_lambda / alpha_pred,
        log_alpha_log_k=K_GEOMETRIC * d_alpha_d_k / alpha_pred,
        d_alpha_d_theta=d_alpha_d_theta,
        log_alpha_log_theta=log_alpha_log_theta,
        log_alpha_log_b11=log_alpha_log_b11,
        r_star=r_star,
        phi_dyn_star=phi_dyn_star,
        phi_dyn_prime_star=phi_dyn_prime_star,
        radicand_prime_star=radicand_prime_star,
    )

def solve_alpha_for_scenario(dc_lock: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    return solve_alpha_from_dc_model(lambda a: dc_of_alpha(a, scenario), dc_lock)

def symmetric_log_sensitivity(
    dc_lock: mp.mpf,
    scenario: ScalarScenario,
    getter: Callable[[ScalarScenario], mp.mpf],
    setter: Callable[[ScalarScenario, mp.mpf], ScalarScenario],
    eps: mp.mpf = SENSITIVITY_EPS,
) -> mp.mpf:
    q0 = getter(scenario)
    alpha_plus = solve_alpha_for_scenario(dc_lock, setter(scenario, q0 * (ONE + eps)))
    alpha_minus = solve_alpha_for_scenario(dc_lock, setter(scenario, q0 * (ONE - eps)))
    return mp.log(alpha_plus / alpha_minus) / mp.log((ONE + eps) / (ONE - eps))

def numerical_sensitivities(dc_lock: mp.mpf, scenario: ScalarScenario) -> Dict[str, mp.mpf]:
    return {
        "Theta_eff": symmetric_log_sensitivity(
            dc_lock,
            scenario,
            getter=lambda s: s.theta_eff,
            setter=lambda s, v: replace(s, theta_eff=v),
        ),
        "A_uv": symmetric_log_sensitivity(
            dc_lock,
            scenario,
            getter=lambda s: s.a_uv,
            setter=lambda s, v: replace(s, a_uv=v),
        ),
        "A_ir": symmetric_log_sensitivity(
            dc_lock,
            scenario,
            getter=lambda s: s.a_ir,
            setter=lambda s, v: replace(s, a_ir=v),
        ),
        "chi": symmetric_log_sensitivity(
            dc_lock,
            scenario,
            getter=lambda s: s.chi,
            setter=lambda s, v: update_dynamic_diagnostics(replace(s, chi=v)),
        ),
        "s_uv": symmetric_log_sensitivity(
            dc_lock,
            scenario,
            getter=lambda s: s.s_uv,
            setter=lambda s, v: update_dynamic_diagnostics(replace(s, s_uv=v)),
        ),
        "s_ir": symmetric_log_sensitivity(
            dc_lock,
            scenario,
            getter=lambda s: s.s_ir,
            setter=lambda s, v: update_dynamic_diagnostics(replace(s, s_ir=v)),
        ),
    }

# =============================================================================
# Report assembly
# =============================================================================

def contribution_share(part: mp.mpf, total: mp.mpf) -> str:
    return fmt_pct(100 * part / total)

def coefficient_error_rows(
    current_rank5: Dict[int, mp.mpf],
    refined_rank100: Dict[int, mp.mpf],
) -> list[tuple[str, str, str, str, str, str, str, str]]:
    rows = []
    for n in range(1, QED_SERIES_ORDER + 1):
        q = QED_A1_UNIVERSAL[n]
        cur = current_rank5[n]
        ref = refined_rank100[n]
        d_cur = cur - q
        d_ref = ref - q
        rel_cur = d_cur / q if q != ZERO else ZERO
        rel_ref = d_ref / q if q != ZERO else ZERO
        rows.append(
            (
                str(n),
                fmt(q, 24),
                fmt(cur, 24),
                fmt(ref, 24),
                fmt_sci(d_cur, 14),
                fmt_sci(d_ref, 14),
                fmt_sci(rel_cur, 14),
                fmt_sci(rel_ref, 14),
            )
        )
    return rows

def build_report_text(
    max_odd_mode: int,
    workers: int,
    dps: int,
    gl_nodes: int,
    print_digits: int,
    output_txt: str,
) -> str:
    global PRINT_DIGITS
    PRINT_DIGITS = print_digits

    out = io.StringIO()
    def write(s: str = "") -> None:
        out.write(s + "\n")

    # -------------------------------------------------------------------------
    # Core computations
    # -------------------------------------------------------------------------
    lambda_data = evaluate_lambda_geom(max_odd_mode, workers, dps, gl_nodes)
    base = build_base_shell_geometry()
    static_data = build_static_completion(base)

    current_scenarios = {n: build_current_scalar_scenario(n, base) for n in RANK_AUDIT_LIST}
    refined_realized_scenarios = {
        n: build_refined_scalar_scenario(n, base, static_data, REFINED_REALIZED_KERNEL_POWER)
        for n in RANK_AUDIT_LIST
    }
    refined_scalar_coeff_scenario = build_refined_scalar_scenario(
        100, base, static_data, REFINED_SCALAR_TABLE_KERNEL_POWER
    )

    article_scenario = current_scenarios[ARTICLE_RANK]
    dc_lock = dc_lock_from_lambda(lambda_data.lambda_geom)
    alpha_article = solve_alpha_for_scenario(dc_lock, article_scenario)
    alpha_article_inv = ONE / alpha_article
    dc_at_alpha_article = dc_of_alpha(alpha_article, article_scenario)
    zeta_star = zeta_from_lambda(lambda_data.lambda_geom)
    d_total_article = d_total_from_dc_and_alpha(dc_at_alpha_article, alpha_article)
    lambda_back = lambda_from_dc_lock(dc_at_alpha_article)
    c_log_star = c_log_from_dc_and_lambda(dc_at_alpha_article, lambda_data.lambda_geom)

    # External QED-induced D_C cross-check
    qdc_universal = build_dc_series_from_qed_ae(QED_A1_UNIVERSAL, QED_SERIES_ORDER)
    alpha_qed_cross = solve_alpha_from_dc_model(lambda a: dc_from_series(a, qdc_universal), dc_lock)
    alpha_qed_cross_inv = ONE / alpha_qed_cross
    dc_at_alpha_qed = dc_from_series(alpha_qed_cross, qdc_universal)

    # Scalar-series data
    dc_coeffs_current = {n: build_dc_series_from_scenario(current_scenarios[n], QED_SERIES_ORDER) for n in (1, 5, 100)}
    dc_coeffs_ref_scalar = build_dc_series_from_scenario(refined_scalar_coeff_scenario, QED_SERIES_ORDER)

    g2_coeffs_current_rank5 = g2_coefficients_from_dc(build_dc_series_from_scenario(current_scenarios[5], PREDICTION_ORDER), PREDICTION_ORDER)
    g2_coeffs_refined_rank100 = g2_coefficients_from_dc(
        build_dc_series_from_scenario(refined_realized_scenarios[100], PREDICTION_ORDER),
        PREDICTION_ORDER,
    )

    # Physical-point tables and rank summary
    current_physical = {n: evaluate_physical_point(ALPHA_REF, current_scenarios[n]) for n in RANK_AUDIT_LIST}
    refined_physical = {n: evaluate_physical_point(ALPHA_REF, refined_realized_scenarios[n]) for n in RANK_AUDIT_LIST}
    max_abs_delta_current = {
        n: max(abs(g2_coefficients_from_dc(build_dc_series_from_scenario(current_scenarios[n], QED_SERIES_ORDER), QED_SERIES_ORDER)[k] - QED_A1_UNIVERSAL[k]) for k in range(1, QED_SERIES_ORDER + 1))
        for n in RANK_AUDIT_LIST
    }
    max_abs_delta_refined = {
        n: max(abs(g2_coefficients_from_dc(build_dc_series_from_scenario(refined_realized_scenarios[n], QED_SERIES_ORDER), QED_SERIES_ORDER)[k] - QED_A1_UNIVERSAL[k]) for k in range(1, QED_SERIES_ORDER + 1))
        for n in RANK_AUDIT_LIST
    }

    # Alpha comparisons
    delta_alpha_article = alpha_article - ALPHA_REF
    delta_alpha_qed = alpha_qed_cross - ALPHA_REF
    delta_alpha_article_qed = alpha_article - alpha_qed_cross

    exact_sens = compute_exact_sensitivities(lambda_data.lambda_geom, dc_lock, alpha_article, article_scenario)
    num_sens = numerical_sensitivities(dc_lock, article_scenario)

    alpha_rank_current = {n: solve_alpha_for_scenario(dc_lock, current_scenarios[n]) for n in RANK_AUDIT_LIST}
    alpha_rank_current_ref = alpha_rank_current[100]

    # -------------------------------------------------------------------------
    # Header
    # -------------------------------------------------------------------------
    write("Λ_geom -> D_C^lock -> α_pred -> pure-photonic g-2 audit")
    write("Primary layer: Λ_geom, D_C^lock, α_pred, and derived Relator g-2 coefficients are computed without benchmark α and without QED coefficients.")
    write("Reference-only layer: benchmark α and universal QED coefficients appear only in explicitly marked comparison tables.")
    write(f"mpmath precision: {dps} decimal digits")
    write(f"odd-mode cutoff: j ≤ {max_odd_mode}")
    write(f"workers used for shell-source modes: {workers}")
    write(f"Gauss--Legendre nodes per shell integral: {gl_nodes}")
    write(f"output text file: {output_txt}")
    write("refined-kernel conventions kept in the file: λ_u for scalar coefficient reproduction; λ_u² for realized refined g-2 tables")

    # -------------------------------------------------------------------------
    # Vector geometry tables
    # -------------------------------------------------------------------------
    source = lambda_data.source
    jac = lambda_data.jacobi

    write(render_table(
        "Table 1. Branch geometry and fixed constants",
        ["Symbol", "Formula", "Value"],
        [
            ("ε₀", "1/√π", fmt(EPSILON_0)),
            ("η₀", "1/π", fmt(ETA_0)),
            ("ℓ₀", "1/(π√π)", fmt(ELL_0)),
            ("γ_E", "Euler–Mascheroni", fmt(GAMMA_E)),
            ("Λ_ind", "ln(8√π) − 2", fmt(LAMBDA_IND)),
            ("C₀^(Gauss)", "(ln2 + γ_E)/2", fmt(C0_GAUSS)),
            ("γ_geom", "½ sinh(η₀)/η₀", fmt(GAMMA_GEOM)),
            ("A", "(ln2)η₀⁴", fmt(A_UV)),
            ("B", "η₀⁴/(8π²)", fmt(B_IR)),
            ("N_ω", "√(2/π)", fmt(N_OMEGA)),
            ("K", "(150π² − 8π⁴ − 315)/(180π⁶)", fmt(K_GEOMETRIC)),
        ],
    ))

    write(render_table(
        "Table 2. TT–χ admission and exact shell source",
        ["Symbol", "Formula", "Value"],
        [
            ("D_IR", "∫₀¹ x² sin²(πx) dx", fmt(lambda_data.p_den)),
            ("N_IR(ℓ₀)", "∫₀¹ x²sin²(πx)[1 − f_swirl/3]e^{-((1−x)/ℓ₀)²} dx", fmt(lambda_data.p_num)),
            ("P_IR^(χ)(ℓ₀)", "N_IR/D_IR", fmt(lambda_data.p_ir)),
            ("ΔΛ^(UV→IR)", "C₀^(Gauss) P_IR^(χ)", fmt(lambda_data.lambda_uv_to_ir)),
            ("‖J_χ‖²", "Σ_{j odd}(J_j^(χ))²", fmt(source.norm_sq)),
            ("Λ_OUT", "−πΣ â_j²/[(j+1)(2j+1)]", fmt(source.lambda_out)),
            ("−½‖J_χ‖²", "quadratic-source identity", fmt(source.lambda_out_from_norm)),
            ("Λ_OUT + ½‖J_χ‖²", "internal check", fmt_sci(source.norm_identity_error)),
        ],
    ))

    low_modes = [mode for mode in source.modes if mode.j <= LOW_MODE_PRINT_MAX]
    write(render_table(
        "Table 3. Low-mode shell-source and return data",
        ["j", "a_j^sh", "â_j", "J_j^(χ)", "w_j", "q_j", "Φ_j", "R_j"],
        [
            (
                str(mode.j),
                fmt(mode.a_sh, 24),
                fmt(mode.a_hat, 24),
                fmt(mode.J, 24),
                fmt(mode.weight, 24),
                fmt(mode.q_j, 24),
                fmt(mode.phi_j, 24),
                fmt(mode.R_j, 24),
            )
            for mode in low_modes
        ],
    ))

    write(render_table(
        "Table 4. Source-weighted mother resolvent and Jacobi audit",
        ["Symbol", "Formula", "Value"],
        [
            ("β", "P_IR^(χ)γ_geom + P_IR^(χ)ΔΛ^(UV→IR)/[2(1+A)]", fmt(lambda_data.beta)),
            ("F_χ", "Σ_j w_j Φ_j", fmt(lambda_data.F_chi)),
            ("R_χ", "1 + βF_χ", fmt(lambda_data.R_chi)),
            ("Λ_OUT R_χ", "exterior-memory contribution", fmt(lambda_data.exterior_memory)),
            ("m₁", "Σ_j w_j g_j", fmt(jac.m1)),
            ("m₂", "Σ_j w_j g_j²", fmt(jac.m2)),
            ("m₃", "Σ_j w_j g_j³", fmt(jac.m3)),
            ("b₁²", "m₂ − m₁²", fmt(jac.b1_sq)),
            ("a₁", "(m₃ − 2m₁m₂ + m₁³)/(m₂ − m₁²)", fmt(jac.a1)),
            ("R_χ^[2]", "1/[1 − a₀ − b₁²/(1 − a₁)]", fmt(jac.R_rank2)),
            ("R_χ − R_χ^[2]", "Jacobi residual", fmt_sci(jac.delta_R)),
            ("δΛ_from_Jacobi", "Λ_OUT(R_χ − R_χ^[2])", fmt_sci(jac.delta_lambda)),
        ],
    ))

    write(render_table(
        "Table 5. Top-level contribution ledger for Λ_geom",
        ["Contribution", "Formula", "Value", "Share"],
        [
            ("Λ_ind", "ln(8√π) − 2", fmt(LAMBDA_IND), contribution_share(LAMBDA_IND, lambda_data.lambda_geom)),
            ("ΔΛ^(UV→IR)", "C₀^(Gauss)P_IR^(χ)", fmt(lambda_data.lambda_uv_to_ir), contribution_share(lambda_data.lambda_uv_to_ir, lambda_data.lambda_geom)),
            ("Λ_OUT R_χ", "Λ_OUT · R_χ", fmt(lambda_data.exterior_memory), contribution_share(lambda_data.exterior_memory, lambda_data.lambda_geom)),
            ("Λ_geom", "Λ_ind + ΔΛ^(UV→IR) + Λ_OUT R_χ", fmt(lambda_data.lambda_geom), "+100.000000%"),
        ],
    ))

    write(render_table(
        "Table 6. Internal decomposition of the exterior-memory contribution",
        ["Symbol", "Formula", "Value"],
        [
            ("Λ_OUT", "exterior source energy", fmt(source.lambda_out)),
            ("R_χ", "source-weighted memory factor", fmt(lambda_data.R_chi)),
            ("Λ_OUT R_χ", "exterior-memory contribution", fmt(lambda_data.exterior_memory)),
            ("weight sum", "Σ_j w_j", fmt(source.weight_sum)),
        ],
    ))

    write(render_table(
        "Table 7. ALP lock data extracted from Λ_geom",
        ["Symbol", "Formula", "Value"],
        [
            ("Λ_geom", "vector-shell input", fmt(lambda_data.lambda_geom)),
            ("ζ_*", "KΛ/(2π²)", fmt(zeta_star)),
            ("D_C^lock", "(3/2)KΛ[1 + KΛ/(2π²)]", fmt(dc_lock)),
            ("C_log", "π² ζ(1+ζ)/D_C", fmt(c_log_star)),
            ("Λ recovered", "Λ(D_C^lock)", fmt(lambda_back)),
            ("Λ recovery error", "Λ_recovered − Λ_geom", fmt_sci(lambda_back - lambda_data.lambda_geom)),
        ],
    ))

    # -------------------------------------------------------------------------
    # Scalar mother tables
    # -------------------------------------------------------------------------
    write(render_table(
        "Table 8. Fixed shell geometry and static hidden-mode data of the scalar block",
        ["Symbol", "Definition", "Value"],
        [
            ("η₀", "1/π", fmt(base.eta0)),
            ("ℓ₀", "1/(π√π)", fmt(base.ell0)),
            ("C_uni", "(1/π)(4/3 + 1/(4π²))", fmt(C_UNI)),
            ("s_uv", "ln 2", fmt(base.s_uv)),
            ("s_ir", "1/(8π²)", fmt(base.s_ir)),
            ("γ_cℓ(0)", "static collar–log overlap", fmt(base.gamma_clog_0)),
            ("N_c", "√(π/8) ℓ₀", fmt(static_data.n_c)),
            ("ρ_stat", "γ_cℓ(0)/√(N_c N_ℓ)", fmt(static_data.rho_stat)),
            ("Λ_h", "(1 − ρ_stat²)^(-1)", fmt(static_data.lambda_h)),
            ("Θ_{1,core}", "2π C_uni", fmt(base.theta1_core)),
            ("Θ_1^(coll)", "(γ_E/π²) C_11(ℓ₀)", fmt(base.theta1_coll)),
            ("Θ_1^(log)", "η₀²/[8(1 − η₀²/4)]", fmt(base.theta1_log)),
            ("B_11^glue", "visible 2-channel static load", fmt(base.b11_glue)),
            ("B_11^glue(recon)", "t_visᵀ(K_vis^(2))^(-1)t_vis", fmt(base.b11_glue_reconstructed)),
            ("ΔB_11^glue", "reconstructed − stored", fmt_sci(base.b11_glue_delta)),
            ("B_11^(3,ρ)", "canonical 3-channel static completion", fmt(static_data.b11_3rho)),
            ("Θ_1^(3,ρ)", "Θ_{1,core} + B_11^(3,ρ)", fmt(static_data.theta1_3rho)),
        ],
    ))

    write(render_table(
        "Table 9. Rank-dependent visible dynamic diagnostics (current evaluator)",
        ["N", "ΔA_mix^(N)", "δA_0^(N)", "χ_cross^(N)"],
        [
            (
                str(n),
                fmt(current_scenarios[n].delta_a_mix, 24),
                fmt(current_scenarios[n].delta_a0, 24),
                fmt(current_scenarios[n].chi, 24),
            )
            for n in RANK_AUDIT_LIST
        ],
    ))

    write(render_table(
        "Table 10. Visible dynamic two-channel parameters",
        ["N", "A_uv^(N)", "A_ir^(N)", "ρ_dyn^(N)", "λ_u^(N)=(1−ρ_dyn²)^(-1)"],
        [
            (
                str(n),
                fmt(current_scenarios[n].a_uv, 24),
                fmt(current_scenarios[n].a_ir, 24),
                fmt(current_scenarios[n].rho_dyn, 24),
                fmt(current_scenarios[n].lambda_u_diag, 24),
            )
            for n in RANK_AUDIT_LIST
        ],
    ))

    write(render_table(
        "Table 11. Scalar roots at the external orientation point x_ref = α_ref/π",
        ["N", "D_C^(current,N)(x_ref)", "D_C^(refined,N)(x_ref)", "ΔD_C^(N)=refined-current"],
        [
            (
                str(n),
                fmt(current_physical[n].dc, 24),
                fmt(refined_physical[n].dc, 24),
                fmt_sci(refined_physical[n].dc - current_physical[n].dc, 14),
            )
            for n in RANK_AUDIT_LIST
        ],
    ))

    write(render_table(
        "Table 12. First coefficients of the scalar branch D_C(x)=Σ c_n x^n",
        ["Coefficient", "N=1 current", "N=5 current", "N=100 current", "N=100 refined (λ_u)", "QED-induced"],
        [
            (
                f"c_{n}",
                fmt(dc_coeffs_current[1][n], 24),
                fmt(dc_coeffs_current[5][n], 24),
                fmt(dc_coeffs_current[100][n], 24),
                fmt(dc_coeffs_ref_scalar[n], 24),
                fmt(qdc_universal[n], 24),
            )
            for n in range(1, QED_SERIES_ORDER + 1)
        ],
    ))

    write(render_table(
        "Table 13. Article lock components at D_* = D_C^lock",
        ["Quantity", "Formula", "Value"],
        [
            ("Θ_eff", "effective scalar stiffness", fmt(article_scenario.theta_eff)),
            ("A_uv", "visible UV dynamic amplitude", fmt(article_scenario.a_uv)),
            ("A_ir", "visible IR dynamic amplitude", fmt(article_scenario.a_ir)),
            ("χ_N", "(ΔA_mix)^2", fmt(article_scenario.chi)),
            ("ρ_dyn", "χ_N/√(s_uv s_ir)", fmt(article_scenario.rho_dyn)),
            ("λ_u diag", "(1 − ρ_dyn²)^(-1)", fmt(article_scenario.lambda_u_diag)),
            ("N(D_*)", "dynamic numerator", fmt(scalar_dynamic_numerator(dc_lock, article_scenario))),
            ("Q(D_*)", "dynamic denominator", fmt(scalar_dynamic_denominator(dc_lock, article_scenario))),
            ("Φ_dyn(D_*)", "kernel contribution", fmt(exact_sens.phi_dyn_star)),
            ("R_moth(D_*)", "scalar radicand", fmt(exact_sens.r_star)),
            ("R_moth'(D_*)", "radicand derivative", fmt(exact_sens.radicand_prime_star)),
        ],
    ))

    # -------------------------------------------------------------------------
    # α lock tables
    # -------------------------------------------------------------------------
    write(render_table(
        "Table 14. Internal values of the current article geometric-lock output",
        ["Item", "Symbol", "Value"],
        [
            ("Current article lock output", "α_lock^(art)", fmt(alpha_article)),
            ("Inverse lock output", "(α_lock^(art))^(-1)", fmt(alpha_article_inv)),
            ("Geometric vector input", "Λ_geom", fmt(lambda_data.lambda_geom)),
            ("Vector overlap strength", "ζ_*", fmt(zeta_star)),
            ("Locked scalar value", "D_C^lock", fmt(dc_lock)),
            ("Scalar block at lock output", "D_C(α_lock^(art))", fmt(dc_at_alpha_article)),
            ("Bridge value", "D_total(α_lock^(art))", fmt(d_total_article)),
            ("Scalar lock residual", "D_C(α_lock^(art)) − D_C^lock", fmt_sci(dc_at_alpha_article - dc_lock, 14)),
        ],
    ))

    write(render_table(
        "Table 15. External comparison against the adopted benchmark α_ref",
        ["Realization", "α", "α⁻¹", "Δα", "Δ(α⁻¹)", "relative ppt", "n_sigma"],
        [
            (
                "Current article lock output",
                fmt(alpha_article, 24),
                fmt(alpha_article_inv, 24),
                fmt_sci(delta_alpha_article, 14),
                fmt_sci(alpha_article_inv - ALPHA_INV_REF, 14),
                fmt_ppt(delta_alpha_article, ALPHA_REF),
                fmt(delta_alpha_article / ALPHA_REF_SIGMA, 16),
            ),
            (
                "Universal-QED cross-check",
                fmt(alpha_qed_cross, 24),
                fmt(alpha_qed_cross_inv, 24),
                fmt_sci(delta_alpha_qed, 14),
                fmt_sci(alpha_qed_cross_inv - ALPHA_INV_REF, 14),
                fmt_ppt(delta_alpha_qed, ALPHA_REF),
                fmt(delta_alpha_qed / ALPHA_REF_SIGMA, 16),
            ),
            (
                "Article minus QED cross-check",
                fmt(delta_alpha_article_qed, 24),
                "—",
                fmt_sci(delta_alpha_article_qed, 14),
                "—",
                fmt_ppt(delta_alpha_article_qed, alpha_article),
                fmt(delta_alpha_article_qed / ALPHA_REF_SIGMA, 16),
            ),
        ],
    ))

    write(render_table(
        "Table 16. QED-induced scalar coefficients c_n^(QED) from the universal pure-photonic series",
        ["Coefficient", "Value"],
        [(f"c_{n}^(QED)", fmt(qdc_universal[n], 30)) for n in range(1, QED_SERIES_ORDER + 1)],
    ))

    write(render_table(
        "Table 17. Universal-QED inverse-α cross-check over the same geometric lock",
        ["Quantity", "Definition", "Value"],
        [
            ("D_C^lock", "same geometric lock target", fmt(dc_lock)),
            ("α_QED→D_C^[5]", "root of D_C^(QED,[5])(α/π)=D_C^lock", fmt(alpha_qed_cross)),
            ("(α_QED→D_C^[5])^(-1)", "inverse form", fmt(alpha_qed_cross_inv)),
            ("D_C^(QED,[5])(α_QED/π)", "QED-induced scalar at the root", fmt(dc_at_alpha_qed)),
            ("Residual", "D_C^(QED,[5]) − D_C^lock", fmt_sci(dc_at_alpha_qed - dc_lock, 14)),
        ],
    ))

    # -------------------------------------------------------------------------
    # Derived pure-photonic g-2 tables
    # -------------------------------------------------------------------------
    rows_rank_summary = []
    for n in RANK_AUDIT_LIST:
        rows_rank_summary.append(
            (
                f"Current rank-{n}",
                fmt(current_physical[n].dc, 24),
                fmt(current_physical[n].g_e_over_2, 24),
                fmt_sci(current_physical[n].g_e_over_2 - refined_physical[100].g_e_over_2, 14),
                fmt_sci(max_abs_delta_current[n], 14),
            )
        )
        rows_rank_summary.append(
            (
                f"Refined rank-{n}",
                fmt(refined_physical[n].dc, 24),
                fmt(refined_physical[n].g_e_over_2, 24),
                fmt_sci(refined_physical[n].g_e_over_2 - refined_physical[100].g_e_over_2, 14),
                fmt_sci(max_abs_delta_refined[n], 14),
            )
        )
    write(render_table(
        "Table 18. Pure-photonic rank summary for the current and refined evaluators",
        ["Model", "D_C(x_ref)", "g_{e,2}(x_ref)", "Δg_{e,2}(x_ref) vs refined-100", "max_{1≤n≤5}|ΔA_1^(2n)|"],
        rows_rank_summary,
    ))

    refined100 = refined_physical[100]
    write(render_table(
        "Table 19. Physical-point values entering the pure-photonic Relator cross-check (refined rank-100 realized)",
        ["Quantity", "Value"],
        [
            ("D_C(x_ref)", fmt(refined100.dc, 24)),
            ("ζ(D_C(x_ref))", fmt(refined100.zeta, 24)),
            ("Λ_π(x_ref)", fmt(refined100.lambda_pi, 24)),
            ("D_total(x_ref)", fmt(refined100.d_total, 24)),
            ("g_{e,2}^Rel(x_ref)", fmt(refined100.g_e_over_2, 24)),
            ("a_e^Rel(x_ref)", fmt(refined100.a_e, 24)),
        ],
    ))

    write(render_table(
        "Table 20. Comparison of the first five pure-photonic coefficients with universal QED",
        ["n", "A_1^(2n) QED", "Â_1^(2n) current-rank-5", "Â_1^(2n) refined-rank-100", "Δ current", "Δ refined", "Δ/QED current", "Δ/QED refined"],
        coefficient_error_rows(g2_coeffs_current_rank5, g2_coeffs_refined_rank100),
    ))

    write(render_table(
        "Table 21. First ten coefficients of the refined rank-100 realized pure-photonic series",
        ["n", "Â_1^(2n) Relator", "A_1^(2n) pure QED", "Status"],
        [
            (
                str(n),
                fmt(g2_coeffs_refined_rank100[n], 24),
                fmt(QED_A1_UNIVERSAL[n], 24) if n in QED_A1_UNIVERSAL else "---",
                "benchmarked" if n in QED_A1_UNIVERSAL else "Relator prediction",
            )
            for n in range(1, PREDICTION_ORDER + 1)
        ],
    ))

    # -------------------------------------------------------------------------
    # Sensitivity and stability tables
    # -------------------------------------------------------------------------
    write(render_table(
        "Table 22. Exact sensitivities at the current article lock",
        ["Quantity", "Value"],
        [
            ("D_* = D_C^lock(Λ_geom)", fmt(dc_lock)),
            ("R_* = R_moth(D_*)", fmt(exact_sens.r_star)),
            ("d D_C^lock / dΛ", fmt(exact_sens.d_lock_d_lambda)),
            ("d D_C^lock / dK", fmt(exact_sens.d_lock_d_k)),
            ("dα / dD at D_*", fmt(exact_sens.d_alpha_d_d)),
            ("dα_pred / dΛ_geom", fmt(exact_sens.d_alpha_d_lambda)),
            ("dα_pred / dK", fmt(exact_sens.d_alpha_d_k)),
            ("∂lnα_pred/∂lnΛ_geom", fmt(exact_sens.log_alpha_log_lambda)),
            ("∂lnα_pred/∂lnK", fmt(exact_sens.log_alpha_log_k)),
            ("∂α_pred/∂Θ_eff", fmt(exact_sens.d_alpha_d_theta)),
            ("∂lnα_pred/∂lnΘ_eff", fmt(exact_sens.log_alpha_log_theta)),
            ("∂lnα_pred/∂lnB_11^glue", fmt(exact_sens.log_alpha_log_b11)),
        ],
    ))

    write(render_table(
        "Table 23. Numerical logarithmic sensitivities of the current article scalar truncation",
        ["Parameter", f"S_q^({fmt_sci(SENSITIVITY_EPS, 1)})"],
        [
            ("Θ_eff", fmt(num_sens["Theta_eff"], 24)),
            ("A_uv^(5)", fmt(num_sens["A_uv"], 24)),
            ("A_ir^(5)", fmt(num_sens["A_ir"], 24)),
            ("χ_5", fmt(num_sens["chi"], 24)),
            ("s_uv", fmt(num_sens["s_uv"], 24)),
            ("s_ir", fmt(num_sens["s_ir"], 24)),
        ],
    ))

    write(render_table(
        "Table 24. Rank-stability of α_pred under the current scalar truncation",
        ["Rank N", "α_pred^(N)", "α_pred^(N) − α_pred^(100)"],
        [
            (
                str(n),
                fmt(alpha_rank_current[n], 34),
                fmt_sci(alpha_rank_current[n] - alpha_rank_current_ref, 14),
            )
            for n in RANK_AUDIT_LIST
        ],
    ))

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    write("Final primary result")
    write("────────────────────")
    write(f"Λ_geom = {fmt(lambda_data.lambda_geom)}")
    write(f"D_C^lock = {fmt(dc_lock)}")
    write(f"α_pred   = {fmt(alpha_article)}")
    write(f"α_pred⁻¹ = {fmt(alpha_article_inv)}")
    write()
    write("Reference-only note")
    write("───────────────────")
    write("Benchmark α and universal QED coefficients are used only in the external comparison tables, in the QED-induced D_C cross-check, and in the coefficient-error ledgers.")
    write()

    return out.getvalue()

# =============================================================================
# CLI
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute Λ_geom, D_C^lock, α_pred, pure-photonic Relator g-2 coefficients, sensitivity audits, and QED comparison tables."
    )
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS, help="mpmath decimal precision")
    parser.add_argument("--max-odd-mode", type=int, default=DEFAULT_MAX_ODD_MODE, help="largest odd shell mode j included")
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1, help="number of worker processes for shell-source modes")
    parser.add_argument("--gl-nodes", type=int, default=DEFAULT_GL_NODES, help="Gauss-Legendre nodes per shell integral")
    parser.add_argument("--print-digits", type=int, default=DEFAULT_PRINT_DIGITS, help="digits printed in tables")
    parser.add_argument("--output-txt", type=str, default="relator_alpha_full_audit_output.txt", help="path for the text report")
    args = parser.parse_args()

    if args.max_odd_mode < 1 or args.max_odd_mode % 2 == 0:
        raise ValueError("--max-odd-mode must be a positive odd integer.")
    if args.workers < 1:
        raise ValueError("--workers must be at least 1.")

    mp.mp.dps = args.dps
    report_text = build_report_text(
        max_odd_mode=args.max_odd_mode,
        workers=args.workers,
        dps=args.dps,
        gl_nodes=args.gl_nodes,
        print_digits=args.print_digits,
        output_txt=args.output_txt,
    )
    output_path = Path(args.output_txt)
    output_path.write_text(report_text, encoding="utf-8")
    print(report_text, end="")

if __name__ == "__main__":
    main()
