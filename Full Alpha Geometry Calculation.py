from __future__ import annotations

"""
Closed numerical reproducer for the Relator alpha workflow.

Core principle
--------------
The geometric lock is built independently from alpha and independently from any
QED coefficient input:

    Lambda_geom  ->  zeta_lock  ->  D_C^lock  ->  alpha_article

Optional external diagnostics
-----------------------------
The experimental alpha benchmark and the first five pure-photonic QED
coefficients are kept in a single structural mode bundle. They are never used
inside the geometric lock. They are used only when the user explicitly switches
from ARTICLE_ONLY_MODE to ARTICLE_PLUS_QED_DIAGNOSTICS.

How to switch modes
-------------------
Comment exactly one of the two lines in the block marked
"ACTIVE_STRUCTURAL_MODE" below.

- ARTICLE_ONLY_MODE:
    * no external alpha comparison
    * no QED-shock branch in the summary values
    * all optional cells are printed as "×"

- ARTICLE_PLUS_QED_DIAGNOSTICS:
    * enables alpha_exp / sigma(alpha) comparisons
    * enables the pure-QED g2 / a_e shock branch
    * enables the regressed D_C^QED(alpha) branch
    * enables benchmark deltas in the coefficient tables

Vector-channel status used for the lock
---------------------------------------
The lock is driven by the current article mean-shell diagnostic representative
Lambda_geom^(mean). This is independent of alpha/QED, but it is still a vector
representative because the exact operator evaluations P_IR^(chi)(ell_0) and
R_chi are not yet closed in the manuscript as theorem-level one-line formulas.

The script also prints a raw Maxwell shell audit computed directly from the
shell-source formulas. That raw audit is informative only; it is not fed into
Lambda_geom^(mean) unless the user rewrites the vector representative block on
purpose.

Important output emphasis
-------------------------
The two most important final values are printed twice:
1) inside the full lock summary table, and
2) inside a dedicated highlighted box near the top of the report.

That dedicated box always prints alpha_article. It prints alpha_QED,regressed
numerically only when ARTICLE_PLUS_QED_DIAGNOSTICS is active; otherwise it
prints the structural placeholder "×".
"""

from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import mpmath as mp
import numpy as np
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from scipy import special


# =============================================================================
# User-facing settings
# =============================================================================
mp.mp.dps = 110

PRINT_DIGITS = 28
SCI_DIGITS = 12
KEY_BOX_DIGITS = 34
PREDICTION_ORDER = 10
CURRENT_ARTICLE_RANK = 5
REFINED_REFERENCE_RANK = 100
VECTOR_SOURCE_MAX_J = 31
VECTOR_SOURCE_GL_NODES = 400
RICH_CONSOLE_WIDTH = 180


@dataclass(frozen=True)
class ExternalDiagnosticInputs:
    alpha_inv_exp: mp.mpf
    alpha_inv_exp_sigma: mp.mpf
    pure_qed_a1: Dict[int, mp.mpf]
    qed_order: int = 5


@dataclass(frozen=True)
class StructuralMode:
    name: str
    external_diagnostics: Optional[ExternalDiagnosticInputs]
    placeholder: str = "×"

    @property
    def diagnostics_enabled(self) -> bool:
        return self.external_diagnostics is not None


OPTIONAL_EXTERNAL_DIAGNOSTICS = ExternalDiagnosticInputs(
    alpha_inv_exp=mp.mpf("137.035999177"),
    alpha_inv_exp_sigma=mp.mpf("0.000000021"),
    pure_qed_a1={
        1: mp.mpf("0.5"),
        2: mp.mpf("-0.328478965579193"),
        3: mp.mpf("1.181241456587"),
        4: mp.mpf("-1.9122457649264455741526471531265"),
        5: mp.mpf("5.891"),
    },
    qed_order=5,
)

ARTICLE_ONLY_MODE = StructuralMode(
    name="ARTICLE_ONLY_MODE",
    external_diagnostics=None,
    placeholder="×",
)

ARTICLE_PLUS_QED_DIAGNOSTICS = StructuralMode(
    name="ARTICLE_PLUS_QED_DIAGNOSTICS",
    external_diagnostics=OPTIONAL_EXTERNAL_DIAGNOSTICS,
    placeholder="×",
)

# =============================================================================
# ACTIVE_STRUCTURAL_MODE
# Comment exactly one line and keep the other one active.
# =============================================================================
# ACTIVE_STRUCTURAL_MODE = ARTICLE_ONLY_MODE
ACTIVE_STRUCTURAL_MODE = ARTICLE_PLUS_QED_DIAGNOSTICS


# =============================================================================
# Exact constants and common helpers
# =============================================================================
ONE = mp.mpf("1")
TWO = mp.mpf("2")
THREE = mp.mpf("3")
FOUR = mp.mpf("4")
HALF = mp.mpf("0.5")
QUARTER = mp.mpf("0.25")

PI = mp.pi
SQRT_PI = mp.sqrt(PI)
GAMMA_E = mp.euler


def geometric_shell_constant_K() -> mp.mpf:
    return (150 * PI**2 - 8 * PI**4 - 315) / (180 * PI**6)


K_GEOMETRIC = geometric_shell_constant_K()


def fmt(value: mp.mpf | float, digits: int = PRINT_DIGITS) -> str:
    return mp.nstr(mp.mpf(value), n=digits)


def fmt_sci(value: mp.mpf | float, digits: int = SCI_DIGITS) -> str:
    return mp.nstr(mp.mpf(value), n=digits, min_fixed=0, max_fixed=0)


def placeholder() -> str:
    return ACTIVE_STRUCTURAL_MODE.placeholder


def fmt_optional(value: Optional[mp.mpf | float], digits: int = PRINT_DIGITS, sci: bool = False) -> str:
    if value is None:
        return placeholder()
    return fmt_sci(value, digits) if sci else fmt(value, digits)


# =============================================================================
# Rich printing helpers
# =============================================================================
console = Console(width=RICH_CONSOLE_WIDTH, record=True)


def make_table(title: str, columns: Sequence[Tuple[str, str]]) -> Table:
    table = Table(title=title, box=box.ROUNDED, title_style="bold")
    for header, justify in columns:
        table.add_column(header, justify=justify, overflow="fold")
    return table


def print_table(title: str, columns: Sequence[Tuple[str, str]], rows: Iterable[Sequence[str]]) -> None:
    table = make_table(title, columns)
    for row in rows:
        table.add_row(*[str(cell) for cell in row])
    console.print(table)


# =============================================================================
# Section A — Exact theorem-path shell source and raw Lambda_OUT(eta0)
# =============================================================================

def minimal_branch_ratio() -> mp.mpf:
    return ONE / PI



def minimal_shell_width() -> mp.mpf:
    return ONE / (PI * SQRT_PI)



def filamentary_reference_core() -> mp.mpf:
    return mp.log(8 * mp.sqrt(mp.pi)) - 2



def gaussian_uv_constant() -> mp.mpf:
    return HALF * (mp.log(2) + mp.euler)


ETA_0 = minimal_branch_ratio()
ELL_0 = minimal_shell_width()
LAMBDA_IND = filamentary_reference_core()
C_UV_GAUSS = gaussian_uv_constant()



def flux_mode_norm_I(j: int) -> float:
    return 2.0 * j * (j + 1) / (2 * j + 1)



def legendre_derivative_np(j: int, x: np.ndarray) -> np.ndarray:
    p_j = special.eval_legendre(j, x)
    p_jm1 = special.eval_legendre(j - 1, x)
    return j * (p_jm1 - x * p_j) / (1.0 - x * x)



def shell_k2(eta: float, v: np.ndarray) -> np.ndarray:
    rho = np.sqrt(np.clip(1.0 - v * v, 0.0, None))
    z = v
    return 4.0 * eta * rho / ((eta + rho) ** 2 + z * z)



def shell_Brho_shape(eta: float, v: np.ndarray) -> np.ndarray:
    rho = np.sqrt(np.clip(1.0 - v * v, 0.0, None))
    z = v
    m = shell_k2(eta, v)
    K = special.ellipk(m)
    E = special.ellipe(m)
    den = np.sqrt((eta + rho) ** 2 + z * z)
    rho_safe = np.where(rho < 1.0e-15, 1.0e-15, rho)
    ratio = (eta * eta + rho * rho + z * z) / (((eta - rho) ** 2) + z * z)
    return z / (2.0 * np.pi * rho_safe * den) * (-K + ratio * E)



def shell_Bz_shape(eta: float, v: np.ndarray) -> np.ndarray:
    rho = np.sqrt(np.clip(1.0 - v * v, 0.0, None))
    z = v
    m = shell_k2(eta, v)
    K = special.ellipk(m)
    E = special.ellipe(m)
    den = np.sqrt((eta + rho) ** 2 + z * z)
    ratio = (eta * eta - rho * rho - z * z) / (((eta - rho) ** 2) + z * z)
    return 1.0 / (2.0 * np.pi * den) * (K + ratio * E)



def shell_btheta_shape(eta: float, v: np.ndarray) -> np.ndarray:
    rho = np.sqrt(np.clip(1.0 - v * v, 0.0, None))
    z = v
    return z * shell_Brho_shape(eta, v) - rho * shell_Bz_shape(eta, v)



def compute_shell_source_coefficients(eta: float, max_j: int, n_gl: int) -> Dict[int, mp.mpf]:
    x, w = np.polynomial.legendre.leggauss(n_gl)
    b_tilde = shell_btheta_shape(eta, x)
    rho = np.sqrt(1.0 - x * x)

    out: Dict[int, mp.mpf] = {}
    for j in range(1, max_j + 1, 2):
        p_prime = legendre_derivative_np(j, x)
        integrand = b_tilde * rho * p_prime
        a_sh_j = np.sum(w * integrand) / flux_mode_norm_I(j)
        out[j] = mp.mpf(str(a_sh_j))
    return out



def compute_shell_source_Jchi(a_sh: Dict[int, mp.mpf]) -> Dict[int, mp.mpf]:
    out: Dict[int, mp.mpf] = {}
    for j, value in a_sh.items():
        a_hat = 2 * (j + 1) * value
        factor = mp.sqrt(2 * mp.pi / ((j + 1) * (2 * j + 1)))
        out[j] = factor * a_hat
    return out



def lambda_out_from_Jchi(jchi: Dict[int, mp.mpf]) -> mp.mpf:
    return -HALF * mp.fsum(value * value for value in jchi.values())



def normalized_source_weights(jchi: Dict[int, mp.mpf]) -> Dict[int, mp.mpf]:
    denom = mp.fsum(value * value for value in jchi.values())
    return {j: (value * value) / denom for j, value in jchi.items()}


RAW_VECTOR_A_SH = compute_shell_source_coefficients(float(ETA_0), VECTOR_SOURCE_MAX_J, VECTOR_SOURCE_GL_NODES)
RAW_VECTOR_JCHI = compute_shell_source_Jchi(RAW_VECTOR_A_SH)
RAW_VECTOR_LAMBDA_OUT = lambda_out_from_Jchi(RAW_VECTOR_JCHI)
RAW_VECTOR_WEIGHTS = normalized_source_weights(RAW_VECTOR_JCHI)


# =============================================================================
# Section B — Article vector representative used for the independent lock
# =============================================================================
@dataclass(frozen=True)
class ArticleVectorRepresentative:
    lambda_uv_to_ir_mean: mp.mpf
    lambda_out_article: mp.mpf
    exterior_memory_total_mean: mp.mpf
    source_weights_article: Dict[int, mp.mpf]

    @property
    def r_chi_mean(self) -> mp.mpf:
        return self.exterior_memory_total_mean / self.lambda_out_article

    @property
    def p_ir_mean(self) -> mp.mpf:
        return self.lambda_uv_to_ir_mean / C_UV_GAUSS

    @property
    def lambda_geom_mean(self) -> mp.mpf:
        return LAMBDA_IND + self.lambda_uv_to_ir_mean + self.exterior_memory_total_mean


ARTICLE_VECTOR_REPRESENTATIVE = ArticleVectorRepresentative(
    lambda_uv_to_ir_mean=mp.mpf("0.0544853495865520906532491466693"),
    lambda_out_article=mp.mpf("-0.0139671580625860254655225000000"),
    exterior_memory_total_mean=mp.mpf("-0.0146086880840972703041996518997"),
    source_weights_article={
        1: mp.mpf("0.9796363359550031216197956500"),
        3: mp.mpf("0.02031062626884619091016822831"),
        5: mp.mpf("4.849145355825759038791211476e-5"),
        7: mp.mpf("4.354340825127690409129907468e-6"),
    },
)

LAMBDA_GEOM = ARTICLE_VECTOR_REPRESENTATIVE.lambda_geom_mean


# =============================================================================
# Section C — Scalar channel: current evaluator and refined no-free-parameter evaluator
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


@dataclass(frozen=True)
class ScalarScenario:
    label: str
    rank: int
    theta1_total: mp.mpf
    a_uv: mp.mpf
    a_ir: mp.mpf
    chi_cross: mp.mpf
    delta_a_mix: mp.mpf
    delta_a0: mp.mpf
    delta_a_diag: mp.mpf
    kernel_scale: mp.mpf


@lru_cache(maxsize=None)
def F_p_cached(p: int, ell: mp.mpf) -> mp.mpf:
    p_mpf = mp.mpf(p)
    phase = 1j * p_mpf * PI * ell / 2
    return (
        (SQRT_PI * ell / 2)
        * mp.e ** (-(p_mpf * PI * ell) ** 2 / 4)
        * mp.re(mp.e ** (1j * p_mpf * PI) * (mp.erf(ONE / ell + phase) - mp.erf(phase)))
    )



def c11_closed_form(ell: mp.mpf) -> mp.mpf:
    term1 = (SQRT_PI * ell / 2) * mp.erf(ONE / ell)
    term2 = (SQRT_PI * ell / 2) * mp.e ** (-(PI**2) * (ell**2))
    term2 *= mp.re(mp.erf(ONE / ell - 1j * PI * ell) - mp.erf(-1j * PI * ell))
    return term1 + term2



def C_mn(m: int, n: int, ell: mp.mpf) -> mp.mpf:
    return m * n * (F_p_cached(abs(m - n), ell) + F_p_cached(m + n, ell))



def gamma_clog(D: mp.mpf, ell: mp.mpf) -> mp.mpf:
    root = mp.sqrt(D + QUARTER)
    return (
        (SQRT_PI * ell / 2)
        * mp.e ** ((ell**2 / 4) * (D + QUARTER))
        * mp.erfc((ell / 2) * root)
    )



def build_base_shell_geometry() -> BaseShellGeometry:
    eta0 = ONE / PI
    ell0 = ONE / (PI * SQRT_PI)
    c_uni = (ONE / PI) * (mp.mpf("4") / 3 + ONE / (4 * PI**2))
    c11 = c11_closed_form(ell0)

    theta1_core = TWO * PI * c_uni
    theta1_coll = (GAMMA_E / PI**2) * c11
    theta1_log = eta0**2 / (8 * (ONE - eta0**2 / 4))
    a_shell = PI**2 * c11 + GAMMA_E / 4
    s_uv = mp.log(2)
    s_ir = ONE / (8 * PI**2)

    gamma0 = gamma_clog(mp.mpf("0"), ell0)
    b11_glue = (
        theta1_coll + theta1_log - gamma0 * theta1_coll * theta1_log
    ) / (ONE - (gamma0**2 / 4) * theta1_coll * theta1_log)

    return BaseShellGeometry(
        c11=c11,
        theta1_core=theta1_core,
        theta1_coll=theta1_coll,
        theta1_log=theta1_log,
        a_shell=a_shell,
        s_uv=s_uv,
        s_ir=s_ir,
        eta0=eta0,
        ell0=ell0,
        gamma_clog_0=gamma0,
        b11_glue=b11_glue,
    )


BASE_GEOMETRY = build_base_shell_geometry()



def higher_mode_sums(rank: int, ell: mp.mpf) -> Tuple[mp.mpf, mp.mpf]:
    delta_a_mix = mp.mpf("0")
    delta_a0 = mp.mpf("0")
    if rank < 2:
        return delta_a_mix, delta_a0
    for n in range(2, rank + 1):
        c1n = C_mn(1, n, ell)
        cnn = C_mn(n, n, ell)
        gap = n**2 - 1
        delta_a_mix += (GAMMA_E / PI**2) * (c1n**2 / gap)
        delta_a0 += -(GAMMA_E / PI**2) * (c1n**2 * cnn / gap**3)
    return delta_a_mix, delta_a0



def build_current_scenario(rank: int) -> ScalarScenario:
    delta_a_mix, delta_a0 = higher_mode_sums(rank, BASE_GEOMETRY.ell0)
    theta1_total = BASE_GEOMETRY.theta1_core + BASE_GEOMETRY.b11_glue
    delta_a_diag = -(BASE_GEOMETRY.eta0**2) * delta_a0
    chi_cross = delta_a_mix**2
    a_uv = BASE_GEOMETRY.a_shell + delta_a_mix + delta_a0 + delta_a_diag
    a_ir = BASE_GEOMETRY.a_shell - delta_a_mix + delta_a0 + delta_a_diag
    return ScalarScenario(
        label=f"current_rank_{rank}",
        rank=rank,
        theta1_total=theta1_total,
        a_uv=a_uv,
        a_ir=a_ir,
        chi_cross=chi_cross,
        delta_a_mix=delta_a_mix,
        delta_a0=delta_a0,
        delta_a_diag=delta_a_diag,
        kernel_scale=ONE,
    )



def build_static_hidden_completion() -> Tuple[mp.mpf, mp.mpf]:
    a = BASE_GEOMETRY.theta1_coll
    b = BASE_GEOMETRY.theta1_log
    gamma0 = BASE_GEOMETRY.gamma_clog_0

    N_c = mp.sqrt(PI / 8) * BASE_GEOMETRY.ell0
    rho_static = gamma0 / mp.sqrt(N_c)
    lambda_hidden = ONE / (ONE - rho_static**2)

    c0 = gamma0 * a * b / TWO
    a_perp = a - c0**2 / b
    b_perp = b - c0**2 / a
    t1 = mp.sqrt(a_perp / a)
    t2 = mp.sqrt(b_perp / b)

    a11 = ONE / a_perp - c0 / lambda_hidden
    a12 = gamma0 / TWO - c0 / lambda_hidden
    a22 = ONE / b_perp - c0 / lambda_hidden
    det2 = a11 * a22 - a12**2
    b11_refined = (t1**2 * a22 - TWO * t1 * t2 * a12 + t2**2 * a11) / det2
    theta1_refined = BASE_GEOMETRY.theta1_core + b11_refined
    return b11_refined, theta1_refined


B11_REFINED, THETA1_REFINED = build_static_hidden_completion()



def build_dynamic_kernel_scale(current: ScalarScenario) -> mp.mpf:
    rho_dyn = current.chi_cross / mp.sqrt(BASE_GEOMETRY.s_uv * BASE_GEOMETRY.s_ir)
    lambda_u = ONE / (ONE - rho_dyn**2)
    return lambda_u**2



def build_refined_scenario(rank: int) -> ScalarScenario:
    current = build_current_scenario(rank)
    return ScalarScenario(
        label=f"refined_rank_{rank}",
        rank=rank,
        theta1_total=THETA1_REFINED,
        a_uv=current.a_uv,
        a_ir=current.a_ir,
        chi_cross=current.chi_cross,
        delta_a_mix=current.delta_a_mix,
        delta_a0=current.delta_a0,
        delta_a_diag=current.delta_a_diag,
        kernel_scale=build_dynamic_kernel_scale(current),
    )


CURRENT_RANK5 = build_current_scenario(CURRENT_ARTICLE_RANK)
REFINED_RANK100 = build_refined_scenario(REFINED_REFERENCE_RANK)



def model_memory_value(D: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    num = scenario.a_uv * (ONE + BASE_GEOMETRY.s_ir * D) + scenario.a_ir * (ONE + BASE_GEOMETRY.s_uv * D)
    num -= TWO * scenario.chi_cross * D * mp.sqrt(scenario.a_uv * scenario.a_ir)
    den = (ONE + BASE_GEOMETRY.s_uv * D) * (ONE + BASE_GEOMETRY.s_ir * D) - (scenario.chi_cross**2) * (D**2)
    return scenario.kernel_scale * (num / den)



def mother_radicand(D: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    return ONE - scenario.theta1_total * D + D**2 * model_memory_value(D, scenario)



def solve_D_of_alpha(alpha: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    x = alpha / PI

    def residual(D: mp.mpf) -> mp.mpf:
        rad = mother_radicand(D, scenario)
        if rad <= 0:
            return mp.mpf("1e100")
        return D - x * mp.sqrt(rad)

    left = mp.mpf("0")
    right = max(mp.mpf("0.003"), TWO * x)
    f_left = residual(left)
    f_right = residual(right)
    if not (f_left < 0 < f_right):
        raise RuntimeError(f"Failed to bracket D_C(alpha) for {scenario.label}.")

    for _ in range(1500):
        mid = (left + right) / TWO
        f_mid = residual(mid)
        if abs(f_mid) < mp.mpf("1e-95") or abs(right - left) < mp.mpf("1e-95"):
            return mid
        if f_mid > 0:
            right = mid
        else:
            left = mid
    return (left + right) / TWO



def alpha_from_locked_D(D_lock: mp.mpf, scenario: ScalarScenario) -> mp.mpf:
    rad = mother_radicand(D_lock, scenario)
    if rad <= 0:
        raise RuntimeError(f"Locked D produces non-positive mother radicand for {scenario.label}.")
    return PI * D_lock / mp.sqrt(rad)


# =============================================================================
# Section D — ALP bridge and geometric lock maps
# =============================================================================

def zeta_from_D(D: mp.mpf) -> mp.mpf:
    return HALF * (mp.sqrt(ONE + FOUR * D / (THREE * PI**2)) - ONE)



def zeta_from_lambda(lambda_value: mp.mpf) -> mp.mpf:
    return (K_GEOMETRIC / (TWO * PI**2)) * lambda_value



def lambda_from_D(D: mp.mpf) -> mp.mpf:
    return (TWO * PI**2 / K_GEOMETRIC) * zeta_from_D(D)



def D_lock_from_lambda(lambda_value: mp.mpf) -> mp.mpf:
    return (THREE / TWO) * K_GEOMETRIC * lambda_value * (ONE + K_GEOMETRIC * lambda_value / (TWO * PI**2))



def D_total_from_D_and_alpha(D: mp.mpf, alpha: mp.mpf) -> mp.mpf:
    x = alpha / PI
    zeta = zeta_from_D(D)
    return D - x * zeta + (x**2) * (zeta**2) / (FOUR * D)



def g2_from_D_and_alpha(D: mp.mpf, alpha: mp.mpf) -> mp.mpf:
    D_total = D_total_from_D_and_alpha(D, alpha)
    return ONE / mp.sqrt(ONE - D_total)


D_LOCK = D_lock_from_lambda(LAMBDA_GEOM)
ZETA_LOCK = zeta_from_lambda(LAMBDA_GEOM)
ALPHA_ARTICLE = alpha_from_locked_D(D_LOCK, CURRENT_RANK5)
ALPHA_ARTICLE_INV = ONE / ALPHA_ARTICLE
G2_ARTICLE = g2_from_D_and_alpha(D_LOCK, ALPHA_ARTICLE)
A_E_ARTICLE = G2_ARTICLE - ONE


# =============================================================================
# Section E — Formal power series tools
# =============================================================================

def zero_series(order: int) -> List[mp.mpf]:
    return [mp.mpf("0") for _ in range(order + 1)]



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
        raise ValueError("Series square root requires positive constant term.")
    out = zero_series(order)
    out[0] = mp.sqrt(a[0])
    for n in range(1, order + 1):
        s = mp.fsum(out[k] * out[n - k] for k in range(1, n))
        out[n] = (a[n] - s) / (TWO * out[0])
    return out



def series_shift_up(a: List[mp.mpf], m: int, order: int) -> List[mp.mpf]:
    out = zero_series(order)
    for n in range(order + 1 - m):
        out[n + m] = a[n]
    return out



def series_shift_down(a: List[mp.mpf], m: int, order: int) -> List[mp.mpf]:
    out = zero_series(order)
    for n in range(order + 1 - m):
        out[n] = a[n + m]
    return out



def series_eval(coeffs: List[mp.mpf], x_value: mp.mpf) -> mp.mpf:
    total = mp.mpf("0")
    power = mp.mpf("1")
    for coefficient in coeffs:
        total += coefficient * power
        power *= x_value
    return total



def dict_to_series(coeffs: Dict[int, mp.mpf], order: int) -> List[mp.mpf]:
    out = zero_series(order)
    for n, value in coeffs.items():
        if 0 <= n <= order:
            out[n] = mp.mpf(value)
    return out



def current_memory_series(D: List[mp.mpf], order: int, scenario: ScalarScenario) -> List[mp.mpf]:
    one = one_series(order)
    den = series_sub(
        series_mul(
            series_add(one, series_scale(D, BASE_GEOMETRY.s_uv)),
            series_add(one, series_scale(D, BASE_GEOMETRY.s_ir)),
            order,
        ),
        series_scale(series_mul(D, D, order), scenario.chi_cross**2),
    )
    den_inv = series_inv(den, order)
    num = series_add(
        series_scale(series_add(one, series_scale(D, BASE_GEOMETRY.s_ir)), scenario.a_uv),
        series_scale(series_add(one, series_scale(D, BASE_GEOMETRY.s_uv)), scenario.a_ir),
    )
    num = series_sub(num, series_scale(D, TWO * scenario.chi_cross * mp.sqrt(scenario.a_uv * scenario.a_ir)))
    return series_scale(series_mul(num, den_inv, order), scenario.kernel_scale)



def build_D_series(order: int, scenario: ScalarScenario) -> List[mp.mpf]:
    D = zero_series(order)
    D[1] = ONE
    one = one_series(order)
    for _ in range(10 * order):
        D2 = series_mul(D, D, order)
        G = current_memory_series(D, order, scenario)
        R = series_sub(one, series_scale(D, scenario.theta1_total))
        R = series_add(R, series_mul(D2, G, order))
        D_new = series_shift_up(series_sqrt(R, order), 1, order)
        if all(abs(D_new[n] - D[n]) < mp.mpf("1e-95") for n in range(order + 1)):
            return D_new
        D = D_new
    return D



def zeta_series_from_D(D: List[mp.mpf], order: int) -> List[mp.mpf]:
    inside = series_add(one_series(order), series_scale(D, FOUR / (THREE * PI**2)))
    return series_scale(series_sub(series_sqrt(inside, order), one_series(order)), HALF)



def D_total_series_from_D(D: List[mp.mpf], order: int) -> Tuple[List[mp.mpf], List[mp.mpf], List[mp.mpf]]:
    zeta = zeta_series_from_D(D, order)
    D_cross = series_shift_up(zeta, 1, order)
    zeta_sq = series_mul(zeta, zeta, order)
    D_hat = series_shift_down(D, 1, order)
    D_hat_inv = series_inv(D_hat, order)
    D_vector_core = series_mul(zeta_sq, D_hat_inv, order)
    D_vector = series_shift_up(series_scale(D_vector_core, QUARTER), 1, order)
    D_total = series_add(series_sub(D, D_cross), D_vector)
    return zeta, D_vector, D_total



def g2_series_from_D_total(D_total: List[mp.mpf], order: int) -> List[mp.mpf]:
    return series_inv(series_sqrt(series_sub(one_series(order), D_total), order), order)



def a_series_from_D_total(D_total: List[mp.mpf], order: int) -> List[mp.mpf]:
    g2 = g2_series_from_D_total(D_total, order)
    out = g2[:]
    out[0] -= ONE
    return out


# =============================================================================
# Section F — Closed Relator coefficient ledgers
# =============================================================================
CURRENT_D_SERIES = build_D_series(PREDICTION_ORDER, CURRENT_RANK5)
REFINED_D_SERIES = build_D_series(PREDICTION_ORDER, REFINED_RANK100)

_, _, CURRENT_D_TOTAL_SERIES = D_total_series_from_D(CURRENT_D_SERIES, PREDICTION_ORDER)
_, _, REFINED_D_TOTAL_SERIES = D_total_series_from_D(REFINED_D_SERIES, PREDICTION_ORDER)
CURRENT_G2_SERIES = g2_series_from_D_total(CURRENT_D_TOTAL_SERIES, PREDICTION_ORDER)
REFINED_G2_SERIES = g2_series_from_D_total(REFINED_D_TOTAL_SERIES, PREDICTION_ORDER)
CURRENT_A_SERIES = a_series_from_D_total(CURRENT_D_TOTAL_SERIES, PREDICTION_ORDER)
REFINED_A_SERIES = a_series_from_D_total(REFINED_D_TOTAL_SERIES, PREDICTION_ORDER)


# =============================================================================
# Section G — Optional pure-QED diagnostics (external only; never used in the lock)
# =============================================================================
@dataclass(frozen=True)
class OptionalDiagnosticsResults:
    alpha_exp: mp.mpf
    alpha_exp_sigma: mp.mpf
    alpha_inv_exp: mp.mpf
    alpha_inv_exp_sigma: mp.mpf
    a_qed_series: List[mp.mpf]
    dc_qed_series: List[mp.mpf]
    g2_qed_at_alpha_exp: mp.mpf
    a_e_qed_at_alpha_exp: mp.mpf
    d_total_qed_at_alpha_exp: mp.mpf
    d_c_qed_at_alpha_exp: mp.mpf
    alpha_qed_regressed: mp.mpf
    alpha_qed_regressed_inv: mp.mpf
    g2_qed_at_regressed_alpha: mp.mpf
    a_e_qed_at_regressed_alpha: mp.mpf
    d_total_qed_at_regressed_alpha: mp.mpf
    d_c_qed_at_regressed_alpha: mp.mpf
    delta_alpha_article_vs_exp: mp.mpf
    delta_alpha_qed_vs_exp: mp.mpf
    delta_alpha_qed_vs_article: mp.mpf
    z_sigma_article: mp.mpf
    z_sigma_qed: mp.mpf



def D_total_series_from_qed_ae(qed_coeffs: Dict[int, mp.mpf], order: int) -> List[mp.mpf]:
    a_e = dict_to_series(qed_coeffs, order)
    g2 = series_add(one_series(order), a_e)
    inv_g2 = series_inv(g2, order)
    inv_sq = series_mul(inv_g2, inv_g2, order)
    return series_sub(one_series(order), inv_sq)



def solve_D_series_from_D_total_target(D_total_target: List[mp.mpf], order: int) -> List[mp.mpf]:
    D = D_total_target[:]
    D[0] = mp.mpf("0")
    D[1] = ONE
    for _ in range(10 * order):
        zeta = zeta_series_from_D(D, order)
        D_cross = series_shift_up(zeta, 1, order)
        zeta_sq = series_mul(zeta, zeta, order)
        D_hat = series_shift_down(D, 1, order)
        D_hat_inv = series_inv(D_hat, order)
        D_vector_core = series_mul(zeta_sq, D_hat_inv, order)
        D_vector = series_shift_up(series_scale(D_vector_core, QUARTER), 1, order)
        D_new = series_sub(series_add(D_total_target, D_cross), D_vector)
        D_new[0] = mp.mpf("0")
        D_new[1] = ONE
        if all(abs(D_new[n] - D[n]) < mp.mpf("1e-95") for n in range(order + 1)):
            return D_new
        D = D_new
    return D



def build_qed_dc_series(qed_coeffs: Dict[int, mp.mpf], order: int) -> List[mp.mpf]:
    D_total_target = D_total_series_from_qed_ae(qed_coeffs, order)
    return solve_D_series_from_D_total_target(D_total_target, order)



def solve_alpha_from_dc_model(dc_model, D_target: mp.mpf) -> mp.mpf:
    def residual(alpha_value: mp.mpf) -> mp.mpf:
        return dc_model(alpha_value) - D_target

    left = mp.mpf("0.001")
    right = mp.mpf("0.02")
    f_left = residual(left)
    f_right = residual(right)
    if not (f_left < 0 < f_right):
        raise RuntimeError("Failed to bracket the alpha root for the QED-induced D_C model.")

    for _ in range(1500):
        mid = (left + right) / TWO
        f_mid = residual(mid)
        if abs(f_mid) < mp.mpf("1e-95") or abs(right - left) < mp.mpf("1e-95"):
            return mid
        if f_mid > 0:
            right = mid
        else:
            left = mid
    return (left + right) / TWO



def build_optional_diagnostics(mode: StructuralMode, d_lock: mp.mpf, alpha_article: mp.mpf) -> Optional[OptionalDiagnosticsResults]:
    if mode.external_diagnostics is None:
        return None

    ext = mode.external_diagnostics
    alpha_inv_exp = ext.alpha_inv_exp
    alpha_inv_exp_sigma = ext.alpha_inv_exp_sigma
    alpha_exp = ONE / alpha_inv_exp
    alpha_exp_sigma = alpha_inv_exp_sigma / (alpha_inv_exp ** 2)

    a_qed_series = dict_to_series(ext.pure_qed_a1, ext.qed_order)
    dc_qed_series = build_qed_dc_series(ext.pure_qed_a1, ext.qed_order)

    def a_e_qed_of_alpha(alpha: mp.mpf) -> mp.mpf:
        return series_eval(a_qed_series, alpha / PI)

    def g2_qed_of_alpha(alpha: mp.mpf) -> mp.mpf:
        return ONE + a_e_qed_of_alpha(alpha)

    def d_c_qed_of_alpha(alpha: mp.mpf) -> mp.mpf:
        return series_eval(dc_qed_series, alpha / PI)

    g2_qed_at_alpha_exp = g2_qed_of_alpha(alpha_exp)
    a_e_qed_at_alpha_exp = g2_qed_at_alpha_exp - ONE
    d_total_qed_at_alpha_exp = ONE - ONE / (g2_qed_at_alpha_exp**2)
    d_c_qed_at_alpha_exp = d_c_qed_of_alpha(alpha_exp)

    alpha_qed_regressed = solve_alpha_from_dc_model(d_c_qed_of_alpha, d_lock)
    alpha_qed_regressed_inv = ONE / alpha_qed_regressed
    g2_qed_at_regressed_alpha = g2_qed_of_alpha(alpha_qed_regressed)
    a_e_qed_at_regressed_alpha = g2_qed_at_regressed_alpha - ONE
    d_total_qed_at_regressed_alpha = ONE - ONE / (g2_qed_at_regressed_alpha**2)
    d_c_qed_at_regressed_alpha = d_c_qed_of_alpha(alpha_qed_regressed)

    delta_alpha_article_vs_exp = alpha_article - alpha_exp
    delta_alpha_qed_vs_exp = alpha_qed_regressed - alpha_exp
    delta_alpha_qed_vs_article = alpha_qed_regressed - alpha_article
    z_sigma_article = delta_alpha_article_vs_exp / alpha_exp_sigma
    z_sigma_qed = delta_alpha_qed_vs_exp / alpha_exp_sigma

    return OptionalDiagnosticsResults(
        alpha_exp=alpha_exp,
        alpha_exp_sigma=alpha_exp_sigma,
        alpha_inv_exp=alpha_inv_exp,
        alpha_inv_exp_sigma=alpha_inv_exp_sigma,
        a_qed_series=a_qed_series,
        dc_qed_series=dc_qed_series,
        g2_qed_at_alpha_exp=g2_qed_at_alpha_exp,
        a_e_qed_at_alpha_exp=a_e_qed_at_alpha_exp,
        d_total_qed_at_alpha_exp=d_total_qed_at_alpha_exp,
        d_c_qed_at_alpha_exp=d_c_qed_at_alpha_exp,
        alpha_qed_regressed=alpha_qed_regressed,
        alpha_qed_regressed_inv=alpha_qed_regressed_inv,
        g2_qed_at_regressed_alpha=g2_qed_at_regressed_alpha,
        a_e_qed_at_regressed_alpha=a_e_qed_at_regressed_alpha,
        d_total_qed_at_regressed_alpha=d_total_qed_at_regressed_alpha,
        d_c_qed_at_regressed_alpha=d_c_qed_at_regressed_alpha,
        delta_alpha_article_vs_exp=delta_alpha_article_vs_exp,
        delta_alpha_qed_vs_exp=delta_alpha_qed_vs_exp,
        delta_alpha_qed_vs_article=delta_alpha_qed_vs_article,
        z_sigma_article=z_sigma_article,
        z_sigma_qed=z_sigma_qed,
    )


OPTIONAL_RESULTS = build_optional_diagnostics(ACTIVE_STRUCTURAL_MODE, D_LOCK, ALPHA_ARTICLE)


# =============================================================================
# Report builders
# =============================================================================

def print_header() -> None:
    mode_explanation = (
        "Switch mode by commenting one line in the ACTIVE_STRUCTURAL_MODE block.\n"
        f"Current mode: {ACTIVE_STRUCTURAL_MODE.name}\n"
        "Main lock path: Lambda_geom (independent vector representative) -> ALP -> D_C^lock -> alpha_article\n"
        "Optional QED/benchmark data never feed Lambda_geom and never modify D_C^lock."
    )
    console.print(Panel(mode_explanation, title="Relator alpha closed / regressive reproducer", expand=False))



def key_alpha_interpretation(optional_results: Optional[OptionalDiagnosticsResults]) -> str:
    base = (
        "alpha_article is the primary closed geometric output of the article branch. "
        "alpha_QED,regressed comes from regressing the same locked D_C through the optional pure-QED g2/a_e shock chain. "
        "Because that branch directly encodes the QED shock information, it is often expected to be benchmark-wise tighter or more precise numerically. "
        "That usually makes alpha_QED,regressed a likely more accurate comparison value against alpha_exp, although it remains an external diagnostic refinement rather than the primary geometric lock itself."
    )
    if optional_results is None:
        return (
            base
            + " In ARTICLE_ONLY_MODE the regressed QED alpha is intentionally suppressed and printed as ×. "
            "Uncomment ARTICLE_PLUS_QED_DIAGNOSTICS to display it numerically."
        )

    article_gap = abs(optional_results.delta_alpha_article_vs_exp)
    qed_gap = abs(optional_results.delta_alpha_qed_vs_exp)
    if qed_gap < article_gap:
        closeness = (
            "In this run the QED-regressed alpha is numerically closer to alpha_exp than alpha_article."
        )
    elif qed_gap > article_gap:
        closeness = (
            "In this run alpha_article is numerically closer to alpha_exp than the QED-regressed alpha."
        )
    else:
        closeness = (
            "In this run both alpha outputs are equally close to alpha_exp within the printed precision."
        )
    return f"{closeness} {base}"



def print_key_alpha_box(optional_results: Optional[OptionalDiagnosticsResults]) -> None:
    key_table = Table(box=box.SIMPLE_HEAVY, show_header=False, expand=False)
    key_table.add_column("Label", style="bold cyan", no_wrap=True)
    key_table.add_column("Value", justify="right")

    key_table.add_row("alpha_article", f"[bold yellow]{fmt(ALPHA_ARTICLE, KEY_BOX_DIGITS)}[/bold yellow]")
    key_table.add_row("alpha_article^-1", f"[bold yellow]{fmt(ALPHA_ARTICLE_INV, KEY_BOX_DIGITS)}[/bold yellow]")
    key_table.add_row(
        "alpha_QED,regressed",
        f"[bold green]{fmt_optional(optional_results.alpha_qed_regressed if optional_results else None, KEY_BOX_DIGITS)}[/bold green]",
    )
    key_table.add_row(
        "alpha_QED,regressed^-1",
        f"[bold green]{fmt_optional(optional_results.alpha_qed_regressed_inv if optional_results else None, KEY_BOX_DIGITS)}[/bold green]",
    )
    key_table.add_row(
        "|alpha_article - alpha_exp|",
        fmt_optional(abs(optional_results.delta_alpha_article_vs_exp) if optional_results else None, 16, sci=True),
    )
    key_table.add_row(
        "|alpha_QED,regressed - alpha_exp|",
        fmt_optional(abs(optional_results.delta_alpha_qed_vs_exp) if optional_results else None, 16, sci=True),
    )
    key_table.add_row("Interpretation", key_alpha_interpretation(optional_results))

    subtitle = (
        "alpha_QED,regressed is numeric only in ARTICLE_PLUS_QED_DIAGNOSTICS"
        if optional_results
        else "ARTICLE_ONLY_MODE active: optional QED regressed output is suppressed"
    )
    console.print(
        Panel(
            key_table,
            title="Key final alpha outputs",
            subtitle=subtitle,
            border_style="bright_magenta",
            expand=False,
        )
    )



def print_structural_mode_table() -> None:
    rows = [
        ("ACTIVE_STRUCTURAL_MODE", ACTIVE_STRUCTURAL_MODE.name, "comment/uncomment exactly one assignment line"),
        (
            "External alpha benchmark",
            "enabled" if ACTIVE_STRUCTURAL_MODE.diagnostics_enabled else "disabled",
            "when disabled, all comparison cells are printed as ×",
        ),
        (
            "Pure-QED g2 / a_e shock branch",
            "enabled" if ACTIVE_STRUCTURAL_MODE.diagnostics_enabled else "disabled",
            "when disabled, D_C^QED and alpha_QED,regressed cells are printed as ×",
        ),
        ("Placeholder symbol", placeholder(), "printed in every suppressed optional field"),
    ]
    print_table(
        "1) Structural output mode",
        (("Item", "left"), ("Value", "left"), ("Meaning", "left")),
        rows,
    )



def print_exact_constants_table() -> None:
    rows = [
        ("eta_0", fmt(ETA_0), "minimal branch ratio = 1/pi"),
        ("ell_0", fmt(ELL_0), "minimal shell width = 1/(pi sqrt(pi))"),
        ("Lambda_ind", fmt(LAMBDA_IND), "filamentary reference core"),
        ("C_UV^Gauss", fmt(C_UV_GAUSS), "Gaussian UV constant"),
        ("K", fmt(K_GEOMETRIC), "exact shell constant entering the ALP map"),
    ]
    print_table(
        "2) Exact closed constants entering the closure",
        (("Quantity", "left"), ("Value", "right"), ("Meaning", "left")),
        rows,
    )



def print_vector_representative_table() -> None:
    rep = ARTICLE_VECTOR_REPRESENTATIVE
    rows = [
        ("P_IR^(mean)(ell_0)", fmt(rep.p_ir_mean), "derived from Lambda_UV->IR^(mean) / C_UV^Gauss"),
        ("Lambda_UV->IR^(mean)(ell_0)", fmt(rep.lambda_uv_to_ir_mean), "article mean-shell diagnostic input"),
        ("Lambda_OUT(eta_0)", fmt(rep.lambda_out_article), "article shell-source diagnostic input"),
        ("R_chi^(mean)", fmt(rep.r_chi_mean), "derived from [Lambda_OUT * R_chi] / Lambda_OUT"),
        ("Lambda_OUT(eta_0) * R_chi^(mean)", fmt(rep.exterior_memory_total_mean), "representative exterior-memory contribution"),
        ("Lambda_geom^(mean)", fmt(rep.lambda_geom_mean), "independent vector representative used in the lock"),
    ]
    print_table(
        "3) Independent vector representative used to build Lambda_geom",
        (("Quantity", "left"), ("Value", "right"), ("Meaning", "left")),
        rows,
    )



def print_raw_maxwell_audit_table() -> None:
    rows = [
        ("Lambda_OUT^raw(eta_0)", fmt(RAW_VECTOR_LAMBDA_OUT), "direct Maxwell shell audit; informative only"),
        ("w_1^raw", fmt(RAW_VECTOR_WEIGHTS[1]), "normalized source weight from direct projection"),
        ("w_3^raw", fmt(RAW_VECTOR_WEIGHTS[3]), "normalized source weight from direct projection"),
        ("w_5^raw", fmt(RAW_VECTOR_WEIGHTS[5]), "normalized source weight from direct projection"),
        ("w_7^raw", fmt(RAW_VECTOR_WEIGHTS[7]), "normalized source weight from direct projection"),
        (
            "tail^raw above j=7",
            fmt(ONE - mp.fsum(RAW_VECTOR_WEIGHTS[j] for j in [1, 3, 5, 7])),
            "direct projection residual tail",
        ),
    ]
    print_table(
        "4) Direct Maxwell shell audit (not fed into Lambda_geom unless the user edits the representative block)",
        (("Quantity", "left"), ("Value", "right"), ("Meaning", "left")),
        rows,
    )



def print_lock_summary_table(optional_results: Optional[OptionalDiagnosticsResults]) -> None:
    rows = [
        ("Lambda_geom", fmt(LAMBDA_GEOM), "independent vector-shell representative; no alpha/QED inserted"),
        ("zeta_lock", fmt(ZETA_LOCK), "zeta = (K / 2pi^2) * Lambda_geom"),
        ("D_C^lock", fmt(D_LOCK), "ALP image of Lambda_geom"),
        ("alpha_article", fmt(ALPHA_ARTICLE), "article branch alpha = pi D_C^lock / sqrt(R_moth(D_C^lock))"),
        ("alpha_article^-1", fmt(ALPHA_ARTICLE_INV), "inverse article alpha"),
        ("g_{e,2}^{article}", fmt(G2_ARTICLE), "computed from D_total(D_C^lock, alpha_article)"),
        ("a_e^{article}", fmt(A_E_ARTICLE), "a_e = g_{e,2} - 1"),
        (
            "alpha_QED,regressed",
            fmt_optional(optional_results.alpha_qed_regressed if optional_results else None),
            "solve D_C^QED(alpha) = D_C^lock from the optional pure-QED g2 shock",
        ),
        (
            "alpha_QED,regressed^-1",
            fmt_optional(optional_results.alpha_qed_regressed_inv if optional_results else None),
            "inverse of the optional regressed QED alpha",
        ),
        (
            "g_{e,2}^{QED}(alpha_QED,regressed)",
            fmt_optional(optional_results.g2_qed_at_regressed_alpha if optional_results else None),
            "optional external QED branch only",
        ),
        (
            "D_C^{QED}(alpha_exp)",
            fmt_optional(optional_results.d_c_qed_at_alpha_exp if optional_results else None),
            "optional g2-shock D_C built from the external pure-QED coefficients",
        ),
        (
            "Delta alpha_article vs alpha_exp",
            fmt_optional(optional_results.delta_alpha_article_vs_exp if optional_results else None, sci=True),
            "external comparison only",
        ),
        (
            "Delta alpha_QED,regressed vs alpha_exp",
            fmt_optional(optional_results.delta_alpha_qed_vs_exp if optional_results else None, sci=True),
            "external comparison only",
        ),
        (
            "Delta alpha_QED,regressed vs alpha_article",
            fmt_optional(optional_results.delta_alpha_qed_vs_article if optional_results else None, sci=True),
            "optional article/QED separation diagnostic",
        ),
        (
            "z_sigma(article)",
            fmt_optional(optional_results.z_sigma_article if optional_results else None),
            "external sigma distance only",
        ),
        (
            "z_sigma(QED regressed)",
            fmt_optional(optional_results.z_sigma_qed if optional_results else None),
            "external sigma distance only",
        ),
    ]
    print_table(
        "5) Independent lock, article alpha, and optional QED / benchmark diagnostics",
        (("Quantity", "left"), ("Value", "right"), ("Meaning", "left")),
        rows,
    )



def print_optional_external_inputs_table(optional_results: Optional[OptionalDiagnosticsResults]) -> None:
    rows = [
        ("alpha_exp^-1", fmt_optional(optional_results.alpha_inv_exp if optional_results else None), "optional external benchmark"),
        ("sigma(alpha^-1)", fmt_optional(optional_results.alpha_inv_exp_sigma if optional_results else None), "optional external benchmark uncertainty"),
        ("alpha_exp", fmt_optional(optional_results.alpha_exp if optional_results else None), "optional external benchmark"),
        ("sigma(alpha)", fmt_optional(optional_results.alpha_exp_sigma if optional_results else None), "propagated optional uncertainty"),
        (
            "g_{e,2}^{QED}(alpha_exp)",
            fmt_optional(optional_results.g2_qed_at_alpha_exp if optional_results else None),
            "optional pure-QED g2 shock evaluated at alpha_exp",
        ),
        (
            "a_e^{QED}(alpha_exp)",
            fmt_optional(optional_results.a_e_qed_at_alpha_exp if optional_results else None),
            "optional pure-QED anomalous moment shock",
        ),
        (
            "D_total^{QED}(alpha_exp)",
            fmt_optional(optional_results.d_total_qed_at_alpha_exp if optional_results else None),
            "optional D_total reconstructed from g2 shock",
        ),
        (
            "D_C^{QED}(alpha_exp)",
            fmt_optional(optional_results.d_c_qed_at_alpha_exp if optional_results else None),
            "optional regressed D_C from the QED g2 shock chain",
        ),
    ]
    print_table(
        "6) Optional external benchmark and g2-shock inputs",
        (("Quantity", "left"), ("Value", "right"), ("Meaning", "left")),
        rows,
    )



def coefficient_rows(
    d_series: List[mp.mpf],
    d_total_series: List[mp.mpf],
    g2_series: List[mp.mpf],
    a_series: List[mp.mpf],
    optional_results: Optional[OptionalDiagnosticsResults],
) -> List[Tuple[str, str, str, str, str, str, str]]:
    rows: List[Tuple[str, str, str, str, str, str, str]] = []
    for n in range(1, PREDICTION_ORDER + 1):
        qed_coeff = None
        delta = None
        if optional_results and n < len(optional_results.a_qed_series):
            qed_coeff = optional_results.a_qed_series[n]
            delta = a_series[n] - qed_coeff
        rows.append(
            (
                str(n),
                fmt(d_series[n], 22),
                fmt(d_total_series[n], 22),
                fmt(g2_series[n], 22),
                fmt(a_series[n], 22),
                fmt_optional(qed_coeff, 22),
                fmt_optional(delta, 12, sci=True),
            )
        )
    return rows



def print_coefficient_tables(optional_results: Optional[OptionalDiagnosticsResults]) -> None:
    columns = (
        ("n", "right"),
        ("D_C coeff", "right"),
        ("D_total coeff", "right"),
        ("g2 coeff", "right"),
        ("a_e coeff", "right"),
        ("pure-QED a_e coeff", "right"),
        ("Delta a_e", "right"),
    )
    print_table(
        "7) Current article rank-5 closed coefficient ledger",
        columns,
        coefficient_rows(CURRENT_D_SERIES, CURRENT_D_TOTAL_SERIES, CURRENT_G2_SERIES, CURRENT_A_SERIES, optional_results),
    )
    print_table(
        "8) Refined rank-100 closed coefficient ledger",
        columns,
        coefficient_rows(REFINED_D_SERIES, REFINED_D_TOTAL_SERIES, REFINED_G2_SERIES, REFINED_A_SERIES, optional_results),
    )



def print_optional_qed_dc_series_table(optional_results: Optional[OptionalDiagnosticsResults]) -> None:
    rows = []
    for n in range(1, 6):
        coeff = optional_results.dc_qed_series[n] if optional_results else None
        rows.append((str(n), fmt_optional(coeff, 24), "coefficient c_n in D_C^QED(x) = sum c_n x^n"))
    print_table(
        "9) Optional regressed D_C^QED(x) series from the pure-QED g2 shock",
        (("n", "right"), ("c_n^{QED}", "right"), ("Meaning", "left")),
        rows,
    )



def print_audit_summary(optional_results: Optional[OptionalDiagnosticsResults]) -> None:
    rows = [
        (
            "Lock path actually used",
            "Lambda_geom^(mean) -> zeta_lock -> D_C^lock -> alpha_article",
            "independent of alpha_exp and independent of pure-QED coefficients",
        ),
        (
            "Optional branch status",
            "enabled" if optional_results else "disabled",
            "all optional cells are printed as × when disabled",
        ),
        (
            "g2 coefficients",
            "printed in tables 7 and 8",
            "same expansion order as the a_e ledger",
        ),
        (
            "Key alpha box",
            "printed near the top of the report",
            "alpha_QED,regressed is numeric only when ARTICLE_PLUS_QED_DIAGNOSTICS is active",
        ),
        (
            "QED regressed alpha row",
            "printed numerically" if optional_results else "printed as ×",
            "controlled only by ACTIVE_STRUCTURAL_MODE",
        ),
    ]
    print_table(
        "10) Audit summary",
        (("Category", "left"), ("Value", "left"), ("Meaning", "left")),
        rows,
    )


# =============================================================================
# Main report
# =============================================================================

def main() -> None:
    print_header()
    print_key_alpha_box(OPTIONAL_RESULTS)
    print_structural_mode_table()
    print_exact_constants_table()
    print_vector_representative_table()
    print_raw_maxwell_audit_table()
    print_lock_summary_table(OPTIONAL_RESULTS)
    print_optional_external_inputs_table(OPTIONAL_RESULTS)
    print_coefficient_tables(OPTIONAL_RESULTS)
    print_optional_qed_dc_series_table(OPTIONAL_RESULTS)
    print_audit_summary(OPTIONAL_RESULTS)


if __name__ == "__main__":
    main()
