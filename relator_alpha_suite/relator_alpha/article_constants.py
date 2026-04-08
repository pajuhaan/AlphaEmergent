from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable

import mpmath as mp

from .common import configure_precision

configure_precision()


# -----------------------------------------------------------------------------
# Core exact or article-fixed constants
# -----------------------------------------------------------------------------

PI = mp.pi
EULER_GAMMA = mp.euler

ETA_0 = 1 / PI
ELL_0 = 1 / (PI * mp.sqrt(PI))

# Eq. (app-closed-geom)
K_ALP = (150 * PI**2 - 8 * PI**4 - 315) / (180 * PI**6)

# Scalar-side exact geometric constants
C_UNI = (1 / PI) * (mp.mpf("4") / 3 + 1 / (4 * PI**2))
THETA_1_CORE = 2 * PI * C_UNI  # = 8/3 + 1/(2π²)
THETA_1_LOG = ETA_0**2 / (8 * (1 - ETA_0**2 / 4))

# Article-fixed closed shell value (Table "Fixed shell geometry and static hidden-mode data")
THETA_1_COLL = mp.mpf("0.01607850686959103319")

# Visible dynamic scales
S_UV = mp.log(2)
S_IR = 1 / (8 * PI**2)

# Vector-channel exact or article-fixed constants
LAMBDA_IND = mp.log(8 * mp.sqrt(PI)) - 2
C_UV_GAUSS = mp.mpf("0.5") * (mp.log(2) + EULER_GAMMA)

DELTA_LAMBDA_OUT_EXACT = mp.mpf("-0.0139671580625860254655225000000")
DELTA_LAMBDA_UV_TO_IR_MEAN = mp.mpf("0.0544853495865520906532491466693")
EXTERIOR_MEMORY_CONTRIBUTION_MEAN = mp.mpf("-0.0146086880840972703041996518997")
LAMBDA_GEOM_MEAN = mp.mpf("0.6916831461069908356724595349")


# -----------------------------------------------------------------------------
# Closed auxiliary formulas
# -----------------------------------------------------------------------------


def gamma_cl(d_value: mp.mpf) -> mp.mpf:
    """Exact collar-log overlap in the closed erfc representation."""
    d_value = mp.mpf(d_value)
    return (
        mp.sqrt(PI)
        * ELL_0
        / 2
        * mp.exp((ELL_0**2 / 4) * (d_value + mp.mpf("0.25")))
        * mp.erfc((ELL_0 / 2) * mp.sqrt(d_value + mp.mpf("0.25")))
    )


GAMMA_CL_0 = gamma_cl(mp.mpf("0"))

# Eq. (B11-glue), stored with the article's extended digits.
B11_GLUE = mp.mpf("0.02904120666648408568405788097756078")



def compute_b11_glue_from_formula() -> mp.mpf:
    """Reconstruct B11^glue from the closed visible two-channel formula."""
    numerator = THETA_1_COLL + THETA_1_LOG - GAMMA_CL_0 * THETA_1_COLL * THETA_1_LOG
    denominator = 1 - (GAMMA_CL_0**2 / 4) * THETA_1_COLL * THETA_1_LOG
    return numerator / denominator


# Canonical three-channel static completion.
N_C = mp.sqrt(PI / 8) * ELL_0
RHO_STAT = GAMMA_CL_0 / mp.sqrt(N_C)
LAMBDA_H = 1 / (1 - RHO_STAT**2)



def compute_b11_3rho_from_formula() -> mp.mpf:
    """Reconstruct B11^(3,ρ) from the canonical no-free-parameter completion."""
    a_value = THETA_1_COLL
    b_value = THETA_1_LOG
    c0 = GAMMA_CL_0 * a_value * b_value / 2
    a_perp = a_value * (1 - GAMMA_CL_0**2 * a_value * b_value / 4)
    b_perp = b_value * (1 - GAMMA_CL_0**2 * a_value * b_value / 4)
    t_perp = mp.matrix(
        [
            mp.sqrt(a_perp / a_value),
            mp.sqrt(b_perp / b_value),
            mp.mpf("0"),
        ]
    )
    k_stat = mp.matrix(
        [
            [1 / a_perp, GAMMA_CL_0 / 2, mp.sqrt(c0)],
            [GAMMA_CL_0 / 2, 1 / b_perp, mp.sqrt(c0)],
            [mp.sqrt(c0), mp.sqrt(c0), LAMBDA_H],
        ]
    )
    return (t_perp.T * k_stat**-1 * t_perp)[0]


B11_3RHO = mp.mpf("0.02904114778429781136756907742390621")
THETA_1_3RHO = THETA_1_CORE + B11_3RHO

# Representative unresolved operator inputs.
R_CHI_MEAN = EXTERIOR_MEMORY_CONTRIBUTION_MEAN / DELTA_LAMBDA_OUT_EXACT
P_IR_CHI_MEAN = DELTA_LAMBDA_UV_TO_IR_MEAN / C_UV_GAUSS


# -----------------------------------------------------------------------------
# Structured scalar-shell data extracted from the article's finite-rank tables
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class ScalarRankData:
    """Finite-rank scalar-shell data reported in the manuscript tables."""

    rank_n: int
    delta_a_mix: mp.mpf
    delta_a0: mp.mpf
    chi_cross: mp.mpf
    a_uv: mp.mpf
    a_ir: mp.mpf
    chi: mp.mpf
    rho_dyn_published: mp.mpf
    hidden_source_norm: mp.mpf | None = None
    jacobi_alpha0: mp.mpf | None = None
    jacobi_beta1: mp.mpf | None = None
    jacobi_alpha1: mp.mpf | None = None
    jacobi_beta2: mp.mpf | None = None
    jacobi_alpha2: mp.mpf | None = None

    @property
    def rho_dyn(self) -> mp.mpf:
        return self.rho_dyn_published

    @property
    def lambda_u_display(self) -> mp.mpf:
        # This is the λ_u shown in Table tab:dc-visible-dynamic-block.
        return 1 / (1 - self.rho_dyn**2)


SCALAR_RANK_DATA: Dict[int, ScalarRankData] = {
    1: ScalarRankData(
        rank_n=1,
        delta_a_mix=mp.mpf("0"),
        delta_a0=mp.mpf("0"),
        chi_cross=mp.mpf("0"),
        a_uv=mp.mpf("2.85766191148700367940"),
        a_ir=mp.mpf("2.85766191148700367940"),
        chi=mp.mpf("0"),
        rho_dyn_published=mp.mpf("0"),
    ),
    3: ScalarRankData(
        rank_n=3,
        delta_a_mix=mp.mpf("5.6297544557207391e-3"),
        delta_a0=mp.mpf("-3.9660140058889997e-4"),
        chi_cross=mp.mpf("3.1694135231707515e-5"),
        a_uv=mp.mpf("2.86293524866547739480"),
        a_ir=mp.mpf("2.85167573975403591666"),
        chi=mp.mpf("3.1694135231707515e-5"),
        rho_dyn_published=mp.mpf("3.3826809862918462e-4"),
    ),
    5: ScalarRankData(
        rank_n=5,
        delta_a_mix=mp.mpf("6.4229315823120090e-3"),
        delta_a0=mp.mpf("-4.0484921668976660e-4"),
        chi_cross=mp.mpf("4.1254050111061048e-5"),
        a_uv=mp.mpf("2.86372101365445760227"),
        a_ir=mp.mpf("2.85087515048983358419"),
        chi=mp.mpf("4.1254050111061048e-5"),
        rho_dyn_published=mp.mpf("4.4030004257257359e-4"),
    ),
    100: ScalarRankData(
        rank_n=100,
        delta_a_mix=mp.mpf("6.4673176500784765e-3"),
        delta_a0=mp.mpf("-4.0504733367299593e-4"),
        chi_cross=mp.mpf("4.1826197587016588e-5"),
        a_uv=mp.mpf("2.86376522167868808086"),
        a_ir=mp.mpf("2.85083058637853112785"),
        chi=mp.mpf("4.1826197587016588e-5"),
        rho_dyn_published=mp.mpf("4.4640651108518807e-4"),
        hidden_source_norm=mp.mpf("0.11058235357030053635"),
        jacobi_alpha0=mp.mpf("0.11874476222472424339"),
        jacobi_beta1=mp.mpf("0.01434871678338358206"),
        jacobi_alpha1=mp.mpf("0.00937286363486739350"),
        jacobi_beta2=mp.mpf("0.00670541329129441449"),
        jacobi_alpha2=mp.mpf("0.01271810767253984155"),
    ),
}


@dataclass(frozen=True)
class ScalarModelPreset:
    """Load-bearing scalar evaluator preset for a specific article table family."""

    key: str
    description: str
    family: str
    rank_data: ScalarRankData
    b11_effective: mp.mpf
    source_renorm_factor: mp.mpf
    source_renorm_label: str

    @property
    def rank_n(self) -> int:
        return self.rank_data.rank_n

    @property
    def a_uv(self) -> mp.mpf:
        return self.rank_data.a_uv

    @property
    def a_ir(self) -> mp.mpf:
        return self.rank_data.a_ir

    @property
    def chi(self) -> mp.mpf:
        return self.rank_data.chi

    @property
    def rho_dyn(self) -> mp.mpf:
        return self.rank_data.rho_dyn

    @property
    def lambda_u_display(self) -> mp.mpf:
        return self.rank_data.lambda_u_display

    @property
    def theta1_effective(self) -> mp.mpf:
        return THETA_1_CORE + self.b11_effective



def build_current_visible_model(rank_n: int) -> ScalarModelPreset:
    rank = SCALAR_RANK_DATA[rank_n]
    return ScalarModelPreset(
        key=f"current_rank{rank_n}",
        description="Current article-consistent glued-shell scalar evaluator",
        family="current_visible",
        rank_data=rank,
        b11_effective=B11_GLUE,
        source_renorm_factor=mp.mpf("1"),
        source_renorm_label="1",
    )



def build_refined_physical_model(rank_n: int) -> ScalarModelPreset:
    rank = SCALAR_RANK_DATA[rank_n]
    lambda_u = rank.lambda_u_display
    return ScalarModelPreset(
        key=f"refined_rank{rank_n}",
        description=(
            "Refined no-free-parameter mother truncation "
            "(main-text physical/g2 convention)"
        ),
        family="refined_physical",
        rank_data=rank,
        b11_effective=B11_3RHO,
        source_renorm_factor=lambda_u**2,
        source_renorm_label="(λ_u^(main text))²",
    )



def build_realized_lock_model(rank_n: int) -> ScalarModelPreset:
    rank = SCALAR_RANK_DATA[rank_n]
    lambda_u = rank.lambda_u_display
    return ScalarModelPreset(
        key=f"realized_rank{rank_n}",
        description=(
            "Realized appendix lock convention "
            "(source prefactor = λ_u^(main text))"
        ),
        family="realized_lock",
        rank_data=rank,
        b11_effective=B11_3RHO,
        source_renorm_factor=lambda_u,
        source_renorm_label="λ_u^(main text)",
    )


CURRENT_MODELS = {rank: build_current_visible_model(rank) for rank in SCALAR_RANK_DATA}
REFINED_MODELS = {rank: build_refined_physical_model(rank) for rank in SCALAR_RANK_DATA}
REALIZED_MODELS = {rank: build_realized_lock_model(rank) for rank in SCALAR_RANK_DATA}

SCALAR_MODEL_PRESETS: Dict[str, ScalarModelPreset] = {
    model.key: model
    for model in [
        *CURRENT_MODELS.values(),
        *REFINED_MODELS.values(),
        *REALIZED_MODELS.values(),
    ]
}


# -----------------------------------------------------------------------------
# Vector-channel shell-source weights (article-fixed operator outputs)
# -----------------------------------------------------------------------------

VECTOR_SHELL_SOURCE_WEIGHTS = {
    1: mp.mpf("0.9796363359550031216197956500"),
    3: mp.mpf("0.02031062626884619091016822831"),
    5: mp.mpf("4.849145355825759038791211476e-5"),
    7: mp.mpf("4.354340825127690409129907468e-6"),
}


# -----------------------------------------------------------------------------
# External orientation and benchmark data quoted in the article
# -----------------------------------------------------------------------------

ALPHA_BENCHMARK_INVERSE = mp.mpf("137.035999177")
ALPHA_BENCHMARK_INVERSE_U = mp.mpf("2.1e-8")
ALPHA_BENCHMARK_DIRECT = 1 / ALPHA_BENCHMARK_INVERSE
ALPHA_BENCHMARK_DIRECT_U = ALPHA_BENCHMARK_INVERSE_U / ALPHA_BENCHMARK_INVERSE**2

ALPHA_CODATA22_INVERSE = mp.mpf("137.035999177")
ALPHA_CODATA22_DIRECT = 1 / ALPHA_CODATA22_INVERSE
X_REFERENCE = ALPHA_CODATA22_DIRECT / PI

QED_G_E2_BENCHMARK_XSTAR = mp.mpf("1.00115965217597119498265100941")
QED_A_E_BENCHMARK_XSTAR = mp.mpf("1.1596521759711949826510094063e-3")


@dataclass(frozen=True)
class ExperimentalAlphaReference:
    year: str
    method: str
    reference_label: str
    alpha_ref: mp.mpf
    u_ref: mp.mpf


EXPERIMENTAL_ALPHA_REFERENCES = [
    ExperimentalAlphaReference(
        year="2006",
        method="Rb recoil (Bloch)",
        reference_label="Clade2006",
        alpha_ref=mp.mpf("0.00729735259"),
        u_ref=mp.mpf("4.8e-11"),
    ),
    ExperimentalAlphaReference(
        year="2008",
        method="Rb recoil + Ramsey–Bordé",
        reference_label="Cadoret2008",
        alpha_ref=mp.mpf("0.00729735255"),
        u_ref=mp.mpf("3.3e-11"),
    ),
    ExperimentalAlphaReference(
        year="2011",
        method="Rb recoil",
        reference_label="Bouchendira2011",
        alpha_ref=mp.mpf("0.007297352572"),
        u_ref=mp.mpf("4.85e-12"),
    ),
    ExperimentalAlphaReference(
        year="2018",
        method="Cs recoil",
        reference_label="Parker2018",
        alpha_ref=mp.mpf("0.007297352571"),
        u_ref=mp.mpf("1.44e-12"),
    ),
    ExperimentalAlphaReference(
        year="2020",
        method="Rb recoil",
        reference_label="Morel2020",
        alpha_ref=mp.mpf("0.0072973525628"),
        u_ref=mp.mpf("5.86e-13"),
    ),
    ExperimentalAlphaReference(
        year="2023",
        method="e⁻ moment + QED (derived)",
        reference_label="Fan2023",
        alpha_ref=mp.mpf("0.0072973525649"),
        u_ref=mp.mpf("8.00e-13"),
    ),
    ExperimentalAlphaReference(
        year="2025",
        method="CODATA (2022 set, publ. 2025)",
        reference_label="CODATA2022",
        alpha_ref=mp.mpf("0.007297352564331425"),
        u_ref=mp.mpf("1.11828e-12"),
    ),
]


# -----------------------------------------------------------------------------
# Universal pure-photonic QED coefficients A_1^(2n)
# -----------------------------------------------------------------------------

QED_A1_COEFFICIENTS = {
    1: mp.mpf("0.5"),
    2: mp.mpf("-0.328478965579193"),
    3: mp.mpf("1.181241456587"),
    4: mp.mpf("-1.9122457649264456"),
    5: mp.mpf("5.891"),
}


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------


def iter_rank_data(ranks: Iterable[int] = (1, 3, 5, 100)) -> list[ScalarRankData]:
    return [SCALAR_RANK_DATA[rank] for rank in ranks]
