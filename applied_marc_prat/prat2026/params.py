"""Base parameters for the Prat (2026) porous salt structure growth law.

Reference
---------
M. Prat, "The Evaporation-Driven Growth Law of Porous Salt Structures",
Langmuir 2026, DOI 10.1021/acs.langmuir.6c02211.

Every value below carries its provenance. Values marked BACK-CALIBRATED are NOT
stated in the paper; they were recovered by requiring the paper's own quoted
results to be reproduced, and any result depending on them must say so.
"""

# --- solution / crystal properties ------------------------------------------
C_EQ = 0.264
"""NaCl solubility mass fraction at normal pressure. Paper, 'Materials and Methods'."""

RHO_CR = 2170.0
"""Halite crystal density [kg/m^3]. Paper, below eq (2)."""

RHO = 1200.0
"""Saturated NaCl solution density [kg/m^3].

BACK-CALIBRATED -- the paper never states it, yet eq (32) needs it.
Recovered by matching two of the paper's own quoted numbers:
  * d_sub = 25 um  ->  (dh/dt)/U_evap = 0.254   (paper text: 0.25 cm/day at
    U_evap = 1 cm/day)
  * d_sub = 10 um  ->  (dh/dt)/U_evap = 0.219   (paper Figure 6 reads ~0.21-0.22)
1200 kg/m^3 is also the textbook density of saturated NaCl brine, so the
calibration is physically consistent rather than merely a fit.
See tests/test_anchors.py::test_growth_rate_anchor_25um -- that test is the sole
support for this value; rerun it after any change here.
"""

# --- salt structure stereology ----------------------------------------------
BETA = 2.0
"""Geometric up-scaling factor in eq (20). Paper text + Appendix (cubic crystal
arrangement of Figure 3)."""

EPS_S = 0.43
"""Salt structure porosity, base value. Paper: measured for efflorescence by
Dai et al. (ref 19); the paper writes 'eps_s ~ 0.43'."""

DC_OVER_DSUB = 0.14
"""Crystal size / substrate particle size ratio. Paper, from Wijnhorst et al.
(ref 20), consistent with SEM images in Norouzi Rad et al. (ref 18).
Used only to derive the 36 in eq (17): a_v = 5(1-eps)/d_c = (5/0.14)(1-eps)/d_sub
                                            ~= 36(1-eps)/d_sub."""

# --- kinetics and diffusion --------------------------------------------------
KR_T21 = 2.3e-3
"""NaCl precipitation reaction coefficient at ~21 C [m/s]. Paper, from the
microfluidic measurement of Naillon et al. (ref 23)."""

DS_T21 = 1.3e-9
"""NaCl molecular diffusion coefficient in water at ~21 C [m^2/s]. Paper text."""

DS_T45 = 2.5e-9
"""NaCl molecular diffusion coefficient at 45 C [m^2/s]. Paper, from Fell &
Hutchison (ref 42), used for the surface-heave experiment."""

EA = 49.1e3
"""Arrhenius activation energy for k_r [J/mol]. Paper, from Linnikov (ref 43)."""

R_GAS = 8.314
"""Universal gas constant [J/(mol K)]."""

T21_K = 294.15
"""Reference temperature for KR_T21 and DS_T21 [K] (21 C)."""

T45_K = 318.15
"""Temperature of the Diouf surface-heave experiment [K] (45 C)."""

# --- unit helpers -------------------------------------------------------------
DAY = 86400.0
"""Seconds per day -- the paper reports velocities in cm/day and mm/day."""


# =============================================================================
# Experimental scenarios used for the comparisons in Figures 8 and 9
# =============================================================================


class NorouziRad:
    """Sand column of Norouzi Rad, Shokri & Sahimi, Phys. Rev. E 88, 032404
    (2013) -- ref 18. Compared in paper Figure 8.

    The paper reads U_evap off their Figure 2 (6 Molal) and compares against
    their Figure 4 (6 Molal).
    """

    U_EVAP = 7.7e-3 / DAY
    """Evaporation velocity [m/s] = 7.7 mm/day. Paper text."""

    D_SUB = 0.58e-3
    """Average sand particle size [m]. Paper text (range 0.17-1 mm)."""

    ETA = 0.76
    """Surface coverage fraction giving best agreement in paper Figure 8a.
    Consistent with Norouzi Rad & Shokri (ref 38) for similar sand."""

    EPS_S = EPS_S
    KR = KR_T21
    DS = DS_T21

    DA_STAR_PAPER = 156.0
    """Da* quoted in the paper text for this scenario. NOTE: our formula gives
    ~170 -- see the Da* discrepancy discussion in the top-level README."""

    AREA = None
    """Cylinder cross-section [m^2] -- NOT given by Prat, lives in ref 18 which
    is paywalled. Recovered instead as the product eta*A from the Figure 8a/8b
    digitisations; see fig08_norouzi_rad/run.py."""


class DioufHeave:
    """Surface-heave experiment of B. Diouf, PhD thesis, Univ. Toulouse (2018)
    -- ref 41. Compared in paper Figure 9.

    Hydrophilic glass beads saturated with C=0.25 NaCl solution below a dry
    hydrophobic bead packing; subflorescence growth pushes the top surface up.
    """

    U_EVAP = 5.1e-3 / DAY
    """Evaporation velocity [m/s] = 5.1 mm/day, from a linear fit to the
    evaporated mass of paper Figure 9a."""

    D_SUB = 90e-6
    """Mean hydrophilic bead diameter [m] (paper takes the mean of 70-110 um)."""

    CYLINDER_DIAMETER = 7e-3
    """Vessel inner diameter [m]. Paper text (not needed for Figure 9b, which
    plots a height, but recorded for completeness)."""

    EPS_S_BASE = 0.43
    """Base-value porosity -> paper quotes Da* = 57."""

    EPS_S_FITTED = 0.34
    """Porosity giving very good agreement in paper Figure 9b -> Da* = 92."""

    DS = DS_T45

    DA_STAR_PAPER_BASE = 57.0
    DA_STAR_PAPER_FITTED = 92.0
    """Da* values quoted in the paper text. Our formula gives ~62.5 and ~103."""
