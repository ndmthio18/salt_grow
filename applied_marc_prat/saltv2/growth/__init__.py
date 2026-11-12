"""The voxel domain of external porous-salt growth.  Regime: evaporation-wicking.

    grid         voxel states and the periodic lattice
    geometry     occupancy growth with a prescribed normal velocity (no physics)
    regularize   the explicit regularisation length of the growth law
    morphology   cluster / spanning / height metrics of a salt mask
    fv, vapor    steady gas diffusion on a fixed geometry; the sealed crust joins the domain as
                 a porous conductor
    crust        per-voxel salt properties fixed at birth, the reaction band, densification and
                 the WET -> SEALED transition
    supply       the quasi-1-D hierarchical supply check that gates growth (the same stack as
                 crust.BirthClosure)
    prat         the Prat (2026) growth law, from the validated prat2026/core.py

Imports ``saltv2.materials`` only; it knows nothing about the pore network.
"""
