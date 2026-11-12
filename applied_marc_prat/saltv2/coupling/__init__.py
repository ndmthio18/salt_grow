"""Where the pore network and the voxel domain meet.

    interface    PoreOpenings: footprints of the surface pores on the voxel grid, the Dirichlet
                 patches of the openings, per-pore sums of the base flux.  No physics.
    ledger       named water / salt mass channels and registered balances
    driver3d     the 3-D closed loop: gas -> supply -> network transport -> growth -> sealing
    driver2d     the 2-D closed loop the 3-D closures were ported from

The drivers are imported as modules (``from saltv2.coupling import driver3d``); nothing is pulled
in at package import, so ``import saltv2`` stays light.
"""
