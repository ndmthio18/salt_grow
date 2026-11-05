# GO / NO-GO assessment

## Decision

**NO-GO for claiming that evolving porosity invalidates Prat's analytical law,
or for launching a broad dynamic two-phase pore-network campaign on the basis
of this prototype alone.**

**Conditional GO for a narrowly targeted crystal/pore-scale follow-up**, aimed
at testing the reaction-layer scale and the low-`k_r` hydraulic-disconnection
regime.

## Evidence for the decision

1. The analytical law passes independent numerical verification under matching
   assumptions.
2. The conservative dynamic baseline differs from Model A by only about 3% in
   velocity. This is smaller than plausible uncertainty in porosity, specific
   surface area, `beta`, reaction kinetics, and brine properties.
3. The constitutive update itself is negligible at baseline: B2 versus B1
   changes terminal velocity by -0.08% and height by -0.04%.
4. The modeled active layer is 1.24 micrometres for 95% precipitation, smaller
   than the estimated 3.5 micrometre crystal scale. A continuum closure cannot
   be treated as physically validated in that regime.
5. Low-`k_r` cases reach `epsilon = 0.02` rapidly while evaporation remains
   imposed. This identifies a regime where permeability and capillary supply
   must enter before continued simulation is meaningful.

## What would justify the next model

A dynamic two-phase pore-network model becomes scientifically worthwhile if it
is designed to answer one of two falsifiable questions:

- Does discrete crystal-scale transport replace the continuum exponential
  reaction layer and alter the effective interface coefficient by more than
  experimental/parameter uncertainty?
- Does pore blocking cause evaporation or capillary supply to fall before the
  continuum model reaches its porosity limit in the low-`k_r` regime?

Before that investment, obtain or bound crystal size/morphology, evolving
specific surface area, local permeability/connectivity, and the interface
upscaling factor `beta`. A pore-network model that merely embeds the same B2
closures would add complexity without resolving the decisive uncertainty.

## Recommended next step

Run a small crystal-resolved or micro-CT-informed unit-cell study of the upper
5-10 micrometres, coupled to a one-dimensional liquid-supply resistance. Use it
to test the interfacial mass balance and determine whether `beta`, `a_v`, and
permeability evolve independently. Only promote to a full two-phase network if
that study predicts a departure larger than the current approximately 3%
continuum effect or reproduces an independently observed hydraulic transition.
