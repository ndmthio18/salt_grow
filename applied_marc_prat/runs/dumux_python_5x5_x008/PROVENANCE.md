# Provenance

## Geometry

- Source host: `liqin`
- Source path: `/home/zhixin/Documents/Simulation_Data_for_paper/Postprocess/Grids/N5x5_R270um_SD10um_L3.5mmx3.5mm_n1_top_varthroat.dgf`
- Archived size: 3470 bytes
- SHA-256: `E48FF929AEAE9149DF9B575FD2BD612C7C19E111810D023DBF7013EE5341B825`

## DuMux

- Isolated worktree: `/home/zhixin/Documents/Humboldt/salt/_work_dumux_5x5_x008`
- Branch: `codex/dumux-5x5-x008-validation`
- `pnm2pncmin2` baseline: `6dd6d62ae80a0d78d8bae0f652b4d06438427fea`
- Validation-only Float64 state/signed-flux diagnostics:
  `ab303806040028cdf9d137dd5b67238b436f2a6f`; this changes only CSV output in
  `appl/2pncmin_isothermal/main.cc`, not residuals, constitutive laws or boundary conditions.
- `dumux`: `a4446704f440faee12d770ff7a18a0d48bd544aa`
- `dune-common`: `cce9c197680c1a9964d21ab4d0e0835ea36affb0`
- `dune-geometry`: `7d5b1d81ad997f81637ac97f753f80a64ff9cdb0`
- `dune-grid`: `75b66b0ebf0656e21af08798188b3d2848c9574d`
- `dune-istl`: `cfa1859941bac063913c5b94e238cd4a388235b9`
- `dune-localfunctions`: `f2c7cfb96327fbfd29744dccf5eac015a1dfa06f`
- `dune-uggrid`: `e26f81ff7d84f5d7b228edb3313beae592d502f7`
- `dune-foamgrid`: `2bf9fe68e6ccaa389adbdbdeff41557b99ed0375`
- Compiler: GCC 13.3.0; flags include `-O3 -DNDEBUG` (the upstream target also adds `-g`).
- Formal-run optimized executable SHA-256:
  `A064BC5549C6781B246FF5C6191CB2FEF49AE69323511384D287F5EB64C5C5F0`.
- Diagnostic-build optimized executable SHA-256:
  `94865A37D73E04D61D1CB845EB9DBE752F79D84E1E21E3C38F5AD729C34CAF3B`.
- Remote run root:
  `/home/zhixin/Documents/Humboldt/salt/_runs/dumux_python_5x5_x008_20261112`.

The original debug-flag runs were retained remotely as `dumux_dt10_debug_partial` and
`dumux_dt5_debug_partial`.  Before switching executables, six primary fields at 0, 1000, 2000,
3000 and 4000 s were identical at the precision written by DuMux.

Formal runs used `DUMUX_NUM_THREADS=1` and `OMP_NUM_THREADS=1`:

- `dumux_dt5_full`: last converged step at 27340 s, then singular UMFPACK matrix, repeated
  Newton time-step division to `9.79e-100 s`, and abort.  Raw archive SHA-256:
  `C81C0D1B774550420CB868AB5EEDBBFFE5FF9B0A7071B8B411CF16771CF445FB`.
- `dumux_dt10_full`: physical time stalled at 27374 s while accepted steps shrank monotonically
  to `3.08e-173 s`; the non-progressing process was terminated.  Raw archive SHA-256:
  `F2889C9FF3DF2D0D9D2F8056CEAF1BEFA6B9DC8C6EFCF4686030DCA9EB5C369A`.
- `dumux_sat_diag_dt5`: 5 s checkpoints to 1300 s; first non-zero solid at 940 s and interpolated
  solubility crossing at 936.228 s.  Raw archive SHA-256:
  `E0B82A4992ABE9E34C8B14FCE039280FEBCD283C20740FD436CC84DF9DCE43D`.

The `dumux_dt10` hash above is also recorded by the local archive tool; the raw run log is the
authoritative record of the user-terminated non-progressing state.

The diagnostic build was exercised separately at 100 s and at 0–1000 s with 5 s and 1 s
checkpoints.  Each checkpoint contains exactly 30 Float64 pore rows and 45 Float64 throat rows.
At 100 s, saturation, mole fraction, solid fraction, pore volume and flux magnitude agree with
the baseline formal run to its Float32 write precision; the new pressure diagnostic resolves a
`4.48e-6 Pa` network range that the Float32 absolute-pressure field rounds away.  The 1 s run
also writes equation-unit water/NaCl component fluxes and direct component inventories for an
independent ledger reconstruction.

## Python

- Python 3.12.4; NumPy 2.0.0; SciPy 1.14.0; Matplotlib 3.9.2.
- The working tree contains the new multi-reservoir DGF reader and dedicated comparison driver.
- Python state uses NaCl mass fraction internally and converts all reported profiles to the
  DuMux NaCl mole fraction.
- Local resource audit before the three Python reruns observed 22 logical CPU scheduling units;
  memory availability was unknown.  The three 30-pore runs were executed concurrently and are
  small enough that memory was not a limiting resource.
- Production Python states: `dt=10, 5, 2.5 s`, all to 40000 s; DuMux-flux replay: `dt=5 s` to
  the last unambiguously pre-saturation common checkpoint, 900 s.

