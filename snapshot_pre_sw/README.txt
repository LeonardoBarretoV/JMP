SNAPSHOT: pre-sticky-wage reference version
Created: 2026-08-25

WHAT THIS IS
------------
A verbatim physical copy of src/ and scripts/ as they stood immediately before
any sticky-wage / endogenous-inflation work began.

This snapshot INCLUDES the uncommitted [Z-NORM] edits, which are intentional and
are part of the reference. They normalize aggregate efficiency
sum(z_state[1,:] .* z_grid) to exactly 1.0 in the INITIAL STEADY STATE
(it was 0.798 before them). The normalization is imposed at the initial steady
state ONLY: sum(z_state[i,:] .* z_grid) still falls ~5.5% on impact and remains
~0.45% below 1.0 at t=40. That time variation is correct and is preserved.

The [Z-NORM] edits touch three files and are tagged on every added line:
  scripts/main_JMP.jl   (+5 / -4)
  src/calibration.jl    (+10 / -3)
  src/transition.jl     (+3 / -3)

WHY IT EXISTS
-------------
Nothing is under version control yet. This directory is the reference for
reproducing the current paper's numbers at any time, including after the
sticky-wage track is finished.

RULES
-----
- Do NOT add this directory to any Julia include path.
- Do NOT edit anything in here.
- Treat it as read-only.

HOW TO USE IT
-------------
To reproduce the current results, run snapshot_pre_sw/scripts/main_JMP.jl
(adjusting ROOT, which is computed relative to the script location).

Companion reproducible artifact, NOT in this directory and NOT in version
control: output/jld2/baseline.jld2 (~820 MB), the solved result of the live
scripts/main_JMP.jl. The sticky-wage track writes elsewhere and never touches it.
