#!/usr/bin/env bash
# Regenerate deriv/rendered/: run the four derivations with rendering on, then build the PDF.
# The scripts also rewrite their src/Matrix/*_generated.jl (byte-identical when nothing changed)
# and stop with exit 1 if a check fails.  Takes about 5 minutes (derive_two_way.wls).
set -euo pipefail
cd "$(dirname "$0")/.."
export FLUIDUM_DERIV_RENDER=1
for s in derive_2p1d_viscous derive_hq_2p1d derive_hq_consistent_2p1d derive_two_way; do
  echo "== $s"
  wolframscript -file "deriv/$s.wls"
done
cd deriv/rendered
latexmk -pdf -interaction=nonstopmode -halt-on-error derivations.tex
