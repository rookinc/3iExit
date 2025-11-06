# 3iExit
**Sun–Earth–3I Minimal-Area Coupling — Analysis Framework**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17546114.svg)](https://doi.org/10.5281/zenodo.17546114)


This repository contains a preregistered, falsifiable workflow to test whether
heliospheric amplification signatures appear **only** within ±24–72 h of
**Sun–Earth–3I minimal-area** epochs, and **only** under moderate–high solar-wind forcing.

## Quick start
```bash
# after you create the empty repo on GitHub:
git clone <YOUR-REPO-URL>.git 3iExit
cd 3iExit

# add files, commit, push
# (if you downloaded the scaffold zip, unzip into this folder first)
git add .
git commit -m "Initial prereg + analysis framework"
git push
```

## Workflow
1. Use JPL Horizons (Sun-centered, ecliptic J2000, AU, 1-hour cadence) to export vectors
   for Earth and the interstellar object (3I/ATLAS). See `HORIZONS_HOWTO.txt`.
2. Compute triangle geometry and minimal-area windows (see geometry notebook/scripts).
3. Collect synchronized datasets (solar wind, geomagnetic, ELF/Schumann).
4. Run `analysis_template.py` to compute:
   - high‑k spectral slopes (β)
   - closure residual R(t) surrogate
   - windowed pass/fail outcomes
5. Inspect `predictions_summary.csv` and iterate thresholds if needed.

## Files
- `preregistration.md` — one-page hypothesis + pass/fail criteria
- `analysis_template.py` — starter analysis pipeline
- `data_schema.json` — expected input schema
- `.zenodo.json` — Zenodo metadata (auto-detected on release)
- `HORIZONS_HOWTO.txt` — exact settings for vector exports
- `.gitignore` — ignores common build artifacts and large data
- `LICENSE` — MIT

## Cite
If you create a Zenodo release, cite the DOI badge in this README.
Replace the badge placeholders after you publish the first release.

---
