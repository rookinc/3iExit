# Preregistration: Sun–Earth–3I Minimal-Area Amplification Test
Date: 2025-11-06 19:50 UTC
Lead: Scott Allen Cave (ORCID 0009-0000-0888-4232)

Hypothesis
---------
Amplification signatures in geospace occur only within ±24–72 h of Sun–Earth–3I minimal-area epochs A(t),
and only when moderate–high solar-wind forcing overlaps those windows. Signatures include:
- spectral steepening in high‑k tails toward −4,
- stress–release in closure residuals R(t),
- cross‑modal synchronization.

Windows
-------
- Use JPL Horizons vectors (Sun-centered, Ecliptic J2000, AU, 1 h cadence) for Earth and 3I/ATLAS.
- Compute A(t) = 0.5‖r_E×r_I‖ and θ(t). Identify local minima of A_scaled = A/(|r_E||r_I|) with ≥24 h separation.
- Analyze ±24–72 h windows around minima only.

Forcing Gate
------------
Any of:
- |B| ≥ 8 nT for ≥3 h, or V_sw ≥ 500 km/s for ≥3 h, or Bz ≤ −5 nT for ≥3 h.

Predictions (pass/fail)
-----------------------
P1 (Geometry gating): signatures occur only inside minimal-area windows.
P2 (Spectral steepening): β ≤ −3.7 for ≥3 consecutive 1 h windows (FDR q ≤ 0.05, ≥2 stations).
P3 (Stress–release): |R(t)| exceeds 90th percentile of 60 d baseline and relaxes within 24 h after driver abates.
P4 (Cross‑modal sync): ≥2 modalities show simultaneous plateau/joint steepening for ≥2 h (perm p < 0.01).
P5 (Quiet Sun null): if no CME/CIR/HCS overlaps a window, no positive signatures should appear.

Decision
--------
If P2–P4 fail across all geometry + forcing overlaps, the amplification‑actor hypothesis is falsified for this epoch.
