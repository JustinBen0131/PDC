# Reconstructing electromagnetic shower positions in sPHENIX

A photon deposits energy in several neighboring calorimeter towers. Estimating where it struck the detector from those energies is harder near a tower boundary or when the photon enters at an angle. Justin Bennett built this study to reconstruct the impact position independently and compare candidate corrections with the established sPHENIX cluster reconstruction. Single-photon simulation supplies a known impact point. A second test uses pairs of showers from neutral-pion decays, whose reconstructed mass depends on their opening angle.

Two outcomes came from the work. The analysis tested alternative position corrections in simulation and Run‑24 proton–proton data. Separately, three merged sPHENIX reconstruction changes made raw and corrected coordinates, cluster timing and signed incidence angles available for downstream studies. The merged changes establish which information analysts can use. The candidate corrections remain comparisons in this repository.

## How the comparisons work

- Follow the energy shared around the seed tower, including showers crossing the azimuthal boundary between tower 255 and tower 0. Compare reconstructed and true positions as functions of energy and position within a tower. [Event-level implementation](src/PositionDependentCorrection.cc) · [Azimuthal boundary demonstration](macros/CoordinateVisualization.py)
- Use the tower orientation and a direction from the collision vertex to determine the signed angle of incidence on the calorimeter. Look for residual trends with geometry and energy. [Geometry sketch](macros/ProjectiveSketch.py) · [ROOT analysis macros](macros)
- Compare neutral-pion mass peaks using the alternative coordinates in single-pion simulation and collision data. This tests how a position change affects a reconstructed two-photon observable. [Single-pion study](SINGLE_PI0_MC) · [Data-analysis outputs](DataOutput)

## Upstream contributions by Justin Bennett

| Merged change | What analysts gained |
|---|---|
| [Raw and corrected coordinates, PR #3879](https://github.com/sPHENIX-Collaboration/coresoftware/pull/3879) | Access to both tower-space centroids for the same cluster. |
| [Energy-weighted cluster time, PR #3956](https://github.com/sPHENIX-Collaboration/coresoftware/pull/3956) | A cluster time calculated from its constituent towers. |
| [Coordinate, timing and incidence properties, PR #4043](https://github.com/sPHENIX-Collaboration/coresoftware/pull/4043) | A common cluster interface with signed incidence angles and the position and timing properties. |
