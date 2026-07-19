# MATS – AWE conjunction finder

Finds when the Swedish **MATS** satellite (NORAD 54227, dawn-dusk SSO at 585 km)
and the **AWE** instrument on the ISS (NORAD 25544) make near-coincident
observations of the mesopause region (~85–90 km altitude).

## Why this matters

Both instruments image the same atmospheric layer — mesospheric airglow / OH
emission — but with very different geometries:

- **AWE** stares nadir from the ISS in a ~600 km wide swath at ~87 km altitude (OH airglow).
- **MATS** images the limb forward of its orbit; the tangent point sits ~3000 km along-track
  ahead of the spacecraft at 75–110 km tangent altitude (O₂ A-band airglow + NLC).

So "close vicinity" in this context can mean two different things, and the script
reports both separately:

| Criterion | Definition | Use |
|---|---|---|
| **A** | MATS limb tangent point falls inside the AWE nadir swath, within ±half-swath cross-track and ±1500 km along-track of the current ISS position | Best science overlap — same air sampled by both |
| **B** | The two spacecraft themselves are within 2000 km of each other in 3-D and within 30 min in time | Looser geometric proximity |

## Files

```
mats_awe_conjunctions.py    # main scan script
fetch_tles.py               # downloads historical TLEs from Space-Track
data/                       # put TLE files here
out/                        # generated CSVs land here
```

## Step-by-step

### 1. Install dependencies

```bash
pip install numpy sgp4 requests
```

### 2. Get a Space-Track account (free)

<https://www.space-track.org/auth/createAccount>

### 3. Download historical TLEs

```bash
export SPACETRACK_USER='your_email'
export SPACETRACK_PASS='your_password'
iop
```

This writes `data/iss_tles.txt` and `data/mats_tles.txt` in 3LE format.

> **Why historical TLEs?** A single TLE drifts; for a multi-month scan you
> want 50–100+ TLEs per satellite (Space-Track typically has one every few
> days), and the script auto-picks the TLE epoch closest to each evaluation
> time.

### 4. Run the scan

```bash
python mats_awe_conjunctions.py \
    --start 2023-11-15 --end 2024-12-31 \
    --step 60 \
    --proximity-km 2000 --proximity-min 30
```

Tunable knobs (defaults shown):

| Flag | Default | Meaning |
|---|---|---|
| `--step` | 60 s | Sampling cadence; 30 s is safer if you want to catch brief overlaps |
| `--proximity-km` | 2000 | Criterion B distance threshold |
| `--proximity-min` | 30 | Gap (min) used to group consecutive samples into a single "event" |
| `--awe-swath-km` | 600 | AWE total cross-track swath at airglow altitude |
| `--awe-altitude` | 87 | AWE OH airglow altitude (km) |
| `--tangent-altitude` | 85 | MATS limb tangent altitude (km) |
| `--tangent-distance` | 3000 | Along-track distance from MATS to its tangent point (km) |

### 5. Read the output

`out/criterion_A_tangent_in_swath.csv` — one row per event:
- `t_start_utc`, `t_end_utc`, `duration_s`, `n_samples`
- `min_sep_km`, `min_sep_utc` — closest spacecraft approach in this event
- `min_sep_iss_lat`, `min_sep_iss_lon`, `min_sep_mats_lat`, `min_sep_mats_lon` — sub-satellite points at closest approach
- `min_tangent_lat`, `min_tangent_lon`, `min_cross_track_km` — MATS tangent point location and distance from ISS ground track

`out/criterion_B_spacecraft_proximity.csv` — same structure, no tangent fields.

`out/criterion_*_raw_samples.csv` — every sample that satisfied the criterion (not grouped), for diagnostic plotting.

`out/summary.txt` — counts and configuration.

## Expected scale

For Nov 2023 – end of 2024 with 60 s sampling that's ~600k samples per satellite.
Runtime on a laptop: a few minutes. Criterion A typically yields a few hundred
events over a year (the geometries are quite specific); Criterion B yields more
because 2000 km is geometrically generous.

## Caveats and where to improve fidelity

1. **MATS attitude.** The script assumes MATS looks straight along its
   velocity vector. In flight MATS pointed slightly off (a few degrees) and
   the attitude profile changed during the mission. If you have access to the
   MATS L1A geolocation product (which contains the actual tangent point
   lat/lon/altitude per image), substitute those for `compute_mats_tangent_point`
   — the conjunction list will be much more precise.

2. **AWE swath.** Treated as a rectangular strip ±300 km cross-track of the
   ISS ground track, at the airglow layer. AWE actually has four cameras with
   slightly overlapping fields; if you have the AWE L1B geolocation product
   you can replace the rectangle with the actual per-image footprint.

3. **TLE accuracy.** SGP4 with TLEs is good to a few km along-track at LEO,
   degrading with epoch age. For a 30 min / 2000 km criterion this is fine.
   For Criterion A (where you care about ground-track overlap to ~150 km) it
   is also fine, but tighter overlap definitions (e.g. <50 km) would benefit
   from the MATS GPS-based precise orbit if available from OHB Sweden.

4. **TEME → ECEF rotation.** Polar motion is neglected (~m-level error,
   irrelevant here).

## Citing the source data

- ISS state vectors: 18 SDS / Space-Track.org
- MATS state vectors: 18 SDS / Space-Track.org
- MATS mission concept and limb geometry: <https://www.eoportal.org/satellite-missions/mats>
- AWE mission concept and instrument: <https://science.nasa.gov/mission/awe/>, <https://www.eoportal.org/satellite-missions/iss-awe>
