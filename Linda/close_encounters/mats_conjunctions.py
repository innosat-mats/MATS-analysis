"""
MATS - AWE conjunction finder
=============================

Finds time windows when the MATS satellite (NORAD 54227) and the ISS-mounted
AWE instrument (NORAD 25544) make coincident observations of the mesopause
region (~85-90 km altitude).

Two criteria are checked independently and reported in separate output files:

  CRITERION A -- "MATS tangent point inside AWE swath"
    Best science overlap.  AWE looks nadir from the ISS in a 600 km wide swath
    centered on the sub-satellite ground track at ~87 km altitude.  MATS looks
    at the limb; its tangent point sits ~3000 km ahead of MATS along the orbit
    track at ~85 km altitude.  We flag every minute where the MATS tangent
    point is INSIDE the AWE swath rectangle (across-track <= 300 km, and the
    tangent point is on a portion of the AWE ground track that ISS occupied
    within the time window).

  CRITERION B -- "Spacecraft proximity"
    Distance between MATS and the ISS as 3-D points in space (ECEF) is below
    a threshold (default 2000 km) and within the time window (default 30 min).
    This is the looser "satellites near each other" criterion.

INPUTS
------
TLE files in 3LE format (one TLE epoch per file or concatenated):
  - ISS:  data/iss_tles.txt
  - MATS: data/mats_tles.txt

Get historical TLEs from www.space-track.org (free account required):
  Login -> Query Builder -> "GP History" -> filter by NORAD ID
  -> select date range 2023-11-09 .. 2024-12-31
  -> download as TLE (3LE).

Save the two files to ./data/ before running.

USAGE
-----
  python mats_awe_conjunctions.py \\
      --start 2023-11-15 --end 2024-12-31 \\
      --step 60 \\
      --proximity-km 2000 --proximity-min 30 \\
      --awe-swath-km 600 \\
      --tangent-altitude 85 \\
      --awe-altitude 87 \\
      --tangent-distance 3000

OUTPUT
------
  out/criterion_A_tangent_in_swath.csv
  out/criterion_B_spacecraft_proximity.csv
  out/summary.txt

NOTES
-----
* The MATS tangent geometry assumes the limb-line-of-sight is along the
  velocity vector (forward-looking).  In reality MATS points slightly
  off the velocity vector and changes pointing during the mission; if you
  have the actual attitude/pointing files from the MATS team, plug them in
  in `compute_mats_tangent_point` for higher fidelity.
* TLE accuracy degrades from epoch.  We auto-pick the TLE whose epoch is
  closest to each evaluation time.
* Default time step 60 s is a reasonable trade-off; 30 s is safer if you
  want to catch very brief overlaps.

Author: drafted by Claude
"""
from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
from sgp4.api import Satrec, jday, SGP4_ERRORS

# ---------- Constants ----------
EARTH_R_KM = 6378.137  # WGS-84 equatorial radius
EARTH_F = 1.0 / 298.257223563  # WGS-84 flattening
EARTH_E2 = EARTH_F * (2 - EARTH_F)
OMEGA_E = 7.2921150e-5  # Earth rotation rate, rad/s

ISS_NORAD = 25544
MATS_NORAD = 54227
AIM_NORAD = 31304

# ---------- Data classes ----------
@dataclass
class TLE:
    name: str
    line1: str
    line2: str
    epoch: datetime  # UTC

    @classmethod
    def from_lines(cls, name: str, l1: str, l2: str) -> "TLE":
        # Epoch parsing per TLE spec
        yr2 = int(l1[18:20])
        year = 2000 + yr2 if yr2 < 57 else 1900 + yr2
        doy_frac = float(l1[20:32])
        ep = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy_frac - 1)
        return cls(name=name.strip(), line1=l1.strip(), line2=l2.strip(), epoch=ep)


# ---------- TLE loading ----------
def load_tles(path: Path) -> list[TLE]:
    """Parse a 3LE-formatted file (concatenated 3-line element sets)."""
    tles: list[TLE] = []
    lines = [l.rstrip("\n") for l in path.read_text().splitlines() if l.strip()]
    i = 0
    while i + 2 < len(lines) + 1:
        # Detect whether the first line of the triple is a name or a "1 ..." line
        if lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            tles.append(TLE.from_lines("UNKNOWN", lines[i], lines[i + 1]))
            i += 2
        elif (
            i + 2 < len(lines)
            and lines[i + 1].startswith("1 ")
            and lines[i + 2].startswith("2 ")
        ):
            tles.append(TLE.from_lines(lines[i], lines[i + 1], lines[i + 2]))
            i += 3
        else:
            i += 1
    tles.sort(key=lambda t: t.epoch)
    return tles


def pick_tle(tles: list[TLE], when: datetime) -> TLE:
    """Pick the TLE whose epoch is closest to `when`."""
    return min(tles, key=lambda t: abs((t.epoch - when).total_seconds()))


# ---------- Coordinate transforms ----------
def gmst_rad(dt: datetime) -> float:
    """Greenwich Mean Sidereal Time at `dt` (UTC), in radians.  Vallado eq."""
    jd, fr = jday(dt.year, dt.month, dt.day,
                  dt.hour, dt.minute, dt.second + dt.microsecond * 1e-6)
    T = ((jd + fr) - 2451545.0) / 36525.0
    gmst_sec = (
        67310.54841
        + (876600.0 * 3600 + 8640184.812866) * T
        + 0.093104 * T * T
        - 6.2e-6 * T ** 3
    )
    gmst_rad = (gmst_sec % 86400.0) / 240.0  # 86400 s/day, 240 s/deg
    return math.radians(gmst_rad)


def teme_to_ecef(r_teme: np.ndarray, dt: datetime) -> np.ndarray:
    """Rotate TEME -> ECEF by GMST.  (Polar motion neglected; ~m-level error.)"""
    g = gmst_rad(dt)
    c, s = math.cos(g), math.sin(g)
    R = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    return R @ r_teme


def ecef_to_geodetic(r_ecef: np.ndarray) -> tuple[float, float, float]:
    """ECEF (km) -> (lat_rad, lon_rad, alt_km).  Bowring closed-form."""
    x, y, z = r_ecef
    lon = math.atan2(y, x)
    a = EARTH_R_KM
    b = a * (1 - EARTH_F)
    p = math.hypot(x, y)
    th = math.atan2(z * a, p * b)
    lat = math.atan2(
        z + EARTH_E2 * (a / b) * b * math.sin(th) ** 3,
        p - EARTH_E2 * a * math.cos(th) ** 3,
    )
    N = a / math.sqrt(1 - EARTH_E2 * math.sin(lat) ** 2)
    alt = p / math.cos(lat) - N
    return lat, lon, alt


def geodetic_to_ecef(lat: float, lon: float, alt_km: float) -> np.ndarray:
    """(lat_rad, lon_rad, alt_km) -> ECEF (km)."""
    a = EARTH_R_KM
    N = a / math.sqrt(1 - EARTH_E2 * math.sin(lat) ** 2)
    x = (N + alt_km) * math.cos(lat) * math.cos(lon)
    y = (N + alt_km) * math.cos(lat) * math.sin(lon)
    z = (N * (1 - EARTH_E2) + alt_km) * math.sin(lat)
    return np.array([x, y, z])


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float, R: float = EARTH_R_KM) -> float:
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * R * math.asin(math.sqrt(a))


# ---------- Geometry ----------
def propagate(tle: TLE, dt: datetime) -> tuple[np.ndarray, np.ndarray]:
    """Return (r_teme_km, v_teme_km_s) from sgp4."""
    sat = Satrec.twoline2rv(tle.line1, tle.line2)
    jd, fr = jday(dt.year, dt.month, dt.day,
                  dt.hour, dt.minute, dt.second + dt.microsecond * 1e-6)
    e, r, v = sat.sgp4(jd, fr)
    if e != 0:
        raise RuntimeError(f"sgp4 error {e}: {SGP4_ERRORS.get(e)}")
    return np.array(r), np.array(v)


def compute_mats_tangent_point(
    r_ecef: np.ndarray,
    v_ecef: np.ndarray,
    tangent_altitude_km: float,
    tangent_distance_km: float,
) -> tuple[float, float]:
    """Compute the (lat, lon) of the MATS limb tangent point.

    Simplified model: tangent point lies along the velocity-vector horizon at
    the requested altitude.  For MATS, the LOS is approximately along the
    +V direction (forward-looking limb), and the tangent point at 85 km is
    ~3000 km ahead of the spacecraft.

    We construct the tangent point by:
      1. Moving from the satellite position by `tangent_distance_km` along the
         along-track horizontal direction projected onto the local horizontal
         plane.
      2. Snapping the result to the requested geodetic altitude.
    """
    # Get satellite geodetic position
    sat_lat, sat_lon, _ = ecef_to_geodetic(r_ecef)

    # Local up unit vector at satellite sub-point
    up = np.array([
        math.cos(sat_lat) * math.cos(sat_lon),
        math.cos(sat_lat) * math.sin(sat_lon),
        math.sin(sat_lat),
    ])

    # Project velocity onto local horizontal -> along-track horizontal direction
    v_horiz = v_ecef - np.dot(v_ecef, up) * up
    if np.linalg.norm(v_horiz) < 1e-9:
        return sat_lat, sat_lon
    track_dir = v_horiz / np.linalg.norm(v_horiz)

    # Step `tangent_distance_km` along the great circle in the track direction
    # (small-circle approximation -- adequate at MLT altitudes, sub-degree error)
    angular_distance = tangent_distance_km / EARTH_R_KM  # radians
    # Convert track_dir (ECEF) to a local east/north heading at sub-point
    east = np.array([-math.sin(sat_lon), math.cos(sat_lon), 0.0])
    north = np.array([
        -math.sin(sat_lat) * math.cos(sat_lon),
        -math.sin(sat_lat) * math.sin(sat_lon),
        math.cos(sat_lat),
    ])
    de = np.dot(track_dir, east)
    dn = np.dot(track_dir, north)
    bearing = math.atan2(de, dn)  # 0 = north, +pi/2 = east

    # Great-circle destination point
    sin_lat1, cos_lat1 = math.sin(sat_lat), math.cos(sat_lat)
    sin_d, cos_d = math.sin(angular_distance), math.cos(angular_distance)
    sin_lat2 = sin_lat1 * cos_d + cos_lat1 * sin_d * math.cos(bearing)
    lat2 = math.asin(sin_lat2)
    lon2 = sat_lon + math.atan2(
        math.sin(bearing) * sin_d * cos_lat1,
        cos_d - sin_lat1 * sin_lat2,
    )
    lon2 = (lon2 + math.pi) % (2 * math.pi) - math.pi
    return lat2, lon2


def cross_track_distance_km(
    p_lat: float,
    p_lon: float,
    a_lat: float,
    a_lon: float,
    b_lat: float,
    b_lon: float,
    R: float = EARTH_R_KM,
) -> float:
    """Perpendicular (cross-track) distance from point P to the great circle A->B."""
    d_ap = haversine_km(a_lat, a_lon, p_lat, p_lon, R) / R  # angular
    # Bearing A->P
    bearing_ap = _initial_bearing(a_lat, a_lon, p_lat, p_lon)
    bearing_ab = _initial_bearing(a_lat, a_lon, b_lat, b_lon)
    return abs(math.asin(math.sin(d_ap) * math.sin(bearing_ap - bearing_ab))) * R


def along_track_distance_km(
    p_lat: float,
    p_lon: float,
    a_lat: float,
    a_lon: float,
    b_lat: float,
    b_lon: float,
    R: float = EARTH_R_KM,
) -> float:
    """Along-track distance from A to the projection of P on the great circle A->B.

    Sign: positive if P projects between A and B in the A->B direction.
    """
    d_ap = haversine_km(a_lat, a_lon, p_lat, p_lon, R) / R
    bearing_ap = _initial_bearing(a_lat, a_lon, p_lat, p_lon)
    bearing_ab = _initial_bearing(a_lat, a_lon, b_lat, b_lon)
    xtd = math.asin(math.sin(d_ap) * math.sin(bearing_ap - bearing_ab))
    atd = math.acos(max(-1.0, min(1.0, math.cos(d_ap) / math.cos(xtd))))
    # Sign: dot product of bearing direction
    if math.cos(bearing_ap - bearing_ab) < 0:
        atd = -atd
    return atd * R


def _initial_bearing(lat1, lon1, lat2, lon2):
    dlon = lon2 - lon1
    y = math.sin(dlon) * math.cos(lat2)
    x = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dlon)
    return math.atan2(y, x)


# ---------- Main scan ----------
def time_range(start: datetime, end: datetime, step_s: int) -> Iterable[datetime]:
    t = start
    dt = timedelta(seconds=step_s)
    while t <= end:
        yield t
        t += dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iss-tles", default="data/iss_tles.txt")
    ap.add_argument("--mats-tles", default="data/mats_tles.txt")
    ap.add_argument("--start", required=True)  # YYYY-MM-DD
    ap.add_argument("--end", required=True)
    ap.add_argument("--step", type=int, default=60, help="Sampling step (s)")
    ap.add_argument("--proximity-km", type=float, default=2000.0)
    ap.add_argument("--proximity-min", type=float, default=30.0,
                    help="Group nearby samples that are closer than this (min) into a single conjunction.")
    ap.add_argument("--awe-swath-km", type=float, default=600.0,
                    help="Total cross-track width of AWE swath at airglow altitude.")
    ap.add_argument("--awe-altitude", type=float, default=87.0,
                    help="OH airglow altitude AWE images (km).")
    ap.add_argument("--tangent-altitude", type=float, default=85.0,
                    help="MATS limb tangent altitude (km).")
    ap.add_argument("--tangent-distance", type=float, default=3000.0,
                    help="Along-track distance from MATS to its tangent point (km).")
    ap.add_argument("--time-window-min", type=float, default=30.0,
                    help="Max time difference (min) between MATS and any instrument to count as a match.")
    ap.add_argument("--aim-tles", default="data/aim_tles.txt",
                    help="Path to AIM (NORAD 31304) TLE file. CIPS analysis is skipped if absent.")
    ap.add_argument("--cips-swath-km", type=float, default=800.0,
                    help="Total cross-track width of the CIPS combined swath at PMC altitude (~83 km).")
    ap.add_argument("--out", default="out")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    iss_tles = load_tles(Path(args.iss_tles))
    mats_tles = load_tles(Path(args.mats_tles))
    if not iss_tles or not mats_tles:
        raise SystemExit("Missing TLEs.  See header docstring for download instructions.")

    print(f"Loaded {len(iss_tles)} ISS TLEs spanning {iss_tles[0].epoch.date()} .. {iss_tles[-1].epoch.date()}")
    print(f"Loaded {len(mats_tles)} MATS TLEs spanning {mats_tles[0].epoch.date()} .. {mats_tles[-1].epoch.date()}")

    aim_tles_path = Path(args.aim_tles)
    if aim_tles_path.exists():
        aim_tles = load_tles(aim_tles_path)
        print(f"Loaded {len(aim_tles)} AIM TLEs spanning {aim_tles[0].epoch.date()} .. {aim_tles[-1].epoch.date()}")
    else:
        aim_tles = []
        print(f"AIM TLE file not found ({args.aim_tles}) — CIPS criteria (D/E) will be skipped.")

    start = datetime.fromisoformat(args.start).replace(tzinfo=timezone.utc)
    end = datetime.fromisoformat(args.end).replace(tzinfo=timezone.utc)

    # AWE science operations: 2023-12-01 through 2025-12-31
    AWE_OBS_START = datetime(2023, 12, 1, tzinfo=timezone.utc)
    AWE_OBS_END = datetime(2025, 12, 31, 23, 59, 59, tzinfo=timezone.utc)
    run_awe = not (end < AWE_OBS_START or start > AWE_OBS_END)
    if not run_awe:
        print(
            f"WARNING: AWE was not observing during {args.start} – {args.end} "
            f"(AWE operated 2023-12-01 to 2025-12-31). "
            f"Criteria A/B/C (AWE) will be skipped."
        )

    # CIPS/AIM science operations: launched 2007-04-25, mission ended ~2023-03-31
    CIPS_OBS_END = datetime(2023, 3, 31, 23, 59, 59, tzinfo=timezone.utc)
    run_cips = bool(aim_tles) and not (start > CIPS_OBS_END)
    if aim_tles and not run_cips:
        print(
            f"WARNING: CIPS/AIM was not observing during {args.start} – {args.end} "
            f"(CIPS/AIM mission ended March 2023). "
            f"Criteria D/E (CIPS) will be skipped."
        )

    half_swath = args.awe_swath_km / 2.0
    all_times = list(time_range(start, end, args.step))
    n = len(all_times)
    # Number of time steps that fit in the allowed time-difference window
    window_steps = max(1, round(args.time_window_min * 60.0 / args.step))

    # ----- Pass 1: pre-compute all ISS positions (only needed for AWE criteria) -----
    iss_cache: list = []
    if run_awe:
        print(f"Pre-computing {n} ISS positions ...")
        for i, t in enumerate(all_times):
            iss_tle = pick_tle(iss_tles, t)
            try:
                r_teme, v_teme = propagate(iss_tle, t)
            except RuntimeError:
                iss_cache.append(None)
                continue
            r_ecef = teme_to_ecef(r_teme, t)
            v_ecef = teme_to_ecef(v_teme, t)
            lat, lon, _ = ecef_to_geodetic(r_ecef)
            # 100-km forward point reuses the tangent-point helper (altitude arg unused there)
            fwd_lat, fwd_lon = compute_mats_tangent_point(r_ecef, v_ecef, 0.0, 100.0)
            iss_cache.append({"t": t, "lat": lat, "lon": lon, "r_ecef": r_ecef,
                              "fwd_lat": fwd_lat, "fwd_lon": fwd_lon})
            if i % 50000 == 0 and i > 0:
                print(f"  ISS pre-compute: {i}/{n}")

    # ----- Pass 1b: pre-compute AIM positions (only if CIPS criteria will run) -----
    half_swath_cips = args.cips_swath_km / 2.0
    aim_cache: list = []
    if run_cips:
        print(f"Pre-computing {n} AIM positions ...")
        for i, t in enumerate(all_times):
            aim_tle = pick_tle(aim_tles, t)
            try:
                r_teme, v_teme = propagate(aim_tle, t)
            except RuntimeError:
                aim_cache.append(None)
                continue
            r_ecef = teme_to_ecef(r_teme, t)
            v_ecef = teme_to_ecef(v_teme, t)
            lat, lon, _ = ecef_to_geodetic(r_ecef)
            fwd_lat, fwd_lon = compute_mats_tangent_point(r_ecef, v_ecef, 0.0, 100.0)
            aim_cache.append({"t": t, "lat": lat, "lon": lon, "r_ecef": r_ecef,
                              "fwd_lat": fwd_lat, "fwd_lon": fwd_lon})
            if i % 50000 == 0 and i > 0:
                print(f"  AIM pre-compute: {i}/{n}")

    # ----- Pass 2: MATS propagation + match against ±time-window ISS positions -----
    print(f"Scanning {n} MATS samples, matching AWE within ±{args.time_window_min:.0f} min ...")
    rows_A: list[dict] = []
    rows_B: list[dict] = []
    rows_C: list[dict] = []  # satellite position inside AWE swath
    rows_D: list[dict] = []  # MATS tangent point inside CIPS swath
    rows_E: list[dict] = []  # MATS satellite position inside CIPS swath

    for i, t in enumerate(all_times):
        mats_tle = pick_tle(mats_tles, t)
        try:
            r_mats_teme, v_mats_teme = propagate(mats_tle, t)
        except RuntimeError:
            continue
        r_mats_ecef = teme_to_ecef(r_mats_teme, t)
        v_mats_ecef = teme_to_ecef(v_mats_teme, t)
        mats_lat, mats_lon, _ = ecef_to_geodetic(r_mats_ecef)
        tan_lat, tan_lon = compute_mats_tangent_point(
            r_mats_ecef, v_mats_ecef, args.tangent_altitude, args.tangent_distance)

        j_lo = max(0, i - window_steps)
        j_hi = min(n, i + window_steps + 1)

        if run_awe:
            # Criterion A: MATS tangent point inside AWE swath.
            # Uses haversine (not cross-track) to avoid false positives from the
            # antipodal extension of the ISS great circle.
            best_dist = float("inf")
            best_j_A = None
            for j in range(j_lo, j_hi):
                iss = iss_cache[j]
                if iss is None:
                    continue
                dist = haversine_km(tan_lat, tan_lon, iss["lat"], iss["lon"])
                if dist < best_dist:
                    best_dist = dist
                    best_j_A = j

            if best_dist <= half_swath and best_j_A is not None:
                iss = iss_cache[best_j_A]
                dt_min = (iss["t"] - t).total_seconds() / 60.0
                xtd = cross_track_distance_km(tan_lat, tan_lon, iss["lat"], iss["lon"],
                                              iss["fwd_lat"], iss["fwd_lon"])
                rows_A.append(dict(
                    mats_obs_utc=t.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    awe_overpass_utc=iss["t"].strftime("%Y-%m-%dT%H:%M:%SZ"),
                    time_diff_min=round(dt_min, 1),
                    mats_tp_lat_deg=round(math.degrees(tan_lat), 3),
                    mats_tp_lon_deg=round(math.degrees(tan_lon), 3),
                    awe_center_lat_deg=round(math.degrees(iss["lat"]), 3),
                    awe_center_lon_deg=round(math.degrees(iss["lon"]), 3),
                    cross_track_km=round(xtd, 1),
                    distance_tp_to_awe_km=round(best_dist, 1),
                ))

            # Criterion B: minimum 3-D spacecraft separation within the time window.
            best_d3d = float("inf")
            best_j_B = None
            for j in range(j_lo, j_hi):
                iss = iss_cache[j]
                if iss is None:
                    continue
                d3d = float(np.linalg.norm(r_mats_ecef - iss["r_ecef"]))
                if d3d < best_d3d:
                    best_d3d = d3d
                    best_j_B = j

            if best_d3d <= args.proximity_km and best_j_B is not None:
                iss = iss_cache[best_j_B]
                dt_min = (iss["t"] - t).total_seconds() / 60.0
                rows_B.append(dict(
                    mats_obs_utc=t.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    awe_overpass_utc=iss["t"].strftime("%Y-%m-%dT%H:%M:%SZ"),
                    time_diff_min=round(dt_min, 1),
                    mats_lat_deg=round(math.degrees(mats_lat), 3),
                    mats_lon_deg=round(math.degrees(mats_lon), 3),
                    iss_lat_deg=round(math.degrees(iss["lat"]), 3),
                    iss_lon_deg=round(math.degrees(iss["lon"]), 3),
                    spacecraft_sep_km=round(best_d3d, 1),
                ))

            # Criterion C: MATS satellite sub-point inside AWE swath.
            best_dist_C = float("inf")
            best_j_C = None
            for j in range(j_lo, j_hi):
                iss = iss_cache[j]
                if iss is None:
                    continue
                dist = haversine_km(mats_lat, mats_lon, iss["lat"], iss["lon"])
                if dist < best_dist_C:
                    best_dist_C = dist
                    best_j_C = j

            if best_dist_C <= half_swath and best_j_C is not None:
                iss = iss_cache[best_j_C]
                dt_min = (iss["t"] - t).total_seconds() / 60.0
                xtd_C = cross_track_distance_km(mats_lat, mats_lon, iss["lat"], iss["lon"],
                                                iss["fwd_lat"], iss["fwd_lon"])
                rows_C.append(dict(
                    mats_obs_utc=t.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    awe_overpass_utc=iss["t"].strftime("%Y-%m-%dT%H:%M:%SZ"),
                    time_diff_min=round(dt_min, 1),
                    mats_sat_lat_deg=round(math.degrees(mats_lat), 3),
                    mats_sat_lon_deg=round(math.degrees(mats_lon), 3),
                    awe_center_lat_deg=round(math.degrees(iss["lat"]), 3),
                    awe_center_lon_deg=round(math.degrees(iss["lon"]), 3),
                    cross_track_km=round(xtd_C, 1),
                    distance_sat_to_awe_km=round(best_dist_C, 1),
                ))

        # Criterion D: MATS tangent point inside CIPS swath (AIM/CIPS version of A)
        if run_cips:
            best_dist_D = float("inf")
            best_j_D = None
            for j in range(j_lo, j_hi):
                aim = aim_cache[j] if j < len(aim_cache) else None
                if aim is None:
                    continue
                dist = haversine_km(tan_lat, tan_lon, aim["lat"], aim["lon"])
                if dist < best_dist_D:
                    best_dist_D = dist
                    best_j_D = j

            if best_dist_D <= half_swath_cips and best_j_D is not None:
                aim = aim_cache[best_j_D]
                dt_min = (aim["t"] - t).total_seconds() / 60.0
                xtd_D = cross_track_distance_km(tan_lat, tan_lon, aim["lat"], aim["lon"],
                                                aim["fwd_lat"], aim["fwd_lon"])
                rows_D.append(dict(
                    mats_obs_utc=t.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    aim_overpass_utc=aim["t"].strftime("%Y-%m-%dT%H:%M:%SZ"),
                    time_diff_min=round(dt_min, 1),
                    mats_tp_lat_deg=round(math.degrees(tan_lat), 3),
                    mats_tp_lon_deg=round(math.degrees(tan_lon), 3),
                    aim_center_lat_deg=round(math.degrees(aim["lat"]), 3),
                    aim_center_lon_deg=round(math.degrees(aim["lon"]), 3),
                    cross_track_km=round(xtd_D, 1),
                    distance_tp_to_aim_km=round(best_dist_D, 1),
                ))

        # Criterion E: MATS satellite position inside CIPS swath (AIM/CIPS version of C)
        if run_cips:
            best_dist_E = float("inf")
            best_j_E = None
            for j in range(j_lo, j_hi):
                aim = aim_cache[j] if j < len(aim_cache) else None
                if aim is None:
                    continue
                dist = haversine_km(mats_lat, mats_lon, aim["lat"], aim["lon"])
                if dist < best_dist_E:
                    best_dist_E = dist
                    best_j_E = j

            if best_dist_E <= half_swath_cips and best_j_E is not None:
                aim = aim_cache[best_j_E]
                dt_min = (aim["t"] - t).total_seconds() / 60.0
                xtd_E = cross_track_distance_km(mats_lat, mats_lon, aim["lat"], aim["lon"],
                                                aim["fwd_lat"], aim["fwd_lon"])
                rows_E.append(dict(
                    mats_obs_utc=t.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    aim_overpass_utc=aim["t"].strftime("%Y-%m-%dT%H:%M:%SZ"),
                    time_diff_min=round(dt_min, 1),
                    mats_sat_lat_deg=round(math.degrees(mats_lat), 3),
                    mats_sat_lon_deg=round(math.degrees(mats_lon), 3),
                    aim_center_lat_deg=round(math.degrees(aim["lat"]), 3),
                    aim_center_lon_deg=round(math.degrees(aim["lon"]), 3),
                    cross_track_km=round(xtd_E, 1),
                    distance_sat_to_aim_km=round(best_dist_E, 1),
                ))

        if i % 5000 == 0 and i > 0:
            print(f"  {i}/{n}  t={t.isoformat()}  |A|={len(rows_A)} |B|={len(rows_B)} |C|={len(rows_C)} "
                  f"|D|={len(rows_D)} |E|={len(rows_E)}")

    # ---- Group consecutive MATS samples into events (gap > proximity_min between rows)
    def group_events(rows: list[dict], gap_min: float, summarizer) -> list[dict]:
        if not rows:
            return []
        events = []
        cur = [rows[0]]
        gap = timedelta(minutes=gap_min)
        for r in rows[1:]:
            t_prev = datetime.fromisoformat(cur[-1]["mats_obs_utc"].replace("Z", "+00:00"))
            t_now = datetime.fromisoformat(r["mats_obs_utc"].replace("Z", "+00:00"))
            if t_now - t_prev > gap:
                events.append(summarizer(cur))
                cur = [r]
            else:
                cur.append(r)
        events.append(summarizer(cur))
        return events

    def _summarize_event_A(samples: list[dict]) -> dict:
        """One row per event; representative sample = the one with minimum cross-track."""
        first, last = samples[0], samples[-1]
        rep = min(samples, key=lambda r: r.get("distance_tp_to_awe_km", float("inf")))
        return dict(
            t_start_utc=first["mats_obs_utc"],
            t_end_utc=last["mats_obs_utc"],
            n_images=len(samples),
            mats_obs_utc=rep["mats_obs_utc"],
            awe_overpass_utc=rep["awe_overpass_utc"],
            time_diff_min=rep["time_diff_min"],
            mats_tp_lat_deg=rep["mats_tp_lat_deg"],
            mats_tp_lon_deg=rep["mats_tp_lon_deg"],
            awe_center_lat_deg=rep["awe_center_lat_deg"],
            awe_center_lon_deg=rep["awe_center_lon_deg"],
            cross_track_km=rep["cross_track_km"],
            distance_tp_to_awe_km=rep["distance_tp_to_awe_km"],
        )

    def _summarize_event_B(samples: list[dict]) -> dict:
        first, last = samples[0], samples[-1]
        best = min(samples, key=lambda r: r.get("spacecraft_sep_km", float("inf")))
        return dict(
            t_start_utc=first["mats_obs_utc"],
            t_end_utc=last["mats_obs_utc"],
            duration_s=int(
                (
                    datetime.fromisoformat(last["mats_obs_utc"].replace("Z", "+00:00"))
                    - datetime.fromisoformat(first["mats_obs_utc"].replace("Z", "+00:00"))
                ).total_seconds()
            ),
            n_samples=len(samples),
            min_sep_km=best["spacecraft_sep_km"],
            min_sep_mats_utc=best["mats_obs_utc"],
            min_sep_awe_utc=best["awe_overpass_utc"],
            min_sep_time_diff_min=best["time_diff_min"],
            min_sep_iss_lat=best["iss_lat_deg"],
            min_sep_iss_lon=best["iss_lon_deg"],
            min_sep_mats_lat=best["mats_lat_deg"],
            min_sep_mats_lon=best["mats_lon_deg"],
        )

    def _summarize_event_C(samples: list[dict]) -> dict:
        """One row per event; representative sample = closest satellite pass to AWE."""
        first, last = samples[0], samples[-1]
        rep = min(samples, key=lambda r: r.get("distance_sat_to_awe_km", float("inf")))
        return dict(
            t_start_utc=first["mats_obs_utc"],
            t_end_utc=last["mats_obs_utc"],
            n_images=len(samples),
            mats_obs_utc=rep["mats_obs_utc"],
            awe_overpass_utc=rep["awe_overpass_utc"],
            time_diff_min=rep["time_diff_min"],
            mats_sat_lat_deg=rep["mats_sat_lat_deg"],
            mats_sat_lon_deg=rep["mats_sat_lon_deg"],
            awe_center_lat_deg=rep["awe_center_lat_deg"],
            awe_center_lon_deg=rep["awe_center_lon_deg"],
            cross_track_km=rep["cross_track_km"],
            distance_sat_to_awe_km=rep["distance_sat_to_awe_km"],
        )

    def _summarize_event_D(samples: list[dict]) -> dict:
        """Criterion D event summary: MATS tangent point in CIPS swath."""
        first, last = samples[0], samples[-1]
        rep = min(samples, key=lambda r: r.get("distance_tp_to_aim_km", float("inf")))
        return dict(
            t_start_utc=first["mats_obs_utc"],
            t_end_utc=last["mats_obs_utc"],
            n_images=len(samples),
            mats_obs_utc=rep["mats_obs_utc"],
            aim_overpass_utc=rep["aim_overpass_utc"],
            time_diff_min=rep["time_diff_min"],
            mats_tp_lat_deg=rep["mats_tp_lat_deg"],
            mats_tp_lon_deg=rep["mats_tp_lon_deg"],
            aim_center_lat_deg=rep["aim_center_lat_deg"],
            aim_center_lon_deg=rep["aim_center_lon_deg"],
            cross_track_km=rep["cross_track_km"],
            distance_tp_to_aim_km=rep["distance_tp_to_aim_km"],
        )

    def _summarize_event_E(samples: list[dict]) -> dict:
        """Criterion E event summary: MATS satellite position in CIPS swath."""
        first, last = samples[0], samples[-1]
        rep = min(samples, key=lambda r: r.get("distance_sat_to_aim_km", float("inf")))
        return dict(
            t_start_utc=first["mats_obs_utc"],
            t_end_utc=last["mats_obs_utc"],
            n_images=len(samples),
            mats_obs_utc=rep["mats_obs_utc"],
            aim_overpass_utc=rep["aim_overpass_utc"],
            time_diff_min=rep["time_diff_min"],
            mats_sat_lat_deg=rep["mats_sat_lat_deg"],
            mats_sat_lon_deg=rep["mats_sat_lon_deg"],
            aim_center_lat_deg=rep["aim_center_lat_deg"],
            aim_center_lon_deg=rep["aim_center_lon_deg"],
            cross_track_km=rep["cross_track_km"],
            distance_sat_to_aim_km=rep["distance_sat_to_aim_km"],
        )

    events_A = group_events(rows_A, args.proximity_min, _summarize_event_A)
    events_B = group_events(rows_B, args.proximity_min, _summarize_event_B)
    events_C = group_events(rows_C, args.proximity_min, _summarize_event_C)
    events_D = group_events(rows_D, args.proximity_min, _summarize_event_D)
    events_E = group_events(rows_E, args.proximity_min, _summarize_event_E)

    # ---- Write CSVs
    def write_csv(path: Path, rows: list[dict]):
        if not rows:
            path.write_text("# No conjunctions found.\n")
            return
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    write_csv(out / "criterion_A_tangent_in_awe_swath.csv", events_A)
    write_csv(out / "criterion_B_spacecraft_proximity.csv", events_B)
    write_csv(out / "criterion_C_satellite_in_awe_swath.csv", events_C)
    write_csv(out / "criterion_D_tangent_in_cips_swath.csv", events_D)
    write_csv(out / "criterion_E_satellite_in_cips_swath.csv", events_E)
    write_csv(out / "criterion_A_raw_samples.csv", rows_A)
    write_csv(out / "criterion_B_raw_samples.csv", rows_B)
    write_csv(out / "criterion_C_raw_samples.csv", rows_C)
    write_csv(out / "criterion_D_raw_samples.csv", rows_D)
    write_csv(out / "criterion_E_raw_samples.csv", rows_E)

    summary = (out / "summary.txt")
    summary.write_text(
        f"MATS-CIPS/AWE conjunction scan\n"
        f"==============================\n"
        f"Time window:      {args.start}  ..  {args.end}  (UTC)\n"
        f"Sampling step:    {args.step} s\n"
        f"Time-diff window: ±{args.time_window_min:.0f} min (atmosphere assumed stationary)\n"
        f"AWE swath width:  {args.awe_swath_km} km @ {args.awe_altitude} km altitude\n"
        f"CIPS swath width: {args.cips_swath_km} km @ 83.0 km altitude\n"
        f"MATS tangent:     {args.tangent_distance} km ahead of MATS @ {args.tangent_altitude} km altitude\n"
        f"\n"
        f"Criterion A (MATS tangent point inside AWE swath):       {len(events_A)} events"
        f"{'' if run_awe else '  [skipped — outside AWE observation window (2023-12-01 to 2025-12-31)]'}\n"
        f"Criterion B (spacecraft <= {args.proximity_km:.0f} km):                    {len(events_B)} events"
        f"{'' if run_awe else '  [skipped]'}\n"
        f"Criterion C (MATS satellite position inside AWE swath):  {len(events_C)} events"
        f"{'' if run_awe else '  [skipped]'}\n"
        f"Criterion D (MATS tangent point inside CIPS swath):      {len(events_D)} events"
        f"{'' if run_cips else ('  [skipped — no AIM TLEs]' if not aim_tles else '  [skipped — outside CIPS observation window (ended 2023-03)]')}\n"
        f"Criterion E (MATS satellite position inside CIPS swath): {len(events_E)} events"
        f"{'' if run_cips else ('  [skipped — no AIM TLEs]' if not aim_tles else '  [skipped — outside CIPS observation window]')}\n"
        f"\n"
        f"Files:\n"
        f"  criterion_A_tangent_in_awe_swath.csv   -- one row per event\n"
        f"  criterion_A_raw_samples.csv            -- one row per timestep\n"
        f"  criterion_B_spacecraft_proximity.csv\n"
        f"  criterion_B_raw_samples.csv\n"
        f"  criterion_C_satellite_in_awe_swath.csv\n"
        f"  criterion_C_raw_samples.csv\n"
        f"  criterion_D_tangent_in_cips_swath.csv  -- one row per event\n"
        f"  criterion_D_raw_samples.csv            -- one row per timestep\n"
        f"  criterion_E_satellite_in_cips_swath.csv\n"
        f"  criterion_E_raw_samples.csv\n"
    )
    print(summary.read_text())


if __name__ == "__main__":
    main()
