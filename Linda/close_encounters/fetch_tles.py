"""
Fetch historical TLEs from Space-Track for ISS, MATS, and AIM, save in 3LE format.

Requires a free Space-Track account: https://www.space-track.org/auth/createAccount

Usage:
    export SPACETRACK_USER='your_email'
    export SPACETRACK_PASS='your_password'
    python fetch_tles.py --start 2023-11-09 --end 2024-12-31

Notes:
    Space-Track's public TLE history endpoint requires authentication via cookie
    login.  This script uses the documented REST workflow.  Output files are
    written to ./data/iss_tles.txt and ./data/mats_tles.txt in 3LE format,
    ready for use with mats_awe_conjunctions.py.
    Output files: data/iss_tles.txt, data/mats_tles.txt, data/aim_tles.txt
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import requests

LOGIN_URL = "https://www.space-track.org/ajaxauth/login"
QUERY_URL = (
    "https://www.space-track.org/basicspacedata/query/class/gp_history/"
    "NORAD_CAT_ID/{norad}/orderby/EPOCH%20asc/EPOCH/{start}--{end}/"
    "format/3le"
)

ISS_NORAD = 25544
MATS_NORAD = 54227
AIM_NORAD = 31304


def fetch(norad: int, start: str, end: str, user: str, pw: str) -> str:
    s = requests.Session()
    r = s.post(LOGIN_URL, data={"identity": user, "password": pw}, timeout=30)
    r.raise_for_status()
    if "Failed" in r.text or r.status_code != 200:
        raise SystemExit(f"Space-Track login failed: {r.text[:200]}")
    url = QUERY_URL.format(norad=norad, start=start, end=end)
    r = s.get(url, timeout=120)
    r.raise_for_status()
    return r.text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True, help="YYYY-MM-DD")
    ap.add_argument("--end", required=True, help="YYYY-MM-DD")
    ap.add_argument("--out-dir", default="data")
    args = ap.parse_args()

    user = os.environ.get("SPACETRACK_USER")
    pw = os.environ.get("SPACETRACK_PASS")
    if not (user and pw):
        sys.exit("Set SPACETRACK_USER and SPACETRACK_PASS environment variables.")

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    for norad, label, fname in [
        (ISS_NORAD, "ISS", "iss_tles.txt"),
        (MATS_NORAD, "MATS", "mats_tles.txt"),
        (AIM_NORAD, "AIM", "aim_tles.txt"),
    ]:
        print(f"Fetching {label} ({norad}) from {args.start} to {args.end} ...")
        text = fetch(norad, args.start, args.end, user, pw)
        n = sum(1 for l in text.splitlines() if l.startswith("1 "))
        path = out / fname
        path.write_text(text)
        print(f"  -> {path}  ({n} TLEs)")


if __name__ == "__main__":
    main()
