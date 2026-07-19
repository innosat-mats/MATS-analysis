"""
Build one video per CCD channel from the July 3rd imagery in the locally
extracted MATS Level0 payload data (see read_MATS_payload_data_from_disk).

Each frame is a single decoded CCD image, contrast-stretched using the 1st/99th
percentile of that channel's whole day of data (so brightness is comparable
frame-to-frame), upscaled for visibility, and stamped with its exposure time.
"""

import datetime as DT
from pathlib import Path

import imageio.v2 as imageio
import numpy as np
import pandas as pd
from matplotlib import cm
from PIL import Image, ImageDraw

from mats_utils.rawdata.read_data import read_MATS_payload_data_from_disk
from mats_l1_processing.read_parquet_functions import convert_image_data, remove_faulty_rows
from mats_l1_processing.read_in_functions import channel_num_to_str

EXTRACTED_PATH = str(Path(__file__).resolve().parents[2] / "extracted")
DAY_START = DT.datetime(2026, 7, 3)
DAY_STOP = DT.datetime(2026, 7, 4)
OUTDIR = Path(__file__).resolve().parent
FPS = 10
TARGET_LONG_EDGE = 600  # px; images are upscaled (integer factor) to about this size


def load_ccd_images():
    df = read_MATS_payload_data_from_disk(DAY_START, DAY_STOP, data_type="CCD", path=EXTRACTED_PATH)
    convert_image_data(df)
    df = remove_faulty_rows(df)
    df = df[(df["EXPDate"] >= pd.Timestamp(DAY_START, tz="UTC"))
            & (df["EXPDate"] < pd.Timestamp(DAY_STOP, tz="UTC"))]
    return df


def make_video_for_channel(sub, ccdsel):
    sub = sub.sort_values("EXPDate")
    images = np.stack(sub["IMAGE"].to_numpy())
    channel_name = channel_num_to_str(ccdsel)

    vmin, vmax = np.percentile(images, [1, 99])
    if vmax <= vmin:
        vmax = vmin + 1

    scale = max(1, TARGET_LONG_EDGE // max(images.shape[1:]))
    fname = OUTDIR / f"ccd_video_{channel_name}_CCDSEL{ccdsel}.mp4"

    writer = imageio.get_writer(str(fname), fps=FPS)
    try:
        for img, ts in zip(images, sub["EXPDate"]):
            norm = np.clip((img.astype(np.float64) - vmin) / (vmax - vmin), 0, 1)
            rgb = (cm.viridis(norm)[:, :, :3] * 255).astype(np.uint8)
            frame = Image.fromarray(rgb).resize(
                (rgb.shape[1] * scale, rgb.shape[0] * scale), Image.NEAREST
            )
            draw = ImageDraw.Draw(frame)
            draw.text((4, 4), f"{channel_name}  {ts}", fill=(255, 255, 255))
            writer.append_data(np.array(frame))
    finally:
        writer.close()

    print(f"{channel_name} (CCDSEL {ccdsel}): {len(images)} frames -> {fname.name}")


def main():
    df = load_ccd_images()
    print(f"Loaded {len(df)} CCD images for {DAY_START.date()}")

    for ccdsel, sub in df.groupby("CCDSEL"):
        if len(sub) == 0:
            continue
        make_video_for_channel(sub, ccdsel)


if __name__ == "__main__":
    main()
