

import xarray as xr
#from xarray import z3fs
#from xarray import zarr


ds = xr.open_zarr("s3://test-release-v1.0.1/mats-level-1b-limb-cropd-UV2-rev-2.zarr/",storage_options={"profile":"mats"})