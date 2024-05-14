'''
This code file is adapted from Lukas
'''
from mats_utils.rawdata.read_data import read_MATS_data
import netCDF4 as nc
import numpy as np
import argparse
import datetime as DT
import pandas
import pickle
import pdb

REQUIRED = ['START_TIME', 'STOP_TIME', 'VERSION', 'CHANNEL']
DEFAULTS = {'VERSION': '0.6'}
NCDF_TYPES = {np.int8: 'i1', np.int16: 'i2', np.int32: 'i4', np.int64: 'i8',
              np.float32: 'f4', np.float16: 'f8', np.bool_: 'i1', str: 'str',
              pandas._libs.tslibs.timestamps.Timestamp: 'timestamp'}
STR_LEN = 20


def get_filter(channel):
    filters = {"IR1": 1, "IR2": 4, "IR3": 3, "IR4": 2, "UV1": 5, "UV2": 6}
    try:
        filt = filters[channel]
    except Exception:
        raise ValueError(f"Invalid channel: {channel}!")

    return {'CCDSEL': [filt, filt]}


def get_args():
    parser = argparse.ArgumentParser(description="Get MATS data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--start_time", dest="START_TIME", type=int, nargs=5,
                        help="Start time for data set.")
    parser.add_argument("--stop_time", dest="STOP_TIME", type=int, nargs=5,
                        help="Start time for data set.")
    parser.add_argument("--version", dest="VERSION", type=str,
                        help="Data version.")
    parser.add_argument("--channel", dest="CHANNEL", type=str, help="Channel to retrieve",
                        choices=["IR1", "IR2", "IR3", "IR4", "UV1", "UV2"])
    parser.add_argument("--conf", type=str, help="Read configuration from file.")
    parser.add_argument("--ncdf_out", type=str, help="Output file name (ncdf)")
    parser.add_argument("--pickle_out", type=str, help="Output file name (pickle)")
    return parser.parse_args()


def write_ncdf_L1b(pdata, outfile, channel, version):
    with nc.Dataset(outfile, 'w') as nf:
        num_images = len(pdata)
        # Global parameters
        nf.date_created = DT.datetime.now().isoformat()
        nf.date_modified = DT.datetime.now().isoformat()
        nf.channel = channel
        nf.data_version = version

        # Create dimensions and their corresponding variables
        # Dimension variable parameters: (<len of variable>, <long name>, <dimension unlimited>)
        params = {"num_images": (num_images, "Image number", True),
                  "im_col": (pdata["NCOL"][0] + 1, "Column number", False),
                  "im_row": (pdata["NROW"][0], "Row number", False)}
        for dim, param in params.items():
            nf.createDimension(dim, None if param[2] else param[0])
            ncdf_create_var(np.arange(param[0]), nf, dim, (dim, ), np.int16, long_name=param[1])

        # Create variable for ImageCalibrated
        ncdf_create_var(np.stack(pdata["ImageCalibrated"].to_numpy(), axis=0),
                        nf, "ImageCalibrated", ("num_images", "im_row", "im_col"), np.float32)

        # Create variable for CalibrationErrors
        error = np.stack([np.stack(pdata["CalibrationErrors"][i], axis=0)
                          for i in range(pdata["CalibrationErrors"].shape[0])], axis=0)
        ncdf_create_var(error, nf, "CalibrationErrors", ("num_images", "im_row", "im_col"), np.int16)

        # Create variable for BadColumns
        bad_cols = np.zeros((num_images, params["im_col"][0]))
        for i, bcols in enumerate(pdata['BadColumns'].to_numpy()):
            for j in bcols:
                bad_cols[i, j] = 1
        ncdf_create_var(bad_cols, nf, "BadColumns", ("num_images", "im_col"), np.int8)

        # Handle the remaining variables automatically
        for var in pdata.keys():
            if var not in nf.variables.keys():
                data = pdata[var].to_numpy()
                ncdf_create_var(data, nf, var, ("num_images", ), type(data.flat[0]))
                                                                                                                

def ncdf_create_var(data, dset, name, dims, dtype, long_name=None, units=None):
    NCDF_TYPES = {np.int8: 'i1', np.int16: 'i2', np.int32: 'i4', np.int64: 'i8',
                  np.float32: 'f4', np.float64: 'f8', np.bool_: 'i1', str: 'str',
                  pandas._libs.tslibs.timestamps.Timestamp: 'timestamp',
                  np.datetime64: 'timestamp', np.ndarray: 'ndarray'}

    if dtype in NCDF_TYPES.keys():
        type_id = NCDF_TYPES[dtype]
    else:
        raise NotImplementedError(f"Input variable {name} is of type {dtype}, writing this to ncdf is not implemented.")
    if type_id == 'str':
        # dname = f"{name}_string"
        # cdata = nc.stringtochar(np.array(data, dtype=str))
        # dset.createDimension(dname, cdata.shape[-1])
        ncvar = dset.createVariable(name, str, dims)
        ncvar[:] = data
    elif type_id == 'timestamp':
        ncvar = dset.createVariable(name, 'f8', dims)
        ncvar.units = "Seconds since 2000.01.01 00:00 UTC"
        if dtype == np.datetime64:
            ncvar[:] = (data - np.datetime64("2000-01-01 00:00:00.0")) / np.timedelta64(1, "s")
        else:
            ncvar[:] = (data - DT.datetime(2000, 1, 1, 0, 0, tzinfo=DT.timezone.utc)) /\
                DT.timedelta(0, 1)
    elif type_id == 'ndarray':
        try:
            size = data[0].shape[0]
            dname = f"{name}_comp"
            if size > 0:
                ctype = type(data[0].flat[0])
                assert np.issubdtype(ctype, np.integer) or ctype == np.float32 or ctype == np.float64
                assert all([len(data[i].shape) == 1 and data[i].shape[0] == size for i in range(data.shape[0])])
                dset.createDimension(dname, size)
                ncvar = dset.createVariable(name, NCDF_TYPES[ctype], (*dims, dname))
                ncvar[:] = np.stack(data, axis=0)
            else:
                assert all([data[i].shape[0] == 0 for i in range(data.shape[0])])
                dset.createDimension(dname, 1)
                ncvar = dset.createVariable(name, 'i1', (*dims, dname))
        except Exception:
            raise NotImplementedError(f"Input variable {name} of type ndarray cannot be imported automatically.")
    else:
        ncvar = dset.createVariable(name, type_id, dims)
        ncvar[:] = data
    if long_name is not None:
        ncvar.long_name = long_name
    if units is not None:
        ncvar.units = units
    return ncvar


def set_vars(vrs, sources):
    for var in vrs:
        unset = True
        for src in sources:
            if (var in src) and (src[var] is not None):
                globals()[var] = src[var]
                unset = False
                break
        if unset:
            raise ValueError(f"The variable {var} is not set!")


def main():
    args = get_args()
    if args.conf:
        exec(open(args.conf).read())
    set_vars(REQUIRED, (vars(args), dict(globals(), **locals()), DEFAULTS))
    dftop = read_MATS_data(DT.datetime(*START_TIME), DT.datetime(*STOP_TIME),
                           get_filter(CHANNEL), level="1b", version=VERSION)
    if args.ncdf_out:
        write_ncdf_L1b(dftop, args.ncdf_out, CHANNEL, VERSION)
    if args.pickle_out:
        with open(args.pickle_out, 'wb') as handle:
            pickle.dump(dftop, handle)
    elif not args.ncdf_out:
        raise RuntimeError("Please specify output file name for at least one format (ncdf and/or pickle).")


if __name__ == "__main__":
    main()
    
    
    