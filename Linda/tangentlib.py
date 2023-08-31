import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize_scalar
from skyfield.positionlib import ICRF
from skyfield.api import wgs84
from skyfield.units import Distance


def rotate(unitvec, yaw, pitch, roll, deg=False):
    def Rx(v, th):
        s = np.sin(th)
        c = np.cos(th)
        return np.matmul([[1, 0, 0], [0, c, -s], [0, s, c]], v)

    def Ry(v, th):
        s = np.sin(th)
        c = np.cos(th)
        return np.matmul([[c, 0, s], [0, 1, 0], [-s, 0, c]], v)

    def Rz(v, th):
        s = np.sin(th)
        c = np.cos(th)
        return np.matmul([[c, -s, 0], [s, c, 0], [0, 0, 1]], v)
    if deg:
        roll *= (np.pi/180)
        pitch *= (np.pi/180)
        yaw *= (np.pi/180)
    return Rz(Ry(Rx(unitvec, roll), pitch), yaw)


def xyz2radec(vector, deg=False, positivera=False):
    ra = np.arctan2(vector[1,:], vector[0,:])
    if positivera:
        if ra < 0:
            ra += 2*np.pi
    dec = np.arcsin(vector[2,:]/norm(vector,axis=0))
    if deg:
        ra *= 180./np.pi
        dec *= 180./np.pi
    return [ra, dec]


def radec2xyz(ra, dec, deg=True):
    if deg:
        ra *= np.pi/180.
        dec *= np.pi/180.
    z = np.sin(dec)
    x = np.cos(ra)*np.cos(dec)
    y = np.sin(ra)*np.cos(dec)
    return [x, y, z]


def funpitch(pitch, t, th, pos, yaw, rotmatrix):
    # print(pitch*180/np.pi)
    FOV = rotate(np.array([1, 0, 0]), yaw, pitch, 0, deg=False)
    FOV = np.matmul(rotmatrix, FOV)
    tp = findtangent(t, pos, FOV)
    return((tp.fun-th)**2)


def funheight(s, t, pos, FOV):
    newp = pos + s * FOV
    newp = ICRF(Distance(m=newp).au, t=t, center=399)
    return wgs84.subpoint(newp).elevation.m


def findtangent(t, pos, FOV):
    res = minimize_scalar(funheight, args=(t, pos, FOV), bracket=(1e5, 3e5))
    return res


def findpitch(th, t, pos, yaw, rotmatrix):
    res = minimize_scalar(funpitch, args=(t, th, pos, yaw, rotmatrix),
                          method="Bounded", bounds=(np.deg2rad(-30), np.deg2rad(-10)))
    return res.x
