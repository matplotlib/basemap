"""Pyrex code to provide python interfaces to proj.4 functions.
Make changes to this file, not the c-wrappers that Pyrex generates."""

import math

cdef double _rad2dg, _dg2rad
_rad2dg = 180./math.pi
_dg2rad = math.pi/180.

cdef extern from "proj_api.h":
    ctypedef void *projPJ
    ctypedef struct projUV:
        double u
        double v
    projPJ pj_init_plus(char *)
    projUV pj_fwd(projUV, projPJ)
    projUV pj_inv(projUV, projPJ)

cdef class Proj:
    cdef char *pjinitstring
    cdef void *projpj
    def __new__(self, char * pjinitstring):
        cdef void *projpj
        self.pjinitstring = pjinitstring
        projpj = pj_init_plus(pjinitstring)
        self.projpj = projpj
    def fwd(self, lons, lats):
        cdef projUV projxyout, projlonlatin
        cdef long ndim, i
        cdef double u, v
        try: # inputs are lists
            ndim = len(lons)
            x = []; y = []
            for i from 0 <= i < ndim:
                projlonlatin.u = _dg2rad*lons[i]
                projlonlatin.v = _dg2rad*lats[i]
                projxyout = pj_fwd(projlonlatin,self.projpj)
                x.append(projxyout.u)
                y.append(projxyout.v)
        except: # inputs are scalars
            projlonlatin.u = lons*_dg2rad
            projlonlatin.v = lats*_dg2rad
            projxyout = pj_fwd(projlonlatin,self.projpj)
            x = projxyout.u
            y = projxyout.v
        return x,y
    def inv(self, x, y):
        cdef projUV projxyin, projlonlatout
        cdef long ndim, i
        cdef double u, v
        try: # inputs are lists
            ndim = len(x)
            lons = []; lats = []
            for i from 0 <= i < ndim:
                projxyin.u = x[i]
                projxyin.v = y[i]
                projlonlatout = pj_inv(projxyin, self.projpj)
                lons.append(projlonlatout.u*_rad2dg)
                lats.append(projlonlatout.v*_rad2dg)
        except: # inputs are scalars.
            projxyin.u = x
            projxyin.v = y
            projlonlatout = pj_inv(projxyin, self.projpj)
            lons = projlonlatout.u*_rad2dg
            lats = projlonlatout.v*_rad2dg
        return lons,lats

