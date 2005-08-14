# This file was created automatically by SWIG.
import shapelibc
class SHPObject:
    def __init__(self,*args):
        self.this = apply(shapelibc.new_SHPObject,args)
        self.thisown = 1

    def __del__(self,shapelibc=shapelibc):
        if self.thisown == 1 :
            shapelibc.delete_SHPObject(self)
    def extents(*args):
        val = apply(shapelibc.SHPObject_extents,args)
        return val
    def vertices(*args):
        val = apply(shapelibc.SHPObject_vertices,args)
        return val
    __setmethods__ = {
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = SHPObject.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "type" : shapelibc.SHPObject_type_get,
        "id" : shapelibc.SHPObject_id_get,
    }
    def __getattr__(self,name):
        method = SHPObject.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C SHPObject instance at %s>" % (self.this,)
class SHPObjectPtr(SHPObject):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = SHPObject



class ShapeFile:
    def __init__(self,*args):
        self.this = apply(shapelibc.new_ShapeFile,args)
        self.thisown = 1

    def __del__(self,shapelibc=shapelibc):
        if self.thisown == 1 :
            shapelibc.delete_ShapeFile(self)
    def close(*args):
        val = apply(shapelibc.ShapeFile_close,args)
        return val
    def info(*args):
        val = apply(shapelibc.ShapeFile_info,args)
        return val
    def read_object(*args):
        val = apply(shapelibc.ShapeFile_read_object,args)
        if val: val = SHPObjectPtr(val) ; val.thisown = 1
        return val
    def write_object(*args):
        val = apply(shapelibc.ShapeFile_write_object,args)
        return val
    def cobject(*args):
        val = apply(shapelibc.ShapeFile_cobject,args)
        return val
    def __repr__(self):
        return "<C ShapeFile instance at %s>" % (self.this,)
class ShapeFilePtr(ShapeFile):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = ShapeFile





#-------------- FUNCTION WRAPPERS ------------------

def open(*args, **kwargs):
    val = apply(shapelibc.open,args,kwargs)
    if val: val = ShapeFilePtr(val); val.thisown = 1
    return val

def create(*args, **kwargs):
    val = apply(shapelibc.create,args,kwargs)
    if val: val = ShapeFilePtr(val); val.thisown = 1
    return val

c_api = shapelibc.c_api

type_name = shapelibc.type_name

part_type_name = shapelibc.part_type_name



#-------------- VARIABLE WRAPPERS ------------------

SHPT_NULL = shapelibc.SHPT_NULL
SHPT_POINT = shapelibc.SHPT_POINT
SHPT_ARC = shapelibc.SHPT_ARC
SHPT_POLYGON = shapelibc.SHPT_POLYGON
SHPT_MULTIPOINT = shapelibc.SHPT_MULTIPOINT
SHPT_POINTZ = shapelibc.SHPT_POINTZ
SHPT_ARCZ = shapelibc.SHPT_ARCZ
SHPT_POLYGONZ = shapelibc.SHPT_POLYGONZ
SHPT_MULTIPOINTZ = shapelibc.SHPT_MULTIPOINTZ
SHPT_POINTM = shapelibc.SHPT_POINTM
SHPT_ARCM = shapelibc.SHPT_ARCM
SHPT_POLYGONM = shapelibc.SHPT_POLYGONM
SHPT_MULTIPOINTM = shapelibc.SHPT_MULTIPOINTM
SHPT_MULTIPATCH = shapelibc.SHPT_MULTIPATCH
SHPP_TRISTRIP = shapelibc.SHPP_TRISTRIP
SHPP_TRIFAN = shapelibc.SHPP_TRIFAN
SHPP_OUTERRING = shapelibc.SHPP_OUTERRING
SHPP_INNERRING = shapelibc.SHPP_INNERRING
SHPP_FIRSTRING = shapelibc.SHPP_FIRSTRING
SHPP_RING = shapelibc.SHPP_RING
