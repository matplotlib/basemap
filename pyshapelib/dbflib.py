# This file was created automatically by SWIG.
import dbflibc
class DBFFile:
    def __init__(self,*args):
        self.this = apply(dbflibc.new_DBFFile,args)
        self.thisown = 1

    def __del__(self,dbflibc=dbflibc):
        if self.thisown == 1 :
            dbflibc.delete_DBFFile(self)
    def close(*args):
        val = apply(dbflibc.DBFFile_close,args)
        return val
    def field_count(*args):
        val = apply(dbflibc.DBFFile_field_count,args)
        return val
    def record_count(*args):
        val = apply(dbflibc.DBFFile_record_count,args)
        return val
    def field_info(*args):
        val = apply(dbflibc.DBFFile_field_info,args)
        return val
    def read_record(*args):
        val = apply(dbflibc.DBFFile_read_record,args)
        return val
    def read_attribute(*args):
        val = apply(dbflibc.DBFFile_read_attribute,args)
        return val
    def add_field(*args):
        val = apply(dbflibc.DBFFile_add_field,args)
        return val
    def write_record(*args):
        val = apply(dbflibc.DBFFile_write_record,args)
        return val
    def commit(*args):
        val = apply(dbflibc.DBFFile_commit,args)
        return val
    def __repr__(self):
        return "<C DBFFile instance at %s>" % (self.this,)
    if not dbflibc._have_commit: del commit
    
    def __del__(self,dbflibc=dbflibc):
        if getattr(self, 'thisown', 0):
            dbflibc.delete_DBFFile(self)
    
class DBFFilePtr(DBFFile):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = DBFFile





#-------------- FUNCTION WRAPPERS ------------------

def open(*args, **kwargs):
    val = apply(dbflibc.open,args,kwargs)
    if val: val = DBFFilePtr(val); val.thisown = 1
    return val

def create(*args, **kwargs):
    val = apply(dbflibc.create,args,kwargs)
    if val: val = DBFFilePtr(val); val.thisown = 1
    return val



#-------------- VARIABLE WRAPPERS ------------------

FTString = dbflibc.FTString
FTInteger = dbflibc.FTInteger
FTDouble = dbflibc.FTDouble
FTInvalid = dbflibc.FTInvalid
_have_commit = dbflibc._have_commit
