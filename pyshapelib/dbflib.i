/* SWIG (www.swig.org) interface file for the dbf interface of shapelib
 *
 * At the moment (Dec 2000) this file is only useful to generate Python
 * bindings. Invoke swig as follows:
 *
 *	swig -python -shadow dbflib.i
 *
 * to generate dbflib_wrap.c and dbflib.py. dbflib_wrap.c defines a
 * bunch of Python-functions that wrap the appripriate dbflib functions
 * and dbflib.py contains an object oriented wrapper around
 * dbflib_wrap.c.
 *
 * This module defines one object type: DBFFile.
 */

/* this is the dbflib module */
%module dbflib

/* first a %{,%} block. These blocks are copied verbatim to the
 * dbflib_wrap.c file and are not parsed by SWIG. This is the place to
 * import headerfiles and define helper-functions that are needed by the
 * automatically generated wrappers.
 */

%{
#include "shapefil.h"


/* Read one attribute from the dbf handle and return it as a new python object
 *
 * If an error occurs, set the appropriate Python exception and return
 * NULL.
 *
 * Assume that the values of the record and field arguments are valid.
 * The name argument will be passed to DBFGetFieldInfo as is and should
 * thus be either NULL or a pointer to an array of at least 12 chars
 */
static PyObject *
do_read_attribute(DBFInfo * handle, int record, int field, char * name)
{
    int type, width;
    PyObject *value;

    type = DBFGetFieldInfo(handle, field, name, &width, NULL);
    /* For strings NULL and the empty string are indistinguishable
     * in DBF files. We prefer empty strings instead for backwards
     * compatibility reasons because older wrapper versions returned
     * emtpy strings as empty strings.
     */
    if (type != FTString && DBFIsAttributeNULL(handle, record, field))
    {
	value = Py_None;
	Py_INCREF(value);
    }
    else
    {
	switch (type)
	{
	case FTString:
	{
	    const char * temp = DBFReadStringAttribute(handle, record, field);
	    if (temp)
	    {
		value = PyString_FromString(temp);
	    }
	    else
	    {
		PyErr_Format(PyExc_IOError,
			     "Can't read value for row %d column %d",
			     record, field);
		value = NULL;
	    }
	    break;
	}
	case FTInteger:
	    value = PyInt_FromLong(DBFReadIntegerAttribute(handle, record,
							   field));
	    break;
	case FTDouble:
	    value = PyFloat_FromDouble(DBFReadDoubleAttribute(handle, record,
							      field));
	    break;
	default:
	    PyErr_Format(PyExc_TypeError, "Invalid field data type %d",
			 type);
	    value = NULL;
	}
    }
    if (!value)
	return NULL;

    return value;
}    

/* the read_attribute method. Return the value of the given record and
 * field as a python object of the appropriate type.
 * 
 * In case of error, set a python exception and return NULL. Since that
 * value will be returned to the python interpreter as is, the
 * interpreter should recognize the exception.
 */

static PyObject *
DBFInfo_read_attribute(DBFInfo * handle, int record, int field)
{
    if (record < 0 || record >= DBFGetRecordCount(handle))
    {
	PyErr_Format(PyExc_ValueError,
		     "record index %d out of bounds (record count: %d)",
		     record, DBFGetRecordCount(handle));
	return NULL;
    }

    if (field < 0 || field >= DBFGetFieldCount(handle))
    {
	PyErr_Format(PyExc_ValueError,
		     "field index %d out of bounds (field count: %d)",
		     field, DBFGetFieldCount(handle));
	return NULL;
    }

    return do_read_attribute(handle, record, field, NULL);
}
    

/* the read_record method. Return the record record as a dictionary with
 * whose keys are the names of the fields, and their values as the
 * appropriate Python type.
 * 
 * In case of error, set a python exception and return NULL. Since that
 * value will be returned to the python interpreter as is, the
 * interpreter should recognize the exception.
 */

static PyObject *
DBFInfo_read_record(DBFInfo * handle, int record)
{
    int num_fields;
    int i;
    int type, width;
    char name[12];
    PyObject *dict;
    PyObject *value;

    if (record < 0 || record >= DBFGetRecordCount(handle))
    {
	PyErr_Format(PyExc_ValueError,
		     "record index %d out of bounds (record count: %d)",
		     record, DBFGetRecordCount(handle));
	return NULL;
    }

    dict = PyDict_New();
    if (!dict)
	return NULL;
	
    num_fields = DBFGetFieldCount(handle);
    for (i = 0; i < num_fields; i++)
    {
	value = do_read_attribute(handle, record, i, name);
	if (!value)
	    goto fail;

	PyDict_SetItemString(dict, name, value);
	Py_DECREF(value);
    }

    return dict;

 fail:
    Py_XDECREF(dict);
    return NULL;
}

/* the write_record method. Write the record record given wither as a
 * dictionary or a sequence (i.e. a list or a tuple).
 *
 * If it's a dictionary the keys must be the names of the fields and
 * their value must have a suitable type. Only the fields actually
 * contained in the dictionary are written. Fields for which there's no
 * item in the dict are not modified.
 *
 * If it's a sequence, all fields must be present in the right order.
 *
 * In case of error, set a python exception and return NULL. Since that
 * value will be returned to the python interpreter as is, the
 * interpreter should recognize the exception.
 *
 * The method is implemented with two c-functions, write_field to write
 * a single field and DBFInfo_write_record as the front-end.
 */


/* write a single field of a record. */
static int
write_field(DBFHandle handle, int record, int field, int type,
	    PyObject * value)
{
    char * string_value;
    int int_value;
    double double_value;

    if (value == Py_None)
    {
	if (!DBFWriteNULLAttribute(handle, record, field))
	{
	    PyErr_Format(PyExc_IOError,
			 "can't write NULL field %d of record %d",
			 field, record);
	    return 0;
	}
    }
    else
    {
	switch (type)
	{
	case FTString:
	    string_value = PyString_AsString(value);
	    if (!string_value)
		return 0;
	    if (!DBFWriteStringAttribute(handle, record, field, string_value))
	    {
		PyErr_Format(PyExc_IOError,
			     "can't write field %d of record %d",
			     field, record);
		return 0;
	    }
	    break;

	case FTInteger:
	    int_value = PyInt_AsLong(value);
	    if (int_value == -1 && PyErr_Occurred())
		return 0;
	    if (!DBFWriteIntegerAttribute(handle, record, field, int_value))
	    {
		PyErr_Format(PyExc_IOError,
			     "can't write field %d of record %d",
			     field, record);
		return 0;
	    }
	    break;

	case FTDouble:
	    double_value = PyFloat_AsDouble(value);
	    if (double_value == -1 && PyErr_Occurred())
		return 0;
	    if (!DBFWriteDoubleAttribute(handle, record, field, double_value))
	    {
		PyErr_Format(PyExc_IOError,
			     "can't write field %d of record %d",
			     field, record);
		return 0;
	    }
	    break;

	default:
	    PyErr_Format(PyExc_TypeError, "Invalid field data type %d", type);
	    return 0;
	}
    }

    return 1;
}

static
PyObject *
DBFInfo_write_record(DBFHandle handle, int record, PyObject *record_object)
{
    int num_fields;
    int i, length;
    int type, width;
    char name[12];
    PyObject * value = NULL;

    num_fields = DBFGetFieldCount(handle);

    /* We used to use PyMapping_Check to test whether record_object is a
     * dictionary like object instead of PySequence_Check to test
     * whether it's a sequence. Unfortunately in Python 2.3
     * PyMapping_Check returns true for lists and tuples too so the old
     * approach doesn't work anymore.
     */
    if (PySequence_Check(record_object))
    {
	/* It's a sequence object. Iterate through all items in the
	 * sequence and write them to the appropriate field.
	 */
	length = PySequence_Length(record_object);
	if (length != num_fields)
	{
	    PyErr_SetString(PyExc_TypeError,
			    "record must have one item for each field");
	    goto fail;
	}
	for (i = 0; i < length; i++)
	{
	    type = DBFGetFieldInfo(handle, i, name, &width, NULL); 
	    value = PySequence_GetItem(record_object, i);
	    if (value)
	    {
		if (!write_field(handle, record, i, type, value))
		    goto fail;
		Py_DECREF(value);
	    }
	    else
	    {
		goto fail;
	    }
	}
    }
    else
    {
	/* It's a dictionary-like object. Iterate over the names of the
         * known fields and write the corresponding item
	 */
	for (i = 0; i < num_fields; i++)
	{
	    type = DBFGetFieldInfo(handle, i, name, &width, NULL);

	    /* if the dictionary has the key name write that object to
	     * the appropriate field, other wise just clear the python
	     * exception and do nothing.
	     */
	    value = PyMapping_GetItemString(record_object, name);
	    if (value)
	    {
		if (!write_field(handle, record, i, type, value))
		    goto fail;
		Py_DECREF(value);
	    }
	    else
	    {
		PyErr_Clear();
	    }
	}
    }

    Py_INCREF(Py_None);
    return Py_None;

 fail:
    Py_XDECREF(value);
    return NULL;
}
%}


/* The commit method implementation
 *
 * The method relies on the DBFUpdateHeader method which is not
 * available in shapelib <= 1.2.10.  setup.py defines
 * HAVE_UPDATE_HEADER's value depending on whether the function is
 * available in the shapelib version the code is compiled with.
 */
%{
static
void
DBFInfo_commit(DBFHandle handle)
{
#if HAVE_UPDATE_HEADER
    DBFUpdateHeader(handle);
#endif
}
%} 


/*
 * The SWIG Interface definition.
 */

/* include some common SWIG type definitions and standard exception
   handling code */
%include typemaps.i
%include exception.i

/* As for ShapeFile in shapelib.i, We define a new C-struct that holds
 * the DBFHandle. This is mainly done so we can separate the close()
 * method from the destructor but it also helps with exception handling.
 *
 * After the DBFFile has been opened or created the handle is not NULL.
 * The close() method closes the file and sets handle to NULL as an
 * indicator that the file has been closed.
 */

%{
    typedef struct {
	DBFHandle handle;
    } DBFFile;
%}


/* The first argument to the DBFFile methods is a DBFFile pointer.
 * We have to check whether handle is not NULL in most methods but not
 * all. In the destructor and the close method, it's OK for handle to be
 * NULL. We achieve this by checking whether the preprocessor macro
 * NOCHECK_$name is defined. SWIG replaces $name with the name of the
 * function for which the code is inserted. In the %{,%}-block below we
 * define the macros for the destructor and the close() method.
 */

%typemap(python,check) DBFFile *{
%#ifndef NOCHECK_$name
    if (!$target || !$target->handle)
	SWIG_exception(SWIG_TypeError, "dbffile already closed");
%#endif
}

%{
#define NOCHECK_delete_DBFFile
#define NOCHECK_DBFFile_close
%}


/* An exception handle for the constructor and the module level open()
 * and create() functions.
 *
 * Annoyingly, we *have* to put braces around the SWIG_exception()
 * calls, at least in the python case, because of the way the macro is
 * written. Of course, always putting braces around the branches of an
 * if-statement is often considered good practice.
 */
%typemap(python,except) DBFFile * {
    $function;
    if (!$source)
    {
    	SWIG_exception(SWIG_MemoryError, "no memory");
    }
    else if (!$source->handle)
    {
	SWIG_exception(SWIG_IOError, "$name failed");
    }
}

/* Exception handler for the add_field method */
%typemap(python,except) int DBFFile_add_field {
    $function;
    if ($source < 0)
    {
    	SWIG_exception(SWIG_RuntimeError, "add_field failed");
    }
}

/* define and use some typemaps for the field_info() method whose
 * C-implementation has three output parameters that are returned
 * through pointers passed into the function. SWIG already has
 * definitions for common types such as int* and we can use those for
 * the last two parameters:
 */

%apply int * OUTPUT { int * output_width }
%apply int * OUTPUT { int * output_decimals }

/* the fieldname has to be defined manually: */
%typemap(python,ignore) char *fieldname_out(char temp[12]) {
    $target = temp;
}

%typemap(python,argout) char *fieldname_out() {
    PyObject * string = PyString_FromString($source);
    $target = t_output_helper($target,string);
}



/*
 * The SWIG-version of the DBFFile struct 
 */

typedef	struct
{
    %addmethods {
	DBFFile(const char *file, const char * mode = "rb") {
	    DBFFile * self = malloc(sizeof(DBFFile));
	    if (self)
		self->handle = DBFOpen(file, mode);
	    return self;
	}
	    
	~DBFFile() {
	    if (self->handle)
		DBFClose(self->handle);
	    free(self);
	}

	void close() {
	    if (self->handle)
		DBFClose(self->handle);
	    self->handle = NULL;
	}

	int field_count() {
	    return DBFGetFieldCount(self->handle);
	}

	int record_count() {
	    return DBFGetRecordCount(self->handle);
	}

	int field_info(int iField, char * fieldname_out,
		       int * output_width, int * output_decimals) {
	    return DBFGetFieldInfo(self->handle, iField, fieldname_out,
				   output_width, output_decimals);
	}
	    
	PyObject * read_record(int record) {
	    return DBFInfo_read_record(self->handle, record);
	}

	PyObject * read_attribute(int record, int field) {
	    return DBFInfo_read_attribute(self->handle, record, field);
	}

	int add_field(const char * pszFieldName, DBFFieldType eType,
		      int nWidth, int nDecimals) {
	    return DBFAddField(self->handle, pszFieldName, eType, nWidth,
			       nDecimals);
	}

	PyObject *write_record(int record, PyObject *dict_or_sequence) {
	    return DBFInfo_write_record(self->handle, record,
					dict_or_sequence);
	}

	void commit() {
	    DBFInfo_commit(self->handle);
	}
	/* Delete the commit method from the class if it doesn't have a
	 * real implementation.
	 */
	%pragma(python) addtomethod="__class__:if not dbflibc._have_commit: del commit"

	/* The __del__ method generated by the old SWIG version we're
	 * tries to access self.thisown which may not be set at all when
	 * there was an exception during construction.  Therefore we
	 * override it with our own version.
	 * FIXME: It would be better to upgrade to a newer SWIG version
	 * or to get rid of SWIG entirely.
	 */
	%pragma(python) addtoclass = "
    def __del__(self,dbflibc=dbflibc):
        if getattr(self, 'thisown', 0):
            dbflibc.delete_DBFFile(self)
    "


    }
} DBFFile;


/*
 * Two module level functions, open() and create() that correspond to
 * DBFOpen and DBFCreate respectively. open() is equivalent to the
 * DBFFile constructor.
 */


%{
    DBFFile * open_DBFFile(const char * file, const char * mode)
    {
	DBFFile * self = malloc(sizeof(DBFFile));
	if (self)
	    self->handle = DBFOpen(file, mode);
	return self;
    }
%}

%name(open) %new DBFFile * open_DBFFile(const char * file,
					const char * mode = "rb");

%{
    DBFFile * create_DBFFile(const char * file)
    {
	DBFFile * self = malloc(sizeof(DBFFile));
	if (self)
	    self->handle = DBFCreate(file);
	return self;
    }
%}
%name(create) %new DBFFile * create_DBFFile(const char * file);



/* constant definitions copied from shapefil.h */
typedef enum {
  FTString,
  FTInteger,
  FTDouble,
  FTInvalid
} DBFFieldType;


/* Put the value of the HAVE_UPDATE_HEADER preprocessor macro into the
 * wrapper so that the __class__ pragma above knows when to remove the
 * commit method
 */
const int _have_commit = HAVE_UPDATE_HEADER;

