/* SWIG (www.swig.org) interface file for shapelib
 *
 * At the moment (Dec 2000) this file is only useful to generate Python
 * bindings. Invoke swig as follows:
 *
 *	swig -python -shadow shapelib.i
 *
 * to generate shapelib_wrap.c and shapelib.py. shapelib_wrap.c
 * defines a bunch of Python-functions that wrap the appripriate
 * shapelib functions and shapelib.py contains an object oriented
 * wrapper around shapelib_wrap.c.
 *
 * Shapelib, and hence this module too, defines two types of objects,
 * shapes and shapefiles.
 */

%module shapelib

/*
 * First, a %{,%}-Block. These blocks are copied verbatim to the
 * shapelib_wrap.c file and are not parsed by SWIG. This is the place to
 * import headerfiles and define helper-functions that are needed by the
 * automatically generated wrappers.
 */

%{

/* import the shapelib headefile. */
#include "shapefil.h"
#include "pyshapelib_api.h"
    
/*
 * Rename a few shapelib functions that are effectively methods with
 * preprocessor macros so that they have the names that swig expects
 * (e.g. the destructor of SHPObject has to be called delete_SHPObject)
 */

#define delete_SHPObject SHPDestroyObject
    
/*
 * The extents() method of SHPObject.
 *
 * Return the extents as a tuple of two 4-element lists with the min.
 * and max. values of x, y, z, m.
 */
static PyObject *
SHPObject_extents(SHPObject *object)
{
    return Py_BuildValue("[dddd][dddd]",
			 object->dfXMin, object->dfYMin, object->dfZMin,
			 object->dfMMin, 
			 object->dfXMax, object->dfYMax, object->dfZMax,
			 object->dfMMax);
}


/*
 * The vertices() method of SHPObject.
 *
 * Return the x and y coords of the vertices as a list of lists of
 * tuples.
 */

static PyObject* build_vertex_list(SHPObject *object, int index, int length);

static PyObject*
SHPObject_vertices(SHPObject *object)
{
    PyObject *result = NULL;
    PyObject *part = NULL;
    int part_idx, vertex_idx;
    int length = 0;


    if (object->nParts > 0)
    {
	/* A multipart shape. Usual for SHPT_ARC and SHPT_POLYGON */
	
	result = PyList_New(object->nParts);
	if (!result)
	    return NULL;

	for (part_idx = 0, vertex_idx = 0; part_idx < object->nParts;
	     part_idx++)
	{
	    if (part_idx < object->nParts - 1)
		length = (object->panPartStart[part_idx + 1]
			  - object->panPartStart[part_idx]);
	    else
		length = object->nVertices - object->panPartStart[part_idx];
	    
	    part = build_vertex_list(object, vertex_idx, length);
	    if (!part)
		goto fail;

	    if (PyList_SetItem(result, part_idx, part) < 0)
		goto fail;

	    vertex_idx += length;
	}
    }
    else
    {
	/* only one part. usual for SHPT_POINT */
	result = build_vertex_list(object, 0, object->nVertices);
    }

    return result;

 fail:
    Py_XDECREF(part);
    Py_DECREF(result);
    return NULL;
}


/* Return the length coordinates of the shape object starting at vertex
 * index as a Python-list of tuples. Helper function for
 * SHPObject_vertices.
 */
static PyObject*
build_vertex_list(SHPObject *object, int index, int length)
{
    int i;
    PyObject * list;
    PyObject * vertex = NULL;

    list = PyList_New(length);
    if (!list)
	return NULL;

    for (i = 0; i < length; i++, index++)
    {
	vertex = Py_BuildValue("dd", object->padfX[index],
			       object->padfY[index]);
	if (!vertex)
	    goto fail;
	if (PyList_SetItem(list, i, vertex) < 0)
	    goto fail;
    }

    return list;

 fail:
    Py_XDECREF(vertex);
    Py_DECREF(list);
    return NULL;
}





/* The constructor of SHPObject. parts is a list of lists of tuples
 * describing the parts and their vertices just likethe output of the
 * vertices() method. part_type_list is the list of part-types and may
 * be NULL. For the meaning of the part-types and their default value
 * see the Shaplib documentation.
 */
SHPObject * new_SHPObject(int type, int id, PyObject * parts,
			  PyObject * part_type_list)
{
    /* arrays to hold thex and y coordinates of the  vertices */
    double *xs = NULL, *ys = NULL;
    /* number of all vertices of all parts */
    int num_vertices;
    /* number of parts in the list parts */
    int num_parts;
    /* start index of in xs and ys of the part currently worked on */
    int part_start;
    /* array of start indices in xs and ys as expected by shapelib */
    int *part_starts = NULL;

    /* generic counter */
    int i;

    /* array of part types. holds the converted content of
     * part_type_list. Stays NULL of part_type_list is NULL
     */
    int *part_types = NULL;

    /* temporary python objects referring to the the list items being
     * worked on.
     */
    PyObject * part = NULL, *tuple = NULL;

    /* The result object */
    SHPObject *result;

    num_parts = PySequence_Length(parts);
    num_vertices = 0; 

    /* parts and part_types have to have the same lengths */
    if (part_type_list
	&& PySequence_Length(parts) != PySequence_Length(part_type_list))
    {
	PyErr_SetString(PyExc_TypeError,
			"parts and part_types have to have the same lengths");
	return NULL;
    }

    /* determine how many vertices there are altogether */
    for (i = 0; i < num_parts; i++)
    {
	PyObject * part = PySequence_GetItem(parts, i);
	if (!part)
	    return NULL;
	num_vertices += PySequence_Length(part);
	Py_DECREF(part);
    }

    /* allocate the memory for the various arrays and check for memory
       errors */
    xs = malloc(num_vertices * sizeof(double));
    ys = malloc(num_vertices * sizeof(double));
    part_starts = malloc(num_parts * sizeof(int));
    if (part_type_list)
	part_types = malloc(num_parts * sizeof(int));

    if (!xs || !ys || !part_starts || (part_type_list && !part_types))
    {
	PyErr_NoMemory();
	goto fail;
    }

    /* convert the part types */
    if (part_type_list)
    {
	for (i = 0; i < num_parts; i++)
	{
	    PyObject * otype = PySequence_GetItem(part_type_list, i);
	    if (!otype)
		return NULL;
	    part_types[i] = PyInt_AsLong(otype);
	    Py_DECREF(otype);
	}
    }

    /* convert the list of parts */
    part_start = 0;
    for (i = 0; i < num_parts; i++)
    {
	int j, length;

	part = PySequence_GetItem(parts, i);
	length = PySequence_Length(part);
	part_starts[i] = part_start;

	for (j = 0; j < length; j++)
	{
	    tuple = PySequence_GetItem(part, j);
	    if (!tuple)
		goto fail;

	    if (!PyArg_ParseTuple(tuple, "dd", xs + part_start + j,
				  ys + part_start + j))
	    {
		goto fail;
	    }
	    Py_DECREF(tuple);
	    tuple = NULL;
	}
	Py_DECREF(part);
	part = NULL;
	part_start += length;
    }

    result = SHPCreateObject(type, id, num_parts, part_starts, part_types,
			     num_vertices, xs, ys, NULL, NULL);
    free(xs);
    free(ys);
    free(part_starts);
    free(part_types);
    return result;

 fail:
    free(xs);
    free(ys);
    free(part_starts);
    free(part_types);
    Py_XDECREF(part);
    Py_XDECREF(tuple);
    return NULL;
}

%}



/*
 * The SWIG Interface definition.
 */

/* include some common SWIG type definitions and standard exception
   handling code */
%include typemaps.i
%include exception.i


/*
 *  SHPObject -- Represents one shape
 */

/* Exception typemap for the SHPObject constructor. The constructor the
   the wrapper function defined above which returns NULL in case of
   error. */
   
%typemap(python,except) SHPObject*new_SHPObject {
    $function;
    if (PyErr_Occurred())
	return NULL;
}

/* Define the SHPObject struct for SWIG. This has to have the same name
 * as the underlying C-struct in shapfil.h, but we don't have to repeat
 * all the fields here, only those we want to access directly, and we
 * can define methods for the object oriented interface.
 */

typedef struct {

    /* The shape object has two read-only attributes: */

    /* The type of the shape. In the c-struct defined the field is
     * called 'nSHPType' but for the python bindings 'type' is more
     * appropriate.
     */
    %readonly %name(type) int nSHPType;

    /* The id of the shape. Here 'id' is a better name than 'nShapeId'. */
    %readonly %name(id) int nShapeId;

    /* The methods */
    %addmethods {

	/* the constructor */
	SHPObject(int type, int id, PyObject * parts,
		  PyObject * part_types = NULL);

	/* The destructor */
	~SHPObject();

	/* extents and vertices correspond to the SHPObject_extents and
	 * SHPObject_vertices defined above
	 */
	PyObject *extents();
	PyObject *vertices();
    }
} SHPObject;


/*
 * ShapeFile --  Represents the shape file
 */

/* Here we do things a little different. We define a new C-struct that
 * holds the SHPHandle. This is mainly done so we can separate the
 * close() method from the destructor but it also helps with exception
 * handling.
 *
 * After the ShapeFile has been opened or created the handle is not
 * NULL. The close() method closes the file and sets handle to NULL as
 * an indicator that the file has been closed.
 */

/* First, define the C-struct */
%{
    typedef struct {
	SHPHandle handle;
    } ShapeFile;
%}

/* define and use some typemaps for the info() method whose
 * C-implementation has four output parameters that are returned through
 * pointers passed into the function. SWIG already has definitions for
 * common types such as int* and we can use those for the first two
 * parameters:
 */
 
%apply int * OUTPUT { int * output_entities }
%apply int * OUTPUT { int * output_type }

/* for the last two, the 4-element arrays of min- and max-values, we
 * have to define our own typemaps:
 */
%typemap (python,ignore) double * extents(double temp[4]) {
    $target = temp;
}

%typemap (python,argout) double * extents {
    PyObject * list = Py_BuildValue("[dddd]",
				    $source[0], $source[1],
				    $source[2], $source[3]);
    $target = t_output_helper($target,list);
}

%apply double * extents { double * output_min_bounds }
%apply double * extents { double * output_max_bounds }

/* The first argument to the ShapeFile methods is a ShapeFile pointer.
 * We have to check whether handle is not NULL in most methods but not
 * all. In the destructor and the close method, it's OK for handle to be
 * NULL. We achieve this by checking whether the preprocessor macro
 * NOCHECK_$name is defined. SWIG replaces $name with the name of the
 * function for which the code is inserted. In the %{,%}-block below we
 * define the macros for the destructor and the close() method.
 */


%typemap(python,check) ShapeFile *{
    %#ifndef NOCHECK_$name
    if (!$target || !$target->handle)
	SWIG_exception(SWIG_TypeError, "shapefile already closed");
    %#endif
}

%{
#define NOCHECK_delete_ShapeFile
#define NOCHECK_ShapeFile_close
%}

/* An exception handle for the constructor and the module level open()
 * and create() functions.
 *
 * Annoyingly, we *have* to put braces around the SWIG_exception()
 * calls, at least in the python case, because of the way the macro is
 * written. Of course, always putting braces around the branches of an
 * if-statement is often considered good practice.
 */
%typemap(python,except) ShapeFile * {
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


/*
 * The SWIG-version of the ShapeFile struct.
 */

typedef struct
{
    /* Only methods and no attributes here: */ 
    %addmethods {

	/* The constructor. Takes two arguments, the filename and the
	 * optinal mode which are passed through to SHPOpen (due to the
	 * renaming trick)
	 */
	ShapeFile(char *file, char * mode = "rb") {
	    ShapeFile * self = malloc(sizeof(ShapeFile));
	    if (self)
		self->handle = SHPOpen(file, mode);
	    return self;
	}

	/* The destructor. Equivalent to SHPClose */
	~ShapeFile() {
	    if (self->handle)
		SHPClose(self->handle);
	    free(self);
	}

	/* close the shape file and set handle to NULL */
	void close() {
	    if (self->handle)
	    {
		SHPClose(self->handle);
		self->handle = NULL;
	    }
	}

	/* info() -- Return a tuple (NUM_SHAPES, TYPE, MIN, MAX) where
	 * NUM_SHAPES is the number of shapes in the file, TYPE is the
	 * shape type and MIN and MAX are 4-element lists with the min.
	 * and max. values of the data.
	 *
	 * The arguments of the underlying shapelib function SHPGetInfo
	 * are all output parameters. To tell SWIG this, we have defined
	 * some typemaps above
	 */
	void info(int * output_entities, int * output_type,
		  double * output_min_bounds, double *output_max_bounds) {
	    SHPGetInfo(self->handle, output_entities, output_type,
		       output_min_bounds, output_max_bounds);
	}

	/* Return object number i */
	%new SHPObject * read_object(int i) {
	    return SHPReadObject(self->handle, i);
	}

	/* Write an object */
	int write_object(int iShape, SHPObject * psObject) {
	    return SHPWriteObject(self->handle, iShape, psObject);
	}

	/* Return the shapelib SHPHandle as a Python CObject */
	PyObject * cobject() {
	    return PyCObject_FromVoidPtr(self->handle, NULL);
	}
    }

} ShapeFile;


/*
 * Two module level functions, open() and create() that correspond to
 * SHPOpen and SHPCreate respectively. open() is equivalent to the
 * ShapeFile constructor.
 */

%{
    ShapeFile * open_ShapeFile(const char *filename, const char * mode) {
	ShapeFile * self = malloc(sizeof(ShapeFile));
	if (self)
	    self->handle = SHPOpen(filename, mode);
	return self;
    }
%}

%name(open) %new ShapeFile *open_ShapeFile(const char *filename,
					   const char * mode = "rb");


%{
    ShapeFile * create_ShapeFile(const char *filename, int type) {
	ShapeFile * self = malloc(sizeof(ShapeFile));
	if (self)
	    self->handle = SHPCreate(filename, type);
	return self;
    }
%}    

%name(create) %new ShapeFile * create_ShapeFile(const char *filename,
						int type);
    

/* Module level function to expose some of the shapelib functions linked
 * with the shapefile C-module to other Python extension modules. This
 * is a kludge to make a Thuban extension work that reads shapes from
 * shapefiles opened by the shapefile module.
 */

%{
    static PyShapeLibAPI the_api = {
	SHPReadObject,
	SHPDestroyObject,
	SHPCreateTree,
	SHPDestroyTree,
	SHPTreeFindLikelyShapes
    };

    PyObject * c_api() {
	return PyCObject_FromVoidPtr(&the_api, NULL);
    }
%}

PyObject * c_api();


/*
 *  Module Level functions 
 */

/* convert shapefile types to names */
%name(type_name) const char *SHPTypeName(int nSHPType);
%name(part_type_name) const char *SHPPartTypeName(int nPartType);


/*
 * Finally, constants copied from shapefil.h
 */

/* -------------------------------------------------------------------- */
/*      Shape types (nSHPType)                                          */
/* -------------------------------------------------------------------- */
#define SHPT_NULL	0
#define SHPT_POINT	1
#define SHPT_ARC	3
#define SHPT_POLYGON	5
#define SHPT_MULTIPOINT	8
#define SHPT_POINTZ	11
#define SHPT_ARCZ	13
#define SHPT_POLYGONZ	15
#define SHPT_MULTIPOINTZ 18
#define SHPT_POINTM	21
#define SHPT_ARCM	23
#define SHPT_POLYGONM	25
#define SHPT_MULTIPOINTM 28
#define SHPT_MULTIPATCH 31


/* -------------------------------------------------------------------- */
/*      Part types - everything but SHPT_MULTIPATCH just uses           */
/*      SHPP_RING.                                                      */
/* -------------------------------------------------------------------- */

#define SHPP_TRISTRIP	0
#define SHPP_TRIFAN	1
#define SHPP_OUTERRING	2
#define SHPP_INNERRING	3
#define SHPP_FIRSTRING	4
#define SHPP_RING	5


