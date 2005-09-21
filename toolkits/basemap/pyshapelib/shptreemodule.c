/* Copyright (c) 2001, 2002 by Intevation GmbH
 * Authors:
 * Bernhard Herzog <bh@intevation.de>
 *
 * This program is free software under the GPL (>=v2)
 * Read the file COPYING coming with Thuban for details.
 */

/* Python wrapper for the shapelib SHPTree */

#include <Python.h>
#include <shapefil.h>

#include "pyshapelib_api.h"

PyShapeLibAPI * api;

typedef struct {
    PyObject_HEAD
    SHPTree * tree; 
} SHPTreeObject;

extern PyTypeObject SHPTreeType;

#define SHPTree_Check(v)		((v)->ob_type == &SHPTreeType)

/* Create a new python wrapper object from a SHPTree pointer */
static PyObject *
SHPTreeObject_FromSHPTree(SHPTree* tree)
{
    SHPTreeObject * self = PyObject_NEW(SHPTreeObject, &SHPTreeType);
    if (!self)
	return NULL;

    self->tree = tree;

    return (PyObject *)self;
}

/* Deallocate the SHPTree wrapper. */
static void
shptree_dealloc(SHPTreeObject * self)
{
    api->SHPDestroyTree(self->tree);
    PyMem_DEL(self);
}

/* Return the repr of the wrapper */
static PyObject *
shptree_repr(SHPTreeObject * self)
{
    char buf[1000];
    sprintf(buf, "<SHPTree at %xul>", (unsigned long)self);
    return PyString_FromString(buf);
}

static PyObject *
shptree_find_shapes(SHPTreeObject * self, PyObject * args)
{
    double min[4] = {0, 0, 0, 0};
    double max[4] = {0, 0, 0, 0};
    int count, idx;
    int * ids;
    PyObject * list = NULL, *temp = NULL;

    if (!PyArg_ParseTuple(args, "(dd)(dd)", min + 0, min + 1,
			  max + 0, max + 1))
	return NULL;

    ids = api->SHPTreeFindLikelyShapes(self->tree, min, max, &count);

    list = PyList_New(count);
    if (!list)
	goto fail;

    /* Turn the returned array of indices into a python list of ints. */
    for (idx = 0; idx < count; idx++)
    {
	temp = PyInt_FromLong(ids[idx]);
	if (!temp)
	    goto fail;

	if (PyList_SetItem(list, idx, temp) == -1)
	{
	    /* temp's refcount has already be decreased. Set temp to
	     * NULL so that the fail code doesn't do it again
	     */
	    temp = NULL;
	    goto fail;
	}
    }

    free(ids);
    return list;
    
 fail:
    free(ids);
    Py_XDECREF(list);
    Py_XDECREF(temp);
    return NULL;
}


static struct PyMethodDef shptree_methods[] = {
    {"find_shapes",	(PyCFunction)shptree_find_shapes,	METH_VARARGS},
    {NULL,	NULL}
};

static PyObject *
shptree_getattr(PyObject * self, char * name)
{
    return Py_FindMethod(shptree_methods, self, name);
}


PyTypeObject SHPTreeType = {
	PyObject_HEAD_INIT(NULL)
	0,
	"SHPTree",
	sizeof(SHPTreeObject),
	0,
	(destructor)shptree_dealloc,	/*tp_dealloc*/
	(printfunc)NULL,		/*tp_print*/
	shptree_getattr,		/*tp_getattr*/
	0,				/*tp_setattr*/
	0,				/*tp_compare*/
	(reprfunc)shptree_repr,		/*tp_repr*/
	0,				/*tp_as_number*/
	0,				/*tp_as_sequence*/
	0,				/*tp_as_mapping*/
	0,				/*tp_hash*/
        0,				/*tp_call*/
        0,				/*tp_str*/
	0,				/*tp_getattro*/
	0,				/*tp_setattro*/
	0,				/*tp_as_buffer*/
};



static PyObject *
shptree_from_shapefile(PyObject * self, PyObject * args)
{
    SHPTree * tree;
    SHPHandle handle;
    PyObject * cobject;
    int dimension, max_depth;

    if (!PyArg_ParseTuple(args, "O!ii", &PyCObject_Type, &cobject,
			  &dimension, &max_depth))
	return NULL;

    handle = PyCObject_AsVoidPtr(cobject);

    tree = api->SHPCreateTree(handle, dimension, max_depth, NULL, NULL);

    /* apparently SHPCreateTree doesn't do any error checking, so we
     * have to assume that tree is valid at this point. */
    return SHPTreeObject_FromSHPTree(tree);
}


static PyMethodDef module_functions[] = {
    {"SHPTree",		shptree_from_shapefile,		METH_VARARGS},
    { NULL, NULL }
};


void
initshptree()
{
    SHPTreeType.ob_type = &PyType_Type;

    Py_InitModule("shptree", module_functions);
    PYSHAPELIB_IMPORT_API(api);
}
