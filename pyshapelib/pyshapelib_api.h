/* Header file for the PyShapelib API for other Python modules */
/* $Revision$ */

#ifndef PYSHAPELIB_API_H
#define PYSHAPELIB_API_H

typedef struct {
    /* Shapefile functions */
    SHPObject * (*SHPReadObject)(SHPHandle hSHP, int iShape);
    void (*SHPDestroyObject)(SHPObject * psObject);

    /* SHPTree functions */
    SHPTree * (*SHPCreateTree)(SHPHandle hSHP, int nDimension, int nMaxDepth,
			       double *padfBoundsMin, double *padfBoundsMax);
    void (*SHPDestroyTree)(SHPTree * hTree);
    int * (*SHPTreeFindLikelyShapes)(SHPTree * hTree, double * padfBoundsMin,
				     double * padfBoundsMax, int *);
} PyShapeLibAPI;


/* Macro to import the shapelib module, extract the API pointer and
 * assign it to the variable given as argument */
#define PYSHAPELIB_IMPORT_API(apivariable)				   \
{									   \
    PyObject * shapelib = PyImport_ImportModule("shapelibc");		   \
    if (shapelib)							   \
    {									   \
	PyObject * c_api_func = PyObject_GetAttrString(shapelib, "c_api"); \
	if (c_api_func)							   \
	{								   \
	    PyObject * cobj = PyObject_CallObject(c_api_func, NULL);	   \
	    if (cobj)							   \
	    {								   \
		(apivariable) = (PyShapeLibAPI*)PyCObject_AsVoidPtr(cobj); \
	    }								   \
	}								   \
    }									   \
}


#endif /* PYSHAPELIB_API_H */
