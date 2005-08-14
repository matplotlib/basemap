import shapelib, dbflib, shptree

#
#       The the shapefile module
#

def make_shapefile(filename):
    # Create a shapefile with polygons
    outfile = shapelib.create(filename, shapelib.SHPT_POLYGON)

    # Create one very simple polygon and write it to the shapefile.  The
    # vertices should be given in clockwise order to comply with the
    # shapefile specification.
    obj = shapelib.SHPObject(shapelib.SHPT_POLYGON, 1,
                             [[(10, 10), (10, 20), (20, 20), (10, 10)]])
    print obj.extents()
    print obj.vertices()
    outfile.write_object(-1, obj)

    # Create a polygon with a hole.  Note that according to the
    # shapefile specification, the vertices of the outer ring have to be
    # in clockwise order and the inner rings have to be in counter
    # clockwise order.
    #
    # There's an optional fourth parameter which when given must be a
    # list of part types, one for each part of the shape.  For polygons,
    # the part type is always shapelib.SHPP_RING, though.  The part
    # types are only relevant for SHPT_MULTIPATCH shapefiles.
    obj = shapelib.SHPObject(shapelib.SHPT_POLYGON, 1,
                             [[(0, 0), (0, 40), (40, 40), (40, 0), (0, 0)],
                              [(10, 10), (20, 10), (20, 20), (10, 20),(10, 10)],
                              ])
    print obj.extents()
    print obj.vertices()
    outfile.write_object(-1, obj)

    # close the file.
    outfile.close()

def read_shapefile(filename):
    # open the shapefile
    shp = shapelib.ShapeFile(filename)

    # the info method returns a tuple (num_shapes, type, min, max) where
    # num_shapes is the number of shapes, type is the type code (one of
    # the SHPT* constants defined in the shapelib module) and min and
    # max are 4-element lists with the min. and max. values of the
    # vertices.
    print shp.info()

    # read_object reads a shape
    obj = shp.read_object(0)

    # The vertices method returns the shape as a list of lists of tuples.
    print obj.vertices()[0][:10]

    # The extents returns a tuple with two 4-element lists with the min.
    # and max. values of the vertices.
    print obj.extents()

    # The type attribute is the type code (one of the SHPT* constants
    # defined in the shapelib module)
    print obj.type

    # The id attribute is the shape id
    print obj.id

    # the cobject method returns a PyCObject containing the shapelib
    # SHPHandle. This is useful for passing shapefile objects to
    # C-Python extensions.
    print shp.cobject()

    # build a quad tree from the shapefile. The first argument must be
    # the return value of the shape file object's cobject method (this
    # is currently needed to access the shape file at the C-level). The
    # second argument is the dimension and the third the maximum depth.
    # 0 means to guess an appropriate depth
    tree = shptree.SHPTree(shp.cobject(), 2, 0)

    # Retrieve the ids for a region. Here we just use the extents of the
    # object previously read from the shapefile
    minima, maxima = obj.extents()
    print tree.find_shapes(minima[:2], maxima[:2])


make_shapefile("testfile")
read_shapefile("testfile")

#
#       Test the DBF file module.
#

def make_dbf(file):
    # create a new dbf file and add three fields.
    dbf = dbflib.create(file)
    dbf.add_field("NAME", dbflib.FTString, 20, 0)
    dbf.add_field("INT", dbflib.FTInteger, 10, 0)
    dbf.add_field("FLOAT", dbflib.FTDouble, 10, 4)

def add_dbf_records(file):
    # add some records to file
    dbf = dbflib.open(file, "r+b")
    # Records can be added as a dictionary...
    dbf.write_record(0, {'NAME': "Weatherwax", "INT":1, "FLOAT":3.1415926535})
    # ... or as a sequence
    dbf.write_record(1, ("Ogg", 2, -1000.1234))

def list_dbf(file):
    # print the contents of a dbf file to stdout
    dbf = dbflib.DBFFile(file)
    print "%d records, %d fields" % (dbf.record_count(), dbf.field_count())
    format = ""
    for i in range(dbf.field_count()):
        type, name, len, decc = dbf.field_info(i)
        if type == 0:
            format = format + " %%(%s)%ds" % (name, len)
        elif type == 1:
            format = format + " %%(%s)%dd" % (name, len)
        elif type == 2:
            format = format + " %%(%s)%dg" % (name, len)
    print format
    for i in range(dbf.record_count()):
        print format % dbf.read_record(i)


make_dbf("testfile")
add_dbf_records("testfile")
list_dbf("testfile")
