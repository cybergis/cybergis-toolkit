import pysal
import sys
import os
from osgeo import ogr

total_size = 10000

driver = ogr.GetDriverByName("ESRI Shapefile")

output_shapefiles = []
sizes = [100]
for size in sizes:
    output = "{0}x{0}.shp".format(size)
    output_shapefiles.append("{0}x{0}.shp".format(size))
    if os.path.exists(output): 
        driver.DeleteDataSource(output)
    datasource = driver.CreateDataSource(output)

    if datasource is None:
        print "Could not create file."

    layer = datasource.CreateLayer(output, geom_type=ogr.wkbPolygon)

    #Create all the fields
    field_ID = ogr.FieldDefn()
    field_ID.SetName('ID')
    field_ID.SetType(ogr.OFTReal)
    field_ID.SetWidth(7)
    field_ID.SetPrecision(4)
    layer.CreateField(field_ID) 
    
    field_SAR1 = ogr.FieldDefn()
    field_SAR1.SetName('SAR1')
    field_SAR1.SetType(ogr.OFTReal)
    field_SAR1.SetWidth(7)
    field_SAR1.SetPrecision(4)
    layer.CreateField(field_SAR1)  
    
    field_Uniform2 = ogr.FieldDefn()
    field_Uniform2.SetName('Uniform2')
    field_Uniform2.SetType(ogr.OFTReal)
    field_Uniform2.SetWidth(7)
    field_Uniform2.SetPrecision(4)
    layer.CreateField(field_Uniform2)     

    field_Adjacent = ogr.FieldDefn()
    field_Adjacent.SetName('Adjacent')
    field_Adjacent.SetType(ogr.OFTString)
    field_Adjacent.SetWidth(50)
    layer.CreateField(field_Adjacent)
    
    field_Touching = ogr.FieldDefn()
    field_Touching.SetName('Touching')
    field_Touching.SetType(ogr.OFTString)
    field_Touching.SetWidth(50)
    layer.CreateField(field_Touching)  
    
    field_Area = ogr.FieldDefn()
    field_Area.SetName('Area')
    field_Area.SetType(ogr.OFTReal)
    field_Area.SetWidth(13)
    field_Area.SetPrecision(11)
    layer.CreateField(field_Area)     
    
    field_Perimeter = ogr.FieldDefn()
    field_Perimeter.SetName('Perimeter')
    field_Perimeter.SetType(ogr.OFTReal)
    field_Perimeter.SetWidth(13)
    field_Perimeter.SetPrecision(11)
    layer.CreateField(field_Perimeter) 
    
    field_CentroidX = ogr.FieldDefn()
    field_CentroidX.SetName('CentroidX')
    field_CentroidX.SetType(ogr.OFTReal)
    field_CentroidX.SetWidth(13)
    field_CentroidX.SetPrecision(11)
    layer.CreateField(field_CentroidX)  
    
    field_CentroidY = ogr.FieldDefn()
    field_CentroidY.SetName('CentroidY')
    field_CentroidY.SetType(ogr.OFTReal)
    field_CentroidY.SetWidth(13)
    field_CentroidY.SetPrecision(11)
    layer.CreateField(field_CentroidY)  
    
    field_MI = ogr.FieldDefn()
    field_MI.SetName('MI')
    field_MI.SetType(ogr.OFTReal)
    field_MI.SetWidth(13)
    field_MI.SetPrecision(11)
    layer.CreateField(field_MI)    
    
    featureDefn = layer.GetLayerDefn()
    step = total_size / size
    counter = 0
    for x1 in range(0,total_size,step):
        x2 = x1 + step
        for y1 in range(0,total_size,step):
            y2 = y1 + step
            
            ring = ogr.Geometry(ogr.wkbLinearRing)
            
            ring.AddPoint(x1, y1)
            ring.AddPoint(x1, y2)
            ring.AddPoint(x2, y2)
            ring.AddPoint(x2, y1)
            ring.AddPoint(x1, y1)
            
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)  

            feature = ogr.Feature(featureDefn)
            feature.SetGeometry(poly)            
            
            feature.SetField('SAR1', 0)
            feature.SetField('Uniform2',0)
            feature.SetField('Adjacent', 'None')
            feature.SetField('Touching', 'Unneeded')
            feature.SetField('Area', (x2-x1)*(y2-y1))
            feature.SetField('Perimeter', (x2-x1)*4)
            feature.SetField('CentroidX', x2-(step/2))
            feature.SetField('CentroidY', y2-(step/2))
            feature.SetField('MI', 18**4 / 6)
            feature.SetField('ID', counter)
            counter += 1
            layer.CreateFeature(feature)
poly.Destroy()
feature.Destroy()
output=None

for shapefile in output_shapefiles:
    w = pysal.rook_from_shapefile(shapefile)
    print w
    infile = ogr.Open(shapefile, 1)
    inlyr = infile.GetLayerByIndex(0)
    counter = 0
    feat = inlyr.GetNextFeature()
    print feat
    while feat is not None:
        feat.SetField('Adjacent', str(w.neighbors[counter]).split("[")[1].split("]")[0])
        print str(w.neighbors[counter]).split("[")[1].split("]")[0]
        counter += 1
        inlyr.SetFeature(feat)
        feat = inlyr.GetNextFeature()

