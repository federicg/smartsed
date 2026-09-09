import processing
import sys
import os


# To run the script you must set PYTHONPATH and path to QGIS binaries
#
# on Linux: export PYTHONPATH=/qgispath/share/qgis/python
# on Windows: set PYTHONPATH=c:\qgispath\python
# on OSX: export PYTHONPATH=/Applications/QGIS.app/Contents/Resources/python
#         export PATH="/Applications/QGIS.app/Contents/MacOS/bin:$PATH"


#os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = '/usr/local/Cellar/qt/5.13.2/plugins/platforms/'
#os.environ['PATH'] = '/Applications/QGIS3.10.app/Contents/MacOS/bin:/usr/local/Cellar/qt/5.13.2/bin:' + os.environ['PATH']



print("QT_QPA_PLATFORM_PLUGIN_PATH:", os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'])
print("PATH:", os.environ['PATH'])


sys.path.append("/Applications/QGIS3.10.app/Contents/Resources/python")
from qgis.core import *
import urllib.parse

QgsApplication.setPrefixPath("/Applications/QGIS3.10.app/Contents/MacOS", True)

# Initialize QGIS Application
qgs = QgsApplication([], False) # if True enables GUI
#qgs.initQgis()

rlayername = 'clay'
uri = QgsDataSourceUri()
uri.setParam('url', 'https://maps.isric.org/mapserv?map=/map/clay.map')
uri.setParam("identifier", rlayername)
rlayer = QgsRasterLayer(str(uri.encodedUri()), 'clay_0-5cm_mean', 'wcs')


#qgs.exitQgis()
exit()

### Raster layer from WCS
rlayername = 'clay'
uri = QgsDataSourceUri()
uri.setParam('url', 'https://maps.isric.org/mapserv?map=/map/clay.map')
uri.setParam("identifier", rlayername)
rlayer = QgsRasterLayer(str(uri.encodedUri()), 'clay_0-5cm_mean', 'wcs')


exit()


# Make WCS Uri
def makeWCSuri( url, layer ):
    params = {  'dpiMode': 7 ,
                'identifier': layer,
                'url': url.split('?')[0]  }

    uri = urllib.parse.unquote( urllib.parse.urlencode(params)  )
    return uri

### Raster layer from WCS
wcsUri = makeWCSuri('https://maps.isric.org/mapserv?map=/map/clay.map', 'clay_0-5cm_mean' )
print(wcsUri)

rlayer = QgsRasterLayer(wcsUri, 'clay_0-5cm_mean', 'wcs')
if not rlayer.isValid():
  print ("Layer failed to load!")

#QgsProject.instance().addMapLayer(rlayer)

mask_filepath = r'D:\temp\mask2.shp'
mask_layer = QgsVectorLayer(mask_filepath, 'mask', 'ogr')
#QgsProject.instance().addMapLayer(mask_layer)

# Save raster
renderer = rlayer.renderer()
provider = rlayer.dataProvider()
crs = rlayer.crs()



pipe = QgsRasterPipe()
projector = QgsRasterProjector()
projector.setCrs(provider.crs(), provider.crs())

if not pipe.set(provider.clone()):
    print("Cannot set pipe provider")

# Commented for extract raw data
# if not pipe.set(renderer.clone()):
    # print("Cannot set pipe renderer")

if not pipe.insert(2, projector):
    print("Cannot set pipe projector")

out_file = 'D:/temp/temporal.tif'
file_writer = QgsRasterFileWriter(out_file)
file_writer.Mode(1)

print ("Saving")

extent = mask_layer.extent()

opts = ["COMPRESS=LZW"]
file_writer.setCreateOptions(opts)
error = file_writer.writeRaster(
    pipe,
    extent.width (),
    extent.height(),
    extent,
    crs)

if error == QgsRasterFileWriter.NoError:
    print ("Raster was saved successfully!")
    layer = QgsRasterLayer(out_file, "result")
    QgsProject.instance().addMapLayer(layer)
else:
    print ("Raster was not saved!")

