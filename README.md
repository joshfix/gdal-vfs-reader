# GDAL VFS Reader for GeoTools

Early stages of a GridCoverage2D reader for GeoTools 
* Accepts Strings as input, namely Strings representing 
GDAL [Virtual File System (VFS)](https://gdal.org/user/virtual_file_systems.html) paths
* Uses GDAL (not a Java input stream or otherwise) to perform the actual reads
* Supports Cloud Optimized Geotiff range reads for the requested coverage area

This project contains very little custom code.  It consists mostly of snippets taken from various places in imageio-ext 
and GeoTools, reconstructed and rearranged to achieve the goal of using the custom 
[VfsPath](./src/main/java/com/joshfix/gdalvfs/geotools/path/VfsPath.java) objects as the source.

The project has only been tested with GeoTIFFs (so far).  There seem to still be some transparency/nodata issues, 
especially when used with ImageMosaic, which are likely due to missing GeoTools settings.
