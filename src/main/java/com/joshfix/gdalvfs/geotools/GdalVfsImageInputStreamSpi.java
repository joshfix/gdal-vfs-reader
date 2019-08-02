package com.joshfix.gdalvfs.geotools;

import com.joshfix.gdalvfs.geotools.path.VfsPath;

import javax.imageio.spi.ImageInputStreamSpi;
import javax.imageio.stream.ImageInputStream;
import java.io.File;
import java.net.URL;
import java.util.Locale;

/**
 * @author joshfix
 * Created on 2019-08-01
 */
public class GdalVfsImageInputStreamSpi extends ImageInputStreamSpi {

    /**
     * ImageReadCRIF.getImageReader method establishes a GdalVfsImageInputStream instance as the input
     * to ImageIO.createImageInputStream, which then checks "spi.getInputClass().isInstance(input)" to
     * verify it has the right SPI.  Something seems not right about this.
     */
    private static final Class<GdalVfsImageInputStream> inputClass = GdalVfsImageInputStream.class;
    private static final String vendorName = "Josh Fix";
    private static final String version = "1.0";


    public GdalVfsImageInputStreamSpi() {
        super(vendorName, version, inputClass);
    }

    /**
     * Returns an instance of the ImageInputStream implementation associated
     * with this service provider.
     *
     * @param input
     *            an object of the class type returned by getInputClass.
     * @param useCache
     *            a boolean indicating whether a cache eraf should be used, in
     *            cases where it is optional.
     *
     * @param cacheDir
     *            a File indicating where the cache eraf should be created, or
     *            null to use the system directory.
     *
     *
     * @return an ImageInputStream instance.
     *
     * @throws IllegalArgumentException
     *             if input is not an instance of the correct class or is null.
     */
    @Override
    public ImageInputStream createInputStreamInstance(Object input, boolean useCache, File cacheDir) {
        if (input instanceof URL) {
            String url = input.toString();
            // TODO figure out how to handle external overviews
            if (url.endsWith(".ovr")) {
                return null;
            }
            return new GdalVfsImageInputStream(new VfsPath(url));
        } else if (input instanceof VfsPath) {
            return new GdalVfsImageInputStream((VfsPath) input);
        } else if (input instanceof String) {
            return new GdalVfsImageInputStream(new VfsPath((String)input));
        } else if (input instanceof GdalVfsImageInputStream) {
            return (GdalVfsImageInputStream)input;
        }
        throw new UnsupportedOperationException("Input type not supported.");
    }

    @Override
    public String getDescription(Locale locale) {
        return "Placeholder input stream for GDAL VFS reader.";
    }
}
