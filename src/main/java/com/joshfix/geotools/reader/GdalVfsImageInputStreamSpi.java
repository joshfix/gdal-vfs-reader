package com.joshfix.geotools.reader;

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

    @Override
    public ImageInputStream createInputStreamInstance(Object input, boolean useCache, File cacheDir) {
        if (input instanceof URL) {
            String url = input.toString();
            if (url.startsWith("wasb")) {
                url = url.replace("wasb://destination@imageryproducts.blob.core.windows.net/", "/vsiaz/destination/");
            }
            return new GdalVfsImageInputStream(new VfsPath(url));
        } else if (input instanceof VfsPath) {
            return new GdalVfsImageInputStream((VfsPath) input);
        } else if (input instanceof String) {
            return new GdalVfsImageInputStream(new VfsPath((String)input));
        }
        throw new UnsupportedOperationException("Input type not supported.");
    }

    @Override
    public String getDescription(Locale locale) {
        return "Placeholder input stream for GDAL VFS reader.";
    }
}
