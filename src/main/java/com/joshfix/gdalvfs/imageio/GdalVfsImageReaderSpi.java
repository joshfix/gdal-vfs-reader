package com.joshfix.gdalvfs.imageio;

import com.joshfix.gdalvfs.geotools.GdalVfsImageInputStream;
import com.joshfix.gdalvfs.geotools.GdalVfsReader;
import com.joshfix.gdalvfs.geotools.path.VfsPath;
import com.sun.media.imageioimpl.common.PackageUtil;
import it.geosolutions.imageio.gdalframework.GDALImageReaderSpi;
import it.geosolutions.imageio.utilities.ImageIOUtilities;

import javax.imageio.ImageReader;
import javax.imageio.spi.ImageReaderSpi;
import javax.imageio.spi.ImageReaderWriterSpi;
import javax.imageio.spi.ImageWriterSpi;
import javax.imageio.spi.ServiceRegistry;
import java.io.IOException;
import java.net.URL;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author joshfix
 * Created on 2019-07-31
 */
public class GdalVfsImageReaderSpi extends GDALImageReaderSpi {

    private static final String[] names = { "tif", "TIF", "tiff", "TIFF", "btiff", "BTIFF" };

    private static final String[] suffixes = { "tif", "tiff", "tf8", "btf", "TIF", "TIFF", "TF8", "BTF" };

    private static final String[] MIMETypes = { "image/tiff" };

    private static final String readerClassName = GdalVfsReader.class.getCanonicalName();

    static final String version = "1.0";

    static final String vendorName = "Josh Fix";
    // writerSpiNames
    static final String[] wSN = { null };

    // StreamMetadataFormatNames and StreamMetadataFormatClassNames
    static final boolean supportsStandardStreamMetadataFormat = false;

    static final String nativeStreamMetadataFormatName = null;

    static final String nativeStreamMetadataFormatClassName = null;

    static final String[] extraStreamMetadataFormatNames = { null };

    static final String[] extraStreamMetadataFormatClassNames = { null };

    // ImageMetadataFormatNames and ImageMetadataFormatClassNames
    static final boolean supportsStandardImageMetadataFormat = false;

    static final String nativeImageMetadataFormatName = null;

    static final String nativeImageMetadataFormatClassName = null;

    static final String[] extraImageMetadataFormatNames = { null };

    static final String[] extraImageMetadataFormatClassNames = { null };

    static final String readerCN = "com.joshfix.gdalvfs.imageio.GdalVfsImageReader";

    private boolean registered = false;

    private static final Logger LOGGER = Logger
            .getLogger(GdalVfsImageReaderSpi.class.getCanonicalName());

    public GdalVfsImageReaderSpi() {
        super(
                vendorName,
                version,
                names,
                suffixes,
                MIMETypes,
                readerCN, // readerClassName
                new Class[] { URL.class, VfsPath.class, String.class, GdalVfsImageInputStream.class},
                wSN, // writer Spi Names
                supportsStandardStreamMetadataFormat,
                nativeStreamMetadataFormatName,
                nativeStreamMetadataFormatClassName,
                extraStreamMetadataFormatNames,
                extraStreamMetadataFormatClassNames,
                supportsStandardImageMetadataFormat,
                nativeImageMetadataFormatName,
                nativeImageMetadataFormatClassName,
                extraImageMetadataFormatNames,
                extraImageMetadataFormatClassNames, Collections
                        .singletonList("JPEG"));

        if (LOGGER.isLoggable(Level.FINE))
            LOGGER.fine("GdalVfsImageReaderSpi Constructor");
    }

    public GdalVfsImageReaderSpi(String vendorName, String version, String[] names, String[] suffixes, String[] MIMETypes, String readerClassName, Class<?>[] inputTypes, String[] writerSpiNames, boolean supportsStandardStreamMetadataFormat, String nativeStreamMetadataFormatName, String nativeStreamMetadataFormatClassName, String[] extraStreamMetadataFormatNames, String[] extraStreamMetadataFormatClassNames, boolean supportsStandardImageMetadataFormat, String nativeImageMetadataFormatName, String nativeImageMetadataFormatClassName, String[] extraImageMetadataFormatNames, String[] extraImageMetadataFormatClassNames, Collection<String> supportedFormats) {
        super(vendorName, version, names, suffixes, MIMETypes, readerClassName, inputTypes, writerSpiNames, supportsStandardStreamMetadataFormat, nativeStreamMetadataFormatName, nativeStreamMetadataFormatClassName, extraStreamMetadataFormatNames, extraStreamMetadataFormatClassNames, supportsStandardImageMetadataFormat, nativeImageMetadataFormatName, nativeImageMetadataFormatClassName, extraImageMetadataFormatNames, extraImageMetadataFormatClassNames, supportedFormats);
    }

    public String getDescription(Locale locale) {
        String desc = PackageUtil.getSpecificationTitle() +
                "GDAL VFS Image Reader";
        return desc;
    }

    public boolean canDecodeInput(Object input) throws IOException {
        return true;
        /*
        if (!(input instanceof ImageInputStream)) {
            return false;
        }

        ImageInputStream stream = (ImageInputStream)input;
        byte[] b = new byte[4];
        stream.mark();
        stream.readFully(b);
        stream.reset();

        return (
                ((b[0] == (byte)0x49 && b[1] == (byte)0x49 &&
                        b[2] == (byte)0x2a && b[3] == (byte)0x00) ||
                        (b[0] == (byte)0x4d && b[1] == (byte)0x4d &&
                                b[2] == (byte)0x00 && b[3] == (byte)0x2a))||

                        ((b[0] == (byte)0x49 && b[1] == (byte)0x49 &&
                                b[2] == (byte)0x2b && b[3] == (byte)0x00) ||
                                (b[0] == (byte)0x4d && b[1] == (byte)0x4d &&
                                        b[2] == (byte)0x00 && b[3] == (byte)0x2b))
        );
        */
    }

    public ImageReader createReaderInstance(Object extension) {
        return new GdalVfsImageReader(this, 1);
    }

    @Override
    public void onRegistration(ServiceRegistry registry, Class category) {
        super.onRegistration(registry, category);
        if (registered) {
            return;
        }
        registered = true;
        Iterator<ImageReaderWriterSpi> readers = ImageIOUtilities
                .getImageReaderWriterSPI(registry, new TIFFFilter(true), "TIFF", true).iterator();
        while (readers.hasNext()) {
            final ImageReaderSpi spi = (ImageReaderSpi) readers.next();
            if (spi == this) {
                continue;
            }
            registry.deregisterServiceProvider(spi);
            registry.setOrdering(category, this, spi);
        }
    }

    /**
     * Filter which returns <code>true</code> if and only if the provider is an ImageReader/WriterSpi which supports the TIFF format.
     */
    static class TIFFFilter implements ServiceRegistry.Filter {
        boolean isReader;

        TIFFFilter(boolean isReader) {
            this.isReader = isReader;
        }

        public boolean filter(Object provider) {
            boolean isSupportedSpi = isReader ? provider instanceof ImageReaderSpi
                    : provider instanceof ImageWriterSpi;
            if (!isSupportedSpi) {
                return false;
            }

            ImageReaderWriterSpi spi = (ImageReaderWriterSpi) provider;
            String[] formatNames = spi.getFormatNames();
            for (int i = 0; i < formatNames.length; i++) {
                if (formatNames[i].equalsIgnoreCase("TIFF")) {
                    return true;
                }
            }

            return false;
        }
    }
}
