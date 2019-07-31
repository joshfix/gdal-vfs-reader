package com.joshfix.geotools.reader;

import com.sun.media.imageioimpl.common.PackageUtil;
import it.geosolutions.imageio.utilities.ImageIOUtilities;
import it.geosolutions.imageioimpl.plugins.tiff.TIFFImageMetadata;
import it.geosolutions.imageioimpl.plugins.tiff.TIFFImageReader;
import it.geosolutions.imageioimpl.plugins.tiff.TIFFImageReaderSpi;
import it.geosolutions.imageioimpl.plugins.tiff.TIFFStreamMetadata;

import javax.imageio.ImageReader;
import javax.imageio.spi.ImageReaderSpi;
import javax.imageio.spi.ImageReaderWriterSpi;
import javax.imageio.spi.ImageWriterSpi;
import javax.imageio.spi.ServiceRegistry;
import javax.imageio.stream.ImageInputStream;
import java.io.IOException;
import java.util.Iterator;
import java.util.Locale;

/**
 * @author joshfix
 * Created on 2019-07-31
 */
public class GdalVfsImageReaderSpi extends ImageReaderSpi {

    private static final String[] names = { "tif", "TIF", "tiff", "TIFF", "btiff", "BTIFF" };

    private static final String[] suffixes = { "tif", "tiff", "tf8", "btf", "TIF", "TIFF", "TF8", "BTF" };

    private static final String[] MIMETypes = { "image/tiff" };

    private static final String readerClassName = GdalVfsReader.class.getCanonicalName();

    private static final String[] writerSpiNames = {
            "it.geosolutions.imageioimpl.plugins.tiff.TIFFImageWriterSpi"
    };

    static final String nativeMetadataFormatName =
            "com_sun_media_imageio_plugins_tiff_stream_1.0";

    private boolean registered = false;

    public GdalVfsImageReaderSpi() {
        super("ImageIO-Ext",
                "1.0",
                names,
                suffixes,
                MIMETypes,
                readerClassName,
                STANDARD_INPUT_TYPE,
                writerSpiNames,
                false,
                nativeMetadataFormatName,
                "it.geosolutions.imageioimpl.plugins.tiff.TIFFStreamMetadataFormat",
                null, null,
                true,
                TIFFImageMetadata.nativeMetadataFormatName,
                "it.geosolutions.imageioimpl.plugins.tiff.TIFFImageMetadataFormat",
                new String[]{""}, new String[]{""}
        );
    }

    public String getDescription(Locale locale) {
        String desc = PackageUtil.getSpecificationTitle() +
                " GDAL VFS Image Reader";
        return desc;
    }

    public boolean canDecodeInput(Object input) throws IOException {
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
    }

    public ImageReader createReaderInstance(Object extension) {
        return new TIFFImageReader(this);
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
