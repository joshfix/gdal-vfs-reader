package com.joshfix.imageio;

import com.joshfix.geotools.reader.GdalVfsImageInputStream;
import com.joshfix.geotools.reader.VfsPath;
import it.geosolutions.imageio.core.CoreCommonIIOStreamMetadata;
import it.geosolutions.imageio.core.GCP;
import it.geosolutions.imageio.gdalframework.GDALCommonIIOImageMetadata;
import it.geosolutions.imageio.gdalframework.GDALImageReader;
import it.geosolutions.imageio.gdalframework.GDALImageReaderSpi;
import it.geosolutions.imageio.gdalframework.GDALUtilities;
import it.geosolutions.imageio.imageioimpl.EnhancedImageReadParam;
import org.gdal.gdal.Band;
import org.gdal.gdal.Dataset;
import org.gdal.gdal.gdal;
import org.gdal.gdalconst.gdalconst;
import org.gdal.gdalconst.gdalconstConstants;

import javax.imageio.ImageReadParam;
import javax.imageio.ImageReader;
import javax.imageio.ImageTypeSpecifier;
import javax.imageio.metadata.IIOMetadata;
import javax.imageio.stream.ImageInputStream;
import java.awt.*;
import java.awt.image.*;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URL;
import java.nio.*;
import java.util.List;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author joshfix
 * Created on 2019-06-24
 */
public class GdalVfsImageReader extends GDALImageReader {

    /**
     * The LOGGER for this class.
     */
    private static final Logger LOGGER = Logger.getLogger(it.geosolutions.imageio.gdalframework.GDALImageReader.class.toString());

    /**
     * list of childs subdatasets names (if any) contained into the source
     */
    private String datasetNames[];

    /**
     * number of subdatasets
     */
    private int nSubdatasets = -1;

    /**
     * The ImageInputStream
     */
    private ImageInputStream imageInputStream;

    private VfsPath vfsPath = null;

    /**
     * {@link HashMap} containing couples (datasetName,
     * {@link GDALCommonIIOImageMetadata}).
     */
    private ConcurrentHashMap<String, GDALCommonIIOImageMetadata> datasetMetadataMap = new ConcurrentHashMap<String, GDALCommonIIOImageMetadata>();

    private ConcurrentHashMap<String, Dataset> datasetsMap = new ConcurrentHashMap<>();

    /**
     * Constructs a <code>GDALImageReader</code> using a
     * {@link GDALImageReaderSpi}.
     *
     * @param originatingProvider The {@link GDALImageReaderSpi} to use for building this
     *                            <code>GDALImageReader</code>.
     */
    public GdalVfsImageReader(final GDALImageReaderSpi originatingProvider) {
        super(originatingProvider);
    }

    /**
     * Constructs a <code>GDALImageReader</code> using a
     * {@link GDALImageReaderSpi}.
     *
     * @param originatingProvider The {@link GDALImageReaderSpi} to use for building this
     *                            <code>GDALImageReader</code>.
     */
    public GdalVfsImageReader(final GDALImageReaderSpi originatingProvider, final int numSubdatasets) {
        super(originatingProvider);
        if (numSubdatasets < 0)
            throw new IllegalArgumentException("The provided number of sub datasets is invalid");
        this.nSubdatasets = numSubdatasets;
    }

    /**
     * Sets the input for the specialized reader.
     *
     * @throws IllegalArgumentException if the provided input is <code>null</code>
     */
    public void setInput(Object input, boolean seekForwardOnly,
                         boolean ignoreMetadata) {
        if (LOGGER.isLoggable(Level.FINE))
            LOGGER.fine("Setting Input");

        // check input
        if (input == null) {
            throw new IllegalArgumentException("The provided input is null!");
        }

        this.seekForwardOnly = seekForwardOnly;
        this.ignoreMetadata = ignoreMetadata;
        this.minIndex = 0;

        // Prior to set a new input, I need to do a pre-emptive reset in order
        // to clear any value-object which was related to the previous input.
        if (this.imageInputStream != null) {
            reset();
            imageInputStream = null;
        }

        String url = null;
        if (input instanceof GdalVfsImageInputStream) {
            vfsPath = ((GdalVfsImageInputStream)input).getVfsPath();
            imageInputStream = (GdalVfsImageInputStream)input;
        }
        /*
        if (input instanceof VfsPath) {
            vfsPath = (VfsPath) input;
        } else if (input instanceof URL || input instanceof URI) {
            url = input.toString();
        } else if (input instanceof String) {
            url = (String) input;
        }
*/
        if (vfsPath == null && url == null) {
            throw new IllegalArgumentException("Input not supported.");
        }

        if (vfsPath == null) {
            if (url.startsWith("wasb")) {
                url = url.replace("wasb://destination@imageryproducts.blob.core.windows.net/", "/vsiaz/destination/");
                vfsPath = new VfsPath(url);
            }
        }

        // is gdal available
        if (!GDALUtilities.isGDALAvailable())
            throw new IllegalStateException("GDAL native libraries are not available.");

        //
        // Checking if this input is of a supported format.
        // Now, I have an ImageInputStream and I can try to see if the input's
        // format is supported by the specialized reader
        //
        boolean isInputDecodable = false;
        String mainDatasetName = null;
        Dataset mainDataSet = null;


        if (vfsPath != null) {
            mainDatasetName = vfsPath.getPath();
            mainDataSet = GDALUtilities.acquireDataSet(vfsPath.getPath(), gdalconstConstants.GA_ReadOnly);
            if (mainDataSet == null) {
                mainDataSet = gdal.Open(vfsPath.getPath());
            }
        }
        if (mainDataSet != null) {
            //isInputDecodable = ((GDALImageReaderSpi) this.getOriginatingProvider()).isDecodable(mainDataSet);
            isInputDecodable = true;
        } else
            isInputDecodable = false;

        if (isInputDecodable) {
            // cache dataset
            datasetsMap.put(mainDatasetName, mainDataSet);

            // input is decodable
            //super.setInput(imageInputStream, seekForwardOnly, ignoreMetadata);

            // Listing available subdatasets
            final java.util.List<String> subdatasets = mainDataSet.GetMetadata_List(GDALUtilities.GDALMetadataDomain.SUBDATASETS);

            // setting the number of subdatasets
            // It is worth to remind that the subdatasets vector
            // contains both Subdataset's Name and Subdataset's Description
            // Thus we need to divide its size by two.
            nSubdatasets = subdatasets.size() / 2;

            // Some formats supporting subdatasets may have no subdatasets.
            // As an instance, the HDF4ImageReader may read HDF4Images
            // which are single datasets containing no subdatasets.
            // Thus, theDataset is simply the main dataset.
            if (nSubdatasets == 0) {
                nSubdatasets = 1;
                datasetNames = new String[1];
                datasetNames[0] = mainDatasetName;
                datasetMetadataMap.put(datasetNames[0], this.createDatasetMetadata(mainDatasetName));

            } else {
                datasetNames = new String[nSubdatasets + 1];
                for (int i = 0; i < nSubdatasets; i++) {
                    final String subdatasetName = (subdatasets.get(i * 2)).toString();
                    final int nameStartAt = subdatasetName.lastIndexOf("_NAME=") + 6;
                    datasetNames[i] = subdatasetName.substring(nameStartAt);
                }
                datasetNames[nSubdatasets] = mainDatasetName;
                datasetMetadataMap.put(datasetNames[nSubdatasets], createDatasetMetadata(mainDataSet, datasetNames[nSubdatasets]));
            }
            // clean list
            subdatasets.clear();

        } else {
            StringBuilder sb = new StringBuilder();
            if (imageInputStream == null) {
                sb.append("Unable to create a valid ImageInputStream for the provided input:");
                sb.append("\n");
                sb.append(input.toString());
            } else
                sb.append("The Provided input is not supported by this reader");
            throw new RuntimeException(sb.toString());
        }
    }

}
