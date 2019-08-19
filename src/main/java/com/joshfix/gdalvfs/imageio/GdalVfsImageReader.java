package com.joshfix.gdalvfs.imageio;

import com.joshfix.gdalvfs.geotools.GdalVfsImageInputStream;
import com.joshfix.gdalvfs.geotools.path.VfsPath;
import it.geosolutions.imageio.gdalframework.GDALCommonIIOImageMetadata;
import it.geosolutions.imageio.gdalframework.GDALImageReader;
import it.geosolutions.imageio.gdalframework.GDALImageReaderSpi;
import it.geosolutions.imageio.gdalframework.GDALUtilities;
import it.geosolutions.imageio.imageioimpl.EnhancedImageReadParam;
import org.gdal.gdal.Band;
import org.gdal.gdal.Dataset;
import org.gdal.gdal.TranslateOptions;
import org.gdal.gdal.gdal;
import org.gdal.gdalconst.gdalconst;
import org.gdal.gdalconst.gdalconstConstants;

import javax.imageio.ImageReadParam;
import javax.imageio.ImageTypeSpecifier;
import javax.imageio.metadata.IIOMetadata;
import javax.imageio.stream.ImageInputStream;
import java.awt.*;
import java.awt.image.*;
import java.io.IOException;
import java.nio.*;
import java.util.*;
import java.util.List;
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
    private static final Logger LOGGER = Logger.getLogger(GdalVfsImageReader.class.toString());

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

    private VfsPath vfsPath;
/*
    private Map<String, String> translatedNames = new HashMap<>();
    private Map<String, Dataset> translatedDatasets = new HashMap<>();
    private Map<String, GDALCommonIIOImageMetadata> translatedMetadata = new HashMap<>();
*/
    /**
     * {@link HashMap} containing couples (datasetName,
     * {@link GDALCommonIIOImageMetadata}).
     */
    private ConcurrentHashMap<String, GDALCommonIIOImageMetadata> datasetMetadataMap = new ConcurrentHashMap<>();

    private ConcurrentHashMap<String, Dataset> datasetsMap = new ConcurrentHashMap<>();

    /**
     * Constructs a <code>GDALImageReader</code> using a
     * {@link GDALImageReaderSpi}.
     *
     * @param originatingProvider The {@link GDALImageReaderSpi} to use for building this
     *                            <code>GDALImageReader</code>.
     */
    public GdalVfsImageReader(final GDALImageReaderSpi originatingProvider) {
        super(originatingProvider, 1);
    }

    /**
     * Constructs a <code>GDALImageReader</code> using a
     * {@link GDALImageReaderSpi}.
     *
     * @param originatingProvider The {@link GDALImageReaderSpi} to use for building this
     *                            <code>GDALImageReader</code>.
     */
    public GdalVfsImageReader(final GDALImageReaderSpi originatingProvider, final int numSubdatasets) {
        super(originatingProvider, numSubdatasets);
        if (numSubdatasets < 0)
            throw new IllegalArgumentException("The provided number of sub datasets is invalid");
        this.nSubdatasets = numSubdatasets;
    }

    public IIOMetadata getImageMetadata(int imageIndex) throws IOException {
        return this.getDatasetMetadata(imageIndex);
    }

    /**
     * Retrieves a {@link GDALCommonIIOImageMetadata} by index.
     *
     * @param imageIndex
     *                is the index of the required
     *                {@link GDALCommonIIOImageMetadata}.
     * @return a {@link GDALCommonIIOImageMetadata}
     */
    public GDALCommonIIOImageMetadata getDatasetMetadata(final int imageIndex) {
        checkImageIndex(imageIndex);
        // getting dataset name
        final String datasetName = datasetNames[imageIndex];

        //String translatedName = translatedNames.get(datasetName);
        GDALCommonIIOImageMetadata retVal = datasetMetadataMap.get(datasetName);
        //GDALCommonIIOImageMetadata retVal = translatedMetadata.get(translatedName);
        if (retVal == null) {
            // do we need to create a dataset
            Dataset ds = datasetsMap.get(datasetName);
/*
            String translateName = "/vsimem/" + UUID.randomUUID();
            Dataset dst = translate(ds, translateName);
            translatedDatasets.put(translateName, dst);
            translatedNames.put(datasetNames[imageIndex], translateName);
*/
            if (ds == null) {
                ds = GDALUtilities.acquireDataSet(datasetName, gdalconst.GA_ReadOnly);
                Dataset dsOld = datasetsMap.putIfAbsent(datasetName, ds);
                if (dsOld != null) {
                    // abandon the DataSet we created
                    GDALUtilities.closeDataSet(ds);
                    ds = dsOld;
                }
            }

            // Add a new GDALCommonIIOImageMetadata to the HashMap
            final GDALCommonIIOImageMetadata datasetMetadataNew = createDatasetMetadata(datasetName);
            retVal = datasetMetadataMap.put(datasetName, datasetMetadataNew);
/*
            GDALCommonIIOImageMetadata tm = createDatasetMetadata(translateName);
            retVal = translatedMetadata.put(translateName, tm);
*/
            if (retVal == null) {
                retVal = datasetMetadataNew;
                //retVal = tm;
            }
        }
        return retVal;
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
            this.input = input;
        }

        if (vfsPath == null && url == null) {
            throw new IllegalArgumentException("Input not supported.");
        }

        if (vfsPath == null) {
            vfsPath = new VfsPath(url);
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
            try {
                isInputDecodable = this.getOriginatingProvider().canDecodeInput(mainDataSet);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if (isInputDecodable) {
            // cache dataset
            datasetsMap.put(mainDatasetName, mainDataSet);
/*
            String translateName = "/vsimem/" + UUID.randomUUID();
            Dataset dst = translate(mainDataSet, translateName);
            translatedDatasets.put(translateName, dst);
            translatedNames.put(mainDatasetName, translateName);
*/
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
                //datasetMetadataMap.put(datasetNames[0], createDatasetMetadata(mainDataSet, mainDatasetName));
                datasetMetadataMap.put(datasetNames[0], createDatasetMetadata(mainDatasetName));

                //translatedMetadata.put(translateName, createDatasetMetadata(translateName));

            } else {
                datasetNames = new String[nSubdatasets + 1];
                for (int i = 0; i < nSubdatasets; i++) {
                    final String subdatasetName = (subdatasets.get(i * 2)).toString();
                    final int nameStartAt = subdatasetName.lastIndexOf("_NAME=") + 6;
                    datasetNames[i] = subdatasetName.substring(nameStartAt);
                }
                datasetNames[nSubdatasets] = mainDatasetName;
                datasetMetadataMap.put(datasetNames[nSubdatasets], createDatasetMetadata(mainDataSet, mainDatasetName));

                //translatedMetadata.put(translateName, createDatasetMetadata(dst, translateName));
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

    public void dispose() {
        super.dispose();

    }

    /**
     * Reset main values
     */
    public void reset() {
        this.input = null;
        dispose();
        super.reset();
        nSubdatasets = -1;
    }

/*
    protected Dataset translate(Dataset dataset, String name) {
        List<String> options = new ArrayList<>();
        options.add("-of");
        options.add("GTiff");
        options.add("-a_nodata");
        options.add("0.0");
        //options.add("-ot");
        //options.add("Byte");
        TranslateOptions to = new TranslateOptions(new Vector(options));

        return gdal.Translate(name, dataset, to);
    }
*/
    /**
     * Copied from imageio-ext GDALImageReader
     *
     * Read data from the required region of the raster.
     *
     * @param itemMetadata
     *                a <code>GDALCommonIIOImageMetadata</code> related to the
     *                actual dataset
     * @param srcRegion
     *                the source Region to be read
     * @param dstRegion
     *                the destination Region of the image read
     * @param selectedBands
     *                an array specifying the requested bands
     * @return the read <code>Raster</code>
     */
    public Raster readDatasetRaster(
            GDALCommonIIOImageMetadata itemMetadata,
            Rectangle srcRegion,
            Rectangle dstRegion,
            int[] selectedBands,
            SampleModel destSampleModel) throws IOException {

        SampleModel destSm = destSampleModel != null ? destSampleModel : itemMetadata.getSampleModel();

        Dataset dataset = datasetsMap.get(itemMetadata.getDatasetName());
        //Dataset dataset = translatedDatasets.get(itemMetadata.getDatasetName());
        if (dataset == null)
            throw new IOException("Error while acquiring the input dataset " + itemMetadata.getDatasetName());

        SampleModel sampleModel = null;
        DataBuffer imgBuffer = null;
        Band pBand = null;
        try {
            int dstWidth = dstRegion.width;
            int dstHeight = dstRegion.height;
            int srcRegionXOffset = srcRegion.x;
            int srcRegionYOffset = srcRegion.y;
            int srcRegionWidth = srcRegion.width;
            int srcRegionHeight = srcRegion.height;

            if (LOGGER.isLoggable(Level.FINE))
                LOGGER.fine("SourceRegion = " + srcRegion.toString());

            // Getting number of bands
            final int nBands = selectedBands != null ? selectedBands.length
                    : destSm.getNumBands();

            int[] banks = new int[nBands];
            int[] offsets = new int[nBands];

            // setting the number of pixels to read
            final int pixels = dstWidth * dstHeight;
            int bufferType = 0, bufferSize = 0;
            int typeSizeInBytes = 0;

            // ////////////////////////////////////////////////////////////////////
            //
            // -------------------------------------------------------------------
            // Raster Creation >>> Step 2: Data Read
            // -------------------------------------------------------------------
            //
            // ////////////////////////////////////////////////////////////////////

            // NOTE: Bands are not 0-base indexed, so we must add 1
            pBand = dataset.GetRasterBand(1);
//pBand.SetNoDataValue(Double.NaN);
            // setting buffer properties
            bufferType = pBand.getDataType();
            typeSizeInBytes = gdal.GetDataTypeSize(bufferType) / 8;
            bufferSize = nBands * pixels * typeSizeInBytes;

            // splitBands = false -> I read n Bands at once.
            // splitBands = false -> I need to read 1 Band at a time.
            boolean splitBands = false;

            if (bufferSize < 0 || destSm instanceof BandedSampleModel) {
                // The number resulting from the product
                // "numBands*pixels*gdal.GetDataTypeSize(buf_type) / 8"
                // may be negative (A very high number which is not
                // "int representable")
                // In such a case, we will read 1 band at a time.
                bufferSize = pixels * typeSizeInBytes;
                splitBands = true;
            }

            int dataBufferType = -1;
            byte[][] byteBands = new byte[nBands][];
            for (int k = 0; k < nBands; k++) {

                // If I'm reading n Bands at once and I performed the first read,
                // I quit the loop
                if (k > 0 && !splitBands)
                    break;

                final byte[] dataBuffer = new byte[bufferSize];

                final int returnVal;
                if (!splitBands) {
                    // I can read nBands at once.
                    final int[] bandsMap = new int[nBands];
                    if (selectedBands != null) {
                        for (int i = 0; i < nBands; i++)
                            // TODO: original code calls selectedBands[i] +1, but that seems to break it.  works if you don't add 1
                            //bandsMap[i] = selectedBands[i] + 1;
                            bandsMap[i] = selectedBands[i];
                    } else {
                        for (int i = 0; i < nBands; i++)
                            bandsMap[i] = i + 1;
                    }

                    returnVal = dataset.ReadRaster(srcRegionXOffset,
                            srcRegionYOffset, srcRegionWidth, srcRegionHeight,
                            dstWidth, dstHeight, bufferType, dataBuffer, bandsMap,
                            nBands * typeSizeInBytes, dstWidth * nBands
                                    * typeSizeInBytes, typeSizeInBytes);

                    byteBands[k] = dataBuffer;
                } else {
                    // I need to read 1 band at a time.
                    Band rBand = null;
                    try{
                        rBand = dataset.GetRasterBand(k + 1);
                        returnVal = rBand.ReadRaster(
                                srcRegionXOffset, srcRegionYOffset, srcRegionWidth,
                                srcRegionHeight, dstWidth, dstHeight, bufferType,
                                dataBuffer);
                        byteBands[k] = dataBuffer;
                    } finally {
                        if (rBand != null){
                            try{
                                // Closing the band
                                rBand.delete();
                            }catch (Throwable e) {
                                if(LOGGER.isLoggable(Level.FINEST))
                                    LOGGER.log(Level.FINEST,e.getLocalizedMessage(),e);
                            }
                        }
                    }
                }
                if (returnVal == gdalconstConstants.CE_None) {
                    if (!splitBands)
                        for (int band = 0; band < nBands; band++) {
                            banks[band] = band;
                            offsets[band] = band;
                        }
                    else {
                        banks[k] = k;
                        offsets[k] = 0;
                    }
                } else {
                    // The read operation was not successfully computed.
                    // Showing error messages.
                    LOGGER.info(new StringBuilder("Last error: ").append(
                            gdal.GetLastErrorMsg()).toString());
                    LOGGER.info(new StringBuilder("Last error number: ").append(
                            gdal.GetLastErrorNo()).toString());
                    LOGGER.info(new StringBuilder("Last error type: ").append(
                            gdal.GetLastErrorType()).toString());
                    throw new RuntimeException(gdal.GetLastErrorMsg());
                }
            }

            // ////////////////////////////////////////////////////////////////////
            //
            // -------------------------------------------------------------------
            // Raster Creation >>> Step 3: Setting DataBuffer
            // -------------------------------------------------------------------
            //
            // //////       //////////////////////////////////////////////////////////////

            // /////////////////////////////////////////////////////////////////////
            //
            // TYPE BYTE
            //
            // /////////////////////////////////////////////////////////////////////
            if (bufferType == gdalconstConstants.GDT_Byte) {
                if (!splitBands) {
                    //                final byte[] bytes = new byte[nBands * pixels];
                    //                bands[0].get(bytes, 0, nBands * pixels);
                    imgBuffer = new DataBufferByte(byteBands[0], nBands * pixels);
                } else {
                    //                final byte[][] bytes = new byte[nBands][];
                    //                for (int i = 0; i < nBands; i++) {
                    ////                    bytes[i] = new byte[pixels];
                    //                    bands[i].get(bytes[i], 0, pixels);
                    //                }
                    imgBuffer = new DataBufferByte(byteBands, pixels);
                }
                dataBufferType = DataBuffer.TYPE_BYTE;
            }
            else {
                ByteBuffer[] bands = new ByteBuffer[nBands];
                for (int k = 0; (splitBands && k < nBands) || (k < 1 && !splitBands); k++) {
                    bands[k]=ByteBuffer.wrap(byteBands[k],0,byteBands[k].length);
                }

                if (bufferType == gdalconstConstants.GDT_Int16
                        || bufferType == gdalconstConstants.GDT_UInt16) {
                    // ////////////////////////////////////////////////////////////////
                    //
                    // TYPE SHORT
                    //
                    // ////////////////////////////////////////////////////////////////

                    if (!splitBands) {
                        // I get short values from the ByteBuffer using a view
                        // of the ByteBuffer as a ShortBuffer
                        // It is worth to create the view outside the loop.
                        short[] shorts = new short[nBands * pixels];
                        bands[0].order(ByteOrder.nativeOrder());
                        final ShortBuffer buff = bands[0].asShortBuffer();
                        buff.get(shorts, 0, nBands * pixels);
                        if (bufferType == gdalconstConstants.GDT_Int16)
                            imgBuffer = new DataBufferShort(shorts, nBands * pixels);
                        else
                            imgBuffer = new DataBufferUShort(shorts, nBands * pixels);
                    } else {
                        short[][] shorts = new short[nBands][];
                        for (int i = 0; i < nBands; i++) {
                            shorts[i] = new short[pixels];
                            bands[i].order(ByteOrder.nativeOrder());
                            bands[i].asShortBuffer().get(shorts[i], 0, pixels);
                        }
                        if (bufferType == gdalconstConstants.GDT_Int16)
                            imgBuffer = new DataBufferShort(shorts, pixels);
                        else
                            imgBuffer = new DataBufferUShort(shorts, pixels);
                    }
                    if (bufferType == gdalconstConstants.GDT_UInt16)
                        dataBufferType = DataBuffer.TYPE_USHORT;
                    else
                        dataBufferType = DataBuffer.TYPE_SHORT;
                } else if (bufferType == gdalconstConstants.GDT_Int32
                        || bufferType == gdalconstConstants.GDT_UInt32) {
                    // ////////////////////////////////////////////////////////////////
                    //
                    // TYPE INT
                    //
                    // ////////////////////////////////////////////////////////////////

                    if (!splitBands) {
                        // I get int values from the ByteBuffer using a view
                        // of the ByteBuffer as an IntBuffer
                        // It is worth to create the view outside the loop.
                        int[] ints = new int[nBands * pixels];
                        bands[0].order(ByteOrder.nativeOrder());
                        final IntBuffer buff = bands[0].asIntBuffer();
                        buff.get(ints, 0, nBands * pixels);
                        imgBuffer = new DataBufferInt(ints, nBands * pixels);
                    } else {
                        int[][] ints = new int[nBands][];
                        for (int i = 0; i < nBands; i++) {
                            ints[i] = new int[pixels];
                            bands[i].order(ByteOrder.nativeOrder());
                            bands[i].asIntBuffer().get(ints[i], 0, pixels);
                        }
                        imgBuffer = new DataBufferInt(ints, pixels);
                    }
                    dataBufferType = DataBuffer.TYPE_INT;

                } else if (bufferType == gdalconstConstants.GDT_Float32) {
                    // /////////////////////////////////////////////////////////////////////
                    //
                    // TYPE FLOAT
                    //
                    // /////////////////////////////////////////////////////////////////////

                    if (!splitBands) {
                        // I get float values from the ByteBuffer using a view
                        // of the ByteBuffer as a FloatBuffer
                        // It is worth to create the view outside the loop.
                        float[] floats = new float[nBands * pixels];
                        bands[0].order(ByteOrder.nativeOrder());
                        final FloatBuffer buff = bands[0].asFloatBuffer();
                        buff.get(floats, 0, nBands * pixels);
                        imgBuffer = new DataBufferFloat(floats, nBands * pixels);
                    } else {
                        float[][] floats = new float[nBands][];
                        for (int i = 0; i < nBands; i++) {
                            floats[i] = new float[pixels];
                            bands[i].order(ByteOrder.nativeOrder());
                            bands[i].asFloatBuffer().get(floats[i], 0, pixels);
                        }
                        imgBuffer = new DataBufferFloat(floats, pixels);
                    }


                    dataBufferType = DataBuffer.TYPE_FLOAT;
                } else if (bufferType == gdalconstConstants.GDT_Float64) {
                    // /////////////////////////////////////////////////////////////////////
                    //
                    // TYPE DOUBLE
                    //
                    // /////////////////////////////////////////////////////////////////////

                    if (!splitBands) {
                        // I get double values from the ByteBuffer using a view
                        // of the ByteBuffer as a DoubleBuffer
                        // It is worth to create the view outside the loop.
                        double[] doubles = new double[nBands * pixels];
                        bands[0].order(ByteOrder.nativeOrder());
                        final DoubleBuffer buff = bands[0].asDoubleBuffer();
                        buff.get(doubles, 0, nBands * pixels);
                        imgBuffer = new DataBufferDouble(doubles, nBands * pixels);
                    } else {
                        double[][] doubles = new double[nBands][];
                        for (int i = 0; i < nBands; i++) {
                            doubles[i] = new double[pixels];
                            bands[i].order(ByteOrder.nativeOrder());
                            bands[i].asDoubleBuffer().get(doubles[i], 0, pixels);
                        }
                        imgBuffer = new DataBufferDouble(doubles, pixels);
                    }
                    dataBufferType = DataBuffer.TYPE_DOUBLE;

                } else {
                    // TODO: Handle more cases if needed. Show the name of the type
                    // instead of the numeric value.
                    LOGGER.info("The specified data type is actually unsupported: "
                            + bufferType);
                }
            }

            // ////////////////////////////////////////////////////////////////////
            //
            // -------------------------------------------------------------------
            // Raster Creation >>> Step 4: Setting SampleModel
            // -------------------------------------------------------------------
            //
            // ////////////////////////////////////////////////////////////////////
            // TODO: Fix this in compliance with the specified destSampleModel
            if (splitBands)
                sampleModel = new BandedSampleModel(dataBufferType, dstWidth,
                        dstHeight, dstWidth, banks, offsets);
            else
                sampleModel = new PixelInterleavedSampleModel(dataBufferType,
                        dstWidth, dstHeight, nBands, dstWidth * nBands, offsets);
        } finally {
            if (pBand != null){
                try{
                    // Closing the band
                    pBand.delete();
                }catch (Throwable e) {
                    if(LOGGER.isLoggable(Level.FINE))
                        LOGGER.log(Level.FINE,e.getLocalizedMessage(),e);
                }
            }
        }

        // ////////////////////////////////////////////////////////////////////
        //
        // -------------------------------------------------------------------
        // Raster Creation >>> Final Step: Actual Raster Creation
        // -------------------------------------------------------------------
        //
        // ////////////////////////////////////////////////////////////////////

        // return Raster.createWritableRaster(sampleModel, imgBuffer, new Point(
        // dstRegion.x, dstRegion.y));
        return Raster.createWritableRaster(sampleModel, imgBuffer, null);
    }

    /**
     * Read the raster and returns a <code>BufferedImage</code>
     *
     * @param imageIndex
     *                the index of the image to be retrieved.
     * @param param
     *                an <code>ImageReadParam</code> used to control the
     *                reading process, or <code>null</code>. Actually,
     *                setting a destinationType allows to specify the number of
     *                bands in the destination image.
     *
     * @return the desired portion of the image as a <code>BufferedImage</code>
     * @throws IllegalArgumentException
     *                 if <code>param</code> contains an invalid specification
     *                 of a source and/or destination band subset or of a
     *                 destination image.
     * @throws IOException
     *                 if an error occurs when acquiring access to the
     *                 underlying datasource
     */
    public BufferedImage read(final int imageIndex, final ImageReadParam param) throws IOException {

        // //
        //
        // Retrieving the requested dataset
        //
        // //
        final GDALCommonIIOImageMetadata item = getDatasetMetadata(imageIndex);
        final int width = item.getWidth();
        final int height = item.getHeight();
        final SampleModel itemSampleModel = item.getSampleModel();
        int itemNBands = itemSampleModel.getNumBands();
        int nDestBands;

        BufferedImage bi = null;
        final ImageReadParam imageReadParam;
        if (param == null)
            imageReadParam = getDefaultReadParam();
        else
            imageReadParam = param;

        // //
        //
        // First, check for a specified ImageTypeSpecifier
        //
        // //
        ImageTypeSpecifier imageType = imageReadParam.getDestinationType();
        SampleModel destSampleModel = null;
        if (imageType != null) {
            destSampleModel = imageType.getSampleModel();
            nDestBands = destSampleModel.getNumBands();
        } else {
            bi = imageReadParam.getDestination();
            if (bi != null)
                nDestBands = bi.getSampleModel().getNumBands();
            else
                nDestBands = itemNBands;
        }

        // //
        //
        // Second, bands settings check
        //
        // //
        checkReadParamBandSettings(imageReadParam, itemNBands, nDestBands);
        int[] srcBands = imageReadParam.getSourceBands();
//        int[] destBands = imageReadParam.getDestinationBands();
//
        // //
        //
        // Third, destination image check
        //
        // //
        // if (bi != null && imageType == null) {
        // if ((srcBands == null) && (destBands == null)) {
        // SampleModel biSampleModel = bi.getSampleModel();
        // if (!bi.getColorModel().equals(item.getColorModel())
        // || biSampleModel.getDataType() != itemSampleModel
        // .getDataType())
        // throw new IllegalArgumentException(
        // "Provided destination image does not have a valid ColorModel or
        // SampleModel");
        // }
        // }

        // //
        //
        // Computing regions of interest
        //
        // //
        Rectangle srcRegion = new Rectangle(0, 0, 0, 0);
        Rectangle destRegion = new Rectangle(0, 0, 0, 0);
        computeRegions(imageReadParam, width, height, bi, srcRegion, destRegion);
        if (imageReadParam != null){
            if (imageReadParam instanceof EnhancedImageReadParam){
                final EnhancedImageReadParam eparam = (EnhancedImageReadParam) imageReadParam;
                final Rectangle dstRegion = eparam.getDestinationRegion();
                if (dstRegion != null){
                    destRegion.height = dstRegion.height;
                    destRegion.width = dstRegion.width;
                }
            }
        }
        if (LOGGER.isLoggable(Level.FINE)) {
            LOGGER.fine("Source Region = " + srcRegion.toString());
            LOGGER.fine("Destination Region = " + destRegion.toString());
        }

        //
        // Getting data
        //
        if (bi == null) {
            // //
            //
            // No destination image has been specified.
            // Creating a new BufferedImage
            //
            // //
            ColorModel cm;
            if (imageType == null) {
                cm = item.getColorModel();
                bi = new BufferedImage(cm, (WritableRaster) readDatasetRaster(
                        item, srcRegion, destRegion, srcBands,null), false, null);
            } else {
                cm = imageType.getColorModel();
                bi = new BufferedImage(cm,
                        (WritableRaster) readDatasetRaster(item, srcRegion,
                                destRegion, srcBands, destSampleModel), false,
                        null);
            }

        } else {
            // //
            //
            // the destination image has been specified.
            //
            // //
            // Rectangle destSize = (Rectangle) destRegion.clone();
            // destSize.setLocation(0, 0);

            Raster readRaster = readDatasetRaster(item, srcRegion, destRegion,
                    srcBands,null);
            WritableRaster raster = bi.getRaster().createWritableChild(0, 0,
                    bi.getWidth(), bi.getHeight(), 0, 0, null);
            // TODO: Work directly on a Databuffer avoiding setRect?
            raster.setRect(destRegion.x, destRegion.y, readRaster);

            // Raster readRaster = readDatasetRaster(item, srcRegion,
            // destRegion,
            // srcBands);
            // WritableRaster raster = bi.getRaster().createWritableChild(
            // destRegion.x, destRegion.y, destRegion.width,
            // destRegion.height, destRegion.x, destRegion.y, null);
            // //TODO: Work directly on a Databuffer avoiding setRect?
            // raster.setRect(readRaster);
        }
        return bi;
    }

    /**
     * Performs a full read operation.
     *
     * @param imageIndex
     *                the index of the image to be retrieved.
     */
    public BufferedImage read(int imageIndex) throws IOException {
        if (LOGGER.isLoggable(Level.FINE))
            LOGGER.fine("read(imageIndex)");
        return read(imageIndex, null);

    }

    /**
     * Checks if the specified ImageIndex is valid.
     *
     * @param imageIndex
     *                the specified imageIndex
     */
    protected void checkImageIndex(final int imageIndex) {

        // When is an imageIndex not valid? 1) When it is negative 2) When the
        // format does not support subdatasets and imageIndex is > 0 3) When the
        // format supports subdatasets, but there isn't any subdataset and
        // imageIndex is greater than zero. 4) When the format supports
        // subdatasets, there are N subdatasets but imageIndex exceeds the
        // subdatasets count.
        //
        // It is worthwhile to remark that in case of nSubdatasets > 0, the
        // mainDataset is stored in the last position of datasetNames array. In
        // such a case the max valid imageIndex is nSubdatasets.

        if (imageIndex < 0 || imageIndex > nSubdatasets) {
            // The specified imageIndex is not valid.
            // Retrieving the valid image index range.
            final int maxImageIndex = nSubdatasets;
            final StringBuilder sb = new StringBuilder("Illegal imageIndex specified = ")
                    .append(imageIndex)
                    .append(", while the valid imageIndex");
            if (maxImageIndex > 0)
                // There are N Subdatasets.
                sb.append(" range should be (0,").append(maxImageIndex).append( ")!!");
            else
                // Only the imageIndex 0 is valid.
                sb.append(" should be 0!");
            throw new IndexOutOfBoundsException(sb.toString());
        }
    }
}
