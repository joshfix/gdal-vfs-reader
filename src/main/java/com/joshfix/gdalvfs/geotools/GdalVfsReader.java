package com.joshfix.gdalvfs.geotools;

import com.joshfix.gdalvfs.geotools.path.VfsPath;
import com.joshfix.gdalvfs.imageio.GdalVfsImageReaderSpi;
import com.sun.media.imageioimpl.common.BogusColorSpace;
import it.geosolutions.imageio.core.CoreCommonImageMetadata;
import it.geosolutions.imageio.gdalframework.GDALCommonIIOImageMetadata;
import it.geosolutions.imageio.gdalframework.GDALUtilities;
import it.geosolutions.imageio.maskband.DatasetLayout;
import it.geosolutions.imageioimpl.plugins.tiff.TIFFStreamMetadata;
import it.geosolutions.imageioimpl.plugins.tiff.TiffDatasetLayoutImpl;
import it.geosolutions.jaiext.range.NoDataContainer;
import it.geosolutions.jaiext.utilities.ImageLayout2;
import org.gdal.gdal.Dataset;
import org.gdal.gdal.gdal;
import org.geotools.coverage.CoverageFactoryFinder;
import org.geotools.coverage.GridSampleDimension;
import org.geotools.coverage.TypeMap;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridCoverageFactory;
import org.geotools.coverage.grid.GridEnvelope2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.coverage.grid.io.AbstractGridCoverage2DReader;
import org.geotools.coverage.grid.io.AbstractGridFormat;
import org.geotools.coverage.grid.io.GridCoverage2DReader;
import org.geotools.coverage.grid.io.OverviewPolicy;
import org.geotools.coverage.grid.io.imageio.MaskOverviewProvider;
import org.geotools.coverage.grid.io.imageio.geotiff.GeoTiffIIOMetadataDecoder;
import org.geotools.coverage.grid.io.imageio.geotiff.GeoTiffMetadata2CRSAdapter;
import org.geotools.coverage.util.CoverageUtilities;
import org.geotools.data.*;
import org.geotools.geometry.GeneralEnvelope;
import org.geotools.geometry.PixelTranslation;
import org.geotools.image.ImageWorker;
import org.geotools.image.util.ImageUtilities;
import org.geotools.referencing.CRS;
import org.geotools.referencing.operation.builder.GridToEnvelopeMapper;
import org.geotools.referencing.operation.transform.IdentityTransform;
import org.geotools.referencing.operation.transform.ProjectiveTransform;
import org.geotools.util.URLs;
import org.geotools.util.Utilities;
import org.geotools.util.factory.GeoTools;
import org.geotools.util.factory.Hints;

import org.opengis.coverage.ColorInterpretation;
import org.opengis.coverage.grid.Format;
import org.opengis.coverage.grid.GridCoverage;
import org.opengis.coverage.grid.GridEnvelope;
import org.opengis.geometry.Envelope;
import org.opengis.parameter.GeneralParameterValue;
import org.opengis.parameter.ParameterDescriptor;
import org.opengis.parameter.ParameterValue;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.ReferenceIdentifier;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.datum.PixelInCell;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.imageio.ImageReader;
import javax.imageio.ImageTypeSpecifier;
import javax.imageio.metadata.IIOMetadata;
import javax.imageio.spi.ImageReaderSpi;
import javax.media.jai.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.image.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author joshfix
 * Created on 2019-07-29
 */
public class GdalVfsReader extends AbstractGridCoverage2DReader {//BaseGDALGridCoverage2DReader {//AbstractGridCoverage2DReader {//implements GridCoverage2DReader {

    private VfsPath vfsPath;

    /**The coverage factory producing a {@link GridCoverage} from an image */
    //private GridCoverageFactory coverageFactory;

    /** This contains the number of overviews.aaa */
    //protected int numOverviews = 0;

    /** 2DGridToWorld math transform. */
    //protected MathTransform raster2Model = null;

    /** Envelope read from file */
    //protected GeneralEnvelope originalEnvelope = null;

    /** Coverage name */
    //protected String coverageName = "geotools_coverage";

    /** Source to read from */
    //protected Object source = null;

    /** Hints used by the {@link AbstractGridCoverage2DReader} subclasses. */
    //protected Hints hints = GeoTools.getDefaultHints();

    /** Highest resolution available for this reader. */
    //protected double[] highestRes = null;

    /** Temp variable used in many readers. */
    //protected boolean closeMe;

    /** In case we are trying to read from a GZipped file this will be set to true. */
    //protected boolean gzipped;

    /** The original {@link GridEnvelope} for the {@link GridCoverage2D} of this reader. */
    //protected GridEnvelope originalGridRange = null;

    /** crs for this coverage */
    //protected CoordinateReferenceSystem crs = null;

    /** scales and offsets for rescaling */
    protected Double[] scales;

    protected Double[] offsets;

    protected ImageLayout2 layout;

    protected Map<String, ArrayList<Resolution>> resolutionsLevelsMap = new HashMap<>();

    ///** Resolutions avialaible through an overviews based mechanism. */
    //protected double[][] overViewResolutions = null;

    ///** Coverage {@link DatasetLayout} containing information about Overviews and Mask management */
    //protected DatasetLayout dtLayout;

    protected double noData = Double.NaN;

    /** Boolean indicating if {@link MaskOverviewProvider} is present */
    private boolean hasMaskOvrProvider = false;

    /** {@link MaskOverviewProvider} instance used for handling internal/external Overviews */
    private MaskOverviewProvider maskOvrProvider;

    /** Adapter for the GeoTiff crs. */
    private GeoTiffMetadata2CRSAdapter gtcs;

    /** Caches an {@code ImageReaderSpi}. */
    protected ImageReaderSpi imageReaderSpi;

    protected ImageReader imageReader;
    private static final Logger LOGGER =
            org.geotools.util.logging.Logging.getLogger(GdalVfsReader.class);

    static {
        //gdal.AllRegister();
        GDALUtilities.loadGDAL();
        gdal.SetConfigOption("CHECK_DISK_FREE_SPACE", "FALSE");
        // config options per https://trac.osgeo.org/gdal/wiki/CloudOptimizedGeoTIFF#HowtoreaditwithGDAL
        gdal.SetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "YES");
        gdal.SetConfigOption("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif");
        gdal.SetConfigOption("CPL_CURL_VERBOSE", "YES");
        gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO");
    }

    public GdalVfsReader(VfsPath vfsPath) throws DataSourceException {
        this(vfsPath, null);
    }

    public GdalVfsReader(VfsPath vfsPath, Hints hints) throws DataSourceException {
        super(vfsPath, null);

        imageReaderSpi = new GdalVfsImageReaderSpi();
        try {
            imageReader = imageReaderSpi.createReaderInstance();
        } catch (IOException e) {
            e.printStackTrace();
        }
        imageReader.setInput(new GdalVfsImageInputStream(vfsPath));

        // GridCoverageFactory initialization
        if (this.hints.containsKey(Hints.GRID_COVERAGE_FACTORY)) {
            final Object factory = this.hints.get(Hints.GRID_COVERAGE_FACTORY);
            if (factory != null && factory instanceof GridCoverageFactory) {
                this.coverageFactory = (GridCoverageFactory) factory;
            }
        }
        if (this.coverageFactory == null) {
            this.coverageFactory = CoverageFactoryFinder.getGridCoverageFactory(this.hints);
        }

        this.vfsPath = vfsPath;

        try {
            setCoverageProperties(imageReader);
            getResolutionInfo();
            setLayout(imageReader);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (TransformException e) {
            e.printStackTrace();
        }
    }

    /**
     * Extract the ImageLayout from the provided reader for the first available image.
     *
     * @param reader an istance of {@link ImageReader}
     * @throws IOException in case an error occurs
     */
    @Override
    protected void setLayout(ImageReader reader) throws IOException {
        Utilities.ensureNonNull("reader", reader);
        // save ImageLayout
        layout = new ImageLayout2();
        ImageTypeSpecifier its = reader.getImageTypes(0).next();
        layout.setColorModel(its.getColorModel()).setSampleModel(its.getSampleModel());
        layout.setMinX(0).setMinY(0).setWidth(reader.getWidth(0)).setHeight(reader.getHeight(0));
        layout.setTileGridXOffset(0)
                .setTileGridYOffset(0)
                .setTileWidth(reader.getTileWidth(0))
                .setTileHeight(reader.getTileHeight(0));
        setlayout(layout);
    }

    /**
     * Gets resolution information about the coverage itself.
     *
     * @throws IOException
     * @throws TransformException
     * @throws DataSourceException
     */
    private void getResolutionInfo() throws IOException, TransformException {
        // ///
        //
        // setting the higher resolution available for this coverage
        //
        // ///
        highestRes = CoverageUtilities.getResolution((AffineTransform) raster2Model);

        if (LOGGER.isLoggable(Level.FINE))
            LOGGER.fine(
                    new StringBuffer("Highest Resolution = [")
                            .append(highestRes[0])
                            .append(",")
                            .append(highestRes[1])
                            .toString());
    }

    /**
     * Taken from BaseGDALGridCoverage2DReader
     *
     * Given a {@link GDALCommonIIOImageMetadata} metadata object, retrieves several properties to
     * properly set envelope, gridrange and crs.
     *
     * @param metadata a {@link GDALCommonIIOImageMetadata} metadata instance from where to search
     *     needed properties.
     */
    private void parseCommonMetadata(final GDALCommonIIOImageMetadata metadata) {

        final GeoTiffIIOMetadataDecoder geoTiffMetadata = new GeoTiffIIOMetadataDecoder(metadata);
        gtcs = new GeoTiffMetadata2CRSAdapter(hints);
        // ////////////////////////////////////////////////////////////////////
        //
        // setting CRS and Envelope directly from GDAL, if available
        //
        // ////////////////////////////////////////////////////////////////////
        // //
        //
        // 1) CRS
        //
        // //
        final Object tempCRS = this.hints.get(Hints.DEFAULT_COORDINATE_REFERENCE_SYSTEM);
        if (tempCRS != null) {
            this.crs = (CoordinateReferenceSystem) tempCRS;
            LOGGER.log(Level.FINE, "Using default coordinate reference system ");
        } else {

            // //
            //
            // Default CRS override it all
            //
            // //
            // //
            //
            // Let the prj file override the internal representation for the undelrying source of
            // information.
            //
            // //
            //parsePRJFile();
            // if there was not prj or the envelope could not be created easily, let's go with the
            // standard metadata.
            if (this.crs == null) {
                final String wkt = metadata.getProjection();

                if ((wkt != null) && !(wkt.equalsIgnoreCase(""))) {
                    try {
                        this.crs = CRS.parseWKT(wkt);
                        final Integer epsgCode = CRS.lookupEpsgCode(this.crs, true);
                        // Force the creation of the CRS directly from the
                        // retrieved EPSG code in order to prevent weird transformation
                        // between "same" CRSs having slight differences.
                        // TODO: cache epsgCode-CRSs
                        if (epsgCode != null) {
                            this.crs = CRS.decode("EPSG:" + epsgCode);
                        }
                    } catch (FactoryException fe) {
                        // unable to get CRS from WKT
                        if (LOGGER.isLoggable(Level.FINE)) {
                            LOGGER.log(
                                    Level.FINE,
                                    "Unable to get CRS from WKT contained in metadata. Looking for a PRJ.");
                        }
                        // reset crs
                        this.crs = null;
                    }
                }
            }
        }

        if (geoTiffMetadata.hasNoData()) {
            noData = geoTiffMetadata.getNoData();
        }

        // collect scales and offsets is present
        collectScaleOffset(metadata);
        // //
        //
        // 2) Grid
        //
        // //
        if (this.originalGridRange == null)
            this.originalGridRange =
                    new GridEnvelope2D(
                            new Rectangle(0, 0, metadata.getWidth(), metadata.getHeight()));

        // //
        //
        // 3) Envelope
        //
        // //
        //
        // Let's look for a world file first.
        //
        //parseWorldFile();
        if (this.originalEnvelope == null) {
            final double[] geoTransform = metadata.getGeoTransformation();
            if ((geoTransform != null) && (geoTransform.length == 6)) {
                final AffineTransform tempTransform =
                        new AffineTransform(
                                geoTransform[1],
                                geoTransform[4],
                                geoTransform[2],
                                geoTransform[5],
                                geoTransform[0],
                                geoTransform[3]);
                // ATTENTION: Gdal geotransform does not use the pixel is
                // centre convention like world files.
                if (this.originalEnvelope == null) {
                    try {
                        // Envelope setting
                        this.originalEnvelope =
                                CRS.transform(
                                        ProjectiveTransform.create(tempTransform),
                                        new GeneralEnvelope(
                                                ((GridEnvelope2D) this.originalGridRange)));
                    } catch (IllegalStateException e) {
                        if (LOGGER.isLoggable(Level.WARNING)) {
                            LOGGER.log(Level.WARNING, e.getLocalizedMessage(), e);
                        }
                    } catch (TransformException e) {
                        if (LOGGER.isLoggable(Level.WARNING)) {
                            LOGGER.log(Level.WARNING, e.getLocalizedMessage(), e);
                        }
                    }
                }
                // Grid2World Transformation
                final double tr = -PixelTranslation.getPixelTranslation(PixelInCell.CELL_CORNER);
                tempTransform.translate(tr, tr);
                this.raster2Model = ProjectiveTransform.create(tempTransform);
            }
        }
    }

    protected void collectScaleOffset(IIOMetadata iioMetadata) {
        if (iioMetadata instanceof CoreCommonImageMetadata) {
            CoreCommonImageMetadata ccm = (CoreCommonImageMetadata) iioMetadata;
            this.scales = ccm.getScales();
            this.offsets = ccm.getOffsets();
        }
    }

    /**
     * This method is responsible for checking the provided coverage name against the coverage name
     * for this {@link GridCoverage2DReader}.
     *
     * @param coverageName the coverage name to check.
     * @return <code>true</code> if this {@link GridCoverage2DReader} contains the provided coverage
     *     name, <code>false</code> otherwise.
     */
    protected boolean checkName(String coverageName) {
        Utilities.ensureNonNull("coverageName", coverageName);
        return coverageName.equalsIgnoreCase(this.coverageName);
    }

    @Override
    public GeneralEnvelope getOriginalEnvelope() {
        return getOriginalEnvelope(coverageName);
    }

    @Override
    public GeneralEnvelope getOriginalEnvelope(String coverageName) {
        if (!checkName(coverageName)) {
            throw new IllegalArgumentException(
                    "The specified coverageName " + coverageName + "is not supported");
        }
        return new GeneralEnvelope(originalEnvelope);
    }

    @Override
    public CoordinateReferenceSystem getCoordinateReferenceSystem() {
        return getCoordinateReferenceSystem(coverageName);
    }

    /**
     * Retrieves the {@link GeneralEnvelope} for this {@link AbstractGridCoverage2DReader}.
     *
     * @return the {@link GeneralEnvelope} for this {@link AbstractGridCoverage2DReader}.
     */
    @Override
    public CoordinateReferenceSystem getCoordinateReferenceSystem(String coverageName) {
        if (!checkName(coverageName)) {
            throw new IllegalArgumentException(
                    "The specified coverageName " + coverageName + "is not supported");
        }

        return crs;
    }

    @Override
    public GridEnvelope getOriginalGridRange() {
        return getOriginalGridRange(coverageName);
    }

    @Override
    public GridEnvelope getOriginalGridRange(String coverageName) {
        if (!checkName(coverageName)) {
            throw new IllegalArgumentException(
                    "The specified coverageName " + coverageName + "is not supported");
        }
        assert originalGridRange.getDimension() == 2;
        return new GridEnvelope2D(
                originalGridRange.getLow(0),
                originalGridRange.getLow(1),
                originalGridRange.getSpan(0),
                originalGridRange.getSpan(1));
    }

    /**
     * Retrieves the original grid to world transformation for this {@link
     * AbstractGridCoverage2DReader}.
     *
     * @param pixInCell specifies the datum of the transformation we want.
     * @return the original grid to world transformation for this {@link
     *     AbstractGridCoverage2DReader}.
     */
    @Override
    public MathTransform getOriginalGridToWorld(PixelInCell pixInCell) {
        // Default implementation for backwards compatibility
        return getOriginalGridToWorld(coverageName, pixInCell);
    }

    @Override
    public MathTransform getOriginalGridToWorld(String coverageName, PixelInCell pixInCell) {
        if (!checkName(coverageName)) {
            throw new IllegalArgumentException(
                    "The specified coverageName " + coverageName + "is not supported");
        }
        synchronized (this) {
            if (raster2Model == null) {
                final GridToEnvelopeMapper geMapper =
                        new GridToEnvelopeMapper(
                                getOriginalGridRange(coverageName),
                                getOriginalEnvelope(coverageName));
                geMapper.setPixelAnchor(PixelInCell.CELL_CENTER);
                raster2Model = geMapper.createTransform();
            }
        }

        // we do not have to change the pixel datum
        if (pixInCell == PixelInCell.CELL_CENTER) return raster2Model;

        // we do have to change the pixel datum
        if (raster2Model instanceof AffineTransform) {
            final AffineTransform tr = new AffineTransform((AffineTransform) raster2Model);
            tr.concatenate(AffineTransform.getTranslateInstance(-0.5, -0.5));
            return ProjectiveTransform.create(tr);
        }
        if (raster2Model instanceof IdentityTransform) {
            final AffineTransform tr = new AffineTransform(1, 0, 0, 1, 0, 0);
            tr.concatenate(AffineTransform.getTranslateInstance(-0.5, -0.5));
            return ProjectiveTransform.create(tr);
        }
        throw new IllegalStateException("This reader's grid to world transform is invalud!");
    }

    @Override
    public Format getFormat() {
        return new GdalVfsFormat();
    }

    @Override
    public String[] getMetadataNames() {
        if (!checkName(coverageName)) {
            throw new IllegalArgumentException(
                    "The specified coverageName " + coverageName + "is not supported");
        }
        return getMetadataNames();
    }

    @Override
    public String[] getMetadataNames(String coverageName) {
        return new String[0];
    }

    @Override
    public String getMetadataValue(String name) {
        return null;
    }

    @Override
    public String getMetadataValue(String coverageName, String name) {
        return null;
    }

    @Override
    public String[] listSubNames() {
        return new String[0];
    }

    @Override
    public String[] getGridCoverageNames() {
        return new String[]{coverageName};
    }

    @Override
    public int getGridCoverageCount() {
        return 1;
    }

    @Override
    public GridCoverage2D read(GeneralParameterValue[] params) throws IOException {
        GeneralEnvelope requestedEnvelope = null;
        Rectangle dim = null;
        Color inputTransparentColor = null;
        OverviewPolicy overviewPolicy = null;
        int[] suggestedTileSize = null;
        boolean rescalePixels = GdalVfsFormat.RESCALE_PIXELS.getDefaultValue();


        //
        // Checking params
        //
        if (params != null) {
            for (int i = 0; i < params.length; i++) {
                final ParameterValue param = (ParameterValue) params[i];
                final ReferenceIdentifier name = param.getDescriptor().getName();
                if (name.equals(AbstractGridFormat.READ_GRIDGEOMETRY2D.getName())) {
                    final GridGeometry2D gg = (GridGeometry2D) param.getValue();
                    requestedEnvelope = new GeneralEnvelope((Envelope) gg.getEnvelope2D());
                    dim = gg.getGridRange2D().getBounds();
                    continue;
                }
                if (name.equals(AbstractGridFormat.OVERVIEW_POLICY.getName())) {
                    overviewPolicy = (OverviewPolicy) param.getValue();
                    continue;
                }
                if (name.equals(AbstractGridFormat.INPUT_TRANSPARENT_COLOR.getName())) {
                    inputTransparentColor = (Color) param.getValue();
                    continue;
                }
                if (name.equals(AbstractGridFormat.SUGGESTED_TILE_SIZE.getName())) {
                    String suggestedTileSize_ = (String) param.getValue();
                    if (suggestedTileSize_ != null && suggestedTileSize_.length() > 0) {
                        suggestedTileSize_ = suggestedTileSize_.trim();
                        int commaPosition = suggestedTileSize_.indexOf(",");
                        if (commaPosition < 0) {
                            int tileDim = Integer.parseInt(suggestedTileSize_);
                            suggestedTileSize = new int[] {tileDim, tileDim};
                        } else {
                            int tileW =
                                    Integer.parseInt(
                                            suggestedTileSize_.substring(0, commaPosition));
                            int tileH =
                                    Integer.parseInt(
                                            suggestedTileSize_.substring(commaPosition + 1));
                            suggestedTileSize = new int[] {tileW, tileH};
                        }
                    }
                    continue;
                }
                if (name.equals(GdalVfsFormat.RESCALE_PIXELS.getName())) {
                    rescalePixels = Boolean.TRUE.equals(param.getValue());
                }
            }
        }

        Hints newHints = null;
        if (suggestedTileSize != null) {
            newHints = hints.clone();
            final ImageLayout layout = new ImageLayout();
            layout.setTileGridXOffset(0);
            layout.setTileGridYOffset(0);
            layout.setTileHeight(suggestedTileSize[1]);
            layout.setTileWidth(suggestedTileSize[0]);
            newHints.add(new RenderingHints(JAI.KEY_IMAGE_LAYOUT, layout));
        }

        BufferedImage bi  = imageReader.read(0);
        PlanarImage coverageRaster =  PlanarImage.wrapRenderedImage(bi);
// applying rescale if needed
        if (rescalePixels) {
            if (!Double.isNaN(noData)) {
                // Force nodata settings since JAI ImageRead may lost that
                // We have to make sure that noData pixels won't be rescaled
                PlanarImage t = PlanarImage.wrapRenderedImage(coverageRaster);
                t.setProperty(NoDataContainer.GC_NODATA, new NoDataContainer(noData));
                coverageRaster = t;
            }
            coverageRaster =
                    PlanarImage.wrapRenderedImage(
                            applyRescaling(
                                    scales, offsets, coverageRaster, newHints));
        }

        //
        // MASKING INPUT COLOR as indicated
        //
        if (inputTransparentColor != null) {
            coverageRaster =
                    new ImageWorker(coverageRaster)
                            .setRenderingHints(newHints)
                            .makeColorTransparent(inputTransparentColor)
                            .getRenderedOperation();
        }

        //
        // External/Internal Masking
        //
        // ROI definition
        ROI roi = null;
        // Using MaskOvrProvider
        /*
        if (hasMaskOvrProvider) {
            // Parameter definition
            GridEnvelope ogr = getOriginalGridRange();
            Rectangle sourceRegion;
            if (readP.getSourceRegion() != null) {
                sourceRegion = readP.getSourceRegion();
            } else {
                sourceRegion = new Rectangle(ogr.getSpan(0), ogr.getSpan(1));
            }

            MaskOverviewProvider.MaskInfo info = maskOvrProvider.getMaskInfo(imageChoice, sourceRegion, readP);
            if (info != null) {
                // Reading Mask
                RenderedOp roiRaster =
                        readROIRaster(
                                info.streamSpi,
                                URLs.fileToUrl(info.file),
                                info.index,
                                newHints,
                                info.readParameters);
                roi = MaskOverviewProvider.scaleROI(roiRaster, coverageRaster.getBounds());
            }
        }
*/

        return createImageCoverage(coverageName, coverageRaster, null);
    }

    @Override
    public GridCoverage2D read(String coverageName, GeneralParameterValue[] parameters) throws IOException {
        // Default implementation for backwards compatibility
        if (coverageName.equalsIgnoreCase(this.coverageName)) {
            return read(parameters);
        }
        // Subclasses should do more checks on coverageName
        throw new IllegalArgumentException(
                "The specified coverageName " + coverageName + "is not supported");
    }

    /**
     * Applies values rescaling if either the scales or the offsets array is non null, or has any
     * value that is not a default (1 for scales, 0 for offsets)
     *
     * @param scales The scales array
     * @param offsets The offsets array
     * @param image The image to be rescaled
     * @param hints The image processing hints, if any (can be null)
     * @return The original image, or a rescaled image
     */
    public static RenderedImage applyRescaling(
            Double[] scales, Double[] offsets, RenderedImage image, Hints hints) {
        // if there is nothing to do, return immediately
        if (scales == null && offsets == null) {
            return image;
        }

        // convert to primitives, apply defaults to nullable elements
        int numBands = image.getSampleModel().getNumBands();
        double[] pscales = toPrimitiveArray(scales, numBands, 1);
        double[] poffsets = toPrimitiveArray(offsets, numBands, 0);
        boolean hasRescaling = false;
        for (int i = 0; i < numBands && !hasRescaling; i++) {
            hasRescaling = poffsets[i] != 0 || pscales[i] != 1;
        }
        if (!hasRescaling) {
            return image;
        }

        // use the input hints if possible, but create a proper layout to impose the target data
        // type
        RenderingHints localHints =
                hints != null
                        ? hints.clone()
                        : (RenderingHints) JAI.getDefaultInstance().getRenderingHints().clone();
        final ImageLayout layout =
                Optional.ofNullable((ImageLayout) localHints.get(JAI.KEY_IMAGE_LAYOUT))
                        .map(il -> (ImageLayout) il.clone())
                        .orElse(new ImageLayout2(image));
        SampleModel sm =
                RasterFactory.createBandedSampleModel(
                        DataBuffer.TYPE_DOUBLE,
                        image.getTileWidth(),
                        image.getTileHeight(),
                        image.getSampleModel().getNumBands());
        layout.setSampleModel(sm);
        layout.setColorModel(
                new ComponentColorModel(
                        new BogusColorSpace(numBands),
                        false,
                        false,
                        Transparency.OPAQUE,
                        DataBuffer.TYPE_DOUBLE));
        localHints.put(JAI.KEY_IMAGE_LAYOUT, layout);

        // at least one band is getting rescaled, apply the operation
        ImageWorker iw = new ImageWorker(image);
        iw.setRenderingHints(localHints);
        iw.rescale(pscales, poffsets);
        return iw.getRenderedImage();
    }

    private static double[] toPrimitiveArray(Double[] input, int numBands, double defaultValue) {
        double[] result = new double[numBands];
        Arrays.fill(result, defaultValue);

        if (input != null) {
            int loopMax = Math.min(input.length, numBands);
            for (int i = 0; i < loopMax; i++) {
                Double v = input[i];
                if (v != null) {
                    result[i] = v;
                }
            }
        }

        return result;
    }

    @Override
    public void dispose() {
        imageReader.dispose();
    }

    @Override
    public Set<ParameterDescriptor<List>> getDynamicParameters() throws IOException {
        return null;
    }

    @Override
    public Set<ParameterDescriptor<List>> getDynamicParameters(String coverageName) throws IOException {
        return null;
    }

    @Override
    public double[] getReadingResolutions(OverviewPolicy policy, double[] requestedResolution) throws IOException {
        // Default implementation for backwards compatibility
        return getReadingResolutions(coverageName, policy, requestedResolution);
    }

    @Override
    public double[] getReadingResolutions(String coverageName, OverviewPolicy policy, double[] requestedResolution) throws IOException {
        if (!checkName(coverageName)) {
            throw new IllegalArgumentException(
                    "The specified coverageName " + coverageName + "is not supported");
        }
        // find the target resolution level
        double[] result;
        if (numOverviews > 0) {
            int imageIdx = pickOverviewLevel(coverageName, policy, requestedResolution);
            result = imageIdx > 0 ? overViewResolutions[imageIdx - 1] : highestRes;
        } else {
            result = getHighestRes();
        }

        // return via cloning to protect internal state
        double[] clone = new double[result.length];
        System.arraycopy(result, 0, clone, 0, result.length);
        return clone;
    }

    @Override
    public int getNumOverviews(String coverageName) {
        if (!checkName(coverageName)) {
            throw new IllegalArgumentException(
                    "The specified coverageName " + coverageName + "is not supported");
        }
        if (dtLayout == null) {
            // Back to the default
            return numOverviews;
        }
        return dtLayout.getNumInternalOverviews()
                + (dtLayout.getNumExternalOverviews() > 0 ? dtLayout.getNumExternalOverviews() : 0);
    }

    @Override
    public int getNumOverviews() {
        // Default implementation for backwards compatibility
        return getNumOverviews(coverageName);
    }


    @Override
    public DatasetLayout getDatasetLayout(String coverageName) {
        return null;
    }

    @Override
    public ImageLayout getImageLayout() throws IOException {
        return layout;
    }

    public String getCoverageName() {
        return coverageName;
    }

    @Override
    public ImageLayout getImageLayout(String coverageName) throws IOException {
        return layout;
    }

    @Override
    public double[][] getResolutionLevels() throws IOException {
        // Default implementation for backwards compatibility
        return getResolutionLevels(coverageName);
    }

    @Override
    public double[][] getResolutionLevels(String coverageName) throws IOException {
        if (!checkName(coverageName)) {
            throw new IllegalArgumentException(
                    "The specified coverageName " + coverageName + "is not supported");
        }

        final double[][] returnValue = new double[numOverviews + 1][2];
        double[] hres = getHighestRes();
        if (hres == null) {
            return null;
        } else {
            System.arraycopy(hres, 0, returnValue[0], 0, 2);
            for (int i = 1; i < returnValue.length; i++) {
                System.arraycopy(overViewResolutions[i - 1], 0, returnValue[i], 0, 2);
            }
            return returnValue;
        }
    }

    @Override
    public ServiceInfo getInfo() {
        DefaultServiceInfo info = new DefaultServiceInfo();
        info.setDescription(source == null ? null : String.valueOf(source));
        info.setTitle(vfsPath.getPath());
        return info;
    }

    @Override
    public ResourceInfo getInfo(String coverageName) {
        return null;
    }

    double[] getHighestRes() {
        return getHighestRes(coverageName);
    }

    protected double[] getHighestRes(String coverageName) {
        if (!checkName(coverageName)) {
            throw new IllegalArgumentException(
                    "The specified coverageName " + coverageName + "is not supported");
        }

        return highestRes;
    }

    private Integer pickOverviewLevel(
            String coverageName, OverviewPolicy policy, double[] requestedRes) {
        // setup policy
        if (policy == null) policy = extractOverviewPolicy();

        ArrayList<Resolution> resolutionsLevels;

        // sort resolutions from smallest pixels (higher res) to biggest pixels (higher res)
        // keeping a reference to the original image choice
        synchronized (this) {
            resolutionsLevels = resolutionsLevelsMap.get(coverageName);
            if (resolutionsLevels == null) {
                resolutionsLevels = new ArrayList<Resolution>();
                resolutionsLevelsMap.put(coverageName, resolutionsLevels);
                // note that we assume what follows:
                // -highest resolution image is at level 0.
                // -all the overviews share the same envelope
                // -the aspect ratio for the overviews is constant
                // -the provided resolutions are taken directly from the grid
                resolutionsLevels.add(new Resolution(1, getHighestRes()[0], getHighestRes()[1], 0));
                if (numOverviews > 0) {
                    for (int i = 0; i < overViewResolutions.length; i++)
                        resolutionsLevels.add(
                                new Resolution(
                                        overViewResolutions[i][0] / getHighestRes()[0],
                                        overViewResolutions[i][0],
                                        overViewResolutions[i][1],
                                        i + 1));
                    Collections.sort(resolutionsLevels);
                }
            }
        }

        // Now search for the best matching resolution.
        // Check also for the "perfect match"... unlikely in practice unless someone
        // tunes the clients to request exactly the resolution embedded in
        // the overviews, something a perf sensitive person might do in fact

        // the requested resolutions
        final double reqx = requestedRes[0];
        final double reqy = requestedRes[1];

        // requested scale factor for least reduced axis
        final Resolution max = resolutionsLevels.get(0);
        final double requestedScaleFactorX = reqx / max.resolutionX;
        final double requestedScaleFactorY = reqy / max.resolutionY;
        final int leastReduceAxis = requestedScaleFactorX <= requestedScaleFactorY ? 0 : 1;
        final double requestedScaleFactor =
                leastReduceAxis == 0 ? requestedScaleFactorX : requestedScaleFactorY;

        // are we looking for a resolution even higher than the native one?
        if (requestedScaleFactor <= 1) return max.imageChoice;
        // are we looking for a resolution even lower than the smallest overview?
        final Resolution min = resolutionsLevels.get(resolutionsLevels.size() - 1);
        if (requestedScaleFactor >= min.scaleFactor) return min.imageChoice;
        // Ok, so we know the overview is between min and max, skip the first
        // and search for an overview with a resolution lower than the one requested,
        // that one and the one from the previous step will bound the searched resolution
        Resolution prev = max;
        final int size = resolutionsLevels.size();
        for (int i = 1; i < size; i++) {
            final Resolution curr = resolutionsLevels.get(i);
            // perfect match check
            if (curr.scaleFactor == requestedScaleFactor) {
                return curr.imageChoice;
            }

            // middle check. The first part of the condition should be sufficient, but
            // there are cases where the x resolution is satisfied by the lowest resolution,
            // the y by the one before the lowest (so the aspect ratio of the request is
            // different than the one of the overviews), and we would end up going out of the loop
            // since not even the lowest can "top" the request for one axis
            if (curr.scaleFactor > requestedScaleFactor || i == size - 1) {
                if (policy == OverviewPolicy.QUALITY) return prev.imageChoice;
                else if (policy == OverviewPolicy.SPEED) return curr.imageChoice;
                else if (requestedScaleFactor - prev.scaleFactor
                        < curr.scaleFactor - requestedScaleFactor) return prev.imageChoice;
                else return curr.imageChoice;
            }
            prev = curr;
        }
        // fallback
        return max.imageChoice;
    }

    /**
     * This method is responsible for checking the overview policy as defined by the provided {@link
     * Hints}.
     */
    private OverviewPolicy extractOverviewPolicy() {
        OverviewPolicy overviewPolicy = null;
        // check if a policy was provided using hints (check even the
        // deprecated one)
        if (this.hints != null)
            if (this.hints.containsKey(Hints.OVERVIEW_POLICY))
                overviewPolicy = (OverviewPolicy) this.hints.get(Hints.OVERVIEW_POLICY);

        // use default if not provided. Default is nearest
        if (overviewPolicy == null) overviewPolicy = OverviewPolicy.getDefaultPolicy();
        assert overviewPolicy != null;
        return overviewPolicy;
    }

    /**
     * Simple support class for sorting overview resolutions
     *
     * @author Andrea Aime
     * @author Simone Giannecchini, GeoSolutions.
     * @since 2.5
     */
    private static class Resolution implements Comparable<Resolution> {
        double scaleFactor;

        double resolutionX;

        double resolutionY;

        int imageChoice;

        public Resolution(
                final double scaleFactor,
                final double resolutionX,
                final double resolutionY,
                int imageChoice) {
            this.scaleFactor = scaleFactor;
            this.resolutionX = resolutionX;
            this.resolutionY = resolutionY;
            this.imageChoice = imageChoice;
        }

        public int compareTo(Resolution other) {
            if (scaleFactor > other.scaleFactor) return 1;
            else if (scaleFactor < other.scaleFactor) return -1;
            else return 0;
        }

        public String toString() {
            return "Resolution[Choice=" + imageChoice + ",scaleFactor=" + scaleFactor + "]";
        }
    }

    public VfsPath getVfsPath() {
        return vfsPath;
    }

    public void setVfsPath(VfsPath vfsPath) {
        this.vfsPath = vfsPath;
    }


    /**
     * Creates a {@link GridCoverage} for the provided {@link PlanarImage} using the  that was provided for this coverage.
     *
     * @param image contains the data for the coverage to create.
     * @return a {@link GridCoverage}
     * @throws IOException

    protected GridCoverage createCoverageFromImage(PlanarImage image) throws IOException {
        return createCoverageFromImage(image, null);
    }
*/
    /**
     * Setting Envelope, GridRange and CRS from the given {@code ImageReader}
     *
     * @param reader the {@code ImageReader} from which to retrieve metadata (if available) for
     *     setting properties
     * @throws IOException
     */
    protected void setCoverageProperties(ImageReader reader) throws IOException {
        // //
        //
        // Getting common metadata from GDAL
        //
        // //
        final IIOMetadata metadata = reader.getImageMetadata(0);
        if (!(metadata instanceof GDALCommonIIOImageMetadata)) {
            throw new DataSourceException(
                    "Unexpected error! Metadata should be an instance of the expected class: GDALCommonIIOImageMetadata.");
        }
        parseCommonMetadata((GDALCommonIIOImageMetadata) metadata);

        // //
        //
        // Envelope and CRS checks
        //
        // //
        if (this.crs == null) {
            LOGGER.info("crs not found, proceeding with default crs");
            this.crs = AbstractGridFormat.getDefaultCRS();
        }

        if (this.originalEnvelope == null) {
            throw new DataSourceException("Unable to compute the envelope for this coverage");
        }

        // setting the coordinate reference system for the envelope, just to make sure we set it
        this.originalEnvelope.setCoordinateReferenceSystem(this.crs);
    }
}
