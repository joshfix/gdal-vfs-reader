package com.joshfix.geotools.reader;

import it.geosolutions.imageio.gdalframework.GDALCommonIIOImageMetadata;
import it.geosolutions.imageio.gdalframework.GDALUtilities;
import it.geosolutions.imageio.maskband.DatasetLayout;
import org.gdal.gdal.Band;
import org.gdal.gdal.Dataset;
import org.gdal.gdal.gdal;
import org.gdal.gdalconst.gdalconstConstants;
import org.geotools.coverage.CoverageFactoryFinder;
import org.geotools.coverage.GridSampleDimension;
import org.geotools.coverage.TypeMap;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridCoverageFactory;
import org.geotools.coverage.grid.GridEnvelope2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.coverage.grid.io.AbstractGridCoverage2DReader;
import org.geotools.coverage.grid.io.GridCoverage2DReader;
import org.geotools.coverage.grid.io.OverviewPolicy;
import org.geotools.coverage.util.CoverageUtilities;
import org.geotools.data.*;
import org.geotools.geometry.GeneralEnvelope;
import org.geotools.geometry.PixelTranslation;
import org.geotools.referencing.CRS;
import org.geotools.referencing.operation.builder.GridToEnvelopeMapper;
import org.geotools.referencing.operation.transform.IdentityTransform;
import org.geotools.referencing.operation.transform.ProjectiveTransform;
import org.geotools.util.Utilities;
import org.geotools.util.factory.GeoTools;
import org.geotools.util.factory.Hints;
import org.opengis.coverage.ColorInterpretation;
import org.opengis.coverage.grid.Format;
import org.opengis.coverage.grid.GridCoverage;
import org.opengis.coverage.grid.GridEnvelope;
import org.opengis.parameter.GeneralParameterValue;
import org.opengis.parameter.ParameterDescriptor;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.datum.PixelInCell;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import javax.imageio.ImageReader;
import javax.media.jai.ImageLayout;
import javax.media.jai.PlanarImage;
import java.awt.*;
import java.awt.geom.AffineTransform;
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
 * Created on 2019-07-29
 */
public class GdalVfsReader implements GridCoverage2DReader {

    private VfsPath vfsPath;
    private Dataset dataset;
    private ConcurrentHashMap<String, GDALCommonIIOImageMetadata> datasetMetadataMap = new ConcurrentHashMap<>();

    /**The coverage factory producing a {@link GridCoverage} from an image */
    private GridCoverageFactory coverageFactory;

    /** Small number used for double comparisons */
    protected static double EPS = 1e-6;

    /** This contains the number of overviews.aaa */
    protected int numOverviews = 0;

    /** 2DGridToWorld math transform. */
    protected MathTransform raster2Model = null;

    /** Envelope read from file */
    protected GeneralEnvelope originalEnvelope = null;

    /** Coverage name */
    protected String coverageName = "geotools_coverage";

    /** Source to read from */
    protected Object source = null;

    /** Hints used by the {@link AbstractGridCoverage2DReader} subclasses. */
    protected Hints hints = GeoTools.getDefaultHints();

    /** Highest resolution available for this reader. */
    protected double[] highestRes = null;

    /** Temp variable used in many readers. */
    protected boolean closeMe;

    /** In case we are trying to read from a GZipped file this will be set to true. */
    protected boolean gzipped;

    /** The original {@link GridEnvelope} for the {@link GridCoverage2D} of this reader. */
    protected GridEnvelope originalGridRange = null;

    /** crs for this coverage */
    protected CoordinateReferenceSystem crs = null;

    private ConcurrentHashMap<String, Dataset> datasetsMap = new ConcurrentHashMap<>();

    /** number of subdatasets */
    private int subDatasetSize = -1;

    /**
     * list of childs subdatasets names (if any) contained into the source
     */
    private String[] datasetNames;

    private Map<String, ArrayList<Resolution>> resolutionsLevelsMap = new HashMap<>();

    /** Resolutions avialaible through an overviews based mechanism. */
    protected double[][] overViewResolutions = null;

    /** Coverage {@link DatasetLayout} containing information about Overviews and Mask management */
    protected DatasetLayout dtLayout;

    private static final Logger LOGGER =
            org.geotools.util.logging.Logging.getLogger(GdalVfsReader.class);

    static {
        gdal.AllRegister();

        gdal.SetConfigOption("CHECK_DISK_FREE_SPACE", "FALSE");
        // config options per https://trac.osgeo.org/gdal/wiki/CloudOptimizedGeoTIFF#HowtoreaditwithGDAL
        gdal.SetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "YES");
        gdal.SetConfigOption("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif");
        gdal.SetConfigOption("CPL_CURL_VERBOSE", "YES");
        gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO");
    }

    public GdalVfsReader(VfsPath vfsPath) {
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
        populateMetadata();
        try {
            getResolutionInfo();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (TransformException e) {
            e.printStackTrace();
        }
    }

    protected void populateMetadata() {
        dataset = GDALUtilities.acquireDataSet(vfsPath.getPath(), gdalconstConstants.GA_ReadOnly);
        /*
        List<String> translateOptions = new ArrayList<>();
        translateOptions.add("-ot");
        translateOptions.add("Byte");

        String name = "/vsimem/" + UUID.randomUUID();

        dataset = gdal.Translate(name, dataset, new TranslateOptions(new Vector(translateOptions)));
*/
        LOGGER.info(new StringBuilder("Opened dataset... last error: ").append(
                gdal.GetLastErrorMsg()).toString());
        if (dataset == null) {
            throw new RuntimeException("Unable to acquire dataset for path " + vfsPath.getPath() + ". Reason: "
                    +  gdal.GetLastErrorMsg());
        }
        System.out.println("Dataset Y size: " + dataset.getRasterYSize());
        System.out.println("Dataset X size: " + dataset.getRasterXSize());

        datasetsMap.put(vfsPath.getPath(), dataset);

        final java.util.List<String> subdatasets = dataset.GetMetadata_List(GDALUtilities.GDALMetadataDomain.SUBDATASETS);

        subDatasetSize = subdatasets.size() / 2;

        // Some formats supporting subdatasets may have no subdatasets.
        // As an instance, the HDF4ImageReader may read HDF4Images
        // which are single datasets containing no subdatasets.
        // Thus, theDataset is simply the main dataset.
        if (subDatasetSize == 0) {
            subDatasetSize = 1;
            datasetNames = new String[1];
            datasetNames[0] = vfsPath.getPath();
            GDALCommonIIOImageMetadata metadata = createDatasetMetadata(vfsPath.getPath());
            //GDALCommonIIOImageMetadata metadata = createDatasetMetadata(dataset, vfsPath.getPath());
            datasetMetadataMap.put(datasetNames[0], metadata);
            parseCommonMetadata(metadata);
        } else {
            datasetNames = new String[subDatasetSize + 1];
            for (int i = 0; i < subDatasetSize; i++) {
                final String subdatasetName = (subdatasets.get(i * 2));
                final int nameStartAt = subdatasetName.lastIndexOf("_NAME=") + 6;
                datasetNames[i] = subdatasetName.substring(nameStartAt);
            }
            datasetNames[subDatasetSize] = vfsPath.getPath();
            GDALCommonIIOImageMetadata metadata = createDatasetMetadata(dataset, datasetNames[subDatasetSize]);
            datasetMetadataMap.put(
                    datasetNames[subDatasetSize],
                    metadata);
            parseCommonMetadata(metadata);
        }
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
     * Build a proper {@link GDALCommonIIOImageMetadata} given the name of a
     * dataset. The default implementation return a
     * {@link GDALCommonIIOImageMetadata} instance.This method should be
     * overridden by the specialized {@link it.geosolutions.imageio.gdalframework.GDALImageReader} in case you need to
     * obtain a specific {@link GDALCommonIIOImageMetadata}'s subclass
     *
     * @param datasetName
     *                the name of the dataset
     */
    protected GDALCommonIIOImageMetadata createDatasetMetadata(final String datasetName) {
        return new GDALCommonIIOImageMetadata(datasetName);
    }

    /**
     * Build a proper {@link GDALCommonIIOImageMetadata} given an input dataset
     * as well as the file name containing such a dataset.
     */
    protected GDALCommonIIOImageMetadata createDatasetMetadata(final Dataset mainDataset, String mainDatasetFileName) {
        return new GDALCommonIIOImageMetadata(mainDataset, mainDatasetFileName, false);
    }

    /**
     * Given a {@link GDALCommonIIOImageMetadata} metadata object, retrieves several properties to
     * properly set envelope, gridrange and crs.
     *
     * @param metadata a {@link GDALCommonIIOImageMetadata} metadata instance from where to search
     *     needed properties.
     */
    private void parseCommonMetadata(final GDALCommonIIOImageMetadata metadata) {

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
        return null;
    }

    @Override
    public Object getSource() {
        return source;
    }

    @Override
    public String[] getMetadataNames() throws IOException {
        if (!checkName(coverageName)) {
            throw new IllegalArgumentException(
                    "The specified coverageName " + coverageName + "is not supported");
        }
        return getMetadataNames();
    }

    @Override
    public String[] getMetadataNames(String coverageName) throws IOException {
        return new String[0];
    }

    @Override
    public String getMetadataValue(String name) throws IOException {
        return null;
    }

    @Override
    public String getMetadataValue(String coverageName, String name) throws IOException {
        return null;
    }

    @Override
    public String[] listSubNames() throws IOException {
        return new String[0];
    }

    @Override
    public String[] getGridCoverageNames() throws IOException {
        return new String[]{coverageName};
    }

    @Override
    public int getGridCoverageCount() throws IOException {
        return 1;
    }

    @Override
    public String getCurrentSubname() throws IOException {
        return null;
    }

    @Override
    public boolean hasMoreGridCoverages() throws IOException {
        return false;
    }

    @Override
    public GridCoverage2D read(GeneralParameterValue[] parameters) throws IOException {
        GridGeometry2D gridGeometry2D = null;

        for (GeneralParameterValue parameter : parameters) {
            if (!(((org.geotools.parameter.Parameter)parameter).getValue() instanceof GridGeometry2D)) {
                continue;
            }
            gridGeometry2D = (GridGeometry2D)((org.geotools.parameter.Parameter)parameter).getValue();
        }

        if (gridGeometry2D == null) {
            throw new IOException("COG support requires the READ_GRIDGEOMETRY2D parameter.");
        }

        GridEnvelope2D gridEnvelope2D = gridGeometry2D.getGridRange2D();
        int minX = Math.min(gridEnvelope2D.getLow().x, originalGridRange.getLow().getCoordinateValues()[0]);
        int minY = Math.min(gridEnvelope2D.getLow().y, originalGridRange.getLow().getCoordinateValues()[1]);
        int maxX = Math.min(gridEnvelope2D.getHigh().x, originalGridRange.getHigh().getCoordinateValues()[0]);
        int maxY = Math.min(gridEnvelope2D.getHigh().y, originalGridRange.getHigh().getCoordinateValues()[1]);

        System.out.println("minX: " + minX);
        System.out.println("minY: " + minY);
        System.out.println("maxX: " + maxX);
        System.out.println("maxY: " + maxY);

        SampleModel sampleModel = datasetMetadataMap.get(vfsPath.getPath()).getSampleModel();
        ColorModel colorModel = GDALUtilities.buildColorModel(sampleModel);


        /*
        BufferedImage bi = getBufferedImageFromDataset(result);
        GridCoverage coverage = createCoverageFromImage(PlanarImage.wrapRenderedImage(bi));
        */

        int[] bands = new int[sampleModel.getNumBands()];
        for (int i = 0; i < sampleModel.getNumBands(); i++) {
            bands[i] = i + 1;
        }

        GDALCommonIIOImageMetadata itemMetadata = datasetMetadataMap.get(vfsPath.getPath());

        Raster raster = readDatasetRaster(
                itemMetadata,
                new Rectangle(minX, minY, maxX, maxY),
                new Rectangle(minX, minY, maxX, maxY),
                //new int[]{1},
                bands,
                sampleModel
                );
        WritableRaster writableRaster = raster.createCompatibleWritableRaster();
        writableRaster.setDataElements(0, 0, raster);
        BufferedImage image = new BufferedImage(colorModel, writableRaster, colorModel.isAlphaPremultiplied(), null);

        GridCoverage coverage = createCoverageFromImage(PlanarImage.wrapRenderedImage(image));


        return (GridCoverage2D) coverage;
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

    @Override
    public void skip() throws IOException {

    }

    @Override
    public void dispose() throws IOException {

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
    public DatasetLayout getDatasetLayout() {
        return null;
    }

    @Override
    public DatasetLayout getDatasetLayout(String coverageName) {
        return null;
    }

    @Override
    public ImageLayout getImageLayout() throws IOException {
        return null;
    }

    @Override
    public ImageLayout getImageLayout(String coverageName) throws IOException {
        return null;
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
    private Raster readDatasetRaster(
            GDALCommonIIOImageMetadata itemMetadata,
            Rectangle srcRegion,
            Rectangle dstRegion,
            int[] selectedBands,
            SampleModel destSampleModel) throws IOException {

        SampleModel destSm = destSampleModel != null ? destSampleModel : itemMetadata.getSampleModel();

        final Dataset dataset = datasetsMap.get(itemMetadata.getDatasetName());
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
                        int length = byteBands[0].length;
                        // I get float values from the ByteBuffer using a view
                        // of the ByteBuffer as a FloatBuffer
                        // It is worth to create the view outside the loop.
                        float[] floats = new float[nBands * pixels];
                        bands[0].order(ByteOrder.nativeOrder());
                        //bands[0].order(ByteOrder.LITTLE_ENDIAN);
                        //bands[0].order(ByteOrder.BIG_ENDIAN);
                        final FloatBuffer buff = bands[0].asFloatBuffer();
/*
                        byte[] bytes = bands[0].array();
                        for (int i = 0; i + 4 < bytes.length; i = i + 4) {

                            int asInt = (bytes[i] & 0xFF)
                                    | ((bytes[i + 1] & 0xFF) << 8)
                                    | ((bytes[i + 2] & 0xFF) << 16)
                                    | ((bytes[i + 3] & 0xFF) << 24);
                            float asFloat = Float.intBitsToFloat(asInt);
                            floats[i] = asFloat;
                        }
                        imgBuffer = new DataBufferFloat(floats, nBands * pixels);
*/


                        buff.get(floats, 0, nBands * pixels);
                        imgBuffer = new DataBufferFloat(floats, nBands * pixels);
                        System.out.println("float array length: " + floats.length);

                        int nanCount = 0;
                        int nonZeroCount = 0;
                        int zeroCount = 0;
                        int unknownCount = 0;
                        for (int i = 0; i < floats.length; i++) {
                            if (Float.isNaN(floats[i])) {
                                nanCount++;
                            } else if (floats[i] != 0) {
                                nonZeroCount++;
                            } else if (floats[i] == 0) {
                                zeroCount++;
                            } else {
                                unknownCount++;
                            }
                        }

                        System.out.println("Float counts:");
                        System.out.println("NaN count: " + nanCount);
                        System.out.println("non-zero count: " + nonZeroCount);
                        System.out.println("zero count: " + zeroCount);
                        System.out.println("unknown count: " + unknownCount);
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
     * Creates a {@link GridCoverage} for the provided {@link PlanarImage} using the {@link
     * #raster2Model} that was provided for this coverage.
     *
     * <p>This method is vital when working with coverages that have a raster to model
     * transformation that is not a simple scale and translate.
     *
     * @param image contains the data for the coverage to create.
     * @param raster2Model is the {@link MathTransform} that maps from the raster space to the model
     *     space.
     * @return a {@link GridCoverage}
     * @throws IOException
     */
    protected GridCoverage createCoverageFromImage(PlanarImage image, MathTransform raster2Model)
            throws IOException {
        // creating bands
        final SampleModel sm = image.getSampleModel();
        final ColorModel cm = image.getColorModel();
        final int numBands = sm.getNumBands();
        final GridSampleDimension[] bands = new GridSampleDimension[numBands];
        // setting bands names.
        Set<String> bandNames = new HashSet<String>();
        for (int i = 0; i < numBands; i++) {
            final ColorInterpretation colorInterpretation = TypeMap.getColorInterpretation(cm, i);
            // make sure we create no duplicate band names
            String bandName;
            if (colorInterpretation == null
                    || colorInterpretation == ColorInterpretation.UNDEFINED
                    || bandNames.contains(colorInterpretation.name())) {
                bandName = "Band" + (i + 1);
            } else {
                bandName = colorInterpretation.name();
            }
            bands[i] = new GridSampleDimension(bandName);
        }

        // creating coverage
        if (raster2Model != null) {
            return coverageFactory.create(
                    coverageName, image, crs, raster2Model, bands, null, null);
        }

        return coverageFactory.create(
                // TODO what envelope???
                //coverageName, image, new GeneralEnvelope(coverageEnvelope), bands, null, null);
                coverageName, image, new GeneralEnvelope(originalEnvelope), bands, null, null);
    }

    /**
     * Creates a {@link GridCoverage} for the provided {@link PlanarImage} using the  that was provided for this coverage.
     *
     * @param image contains the data for the coverage to create.
     * @return a {@link GridCoverage}
     * @throws IOException
     */
    protected GridCoverage createCoverageFromImage(PlanarImage image) throws IOException {
        return createCoverageFromImage(image, null);
    }

}
