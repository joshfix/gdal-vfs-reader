package com.joshfix.gdalvfs.geotools.path;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * Need to determine best way to handle/parse URLs to the appropriate virtual file system prefix
 *
 * @author joshfix Created on 2019-07-29 */
public class VfsPath {

    private String path;
    private URL url;
    private Prefix prefix;
    private static Map<Prefix, VfsPathConfigurer> configurers = new HashMap<>();

    static {
        configurers.put(Prefix.AZURE, new AzureVfsPathConfigurer());
    }

    private VfsPath(){}

    public VfsPath(String urlString) {
        try {
            url = new URL(urlString);
            determinePrefix();
            buildPath();
        } catch (MalformedURLException mue) {
            throw new IllegalArgumentException("Invalid url: " + urlString);
        }
    }

    public VfsPath(String urlString, Prefix prefix) {
        try {
            url = new URL(urlString);
            this.prefix = prefix;
            buildPath();
        } catch (MalformedURLException mue) {
            throw new IllegalArgumentException("Invalid url: " + urlString);
        }
    }

    public VfsPath(URL url, Prefix prefix) {
        this.url = url;
        this.prefix = prefix;
        buildPath();
    }

    public VfsPath(URL url) {
        this.url = url;
        determinePrefix();
        buildPath();
    }

    protected String buildPath() {
        if (configurers.keySet().contains(prefix)) {
            configurers.get(prefix).configure(this);
        } else {
            path = prefix + url.toString();
        }
        return path;
    }

    public void determinePrefix() {
        switch (url.getProtocol().toLowerCase()) {
            case "wasb":
            case "wasbs":
                prefix = Prefix.AZURE;
                break;
            case "s3":
                prefix = Prefix.S3;
                break;
            case "file":
                prefix = Prefix.FILE;
                break;
            case "http":
            case "https":
            default:
                prefix = Prefix.CURL;
                break;
        }
    }

    protected void setPath(String path) {
        this.path = path;
    }

    public String getPath() {
        return path;
    }

    public URL getUrl() {
        return url;
    }

    public Prefix getPrefix() {
        return prefix;
    }

    @Override
    public String toString() {
        return path;
    }

    @Override
    public boolean equals(Object other) {
        if (other instanceof VfsPath && other.toString().equals(path)) {
            return true;
        }
        return false;
    }

    public enum Prefix {
        FILE(""),
        ZIP("/vsizip/"),
        GZIP("/vsigzip/"),
        TAR("/vsitar/"),
        CURL("/vsicurl/"),
        CURL_STREAMING("/vsicurl_streaming/"),
        S3("/vsis3/"),
        S3_STREAMING("/vsis3_streaming/"),
        GOOGLE("/vsigs/"),
        GOOGLE_STREAMING("/vsigs_streaming/"),
        AZURE("/vsiaz/"),
        AZURE_STREAMING("/vsiaz_streaming/"),
        ALIBABA("/vsioss/"),
        ALIBABA_STREAMING("/vsioss_streaming/"),
        SWIFT("/vsiswift/"),
        SWIFT_STREAMING("/vsiswift_streaming/"),
        HADOOP("/vsihdfs/"),
        WEB_HADOOP("/vsiwebhdfs/"),
        STANDARD_IN("/vsistdin/"),
        STANDARD_OUT("/vsistdout/"),
        MEMORY("/vsimem/"),
        SUBFILE("/vsisubfile/"),
        SPARSE("/vsisparce/"),
        ENCRYPTED("/vsicrypt/");

        private String value;

        Prefix(String value) {
            this.value = value;
        }

        public String getValue() {
            return value;
        }

        @Override
        public String toString() {
            return String.valueOf(value);
        }

        public static Prefix fromValue(String text) {
            for (Prefix b : Prefix.values()) {
                if (String.valueOf(b.value).equals(text)) {
                    return b;
                }
            }
            return null;
        }

    }
}

