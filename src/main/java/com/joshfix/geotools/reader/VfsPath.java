package com.joshfix.geotools.reader;

/**
 * This could do some regex test to ensure it's a valid vfs format, however given that we're just passing strings into
 * gdal.Open methods, any file path is valid.  Only allowing vfs paths would artificially restrict capabilities.
 * What to do?
 *
 * @author joshfix Created on 2019-07-29 */
public class VfsPath {

    private String path;

    private VfsPath(){}

    public VfsPath(String path) {
        this.path = path;
    }

    public String getPath() {
        return path;
    }

    @Override
    public String toString() {
        return path;
    }
}
