package com.joshfix.gdalvfs.geotools;

import com.joshfix.gdalvfs.geotools.path.VfsPath;

import javax.imageio.stream.ImageInputStreamImpl;
import java.io.IOException;

/**
 * @author joshfix
 * Created on 2019-08-01
 */
public class GdalVfsImageInputStream extends ImageInputStreamImpl {

    private VfsPath vfsPath;

    public GdalVfsImageInputStream(VfsPath vfsPath) {
        this.vfsPath = vfsPath;
    }

    public VfsPath getVfsPath() {
        return vfsPath;
    }

    @Override
    public int read() throws IOException {
        return 1;
    }

    @Override
    public int read(byte[] b, int off, int len) throws IOException {
        b = new byte[len];
        return len;
        //return 1;
    }

}
