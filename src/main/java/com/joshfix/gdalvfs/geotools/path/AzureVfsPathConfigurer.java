package com.joshfix.gdalvfs.geotools.path;

import java.net.URL;

/**
 * Parses http://<container>@<storage_account>blob.core.windows.net/SomeDirectory/a.tif to
 * /vsiaz/<container>/SomeDirectory/a.tif
 * @author joshfix
 * Created on 2019-08-02
 */
public class AzureVfsPathConfigurer implements VfsPathConfigurer {

    @Override
    public void configure(VfsPath vfsPath) {
        URL url = vfsPath.getUrl();
        String authority = url.getAuthority().split("@")[0];
        String path = vfsPath.getPrefix() + authority + vfsPath.getUrl().getPath();
        vfsPath.setPath(path);
    }

}
