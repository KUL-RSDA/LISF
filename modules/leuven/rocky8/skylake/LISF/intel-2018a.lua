help([==[

  Description
  ===========
  LISF
  DEPRACATED! compiler too old

  ]==])

-- Define the name and version of the module
whatis("Name: LISF")
whatis("Version: 0.1")
whatis("Description: Land Information System Framework.")

-- Load dependency modules
depends_on("intel/2018a")
depends_on("Python/3.6.4-intel-2018a")
-- depends_on("Python/2.7.14-GCCcore-6.4.0-bare")
depends_on("Perl/5.28.1-GCCcore-6.4.0")
depends_on("HDF/4.2.14-intel-2018a-w-fortran-no-netcdf")
depends_on("HDF5/1.10.1-intel-2018a")
depends_on("netCDF/4.6.0-intel-2018a")
-- depends_on("JasPer/2.0.14-GCCcore-6.4.0")
-- depends_on("grib_api/1.24.0-intel-2018a")
depends_on("ESMF/7.1.0r-intel-2018a")
-- depends_on("zlib/1.2.11-GCCcore-6.4.0")
-- depends_on("Szip/2.1.1-GCCcore-6.4.0")
-- depends_on("libxml2/2.9.7-GCCcore-6.4.0")
depends_on("HDF-EOS2/20.1.00-intel-2018a-HDF4-w-fortran")
depends_on("GDAL/2.4.1-intel-2018a-Python-3.6.4")
-- depends_on("GDAL/2.2.3-intel-2018a-Python-2.7.14")

-- Set environment variables
setenv("LDT_ARCH", "linux_ifc")
setenv("LIS_ARCH", "linux_ifc")
setenv("LDT_SPMD", "parallel")
setenv("LIS_SPMD", "parallel")
setenv("LDT_FC", "mpiifort")
setenv("LIS_FC", "mpiifort")
setenv("LDT_CC", "mpiicc")
setenv("LIS_CC", "mpiicc")
setenv("LDT_MODESMF", os.getenv("EBROOTESMF") .. "/mod")
setenv("LIS_MODESMF", os.getenv("EBROOTESMF") .. "/mod")
setenv("LDT_LIBESMF", os.getenv("EBROOTESMF") .. "/lib")
setenv("LIS_LIBESMF", os.getenv("EBROOTESMF") .. "/lib")
-- setenv("LDT_OPENJPEG", os.getenv("EBROOTOPENJPEG"))
-- setenv("LIS_OPENJPEG", os.getenv("EBROOTOPENJPEG"))
-- setenv("LDT_ECCODES", os.getenv("EBROOTECCODES"))
-- setenv("LIS_ECCODES", os.getenv("EBROOTECCODES"))
setenv("LDT_NETCDF", os.getenv("EBROOTNETCDF"))
setenv("LIS_NETCDF", os.getenv("EBROOTNETCDF"))
setenv("LDT_HDF4", os.getenv("EBROOTHDF"))
setenv("LIS_HDF4", os.getenv("EBROOTHDF"))
setenv("LDT_HDF5", os.getenv("EBROOTHDF5"))
setenv("LIS_HDF5", os.getenv("EBROOTHDF5"))
setenv("LDT_HDFEOS", os.getenv("EBROOTHDFMINEOS2"))
setenv("LIS_HDFEOS", os.getenv("EBROOTHDFMINEOS2"))
setenv("LDT_GDAL", os.getenv("EBROOTGDAL"))
setenv("LIS_GDAL", os.getenv("EBROOTGDAL"))
-- setenv("LDT_LIBGEOTIFF", os.getenv("EBROOTGEOTIFF"))

-- setenv("LDT_JASPER", os.getenv("EBROOTJASPER"))
-- setenv("LIS_JASPER", os.getenv("EBROOTJASPER"))
-- setenv("LDT_GRIBAPI", os.getenv("EBROOTGRIB_API"))
-- setenv("LIS_GRIBAPI", os.getenv("EBROOTGRIB_API"))
