help([==[

  Description
  ===========
  LISF

  ]==])

-- Define the name and version of the module
whatis("Name: LISF")
whatis("Version: 0.1")
whatis("Description: Land Information System Framework.")

-- Load dependency modules
depends_on("intel/2023b")
depends_on("Python/3.11.5-GCCcore-13.2.0")
depends_on("Perl/5.38.0-GCCcore-13.2.0")
depends_on("HDF/4.2.16-2-GCCcore-13.2.0-w-fortran-no-netCDF")
depends_on("HDF5/1.14.3-iimpi-2023b")
depends_on("netCDF/4.9.2-iimpi-2023b")
depends_on("ecCodes/2.31.0-iimpi-2023b")
depends_on("ESMF/8.6.1-intel-2023b")
-- depends_on("HDF-EOS2/3.0-GCCcore-13.2.0")
-- depends_on("GDAL/3.9.0-intel-2023b")

-- Set environment variables
setenv("LDT_ARCH", "linux_ifc")
setenv("LIS_ARCH", "linux_ifc")
setenv("LDT_FC", "mpiifx")
setenv("LIS_FC", "mpiifx")
setenv("LDT_CC", "mpiicx")
setenv("LIS_CC", "mpiicx")
setenv("LDT_MODESMF", os.getenv("EBROOTESMF") .. "/mod")
setenv("LIS_MODESMF", os.getenv("EBROOTESMF") .. "/mod")
setenv("LDT_LIBESMF", os.getenv("EBROOTESMF") .. "/lib")
setenv("LIS_LIBESMF", os.getenv("EBROOTESMF") .. "/lib")
setenv("LDT_JASPER", os.getenv("EBROOTJASPER"))
setenv("LIS_JASPER", os.getenv("EBROOTJASPER"))
-- setenv("LDT_OPENJPEG", os.getenv("EBROOTOPENJPEG"))
-- setenv("LIS_OPENJPEG", os.getenv("EBROOTOPENJPEG"))
setenv("LDT_ECCODES", os.getenv("EBROOTECCODES"))
setenv("LIS_ECCODES", os.getenv("EBROOTECCODES"))
setenv("LDT_NETCDF", os.getenv("EBROOTNETCDF"))
setenv("LIS_NETCDF", os.getenv("EBROOTNETCDF"))
setenv("LDT_HDF4", os.getenv("EBROOTHDF"))
setenv("LIS_HDF4", os.getenv("EBROOTHDF"))
setenv("LDT_HDF5", os.getenv("EBROOTHDF5"))
setenv("LIS_HDF5", os.getenv("EBROOTHDF5"))
-- setenv("LDT_HDFEOS", os.getenv("EBROOTHDFMINEOS2"))
-- setenv("LIS_HDFEOS", os.getenv("EBROOTHDFMINEOS2"))
-- setenv("LDT_GDAL", os.getenv("EBROOTGDAL"))
-- setenv("LIS_GDAL", os.getenv("EBROOTGDAL"))
-- setenv("LDT_LIBGEOTIFF", os.getenv("EBROOTGEOTIFF"))
