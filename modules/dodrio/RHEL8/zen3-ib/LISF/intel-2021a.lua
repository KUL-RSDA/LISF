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
local root = "/dodrio/scratch/projects/2022_200/project_input/rsda/src_code/nu-wrf_ekman_v11.0/baselibs/intel-intelmpi"
prepend_path("MODULEPATH", root)
depends_on("intel/2021a")
depends_on("GCC/10.3.0")
depends_on("libjpeg-turbo/2.0.6-GCCcore-10.3.0")
depends_on("zlib/1.2.11-GCCcore-10.3.0")
depends_on("libtirpc/1.3.2-GCCcore-10.3.0")
depends_on("JasPer/1.900.1-GCCcore-10.3.0")
depends_on("Szip/2.1.1-GCCcore-10.3.0")
depends_on("cURL/7.76.0-GCCcore-10.3.0")
depends_on("HDF5/1.10.7-iimpi-2021a")

-- Set environment variables
setenv("LDT_ARCH", "linux_ifc")
setenv("LIS_ARCH", "linux_ifc")
setenv("LDT_SPMD", "parallel")
setenv("LIS_SPMD", "parallel")
setenv("LDT_FC", "mpiifort")
setenv("LIS_FC", "mpiifort")
setenv("LDT_CC", "mpiicc")
setenv("LIS_CC", "mpiicc")
setenv("LDT_JASPER", os.getenv("EBROOTJASPER"))
setenv("LIS_JASPER", os.getenv("EBROOTJASPER"))
setenv("LDT_GRIBAPI", root .. "/grib_api/")
setenv("LIS_GRIBAPI", root .. "/grib_api/")
setenv("LDT_NETCDF", root .. "/netcdf4/")
setenv("LIS_NETCDF", root .. "/netcdf4/")
setenv("LDT_HDF4", root .. "/hdf4/")
setenv("LIS_HDF4", root .. "/hdf4/")
setenv("LDT_HDFEOS", root .. "/hdfeos/")
setenv("LIS_HDFEOS", root .. "/hdfeos/")
setenv("LDT_HDF5", os.getenv("EBROOTHDF5"))
setenv("LIS_HDF5", root .. "/hdf5/")
setenv("LDT_MODESMF", root .. "/esmf/mod/modO/Linux.intel.64.intelmpi.default/")
setenv("LIS_MODESMF", root .. "/esmf/mod/modO/Linux.intel.64.intelmpi.default/")
setenv("LDT_LIBESMF", root .. "/esmf/lib/libO/Linux.intel.64.intelmpi.default/")
setenv("LIS_LIBESMF", root .. "/esmf/lib/libO/Linux.intel.64.intelmpi.default/")

-- Add paths to the PATH environment variable
prepend_path("PATH", root .. "netcdf4/bin")
