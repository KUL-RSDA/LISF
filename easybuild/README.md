There are several options to compile and run the LISF components LDT and LIS.

## EasyBuild

### Building

Start an interactive job on a (compute) node with an architecture similar to the nodes where you would like to run your software. Request sufficient cores (e.g., 8) using the `-c` flag.

Next, source the `load_easybuild.sh` script. This will ask you where to install and set your EasyBuild environment accordingly. 
```bash
source /data/leuven/314/vsc31497/software/src/HPC/tools/stable/easybuild/load_easybuild.sh
```
A directory in your `$VSC_DATA` is proposed. Please accept it unless you know what you are doing.

Now, you are ready to start building LIS or LDT. First, make sure you are in the top directory of LISF, so not in the `ldt` or the `lis` subdirectory. Then run:
```bash
eb easybuild/LISF-LDT-dev-intel-2023b.eb
```
for LDT, or:
```bash
eb easybuild/LISF-LIS-dev-intel-2023b.eb
```
for LIS. If everything goes right, this should build and install LDT/LIS as currently found in the LISF directory.

### Running

On a (compute) node that is similar in architecture to where you built the software, add the path where the software was installed to your module path. You can do this by running
```bash
module use "${VSC_DATA}/software/apps/${VSC_INSTITUTE_CLUSTER}/${VSC_OS_LOCAL}/${VSC_ARCH_LOCAL}${VSC_ARCH_SUFFIX}/modules/all"
```
in case you accepted the default path that was proposed by the `load_easybuild.sh` script. For your own convenience, you could put this line in your `.bashrc` such that the LISF modules you installed are always readily available to you when you log in to the cluster.

Running
```bash
module avail LISF
```
will now show the LISF components you have installed with the above procedure. Notice that the **git tag** and **git hash** are present in the version suffix of the installed LISF modules. This allows you to easily keep multiple LISF versions next to each other.

## `compile` script

### Building

Run
```bash
source compile
```

This is a convenience script that underlyingly calls: `source load_modules`, `configure` and `./compile -j <N>` in the right directories with the settings that you will have provided after having answered some questions.

### Running

First, load the modules needed by LISF:
```bash
source load_modules
```
Then run the executable `LDT` or `LIS` as found in the `make` directory of the corresponding component. You can move the executables to wherever you like.

## Bare metal

### Building

Start by loading the modules needed by LISF:
```bash
source load_modules
```

Now, proceed with running `configure` and `./compile -j <N>` as found in the respective component's directory.

### Running

See [above](#running-1)
