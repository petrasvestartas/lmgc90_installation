# Installation of LMGC90

### Download LMGC90
https://lmgc90.pages-git-xen.lmgc.univ-montp2.fr/lmgc90_dev/downloads/lmgc90_user_2025.rc1.zip

### Create Environment
```bash
conda env create -f environment.yml
conda activate lmgc
```
### Build LMGC90
```bash
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DVENV_PATH=$CONDA_PREFIX ..
make -j$(nproc)
```
 
### Run Examples
```bash
export PYTHONPATH=/home/pv/brg/code_fortran/lmgc90_user_2025.rc1/build:$PYTHONPATH
cd examples/compas
python command.py
```

## Ubuntu


<img width="1693" height="1281" alt="image" src="https://github.com/user-attachments/assets/9a4c4711-8879-4201-9a69-907ed89e0271" />

Missing HDF5 with Fortran Support:

```bash
sudo apt update
sudo apt install libhdf5-dev libhdf5-fortran-102
cd build
cmake ..
make
```
