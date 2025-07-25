# Installation of LMGC90

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
