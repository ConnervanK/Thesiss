call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
cd C:\Users\Administrator\Thesis\Fast-GPR-FWI\src
nvcc -shared -o lib/pml_updates_e.dll lib/pml_updates_e.cu
nvcc -shared -o lib/pml_updates_h.dll lib/pml_updates_h.cu
nvcc -shared -o lib/back.dll lib/back.cu
nvcc -shared -o lib/uc.dll lib/uc.cu
nvcc -shared -o lib/fields_updates_gpu.dll lib/fields_updates_gpu.cu
nvcc -shared -o lib/sourcereceiver.dll lib/sourcereceiver.cu
