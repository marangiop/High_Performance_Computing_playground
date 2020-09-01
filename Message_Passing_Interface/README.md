####Install MPI
sudo apt install libopenmpi-dev
####Compile code
mpicc mpp_cw_submission.c pgmio.c -o mpi_code_executable
####Run code
./mpi_code_executable

