#ifndef MY_BROADCAST_HPP
#define MY_BROADCAST_HPP

#include <mpi.h>

template <typename T>
void my_broadcast(T* data, int count, int root, MPI_Comm comm){
    MPI_Datatype mpi_type;
    if (typeid(T) == typeid(char))          mpi_type = MPI_CHAR;
    else if (typeid(T) == typeid(int))      mpi_type = MPI_INT;
    else if (typeid(T) == typeid(float))    mpi_type = MPI_FLOAT;
    else if (typeid(T) == typeid(double))   mpi_type = MPI_DOUBLE;
    else {
        mpi_type = MPI_BYTE;
        count *= sizeof(T);
    }

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const int tag = 0;
    if (rank == root){
        for (int i = 0; i < size; ++i) {
            if (i != root) {
                MPI_Send(data, count, mpi_type, i, tag, comm);
            }
        }
    } else{
        MPI_Recv(data, count, mpi_type, root, tag, comm, MPI_STATUS_IGNORE);
    }
}

template void my_broadcast<char>(char*, int, int, MPI_Comm);
template void my_broadcast<int>(int*, int, int, MPI_Comm);
template void my_broadcast<float>(float*, int, int, MPI_Comm);
template void my_broadcast<double>(double*, int, int, MPI_Comm);

#endif