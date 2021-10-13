#include <cstdio>
#include <iostream>

// Device code
__global__ void MyKernel(int* devPtr, size_t pitch)
{
    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    
    int* row = (int*)((char*)devPtr + ty * pitch) + tx;
    *row *= 2;
}


int main()
{
    // Host code
    int* devPtr;
    size_t pitch;
    cudaError_t err; 
    err = cudaMallocPitch(&devPtr, &pitch,
            width * sizeof(float), height);
    std::cout << pitch << std::endl;
    if(err != cudaSuccess){
        fprintf(stderr, "%s", cudaGetErrorString(err));
    }else{
        std::cout << "coppied successfully" <<std::endl;
    }
    err = cudaMemcpy2D(devPtr, pitch, &M, width*sizeof(int), width*sizeof(int), height, cudaMemcpyDefault); 
    if(err != cudaSuccess){
        fprintf(stderr, "%s", cudaGetErrorString(err));
    }else{
        std::cout << "coppied successfully" <<std::endl;
    }
    dim3 block_dim(20,20);
    MyKernel<<<1, block_dim>>>(devPtr, pitch);

    err = cudaMemcpy2D(&C, width*sizeof(int), devPtr, pitch, width*sizeof(int), height, cudaMemcpyDeviceToHost); 
    if(err != cudaSuccess){
        fprintf(stderr, "%s", cudaGetErrorString(err));
    }
    print_arr(C);

}
