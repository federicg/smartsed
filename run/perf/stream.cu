// Minimal STREAM copy/triad to measure achievable DRAM bandwidth (double precision).
#include <cstdio>
#include <cuda_runtime.h>
__global__ void copy(const double* a, double* c, size_t n){ for(size_t i=blockIdx.x*(size_t)blockDim.x+threadIdx.x;i<n;i+=(size_t)gridDim.x*blockDim.x) c[i]=a[i]; }
__global__ void triad(double* a, const double* b, const double* c, double s, size_t n){ for(size_t i=blockIdx.x*(size_t)blockDim.x+threadIdx.x;i<n;i+=(size_t)gridDim.x*blockDim.x) a[i]=b[i]+s*c[i]; }
int main(){
  const size_t n = 1ull<<27; double *a,*b,*c; cudaMalloc(&a,n*8); cudaMalloc(&b,n*8); cudaMalloc(&c,n*8);
  cudaMemset(a,0,n*8); cudaMemset(b,0,n*8); cudaMemset(c,0,n*8);
  int sms; cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, 0);
  dim3 g(sms*32), t(256); cudaEvent_t e0,e1; cudaEventCreate(&e0); cudaEventCreate(&e1);
  float best_copy=1e9f, best_triad=1e9f, ms;
  for(int r=0;r<20;r++){
    cudaEventRecord(e0); copy<<<g,t>>>(a,c,n); cudaEventRecord(e1); cudaEventSynchronize(e1); cudaEventElapsedTime(&ms,e0,e1); if(r&&ms<best_copy) best_copy=ms;
    cudaEventRecord(e0); triad<<<g,t>>>(a,b,c,0.4,n); cudaEventRecord(e1); cudaEventSynchronize(e1); cudaEventElapsedTime(&ms,e0,e1); if(r&&ms<best_triad) best_triad=ms;
  }
  printf("copy  %.1f GB/s\ntriad %.1f GB/s\n", 2.0*n*8/best_copy/1e6, 3.0*n*8/best_triad/1e6);
  return cudaGetLastError();
}
