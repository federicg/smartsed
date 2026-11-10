// Peak double-precision FLOP rate of the GPU: register-only FMA chains (no memory traffic).
// 8 independent chains per thread hide the FMA latency; 1 FMA = 2 flops.
#include <cstdio>
#include <cuda_runtime.h>
__global__ void fma_loop(double* out, int iters){
  double a0=1,a1=2,a2=3,a3=4,a4=5,a5=6,a6=7,a7=8; const double m=1.0000001+threadIdx.x*1e-12, c=1e-9;
  for(int i=0;i<iters;i++){ a0=fma(a0,m,c);a1=fma(a1,m,c);a2=fma(a2,m,c);a3=fma(a3,m,c);a4=fma(a4,m,c);a5=fma(a5,m,c);a6=fma(a6,m,c);a7=fma(a7,m,c); }
  if(threadIdx.x==0) out[blockIdx.x]=a0+a1+a2+a3+a4+a5+a6+a7;
}
int main(){
  int sms; cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, 0);
  const int blocks=sms*8, threads=1024, iters=200000; double* out; cudaMalloc(&out,blocks*8);
  cudaEvent_t e0,e1; cudaEventCreate(&e0); cudaEventCreate(&e1); float ms, best=1e9f;
  for(int r=0;r<6;r++){ cudaEventRecord(e0); fma_loop<<<blocks,threads>>>(out,iters); cudaEventRecord(e1); cudaEventSynchronize(e1); cudaEventElapsedTime(&ms,e0,e1); if(r&&ms<best) best=ms; }
  printf("SMs %d  peak DP %.1f GFLOP/s\n", sms, 16.0*blocks*threads*iters/best/1e6);
  return cudaGetLastError();
}
