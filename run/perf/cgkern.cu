// Effective DRAM bandwidth of the CG building blocks used by smartsed
// (cuSPARSE SpMV + SpSV with an IC0-shaped lower factor, cuBLAS dot/axpy/nrm2)
// on a 5-point matrix of N unknowns. Bytes = minimum traffic; time = CUDA events.
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>
#define CK(x) do{ if((x)!=0){printf("err %s:%d\n",__FILE__,__LINE__); exit(1);} }while(0)
int main(int argc,char**argv){
  int nx = argc>1? atoi(argv[1]) : 1030;           // nx*nx unknowns
  int N = nx*nx; double peak = argc>2? atof(argv[2]) : 551.0;
  // A: 5-point SPD (4 on diag, -1 off), L: lower part incl. diag (IC0 sparsity)
  std::vector<int> ar, ac, lr, lc; std::vector<double> av, lv; ar.push_back(0); lr.push_back(0);
  for(int i=0;i<nx;i++) for(int j=0;j<nx;j++){ int r=i*nx+j;
    auto add=[&](int c,double v){ ac.push_back(c); av.push_back(v); if(c<=r){ lc.push_back(c); lv.push_back(c==r?2.0:-0.5);} };
    if(i>0) add(r-nx,-1); if(j>0) add(r-1,-1); add(r,4.2); if(j<nx-1) add(r+1,-1); if(i<nx-1) add(r+nx,-1);
    ar.push_back(ac.size()); lr.push_back(lc.size()); }
  int nnz=ac.size(), lnnz=lc.size();
  int *dar,*dac,*dlr,*dlc; double *dav,*dlv,*x,*y,*z;
  CK(cudaMalloc(&dar,(N+1)*4)); CK(cudaMalloc(&dac,nnz*4)); CK(cudaMalloc(&dav,nnz*8));
  CK(cudaMalloc(&dlr,(N+1)*4)); CK(cudaMalloc(&dlc,lnnz*4)); CK(cudaMalloc(&dlv,lnnz*8));
  CK(cudaMalloc(&x,N*8)); CK(cudaMalloc(&y,N*8)); CK(cudaMalloc(&z,N*8));
  cudaMemcpy(dar,ar.data(),(N+1)*4,cudaMemcpyHostToDevice); cudaMemcpy(dac,ac.data(),nnz*4,cudaMemcpyHostToDevice); cudaMemcpy(dav,av.data(),nnz*8,cudaMemcpyHostToDevice);
  cudaMemcpy(dlr,lr.data(),(N+1)*4,cudaMemcpyHostToDevice); cudaMemcpy(dlc,lc.data(),lnnz*4,cudaMemcpyHostToDevice); cudaMemcpy(dlv,lv.data(),lnnz*8,cudaMemcpyHostToDevice);
  std::vector<double> h(N,1.0); cudaMemcpy(x,h.data(),N*8,cudaMemcpyHostToDevice); cudaMemcpy(y,h.data(),N*8,cudaMemcpyHostToDevice);
  cublasHandle_t cb; cublasCreate(&cb); cusparseHandle_t cs; cusparseCreate(&cs);
  cusparseSpMatDescr_t A,L; cusparseDnVecDescr_t vx,vy,vz;
  CK(cusparseCreateCsr(&A,N,N,nnz,dar,dac,dav,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO,CUDA_R_64F));
  CK(cusparseCreateCsr(&L,N,N,lnnz,dlr,dlc,dlv,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO,CUDA_R_64F));
  cusparseFillMode_t fm=CUSPARSE_FILL_MODE_LOWER; cusparseDiagType_t dt=CUSPARSE_DIAG_TYPE_NON_UNIT;
  cusparseSpMatSetAttribute(L,CUSPARSE_SPMAT_FILL_MODE,&fm,sizeof(fm)); cusparseSpMatSetAttribute(L,CUSPARSE_SPMAT_DIAG_TYPE,&dt,sizeof(dt));
  cusparseCreateDnVec(&vx,N,x,CUDA_R_64F); cusparseCreateDnVec(&vy,N,y,CUDA_R_64F); cusparseCreateDnVec(&vz,N,z,CUDA_R_64F);
  double one=1,zero=0; size_t bs; void* buf;
  cusparseSpMV_bufferSize(cs,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,A,vx,&zero,vy,CUDA_R_64F,CUSPARSE_SPMV_ALG_DEFAULT,&bs); cudaMalloc(&buf,bs+8);
  cusparseSpSVDescr_t sl,slt; cusparseSpSV_createDescr(&sl); cusparseSpSV_createDescr(&slt); size_t b1,b2; void *buf1,*buf2;
  cusparseSpSV_bufferSize(cs,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,L,vx,vz,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,sl,&b1); cudaMalloc(&buf1,b1+8);
  cusparseSpSV_bufferSize(cs,CUSPARSE_OPERATION_TRANSPOSE,&one,L,vz,vy,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,slt,&b2); cudaMalloc(&buf2,b2+8);
  CK(cusparseSpSV_analysis(cs,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,L,vx,vz,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,sl,buf1));
  CK(cusparseSpSV_analysis(cs,CUSPARSE_OPERATION_TRANSPOSE,&one,L,vz,vy,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,slt,buf2));
  cudaEvent_t e0,e1; cudaEventCreate(&e0); cudaEventCreate(&e1); double res; const int R=50;
  auto T=[&](const char* name,double bytes,auto f){ f(); cudaDeviceSynchronize(); cudaEventRecord(e0); for(int k=0;k<R;k++) f(); cudaEventRecord(e1); cudaEventSynchronize(e1);
    float ms; cudaEventElapsedTime(&ms,e0,e1); double us=ms*1e3/R, gbs=bytes/(us*1e3);
    printf("%-22s %9.1f us  %7.1f GB/s  %5.1f%% of %.0f\n",name,us,gbs,100*gbs/peak,peak); };
  double vN=8.0*N;
  printf("N=%d nnz(A)=%d nnz(L)=%d  vector=%.1f MB\n",N,nnz,lnnz,vN/1e6);
  T("SpMV A*x", 12.0*nnz+4.0*(N+1)+2*vN, [&]{ cusparseSpMV(cs,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,A,vx,&zero,vy,CUDA_R_64F,CUSPARSE_SPMV_ALG_DEFAULT,buf); });
  T("SpSV L (forward)", 12.0*lnnz+4.0*(N+1)+2*vN, [&]{ cusparseSpSV_solve(cs,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,L,vx,vz,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,sl); });
  T("SpSV L^T (backward)", 12.0*lnnz+4.0*(N+1)+2*vN, [&]{ cusparseSpSV_solve(cs,CUSPARSE_OPERATION_TRANSPOSE,&one,L,vz,vy,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,slt); });
  T("dot (host result)", 2*vN, [&]{ cublasDdot(cb,N,x,1,y,1,&res); });
  T("nrm2 (host result)", vN, [&]{ cublasDnrm2(cb,N,x,1,&res); });
  T("axpy", 3*vN, [&]{ cublasDaxpy(cb,N,&one,x,1,y,1); });
  return 0;
}
