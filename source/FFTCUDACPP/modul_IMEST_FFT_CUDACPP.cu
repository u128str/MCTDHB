#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cufft.h>


#include <iostream>

 #include <complex>


#include <cuda_runtime.h>

#define imin(a,b) (a<b?a:b)     

 


__global__ void  prod( cuDoubleComplex *vec2, cuDoubleComplex *a, int NDIM)
 {    
            
          int id = threadIdx.x  + blockIdx.x* blockDim.x;
       cuDoubleComplex  d;
     // make_cuDoubleComplex(1.0/NDIM,0.0); 
         d.x=(1.0/NDIM);
         d.y=0.0;
      
      while(id<NDIM)
      {
                        
            vec2[id]= cuCmul (cuCmul(a[id], d) ,  vec2[id] );

  //        vec2[id]= cuCmul (a[id] / NDIM ,  vec2[id] );

 
         id += blockDim.x*gridDim.x;

      }
  } 


/*--------- function called from main fortran programn ---------------*/

extern "C" void kernel_imestfft_(cuDoubleComplex *wsl,cuDoubleComplex *psi,int *NX,int *NY,int *NZ,  int *Xcase,cuDoubleComplex *vtilde)
{
           cuDoubleComplex  *psi_d, *vtilde_d;

        int NDX = *NX;
        int NDY = *NY;
        int NDZ = *NZ;
        int NDIM= NDX*NDY*NDZ;
        int dcase= *Xcase;

// ---- CUDA  variables: ----------------------------------------------------------------------- 

            const int threads = imin(64,NDIM);   //dim3 (16,1,1);
            const int blocks  = imin(32, NDIM/threads );  // ;  //imin(8, (M+threads)/threads ); 
//----------------------------------------------------------------------------------------------

         cufftHandle plan;
         cudaMalloc( (void **)&psi_d, sizeof(cuDoubleComplex) * NDIM );
         cudaMalloc( (void **)&vtilde_d, sizeof(cuDoubleComplex) *NDIM);

         cudaMemcpy( psi_d,       psi, sizeof(cuDoubleComplex)*NDIM, cudaMemcpyHostToDevice );
         cudaMemcpy( vtilde_d, vtilde, sizeof(cuDoubleComplex)*NDIM, cudaMemcpyHostToDevice );

//                 printf("Kernel GPU   \n");  

      if (dcase == 1)  cufftPlan1d(&plan, NDX, CUFFT_Z2Z, 1); 
     
     if (dcase == 2)   cufftPlan2d(&plan, NDY, NDX,      CUFFT_Z2Z); 
   
     if (dcase == 3)   cufftPlan3d(&plan, NDZ, NDY, NDX, CUFFT_Z2Z); 

          cufftExecZ2Z(plan,  psi_d,  psi_d, CUFFT_FORWARD);

          prod<<<blocks,threads>>>( psi_d, vtilde_d, NDIM); 

          cufftExecZ2Z(plan,  psi_d,  psi_d, CUFFT_INVERSE);


        /*  copy vectors from GPU to CPU   */
           
            cudaMemcpy(wsl, psi_d, sizeof(cuDoubleComplex) * NDIM, cudaMemcpyDeviceToHost);

           cufftDestroy(plan);
           cudaFree(psi_d);
           cudaFree(vtilde_d);
       return;

}
