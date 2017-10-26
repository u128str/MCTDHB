#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cufft.h>


#include <iostream>

// #include <complex>


#include <cuda_runtime.h>

#define imin(a,b) (a<b?a:b)     

 

/*--------- function called from main fortran programn ---------------*/
// -------------  Onli 1D, 2D or 3D  cuFFT transform !!!!!!!!!!!!!!!!!!!


extern "C" void kernel_imestfft_t_(cuDoubleComplex *wtilde,int *NX,int *NY,int *NZ,  int *Xcase)
{
           cuDoubleComplex  *psi_d;

        int NDX = *NX;
        int NDY = *NY;
        int NDZ = *NZ;
        int NDIM= NDX*NDY*NDZ;
        int dcase= *Xcase;


//----------------------------------------------------------------------------------------------
//                 printf("IN T GPU _Transform  \n");  

         cufftHandle plan;
         cudaMalloc( (void **)&psi_d, sizeof(cuDoubleComplex) * NDIM );
          cudaMemcpy( psi_d, wtilde, sizeof(cuDoubleComplex)*NDIM, cudaMemcpyHostToDevice );

 //                printf("Kernel GPU _Transform  \n");  

      if (dcase == 1)  cufftPlan1d(&plan, NDX, CUFFT_Z2Z, 1); 
     
     if (dcase == 2)   cufftPlan2d(&plan, NDY, NDX,      CUFFT_Z2Z); 
   
     if (dcase == 3)   cufftPlan3d(&plan, NDZ, NDY, NDX, CUFFT_Z2Z); 

          cufftExecZ2Z(plan,  psi_d,  psi_d, CUFFT_FORWARD);

 
        /*  copy vectors from GPU to CPU   */
           
            cudaMemcpy(wtilde, psi_d, sizeof(cuDoubleComplex) * NDIM, cudaMemcpyDeviceToHost);

           cufftDestroy(plan);
           cudaFree(psi_d);
        return;

}
