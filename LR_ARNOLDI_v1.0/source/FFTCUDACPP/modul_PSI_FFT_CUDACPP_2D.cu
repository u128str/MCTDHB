#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cufft.h>

#include <iostream>
#include <fstream>
#include <sstream> 

// #include <complex>


#include <cuda_runtime.h>
#include <cublas.h> 

//<cublas_v2.h>   
                //  old version <cublas.h> 
//-------------------------------------------------------------------------------------

 //---//----  CUDA variables: 

            const int threads = 256;  //dim3 (16,1,1);
            const int blocks  = 16;  //imin(8, (M+threads)/threads ); 


//-----------------------------------------------------------------------------------//

__global__ void mltpl_2D( cuDoubleComplex *vec2, cuDoubleComplex *ax,cuDoubleComplex *ay, int NX, int NY)
   {    
            
          int jd = threadIdx.x  + blockIdx.x* blockDim.x;

          //  for (int jd = 0; jd < NDY; ++jd ) 
    
            while(jd<NY)

            {
             for (int id = 0; id < NX; ++id ) 
             {
               int ind =id+NX*jd;
             vec2[ind]= cuCmul (cuCadd(ax[id], ay[jd] ),  vec2[ind] );
             }
            jd += blockDim.x*gridDim.x;

              __syncthreads(); 

            }
//------------------------------------------------------------------------------------//


   }






// function called from main fortran program

extern "C" void   kernel_fft_2d_(cufftDoubleComplex *psi,
                cuDoubleComplex *ortkx, cuDoubleComplex  *ortky, int *NX, int *NY)
           {
           
          cuDoubleComplex   *ortkx_d, *ortky_d;  //, *psiV_d;  // declare GPU vector copies cufftDoubleComplex  *psi_d, *ortkx_d;
    
          cufftDoubleComplex   *psi_d;         

           int NDX = *NX;       
           int NDY = *NY;  
           cufftHandle plan; 
 
  //  int  BATCH = 1;     // batch Number of transforms of size nx

   //         printf("Kernel CPP FFT 2D\n"); 
       // Allocate memory on GPU
       
          cudaMalloc( (void **)&psi_d, sizeof(cufftDoubleComplex) * NDX*NDY );
    
          cudaMalloc( (void **)&ortkx_d, sizeof(cuDoubleComplex) * NDX );
          cudaMalloc( (void **)&ortky_d, sizeof(cuDoubleComplex) * NDY );

   //----------  copy vectors from CPU to GPU
      
//      cudaMemcpy( psi_d, psi, sizeof(cufftDoubleComplex) * NX, cudaMemcpyHostToDevice );
 
             /* Initialize CUBLAS */
   
 //          cublasHandle_t handle;
          
   //        cublasCreate(&handle);    // cublasCreate_v2(&handle); 
          
     //    cublasSetVector(NDX*NDY,sizeof(cuDoubleComplex),psi,1, psi_d,1);
       //  cublasSetVector(NDX,sizeof(cuDoubleComplex),ortkx,1, ortkx_d,1);
         //cublasSetVector(NDY,sizeof(cuDoubleComplex),ortky,1, ortky_d,1);

         
         cudaMemcpy( psi_d, psi, sizeof(cufftDoubleComplex)*NDX*NDY, cudaMemcpyHostToDevice );

         cudaMemcpy( ortkx_d, ortkx, sizeof(cuDoubleComplex) * NDX, cudaMemcpyHostToDevice );
         cudaMemcpy( ortky_d, ortky, sizeof(cuDoubleComplex) * NDY, cudaMemcpyHostToDevice );
   
        /* Create a 2D FFT plan. */
  
             cufftPlan2d(&plan, NDY, NDX, CUFFT_Z2Z);

      /* Use the CUFFT plan to transform the signal in place. */
           
                       
            cufftExecZ2Z(plan,  psi_d,  psi_d, CUFFT_FORWARD);

   // call function on GPU 
               // mltpl_2D<<<blocks,threads>>>(psiV_d, psi_d, ortkx_d, ortky_d,NDX, NDY);

//              printf("Cufft forvard \n"); 
               
                mltpl_2D<<<blocks,threads>>>( psi_d, ortkx_d, ortky_d, NDX, NDY);

//             printf("Kernel 2D\n"); 
 
 //----------------------- kernel on CPU -----------------------------------------
//------------------------PROOF - rabotaet pravil'no -----------------------
/*
             cudaMemcpy(psi, psi_d, sizeof(cuDoubleComplex) * NDX*NDY, cudaMemcpyDeviceToHost);
            for (int jd = 0; jd < NDY; ++jd ) 
            {
             for (int id = 0; id < NDX; ++id ) 
             {
               int ind =id+NDX*jd;
             psi[ind]= cuCmul (cuCadd(ortkx[id], ortky[jd] ),  psi[ind] );
             }
            }
            
            cudaMemcpy( psi_d, psi, sizeof(cuDoubleComplex)*NDX*NDY, cudaMemcpyHostToDevice );

*/
//--------------------------------------------------------------------------------------------

/*           
  Do i=1,NDX
         Do j=1,NDY*NDZ
           ind=i+NDX*(J-1)
            psi(ind)=(ort_kx(i)**2+ ort_ky(j)**2)/(2*Time_mass)/NDX/NDY*
     &            psi(ind)
         EndDo
      EndDo 
 */           
                 
            cufftExecZ2Z(plan,  psi_d,  psi_d, CUFFT_INVERSE);

  //           printf("Cufft-Inverse\n"); 
  
      // copy vector back from GPU to CPU  
       
      //      cublasGetVector(NDX*NDY,sizeof(cuDoubleComplex),psi_d,1, psi,1); 
            
             cudaMemcpy(psi, psi_d, sizeof(cufftDoubleComplex) * NDX*NDY, cudaMemcpyDeviceToHost);

//                 printf("Copy from device psi[NDX*NDY-1]= %8.6f\n", psi[NDX*NDY-1]); 

                     
   // copy vectors back from GPU to CPU
 //  cudaMemcpy( a, a_d, sizeof(float) * N, cudaMemcpyDeviceToHost );
  // cudaMemcpy( b, b_d, sizeof(float) * N, cudaMemcpyDeviceToHost );

   // free GPU memory
     
           cufftDestroy(plan);
  
           cudaFree(psi_d);
       
           cudaFree(ortkx_d);
           cudaFree(ortky_d);

        // cublasDestroy(handle);          //  cublasDestroy_v2(handle);
           
        

       return;

         }
