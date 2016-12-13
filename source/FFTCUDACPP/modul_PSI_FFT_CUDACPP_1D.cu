#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cufft.h>


#include <iostream>

#include <complex>


#include <cuda_runtime.h>
#include <cublas_v2.h>   

 //#include <cublas.h> 
//------------------------------------------------------------------

#define imin(a,b) (a<b?a:b)         


 
 __global__ void kinetic( cuDoubleComplex *vec, cuDoubleComplex *a, int Nx, int Ny, int Nz, int mcase)
{ 
       
      if (mcase==2)
        {
          int idx = threadIdx.x  + blockIdx.x* blockDim.x;

             while(idx<Nx)
             {
              for (int idy = 0; idy < Ny; ++idy ) 
              {
                for (int idz = 0; idz < Nz; ++idz ) 
                {

//            ind=k+(J-1)*NDX+(I-1)*NDX*NDY
//            psi(ind)=psi(ind)*ortky(J)**2/(2*Time_mass)/NDY
               int ind =idx + idy*Nx+idz*Nx*Ny;
              
               vec[ind]=cuCmul(vec[ind],a[idy]);
                 }
              }

            idx += blockDim.x*gridDim.x;
          }
        }

      if (mcase==1)
        {
          int idy = threadIdx.x  + blockIdx.x* blockDim.x;

             while(idy<Ny)
             {
              for (int idz = 0; idz < Nz; ++idz ) 
              {
                for (int idx = 0; idx < Nx; ++idx ) 
                {

               int ind =idx + idy*Nx+idz*Nx*Ny;
              
               vec[ind]=cuCmul(vec[ind],a[idx]);
                 }
              }
//    __syncthreads(); 
            idy += blockDim.x*gridDim.x;
         }
      }
    if (mcase==3)
        {
          int idz = threadIdx.x  + blockIdx.x* blockDim.x;

             while(idz<Nz)
             {
              for (int idx = 0; idx < Nx; ++idx ) 
              {
                for (int idy = 0; idy < Ny; ++idy ) 
                {

               int ind =idx + idy*Nx+idz*Nx*Ny;
              
               vec[ind]=cuCmul(vec[ind],a[idz]);
                 }
              }
//   __syncthreads(); 
            idz += blockDim.x*gridDim.x;
            }
        }
}       

   
  

/*--------- function called from main fortran programn ---------------*/

extern "C" void kernel_fft_1d_(cuDoubleComplex *psi,cuDoubleComplex *ortk,int * Nx,int * Ny,int * Nz, int *Xcase)
{
           cuDoubleComplex  *psi_d, *ortk_d;  // declare GPU vector copies cufftDoubleComplex  *psi_d, *ortkx_d;

     
       int NX = *Nx;
       int NY = *Ny;        
       int NZ = *Nz;
       int NDIM= NX*NY*NZ;
       int ncase= *Xcase;
      
   //    cublasStatus stat ;

//---//----  CUDA variables: 

            const int threads = imin(256,NDIM);   //dim3 (16,1,1);
            const int blocks  = imin(256, NDIM/threads ); 
            cufftHandle plan;
          
  //        cublasInit ( ) ;
 
       // Allocate memory on GPU
       
          cudaMalloc( (void **)&psi_d, sizeof(cuDoubleComplex) * NDIM );

 
  //      cudaMemcpy( psi_d, psi, sizeof(cuDoubleComplex) *NDIM, cudaMemcpyHostToDevice );
        
         //  cublasAlloc (NDIM , sizeof(cuDoubleComplex),(void **)&psi_d ) ;


  /*
         if ( stat != CUBLAS_STATUS_SUCCESS ) 
          {
         printf ( " device memory allocationfailed " ) ;
         cublasShutdown ( ) ;
         return EXIT_FAILURE;
          }
*/

//=====================================================================================          
    if (ncase == 1)
     {
          cudaMalloc( (void **)&ortk_d, sizeof(cuDoubleComplex) * NX );
          /* Initialize CUBLAS */
             cublasHandle_t handle;
             cublasCreate_v2(&handle);    // cublasCreate_v2(&handle); 
             cublasSetVector(NDIM,sizeof(cuDoubleComplex),psi,1, psi_d,1);
             cublasSetVector(NX,sizeof(cuDoubleComplex),ortk,1, ortk_d,1);    
         // cudaMemcpy( ortk_d, ortk, sizeof(cuDoubleComplex) * NX, cudaMemcpyHostToDevice );
         /* Create a 1D FFT plan. */ 
//!cufftPlanMany(cufftHandle *plan, int rank, int *n, int *inembed, int istride, int idist, int *onembed, int ostride, int odist, cufftType type, int batch);
                  int num[1], inembed[1], onembed[1];
                  num[0]=NX;   
                  inembed[0]= 0;  
                  onembed[0]= 0;
              cufftPlanMany(&plan,1, num,inembed,1,NX, onembed, 1,NX,CUFFT_Z2Z,NY*NZ);
//               cufftPlanMany(&plan,1, num ,NULL,1,0,NULL,1,0,CUFFT_Z2Z,NY*NZ);
//         printf ( " 1D FFT CUDA CPP in X " ) ;
              cufftExecZ2Z(plan,  psi_d,  psi_d, CUFFT_FORWARD);
            kinetic<<<blocks,threads>>>(psi_d, ortk_d, NX,NY,NZ, ncase);
               cufftExecZ2Z(plan,  psi_d,  psi_d, CUFFT_INVERSE);
        /*  copy vectors from GPU to CPU   */
          cublasGetVector(NDIM,sizeof(cuDoubleComplex),psi_d,1, psi,1); 
     //         cudaMemcpy(psi, psi_d, sizeof(cuDoubleComplex) * NX, cudaMemcpyDeviceToHost);
           cublasDestroy(handle);          //  cublasDestroy_v2(handle);
           cufftDestroy(plan);
           cudaFree(psi_d);
           cudaFree(ortk_d);
      }
//=====================================================================================          
    if (ncase == 2)
     {
          cudaMalloc( (void **)&ortk_d, sizeof(cuDoubleComplex) * NY );
          /* Initialize CUBLAS */
             cublasHandle_t handle;
             cublasCreate_v2(&handle);    // cublasCreate_v2(&handle); 
             cublasSetVector(NDIM,sizeof(cuDoubleComplex),psi,1, psi_d,1);
             cublasSetVector(NY,sizeof(cuDoubleComplex),ortk,1, ortk_d,1);    
         // cudaMemcpy( ortk_d, ortk, sizeof(cuDoubleComplex) * NX, cudaMemcpyHostToDevice );
         /* Create a 1D FFT plan. */ 


       cufftResult res;





//!cufftPlanMany(cufftHandle *plan, int rank, int *n, int *inembed, int istride, int idist, int *onembed, int ostride, int odist, cufftType type, int batch);
                 // int num[1]={NY}, inembed[1]={0}, onembed[1]={0};
                     int num[1]={NY}, inembed[1]={0}, onembed[1]={0};

/////             Everything 2D WORKS in 1D along Y direction!!!!
////////            NEpravil'noe rabotaet v 3D po Y direction pochemu??? - hren ego znaet!!!!
//X works              cufftPlanMany(&plan,1, num,inembed,1    ,NX, onembed, 1    ,NX,CUFFT_Z2Z,NY*NZ);
               // res =  cufftPlanMany(&plan,1, num,inembed,NX   ,1 , onembed, NX   ,1 ,CUFFT_Z2Z,NX*NZ);
               
                  res =  cufftPlanMany(&plan,1, num,inembed,NX,1 , onembed, NX, 1 ,CUFFT_Z2Z, NX);  // preobrazovanie 1 sloja

//Z works              cufftPlanMany(&plan,1, num,inembed,NX*NY,1 , onembed, NX*NY,1 ,CUFFT_Z2Z,NX*NY);
        //    cufftPlanMany(&plan,1, num,inembed,NX,1, onembed, NX,1,CUFFT_Z2Z,NX*NZ);
//           printf ("1D FFT CUDA CPP in Y res %d \n",res);  
            // cufftPlan1d(&plan, NX, CUFFT_Z2Z, BATCH);
        /* Use the CUFFT plan to transform the signal in place. */

                 for (int idz = 0; idz < NZ; ++idz ) 
                {
     
               cufftExecZ2Z(plan,  &psi_d[idz*NX*NY],  &psi_d[idz*NX*NY], CUFFT_FORWARD);
                }

            kinetic<<<blocks,threads>>>(psi_d, ortk_d,  NX,NY,NZ, ncase);
              

                        for (int idz = 0; idz < NZ; ++idz ) 
                {

                 // cufftExecZ2Z(plan,  psi_d,  psi_d, CUFFT_INVERSE);

             cufftExecZ2Z(plan,  &psi_d[idz*NX*NY],  &psi_d[idz*NX*NY], CUFFT_INVERSE);
                }   

 /*          Do J=1,NDY
            Do I=1,NDZ
            Do K=1,NDX
            ind=k+(J-1)*NDX+(I-1)*NDX*NDY
            psi(ind)=psi(ind)*ortky(J)**2/(2*Time_mass)/NDY
            EndDo
            EndDo
            EndDo 
*/
        /*  copy vectors from GPU to CPU   */
          cublasGetVector(NDIM,sizeof(cuDoubleComplex),psi_d,1, psi,1); 
     //         cudaMemcpy(psi, psi_d, sizeof(cuDoubleComplex) * NX, cudaMemcpyDeviceToHost);
           cublasDestroy(handle);          //  cublasDestroy_v2(handle);
           cufftDestroy(plan);
           cudaFree(psi_d);
           cudaFree(ortk_d);

      }
//======================================================================================== 
   if (ncase == 3)
     {
         cudaMalloc( (void **)&ortk_d, sizeof(cuDoubleComplex) * NZ );
            

   // copy vectors from CPU to GPU
      
          /* Initialize CUBLAS */
   
             cublasHandle_t handle;
             cublasCreate_v2(&handle);    // cublasCreate_v2(&handle); 
          
             cublasSetVector(NDIM,sizeof(cufftDoubleComplex),psi,1, psi_d,1);
             cublasSetVector(NZ,sizeof(cufftDoubleComplex),ortk,1, ortk_d,1);    
         
         // cudaMemcpy( ortk_d, ortk, sizeof(cuDoubleComplex) * NX, cudaMemcpyHostToDevice );

         /* Create a 1D FFT plan. */ 

// ETO NE RABOTAET!!!!               cufftPlanMany(&plan,1, Nz,NULL,NX*NY,1,NULL,NX*NY,1,CUFFT_Z2Z,NY*NX);
         //    cufftPlan1d(&plan, NX, CUFFT_Z2Z, BATCH);
                  int num[1], inembed[1], onembed[1];
                  num[0]=NZ;   
                  inembed[0]= 0;  // 0;  oksa: STAVLU NDIM --NE NADo I TAK RABOTAET
                  onembed[0]= 0;      //0;
///   Wrode rabotaet.... v 3D da FFT XZ rabotaet no esli vzyt' FFT Y to vse letot na JUX
              cufftPlanMany(&plan,1, num,inembed,NX*NY,1, onembed, NX*NY,1,CUFFT_Z2Z,NX*NY);
//         printf ( " 1D FFT CUDA CPP in Z " ) ;

       
        /* Use the CUFFT plan to transform the signal in place. */
           
               cufftExecZ2Z(plan,  psi_d,  psi_d, CUFFT_FORWARD);
        
            kinetic<<<blocks,threads>>>(psi_d, ortk_d,  NX,NY,NZ, ncase);
            
               cufftExecZ2Z(plan,  psi_d,  psi_d, CUFFT_INVERSE);
        
        /*  copy vectors from GPU to CPU   */
           
          cublasGetVector(NDIM,sizeof(cufftDoubleComplex),psi_d,1, psi,1); 
            
     //         cudaMemcpy(psi, psi_d, sizeof(cuDoubleComplex) * NX, cudaMemcpyDeviceToHost);
 
           cublasDestroy(handle);          //  cublasDestroy_v2(handle);
           cufftDestroy(plan);

           cudaFree(psi_d);
           cudaFree(ortk_d);

      }
 


       return;

}
