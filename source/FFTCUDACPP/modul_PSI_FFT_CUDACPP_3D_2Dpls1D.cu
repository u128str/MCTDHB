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
#include <cublas_v2.h>  

//#include <cublas.h> 

//-------------------------------------------------------------------------------------

// ---------   Neverno rabotaet "chistoe 3D _cufft" iz-za raznogo razmeshenija dannux
//----------:             SDELANO 2D+ 1D
//       
//--------------------------------------------------------------------------------------

//   cudaMemcpy rabotaet bustree, chem CUBLAS 
//--------------------------------------------------------------------------------------

#define imin(a,b) (a<b?a:b)     


__global__ void  kinetic_3D2pls1( cuDoubleComplex *vec2, 
                             cuDoubleComplex *ax, 
                             cuDoubleComplex *ay, 
                             cuDoubleComplex *az,
                             int NX, int NY, int NZ)
 {    
            
          int kd = threadIdx.x  + blockIdx.x* blockDim.x;

      
      while(kd<NZ)
      {
          for (int id = 0; id < NX; ++id ) 
          {
             for (int jd = 0; jd < NY; ++jd ) 
             {
               
               int ind =id+NX*jd+kd*NX*NY;
              
           cuDoubleComplex sum =cuCadd(ax[id], ay[jd]); 
                           sum= cuCadd( sum, az[kd] );

             vec2[ind]= cuCmul ( sum ,  vec2[ind] );
 
             }
          }
      __syncthreads(); 
        kd += blockDim.x*gridDim.x;

      }
  }
 
  
     void  kinetic_CPU_3D2pls1(   cuDoubleComplex *vec2, 
                             cuDoubleComplex *ax, 
                             cuDoubleComplex *ay, 
                             cuDoubleComplex *az,
                             int NX, int NY, int NZ)
{

         int ind;
         cuDoubleComplex sum;

         for (int idz = 0; idz < NZ; ++idz ) 
         {
            for (int jd = 0; jd < NY; ++jd ) 
            {
              for (int kd = 0;  kd < NX; ++kd ) 
              {
              
             ind = kd+NX*jd+ idz*NX*NY;                   

             sum = cuCadd(ax[kd], ay[jd]); 
             sum = cuCadd( sum,   az[idz] );

             vec2[ind]= cuCmul ( sum,  vec2[ind] );
              }
            }
          }

         return;
}
 

//------------------------------------------------------------------------------------//


//--------------------- function called from main fortran program -------------------//


extern "C" void   kernel_fft_3d2pls1_(cuDoubleComplex *psi,
                                 cuDoubleComplex *ortkx, 
                                 cuDoubleComplex *ortky,
                                 cuDoubleComplex *ortkz, int *NX, int *NY, int *NZ)
     {
       
  //                 printf("Cufft_3D-as 2D plus 1D\n");
//   printf ( " 3d-----------full -----------oksa" ) ;            

          cuDoubleComplex   *ortkx_d, *ortky_d, *ortkz_d;  //, *psiV_d;  // declare GPU vector copies cufftDoubleComplex  *psi_d, *ortkx_d;
    
          cuDoubleComplex   *psi_d;      

           int NDX = *NX;       
           int NDY = *NY;  
           int NDZ = *NZ; 

           int NDIM= NDX*NDY*NDZ; 
           
           cufftHandle plan2,plan; 

// ---- CUDA  variables: ----------------------------------------------------------------------- 

            const int threads = imin(64,NDIM);   //dim3 (16,1,1);
            const int blocks  = imin(32, NDIM/threads );  // ;  //imin(8, (M+threads)/threads ); 
//----------------------------------------------------------------------------------------------

          cudaMalloc( (void **)&psi_d, sizeof(cuDoubleComplex) * NDIM );
    
          cudaMalloc( (void **)&ortkx_d, sizeof(cuDoubleComplex) * NDX );
          cudaMalloc( (void **)&ortky_d, sizeof(cuDoubleComplex) * NDY );
          cudaMalloc( (void **)&ortkz_d, sizeof(cuDoubleComplex) * NDZ );
   


/* --------------- Initialize CUBLAS --------------------------------------- */
// ------------  cublas  -----------------------
// ------------  vozmojno nado allokirovat' memory  cublasAlloc

/*             cublasHandle_t handle;
             
             cublasCreate_v2(&handle);    // cublasCreate_v2(&handle); 

             cublasSetVector(NDIM,sizeof(cuDoubleComplex),psi,1, psi_d,1);          // 1
            
             cublasSetVector(NDX,sizeof(cuDoubleComplex),ortkx,1, ortkx_d,1); 
             cublasSetVector(NDY,sizeof(cuDoubleComplex),ortky,1, ortky_d,1);    
             cublasSetVector(NDZ,sizeof(cuDoubleComplex),ortkz,1, ortkz_d,1); 
*/  
      
 //-----------------------------------------------------------------------------------------------
            
              cudaMemcpy( psi_d, psi,     sizeof(cuDoubleComplex) * NDIM, cudaMemcpyHostToDevice );
              cudaMemcpy( ortkx_d, ortkx, sizeof(cuDoubleComplex) * NDX,  cudaMemcpyHostToDevice );
              cudaMemcpy( ortky_d, ortky, sizeof(cuDoubleComplex) * NDY,  cudaMemcpyHostToDevice ); 
              cudaMemcpy( ortkz_d, ortkz, sizeof(cuDoubleComplex) * NDZ,  cudaMemcpyHostToDevice ); 
       
  

/* -------------- Create a 2D + 1D FFT plans. ------------------------------------  */
               
                        
                 int num[1], inembed[1], onembed[1];
                   num[0]=NDZ;   
                   inembed[0]= 0; 
                   onembed[0]= 0;  
               
                     int inembed2[2]= {0,0}; 
                     int onembed2[2]= {0,0};      
            
               // cufftPlanMany(&plan,1, num,inembed,1, 0, onembed, 1,0,CUFFT_Z2Z,NDX*NDY); 

                  cufftPlanMany(&plan,1, num,inembed,NDX*NDY,1, onembed, NDX*NDY,1,CUFFT_Z2Z,NDX*NDY);
           
              int num2[2]={ NDY, NDX };  //    !!!!!!!!!!!!!!!! pomenjala X<-> Y mestami
               
               cufftPlanMany(&plan2,2, num2,NULL,1,0,NULL,1,0,CUFFT_Z2Z, NDZ);
             
                //   cufftPlanMany(&plan2,2, num2,inembed2,0,NDX*NDY, onembed2, 0, NDX*NDY,CUFFT_Z2Z, NDZ);

             //1       cufftPlan2d(&plan2, NDY, NDX, CUFFT_Z2Z); 

/* ------    Use the CUFFT plan to transform the signal in place.  ---------  */
           
         //1 for (int idz = 0; idz < NDZ; ++idz ) 
         //1       {
     
            //1   cufftExecZ2Z(plan2,  &psi_d[idz*NDX*NDY],  &psi_d[idz*NDX*NDY], CUFFT_FORWARD);
              //1  }


                     
                cufftExecZ2Z(plan2, psi_d, psi_d, CUFFT_FORWARD);
         
                cufftExecZ2Z(plan,   psi_d,  psi_d, CUFFT_FORWARD);

/*---------- call function on GPU -------------------------------------------  */
                 
       
             kinetic_3D2pls1<<<blocks,threads>>>( psi_d, ortkx_d, ortky_d, ortkz_d, NDX, NDY,NDZ);


//----------------------- kernel on CPU -----------------------------------------
//------------------------  PROOF -  -----------------------
   /*    
             cudaMemcpy(psi, psi_d, sizeof(cuDoubleComplex) *NDIM, cudaMemcpyDeviceToHost);
              
              kinetic_CPU_3D2pls1 ( psi, ortkx, ortky, ortkz, NDX, NDY, NDZ);

           
  //         cublasGetVector(NDIM,sizeof(cuDoubleComplex),psi_d,1, psi,1); 
     
            
  
          cudaMemcpy( psi_d, psi,     sizeof(cuDoubleComplex) * NDIM, cudaMemcpyHostToDevice );

//       cublasSetVector(NDIM,sizeof(cuDoubleComplex),psi,1, psi_d,1);
   
  */


//--------------------------------------------------------------------------------------------
   
                cufftExecZ2Z(plan2, psi_d, psi_d, CUFFT_INVERSE);
         
                cufftExecZ2Z(plan,  psi_d, psi_d, CUFFT_INVERSE);

      //1    for (int idz = 0; idz < NDZ; ++idz ) 
         //1       {
     
            //1   cufftExecZ2Z(plan2,  &psi_d[idz*NDX*NDY],  &psi_d[idz*NDX*NDY], CUFFT_INVERSE);
             //1   }



              // printf("Cufft_3D-Inverse\n");

//---------------- copy vector back from GPU to CPU  --------------------------------------------- 
       
         //       cublasGetVector(NDIM,sizeof(cuDoubleComplex),psi_d,1, psi,1);           // 2
               
              cudaMemcpy(psi, psi_d, sizeof(cuDoubleComplex) *NDIM, cudaMemcpyDeviceToHost); 
            
  
   // free GPU memory
     
           //cublasDestroy_v2(handle);          //  cublasDestroy_v2(handle);
           cufftDestroy(plan);
           cufftDestroy(plan2);
  

           cudaFree(psi_d);
       
           cudaFree(ortkx_d);
           cudaFree(ortky_d);
           cudaFree(ortkz_d);

                   
        

       return;

         }
