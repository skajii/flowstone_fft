#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

// by satoshi kajii Oct/2015

//FRAMES_API 
// Use MinGW
//$ g++  -c gsl_fft_buf4.cpp -I /local/include/              
//$ g++ -static-libgcc -static-libstdc++ -shared -o gsl_fft_buf4.dll gsl_fft_buf4.o -lgsl -L /local/lib


#define GETFLOAT(p) *((float*)&p)
#define GETBOOL(p) *((bool*)&p)
#define GETINT(p) p
#define GETSTRING(p) *((char**)&p)
#define GETFLOATARRAY(p) p ? ((float*)p+1) : 0
#define GETINTARRAY(p) p ? ((int*)p+1) : 0
#define GETSTRINGARRAY(p) p ? ((char**)p+1) : 0
#define GETARRAYSIZE(p) p ? *((int*)p) : 0
#define GETFRAME(p) p ? ((float*)p+1) : 0
#define GETFRAMESIZE(p) p ? *((int*)p) : 0
#define GETBITMAPWIDTH(p) p ? *((int*)p) : 0
#define GETBITMAPHEIGHT(p) p ? *((int*)p+1) : 0
#define GETBITMAPCHANNELS(p) p ? *((int*)p+2) : 0
#define GETBITMAPDATA(p) p ? ((BYTE*)p+12) : 0
#define GETBITMAPBYTES(p) p ? *((int*)p) * *((int*)p+1) * *((int*)p+2) : 0
#define NEWINTARRAY(p,n) if(n>0) { *((int**)&p)=new int[n+1]; ((int*)p)[0]=n; }
#define NEWFLOATARRAY(p,n) if(n>0) { *((float**)&p)=new float[n+1]; ((int*)p)[0]=n; }
#define NEWSTRINGARRAY(p,n) if(n>0) { *((char***)&p)=new char*[n+1]; ((int*)p)[0]=n; }
#define DELETESTRING(p) if(p) { delete *((char**)&p); p=0; }
#define DELETEINTARRAY(p) if(p) { delete *((int**)&p); p=0; }
#define DELETEFLOATARRAY(p) if(p) { delete *((float**)&p); p=0; }
#define DELETESTRINGARRAY(p) if(p) { for( int j=0; j<*((int*)p); j++ ) { if( ((char**)p+1)[j] ) delete ((char**)p+1)[j]; } delete *((char***)&p); p=0; }

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

extern "C" __declspec(dllexport) void fft( int nParams, int* pIn, int* pOut )
{

  if( pIn && pOut && nParams >= 3 )
    {
      if( pIn[0] )
	{
	  float* pData = GETFRAME(pIn[0]);
	  int n = GETFRAMESIZE(pIn[0]);
	  
	  float* pData2 = GETFRAME(pIn[1]);
		  
	  float* pData3 = GETFRAME(pIn[2]);

	  double *data = new double[n*2];
	  //data = (double*) malloc(sizeof(double)*n*2);
	  
	  if(pData2 && pData3)
	    {
	      float curr;
	      float re,im,mag,ang;
	      
	      for(int i = 0; i < n; i++)
		{
		  REAL(data,i) = (double) pData[i];         
		  IMAG(data,i) = 0.0;    
		}
	      
	      gsl_fft_complex_radix2_forward(data, 1, n);

	      for( int i=0; i<n; i++ )
		{
		  re = (float) REAL(data, i);
		  im = (float) IMAG(data, i);
		  pData[i] = re; 
		  pData2[i] = im;     
		  pData3[i] = i ;
		}
	      //free(data);
		      delete[] data;
	    }  
	}
    }
}

extern "C" __declspec(dllexport) void ifft( int nParams, int* pIn, int* pOut )
{

  if( pIn && pOut && nParams >= 1 )
    {
      if( pIn[0] )
	{
	  float* pData = GETFRAME(pIn[0]);
	  int n = GETFRAMESIZE(pIn[0]);
	  float* pDataOut = GETFRAME(pOut[0]);

	  float* pData2 = GETFRAME(pIn[1]);
		  
	  float* pDataOut2 = GETFRAME(pOut[1]);
		  
	  //double data[2*n];
	  double *data = new double[n*2];
	  //data = (double*) malloc(sizeof(double)*n*2);

	  if(pData2)
	    {
	      float curr;
	      float re,im,mag,ang;
	      
	      for(int i = 0; i < n; i++)
		{
		  REAL(data,i) = (double) pData[i];         
		  IMAG(data,i) = (double) pData2[i];    
		}

	      gsl_fft_complex_radix2_inverse(data, 1, n);
	      
	      
	      for( int i=0; i<n; i++ )
		{
		  re = (float) REAL(data, i);
		  im = (float) IMAG(data, i);
		  pDataOut[i] = re;
		  pDataOut2[i] = im;     
		}
	      //     free(data);
		    delete[] data;
	    }
	}
    
    }
}
