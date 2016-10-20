/**********************************************
* Author: Yingying Gu (ying.y.gu@gmail.com)
* version 1.0
* Copyright 2015
* University of Wisconsin-Milwaukee
***********************************************/
#include "mex.h"
#include "stdio.h"
#include "stdint.h"
#include "math.h"
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Declare variables. */
    
    double low_bkg, high_bkg;
    int num_data, num_K;
    double *p_data, *p_mu, *p_sigma, *p_pipi;
    double *p_low_bkg, *p_high_bkg;
    double *r_state;
    
    /* Get the number of elements in the input argument. */
    num_data=mxGetNumberOfElements(prhs[0]);
    num_K=mxGetNumberOfElements(prhs[1]);
    
    /* Get the data. */
    p_data = (double *)mxGetPr(prhs[0]);
    p_mu=(double *)mxGetPr(prhs[1]);
    p_sigma=(double *)mxGetPr(prhs[2]);
    p_pipi=(double *)mxGetPr(prhs[3]);
    
    /*background*/
    p_low_bkg=(double *)mxGetPr(prhs[4]);
    p_high_bkg=(double *)mxGetPr(prhs[5]);
    low_bkg=p_low_bkg[0];
    high_bkg=p_high_bkg[0];
    
    /*get the rows and cols of data*/
    int number_of_dims, c;
    int  *dim_array;
    number_of_dims = mxGetNumberOfDimensions(prhs[0]);
    dim_array = mxGetDimensions(prhs[0]);
    int rows=*dim_array++;
    int cols=*dim_array++;
    
    int i,j,h,k,nn;
    
    double pi=3.1416;
    double rlgs=0, Umf_1=0, Umf_2=0, Umf_3=0, Umf;
    double meanf[num_K];
    int k2;
    
    /* mxSetPr requires dynamic array mxCalloc/mxMalloc/mxRealloc  */
    double* p_state = mxCalloc(rows*cols*num_K, sizeof(double));
    
    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++){
            for (k=0;k<num_K;k++){
                
                rlgs=log(1/(sqrt(2*pi*p_sigma[k])))+(-0.5*(p_data[i+j*rows]-p_mu[k])*(p_data[i+j*rows]-p_mu[k])/(p_sigma[k]));
                Umf_1=rlgs;
                Umf_2=log(p_pipi[k]);
                
                Umf=Umf_1+Umf_2;
                meanf[k]=Umf;
                
            }
            /*find max in meanf*/
            double temp_max=-100;
            for (k=0;k<num_K;k++){
                if (temp_max<meanf[k])
                    temp_max=meanf[k];
            }
            
            /*update the state[i][j][k]*/
            double p_z;
            for (k=0;k<num_K;k++){
                
                p_z=0;
                for (k2=0;k2<num_K;k2++){
                    
                    p_z=p_z+exp(meanf[k2]-temp_max);
                    
                }
                
                p_state[i+j*rows+k*rows*cols]=exp(meanf[k]-temp_max)/p_z;
                
                if(p_data[i+j*rows]<low_bkg | p_data[i+j*rows]>high_bkg) p_state[i+j*rows+k*rows*cols]=0;
                
                if(mxIsNaN(p_state[i+j*rows+k*rows*cols])) p_state[i+j*rows+k*rows*cols]=0;
                
            }
        }
        
    }
    
    /*Allocate the space for the return argument */
    
    int temp=rows*cols*num_K;
    plhs[0] =mxCreateDoubleMatrix(temp, 1, mxREAL);
    r_state = mxGetPr(plhs[0]);
    
    
    for(nn=0;nn<temp;nn++){
        r_state[nn]=p_state[nn];
    }
   
    
}