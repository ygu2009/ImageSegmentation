#include "mex.h"
#include "stdio.h"
#include "stdint.h"
#include "math.h"
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Declare variables. */
    double MRF_beta;   /*MRF parameters*/
    double low_bkg, high_bkg, bkg_num;
    int num_data, num_K;
    double *p_data, *p_mu, *p_sigma, *p_pipi, *p_state, *p_flag, *p_MRFbeta;
    double *p_MAP_alpha, *p_MAP_beta_0, *p_MAP_m0, *p_MAP_S0, *p_MAP_v0;
    double *p_low_bkg, *p_high_bkg, *p_bkg_num;
    double *r_state, *r_mu, *r_sigma, *r_pipi, *r_nc;
    
    /* Get the number of elements in the input argument. */
    num_data=mxGetNumberOfElements(prhs[0]);
    num_K=mxGetNumberOfElements(prhs[1]);
    
    /* Get the data and parameters */
    p_data = (double *)mxGetPr(prhs[0]);
    p_mu=(double *)mxGetPr(prhs[1]);
    p_sigma=(double *)mxGetPr(prhs[2]);
    p_pipi=(double *)mxGetPr(prhs[3]);
    p_state=(double *)mxGetPr(prhs[4]);
    p_flag=(double *)mxGetPr(prhs[5]);
    p_MRFbeta=(double *)mxGetPr(prhs[6]);
    
    /*MAP Priors*/
    p_MAP_alpha=(double *)mxGetPr(prhs[7]);
    p_MAP_beta_0=(double *)mxGetPr(prhs[8]);
    p_MAP_m0=(double *)mxGetPr(prhs[9]);
    p_MAP_S0=(double *)mxGetPr(prhs[10]);
    p_MAP_v0=(double *)mxGetPr(prhs[11]);
    
    /*background*/
    p_low_bkg=(double *)mxGetPr(prhs[12]);
    p_high_bkg=(double *)mxGetPr(prhs[13]);
    p_bkg_num=(double *)mxGetPr(prhs[14]);
    
    
    MRF_beta=p_MRFbeta[0];
    
    low_bkg=p_low_bkg[0];
    high_bkg=p_high_bkg[0];
    bkg_num=p_bkg_num[0];
    
    /*get the rows and cols of data*/
    mwSize rows=mxGetM(prhs[0]);
    mwSize cols=mxGetN(prhs[0]);
    
    int i,j,k,nn;
    
    
    /*clique function*/
    double v2[num_K][num_K];
    for (i=0;i<num_K;i++){
        for (j=0;j<num_K;j++){
            v2[i][j]=0.5;
            if (i==j)
                v2[i][j]=-0.5;
        }
    }
    
    /*E step*/
    double pi=3.1416;
    double rlgs=0, Umf_1=0, Umf_2=0, Umf_3=0, Umf;
    double meanf[num_K];
    int k2;
    double p_z;
    
    for (i=1;i<rows-1;i++){
        for (j=1;j<cols-1;j++){
            for (k=0;k<num_K;k++){
                if (p_flag[k]<0){
                    rlgs=log(1/(sqrt(2*pi*p_sigma[k])))+(-0.5*(p_data[i+j*rows]-p_mu[k])*(p_data[i+j*rows]-p_mu[k])/(p_sigma[k]));
                    Umf_1=-(1/MRF_beta)*rlgs;
                    Umf_2=-(1/MRF_beta)*log(p_pipi[k]);
                    Umf_3=0;
                    for (k2=0;k2<num_K;k2++)
                        Umf_3=Umf_3-v2[k][k2]*p_state[i+j*rows+k2*rows*cols];
                    int s1=0,s2=0;
                    for (s1=0;s1<3;s1++){
                        for (s2=0;s2<3;s2++){
                            for (k2=0;k2<num_K;k2++){
                                Umf_3=Umf_3+v2[k][k2]*p_state[i-1+s1+(j-1+s2)*rows+k2*rows*cols];
                            }
                        }
                    }
                    Umf=Umf_1+Umf_2+Umf_3;
                    meanf[k]=MRF_beta*Umf;
                }
            }
            /*find max in meanf*/
            double temp_max=-100;
            for (k=0;k<num_K;k++){
                if (temp_max<(-1)*meanf[k]){
                    temp_max=(-1)*meanf[k]; 
                }
                
            }
            
            /*update the state[i][j][k]*/
            
            for (k=0;k<num_K;k++){
                if (p_flag[k]<0){
                    p_z=0;
                    for (k2=0;k2<num_K;k2++){
                        if (p_flag[k2]<0){
                            p_z=p_z+exp((-1)*meanf[k2]-temp_max);
                        }
                    }
                    p_state[i+j*rows+k*rows*cols]=exp((-1)*meanf[k]-temp_max)/p_z;
                    
                    if(p_data[i+j*rows]<low_bkg | p_data[i+j*rows]>high_bkg) p_state[i+j*rows+k*rows*cols]=0;
                    
                    if(mxIsNaN(p_state[i+j*rows+k*rows*cols])) p_state[i+j*rows+k*rows*cols]=0;
                }
            }
        }
    }
    
    /*M step*/
    double mu[num_K];
    double n_c[num_K];
    /*estimate mu*/
    for(k=0;k<num_K;k++) {
        if (p_flag[k]<0){
            mu[k]=0;
            n_c[k]=0;
            for (i=1;i<rows-1;i++){
                for (j=1;j<cols-1;j++){
                    mu[k]=mu[k]+p_state[i+j*rows+k*rows*cols]*p_data[i+j*rows];
                    n_c[k]=n_c[k]+p_state[i+j*rows+k*rows*cols];
                }
            }
            /*mu[k]=mu[k]/n_c[k];*/ /*regular EM*/
            mu[k]=(mu[k]+p_MAP_beta_0[k]*p_MAP_m0[k])/(n_c[k]+p_MAP_beta_0[k]); /*MAP EM*/
        }
    }
    
    /*estimate sigma*/
    double sigma[num_K];
    double D=2;
    for(k=0;k<num_K;k++) {
        if (p_flag[k]<0){
            sigma[k]=0;
            for (i=1;i<rows-1;i++){
                for (j=1;j<cols-1;j++){
                    sigma[k]=sigma[k]+p_state[i+j*rows+k*rows*cols]*(p_data[i+j*rows]-mu[k])*(p_data[i+j*rows]-mu[k]);
                }
            }
            /*sigma[k]=sigma[k]/n_c[k];*/ /*regular EM*/
            sigma[k]=(p_MAP_S0[k]+sigma[k]+(p_MAP_beta_0[k]*(mu[k]-p_MAP_m0[k])*(mu[k]-p_MAP_m0[k])))/(n_c[k]+p_MAP_beta_0[k]+p_MAP_v0[k]+D+2); /*MAP EM*/
        }
    }
    
    /*estimate pi*/
    double temp_sum_alpha=0;
    double temp_K=0;
    for(k=0;k<num_K;k++) {
        if (p_flag[k]<0){
            temp_sum_alpha=temp_sum_alpha+p_MAP_alpha[k];
            temp_K=temp_K+1;
        }
    }
    
    double pipi[num_K];
    for(k=0;k<num_K;k++) {
        if (p_flag[k]<0){
            /*pipi[k]=n_c[k]/num_data;*/ /*regular EM*/
            pipi[k]=(n_c[k]+p_MAP_alpha[k]-1)/(num_data-bkg_num+temp_sum_alpha-temp_K); /*MAP EM*/
        }
    }
    
    /*Allocate the space for the return argument */
    
    /*ouput new state*/
    int temp=rows*cols*num_K;
    plhs[0] = mxCreateDoubleMatrix(temp, 1, mxREAL);
    r_state = mxGetPr(plhs[0]);
    
    for(nn=0;nn<temp;nn++){
        r_state[nn]=p_state[nn];
    }
    
    /*output new mu*/
    plhs[1] = mxCreateDoubleMatrix(num_K, 1, mxREAL);
    r_mu = mxGetPr(plhs[1]);
    nn=0;
    for(k=0;k<num_K;k++)
        r_mu[nn++]=mu[k];
    
    /*output new mu*/
    plhs[2] = mxCreateDoubleMatrix(num_K, 1, mxREAL);
    r_sigma = mxGetPr(plhs[2]);
    nn=0;
    for(k=0;k<num_K;k++)
        r_sigma[nn++]=sigma[k];
    
    /*output new pipi*/
    plhs[3] = mxCreateDoubleMatrix(num_K, 1, mxREAL);
    r_pipi = mxGetPr(plhs[3]);
    nn=0;
    for(k=0;k<num_K;k++)
        r_pipi[nn++]=pipi[k];
    
    /*output new n_c*/
    plhs[4] = mxCreateDoubleMatrix(num_K, 1, mxREAL);
    r_nc = mxGetPr(plhs[4]);
    nn=0;
    for(k=0;k<num_K;k++)
        r_nc[nn++]=n_c[k];
    
}