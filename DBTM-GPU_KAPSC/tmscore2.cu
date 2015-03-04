#include"head.h"
__device__ double rms[350];
__device__  double rr[350][3][3]= {0};
__device__  double hm[350][3] = {0};
//temporarily
__device__  double tmscore[350]= {0};
__device__  double tmscore2[350]= {0};

__device__  double rp1[350][200][3]= {0};
__device__  double rp2[350][200][3]= {0};

__device__  double dist[65000]= {0};
__device__  double u[350][200][3]= {0};

//__device__  double tdp[6000]= {0};
__device__  double ialignt[200*350]= {0};

__device__  int nalign[350]= {0};

__device__ int iter[350]= {0};

__device__ int nchange[350]= {0};

__device__ int tn=301;

__device__
//__global__


__device__
//__global__
void
calculate_tm(double d0,double d,int pl,double *p0,double *p1)
{

   
    int tx=threadIdx.x+blockIdx.x*blockDim.x;
    int ty=threadIdx.y;

    if(tx<tn)
    {
        double d2=0,dis=0;
        //int nchange=0;
        int k,k2;
        tmscore2[tx]=0.0;

        for(k=0; k<pl; k++)
        {
			for(k2=0; k2<3; k2++)
				u[tx][k][k2]=hm[tx][k2] +(rr[tx][k2][0]*p0[k]+rr[tx][k2][1]*p0[pl+k]+rr[tx][k2][2]*p0[2*pl+k])-\
							   p1[tx*3*pl+k2*pl+k];
			
            dis=u[tx][k][0]*u[tx][k][0]+u[tx][k][1]*u[tx][k][1]+u[tx][k][2]*u[tx][k][2];
            tmscore2[tx]+=1.0/(1.0+dis/d0/d0);
            dist[tx*pl+k]=dis;

        }


        int ncut=0;
        while(ncut<3)
        {
            d2=d*d;
            ncut=0;
            for(k=0; k<pl; k++)
            {
                if(dist[tx*pl+k]<d2) ncut++;
            }
            d+=0.5;
        }

        nchange[tx]=0;

        ncut=0;
        for(k=0; k<pl; k++)
        {
            //	if(ty==0 && tx%3==0)
            //		cuPrintf(" =======%lf   %lf \n",dist[tx*pl+k],d2);
            if(dist[tx*pl+k]<d2)
            {
                if(ncut<nalign[tx] && ialignt[tx*pl+ncut]==k) ncut++;
                else
                {
                    nchange[tx]=1;
                    ialignt[tx*pl+ncut]=k;
                    ncut++;
                }

            }

            __syncthreads();
        }

        __syncthreads();
        //if(ty==0 && tx%3==0)
        //cuPrintf("%lf \n",tmscore[tx]);

        if(tmscore2[tx]/(double)pl>tmscore[tx])
            tmscore[tx]=tmscore2[tx]/(double)pl;


        int n=0;
        for(k=0; k<ncut; k++)
        {
            int m=ialignt[tx*pl+k];
            
		    rp1[tx][n][0]=p0[0*pl+m];
            rp1[tx][n][1]=p0[1*pl+m];
            rp1[tx][n][2]=p0[2*pl+m];
						 
			rp2[tx][n][0]=p1[tx*pl*3+0*pl+m];
            rp2[tx][n][1]=p1[tx*pl*3+1*pl+m];
            rp2[tx][n][2]=p1[tx*pl*3+2*pl+m];
						
			
            n++;

        }

        __syncthreads();
        //if(ty==0&&tx%3==0)
        //	cuPrintf("dddddddd  %d \n",ncut);
        nalign[tx]=ncut;
    }
}


__device__ double d0[350]= {0};
__device__ double d[350]= {0};
__device__ int seed[350]={0};
__device__ int i[350]={0};
__device__  int istart[350]={0};

__global__
void
tmscore_gpu(int pl,double *p0,double *p1,double *tm2)
{

    

    int tx=threadIdx.x+blockIdx.x*blockDim.x;
    int ty=threadIdx.y;


    if(tx<tn)
    {
            d0[tx]=1.24*powf((pl-15),(1.0/3.0))-1.8;
            if(d0[tx]<0.5) d0[tx]=0.5;
            double d0_search=d0[tx];
            if(d0_search>8)d0_search=8;
            if(d0_search<4.5) d0_search=4.5;

            int n_it=20;
            //int n_init_max=6;
            int n_init=0;
            int L_init_min=4;
            int L_ini[6];

            if(pl<4) L_init_min=pl;
            int len=pl;
            int divisor=1;

            while(len>L_init_min && n_init<5)
            {
                L_ini[n_init++]=len;
                divisor*=2;
                len=pl/divisor;
            }
            L_ini[n_init++]=4;
            if(L_ini[n_init-1]>L_init_min) L_ini[n_init++]=L_init_min;

		threadfence();
        for(seed[tx]=0; seed[tx]<n_init; seed[tx]++)
        {
            int L_init=L_ini[seed[tx]];
            //if(ty==0&&tx%3==0)
            //	cuPrintf("L_init ====================== %d \n",L_init);
            for(istart[tx]=0; istart[tx]<=pl-L_init; istart[tx]++)
            {
                int n;
                nalign[tx]=L_init;
                {
                    n=0;
                    for(i[tx]=0; i[tx]<pl; i[tx]++)
                    {
                        if(i[tx]>=istart[tx] && i[tx]<istart[tx]+L_init)
                        {

                            rp1[tx][n][0]=p0[0*pl+i[tx]];
                            rp1[tx][n][1]=p0[1*pl+i[tx]];
                            rp1[tx][n][2]=p0[2*pl+i[tx]];
						 
							rp2[tx][n][0]=p1[tx*pl*3+0*pl+i[tx]];
                            rp2[tx][n][1]=p1[tx*pl*3+1*pl+i[tx]];
                            rp2[tx][n][2]=p1[tx*pl*3+2*pl+i[tx]];
						
                            n++;
                        }
                    }
					__syncthreads();
                }


				Kabsch(*(rp1+tx),*(rp2+tx),nalign[tx],1,&rms[tx],*(hm+tx),*(rr+tx));


				calculate_tm(d0[tx],d0_search-1,pl,p0,p1);
                
				d[tx]=d0_search+1;


                //if(tx%3>1000)
                {
                    for(iter[tx]=0; iter[tx]<n_it&&nchange[tx]; iter[tx]++)
                    {
                
						Kabsch(*(rp1+tx),*(rp2+tx),nalign[tx],1,&rms[tx],*(hm+tx),*(rr+tx));
						calculate_tm(d0[tx],d[tx],pl,p0,p1);
					 }
                }
            }
        }
        tm2[tx]=tmscore[tx];
        __syncthreads();

    }
}

int main(void)
{
    int length;
	int i;
    int size=350*200*3;
    double  *p=(double  *)malloc(size*sizeof(double ));

    double *tm=(double  *)malloc(350*sizeof(double ));

    getnpdb("l1",p,&length);

    double  *dp1; //pdb info
    cudaMalloc((void**)&dp1,size*sizeof(double));
    cudaMemcpy(dp1,p,size*sizeof(double),cudaMemcpyHostToDevice);

    double  *dp0; //pdb info
    cudaMalloc((void**)&dp0,3*200*sizeof(double));
    cudaMemcpy(dp0,p,3*200*sizeof(double),cudaMemcpyHostToDevice);


    double  *dtm;
    cudaMalloc((void**)&dtm,350*sizeof(double));
    cudaMemset(dtm,0,350*sizeof(double));


    dim3 grid(100,1);
    dim3 block(3,1);

    cudaPrintfInit();
    tmscore_gpu<<<grid,block>>>(140,dp0,dp1,dtm);
    cudaPrintfDisplay(stdout, true);
    cudaPrintfEnd();

    cudaMemcpy(tm,dtm,350*sizeof(double),cudaMemcpyDeviceToHost);
    
	for(i=0;i<300;i++)
	{
		printf("%lf \n",tm[i]);
	}

    cudaFree(dp0);
    cudaFree(dp1);
    cudaFree(dtm);
    free(tm);
    free(p);

    return 0;
}
