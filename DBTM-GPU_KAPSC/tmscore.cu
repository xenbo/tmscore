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

    int bidx=blockIdx.x;
    //int tx=threadIdx.x+blockIdx.x*blockDim.x;
    int ty=threadIdx.y;

    if(bidx<tn)
    {
        double d2=0,dis=0;
        //int nchange=0;
        int k,k2;
        tmscore2[bidx]=0.0;

        for(k=0; k<pl; k++)
        {
			for(k2=0; k2<3; k2++)
				u[bidx][k][k2]=hm[bidx][k2] +(rr[bidx][k2][0]*p0[k]+rr[bidx][k2][1]*p0[pl+k]+rr[bidx][k2][2]*p0[2*pl+k])-\
							   p1[bidx*3*pl+k2*pl+k];
			
            dis=u[bidx][k][0]*u[bidx][k][0]+u[bidx][k][1]*u[bidx][k][1]+u[bidx][k][2]*u[bidx][k][2];
            tmscore2[bidx]+=1.0/(1.0+dis/d0/d0);
            dist[bidx*pl+k]=dis;

        }


        int ncut=0;
        while(ncut<3)
        {
            d2=d*d;
            ncut=0;
            for(k=0; k<pl; k++)
            {
                if(dist[bidx*pl+k]<d2) ncut++;
            }
            d+=0.5;
        }

        nchange[bidx]=0;

        ncut=0;
        for(k=0; k<pl; k++)
        {
            //	if(ty==0 && tx%3==0)
            //		cuPrintf(" =======%lf   %lf \n",dist[bidx*pl+k],d2);
            if(dist[bidx*pl+k]<d2)
            {
                if(ncut<nalign[bidx] && ialignt[bidx*pl+ncut]==k) ncut++;
                else
                {
                    nchange[bidx]=1;
                    ialignt[bidx*pl+ncut]=k;
                    ncut++;
                }

            }

            __syncthreads();
        }

        __syncthreads();
        //if(ty==0 && tx%3==0)
        //cuPrintf("%lf \n",tmscore[bidx]);

        if(tmscore2[bidx]/(double)pl>tmscore[bidx])
            tmscore[bidx]=tmscore2[bidx]/(double)pl;


        int n=0;
        for(k=0; k<ncut; k++)
        {
            int m=ialignt[bidx*pl+k];
            
		    rp1[bidx][n][0]=p0[0*pl+m];
            rp1[bidx][n][1]=p0[1*pl+m];
            rp1[bidx][n][2]=p0[2*pl+m];
						 
			rp2[bidx][n][0]=p1[bidx*pl*3+0*pl+m];
            rp2[bidx][n][1]=p1[bidx*pl*3+1*pl+m];
            rp2[bidx][n][2]=p1[bidx*pl*3+2*pl+m];
						
			
            n++;

        }

        __syncthreads();
        //if(ty==0&&tx%3==0)
        //	cuPrintf("dddddddd  %d \n",ncut);
        nalign[bidx]=ncut;
    }
}


__device__ double d0[35]= {0};
__device__ double d[35]= {0};
__device__ int seed[35]={0};
__device__ int i[35]={0};
__device__  int istart[35]={0};

__global__
void
tmscore_gpu(int pl,double *p0,double *p1,double *tm2)
{

    int bidx=blockIdx.x;

    int tx=threadIdx.x+blockIdx.x*blockDim.x;
    int ty=threadIdx.y;


    if(bidx<tn)
    {
            d0[bidx]=1.24*powf((pl-15),(1.0/3.0))-1.8;
            if(d0[bidx]<0.5) d0[bidx]=0.5;
            double d0_search=d0[bidx];
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
        for(seed[bidx]=0; seed[bidx]<n_init; seed[bidx]++)
        {
            int L_init=L_ini[seed[bidx]];
            //if(ty==0&&tx%3==0)
            //	cuPrintf("L_init ====================== %d \n",L_init);
            for(istart[bidx]=0; istart[bidx]<=pl-L_init; istart[bidx]++)
            {
                int n;
                nalign[bidx]=L_init;
                {
                    n=0;
                    for(i[bidx]=0; i[bidx]<pl; i[bidx]++)
                    {
                        if(i[bidx]>=istart[bidx] && i[bidx]<istart[bidx]+L_init)
                        {

                            rp1[bidx][n][0]=p0[0*pl+i[bidx]];
                            rp1[bidx][n][1]=p0[1*pl+i[bidx]];
                            rp1[bidx][n][2]=p0[2*pl+i[bidx]];
						 
							rp2[bidx][n][0]=p1[bidx*pl*3+0*pl+i[bidx]];
                            rp2[bidx][n][1]=p1[bidx*pl*3+1*pl+i[bidx]];
                            rp2[bidx][n][2]=p1[bidx*pl*3+2*pl+i[bidx]];
						
                            n++;
                        }
                    }
                }


				Kabsch(*(rp1+bidx),*(rp2+bidx),nalign[bidx],1,&rms[bidx],*(hm+bidx),*(rr+bidx));

/*
				int i,j;
				for(i=0;i<3;i++)
				{
					cuPrintf("%lf ",hm[bidx][i]);
					for(j=0;j<3;j++)
					{
						
						cuPrintf("%lf  \n",rr[bidx][i][j]);
					}
				}
*/
				calculate_tm(d0[bidx],d0_search-1,pl,p0,p1);
                
				d[bidx]=d0_search+1;


                //if(tx%3>1000)
                {
                    for(iter[bidx]=0; iter[bidx]<n_it&&nchange[bidx]; iter[bidx]++)
                    {
                
						Kabsch(*(rp1+bidx),*(rp2+bidx),nalign[bidx],1,&rms[bidx],*(hm+bidx),*(rr+bidx));
						calculate_tm(d0[bidx],d[bidx],pl,p0,p1);
					 }
                }
            }
        }
        tm2[bidx]=tmscore[bidx];
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
    dim3 block(1,1);

    cudaPrintfInit();
    tmscore_gpu<<<grid,block>>>(140,dp0,dp1,dtm);
    cudaPrintfDisplay(stdout, true);
    cudaPrintfEnd();

    cudaMemcpy(tm,dtm,350*sizeof(double),cudaMemcpyDeviceToHost);
    
	for(i=0;i<100;i++)
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
