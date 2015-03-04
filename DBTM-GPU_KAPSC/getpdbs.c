void getnpdb(char p[300] ,double  pdblist[],int *l)
{
	 int nline2=0;
	FILE *fp=NULL;
	FILE *fppdb=NULL;
	if((fp=fopen(p,"r"))==NULL)
	{
			printf("open list error \n");
	}

	static	char line[300][100];
	static	char pdbline[60];
	static	int nline=0,i=0,j=0;
	static	char tmpCA[5];
	static	char ac_name[5];
	
	static  char x[10];
	static	char y[10];
	static  char z[10];

	double pdblist2[300];int nline22=0;
	double pdblist3[300];int nline23=0;

	while(fgets(line[nline++],100,fp))
	{
		sscanf(line[nline-1],"%s",line[nline-1]);
		if((fppdb=fopen(line[nline-1],"r"))==NULL)
		{
			printf("open NO. %d pdb error %s \n",nline,line[nline-1]);
		}
		
		while(fgets(pdbline,60,fppdb))
		{
			for(i=12;i<=15;i++)
			tmpCA[j++]=pdbline[i];
			j=0;		
			sscanf(tmpCA,"%s",ac_name);
		
			if(strcmp(ac_name,"CA")==0)
			{

				for(i=30;i<=39;i++)
					x[j++]=pdbline[i];
				j=0;
				sscanf(x,"%lf",&(pdblist[nline2++]));
				
				for(i=38;i<=47;i++)
					y[j++]=pdbline[i];
				j=0;
				sscanf(y,"%lf",&(pdblist2[nline22++]));

				for(i=46;i<=54;i++)
					z[j++]=pdbline[i];
				j=0;
				sscanf(z,"%lf",&(pdblist3[nline23++]));
			}
		}
		fclose(fppdb);

		for(j=0;j<nline22;j++)
		{
			pdblist[nline2++]=pdblist2[j];
		}
		for(j=0;j<nline23;j++)
		{
			pdblist[nline2++]=pdblist3[j];
		}
		nline22=0;nline23=0;
	}
	fclose(fp);
	*l=nline2;		
}
