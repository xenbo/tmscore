/*************************************************************************
	> File Name: getpadbpath.c
	> Author: ma6174
	> Mail: ma6174@163.com 
	> Created Time: 2014年03月26日 星期三 21时29分17秒
 ************************************************************************/


void getpdbpath(char *p,char p2[100][100])
{
	FILE *fp;
	if((fp=fopen(p,"r"))==NULL) 
			printf(" open file  error !!\n");
	int nl=0;		
	while(fgets(p2[nl],100,fp))
	{
		nl++;
	}
	fclose(fp);
	
	for(nl=0;nl<6;nl++)
		printf(p2[nl]);
}


