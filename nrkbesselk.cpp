/* Driver for routine bessk */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define MAXSTR 80

int main(void)
{
	char txt[MAXSTR];
	int i,nval,n;
	float val,x;
	FILE *fp;

	if ((fp = fopen("fncval.dat","r")) == NULL)
		nrerror("Data file fncval.dat not found\n");
	fgets(txt,MAXSTR,fp);
	while (strncmp(txt,"Modified Bessel Function Kn",27)) {
		fgets(txt,MAXSTR,fp);
		if (feof(fp)) nrerror("Data not found in fncval.dat\n");
	}
	fscanf(fp,"%d %*s",&nval);
	printf("\n%s\n",txt);
	printf("%4s %7s %14s %19s\n","n","x","actual","bessk(n,x)");
	for (i=1;i<=nval;i++) {
		fscanf(fp,"%d %f %f",&n,&x,&val);
		printf("%4d %8.2f %18.7e %16.7e \n",n,x,val,bessk(n,x));
	}
	fclose(fp);
	return 0;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 3.9D2. */
