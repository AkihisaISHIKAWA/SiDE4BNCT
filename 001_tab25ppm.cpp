#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <direct.h>

//Physical constant
#define PI 3.141592653589793	//pi
#define MP 1.6605655e-27	// MP=1.6605655e-27[kg] Weight of proton
#define NAV  6.022045e23	// NAV=6.022045e23[1/mol] Avogadro number */
#define K   1.380662e-23	// K=1.380662e-23[J/K] Boltzmann's constant */
#define R 8.31441		// R=8.31441[J/K/mol] gas constant */
#define E 1.6021892e-19		//  E=1.6021892e-19[C] Elementary electric charge */
#define C 2.9979246e8		//  C = 2.9979246e8 [m/s] light velocity  */
#define H 6.626176e-34		//  H = 6.626176e-34 [Js] Planck's constant  */
#define Cbn 25.0       //Cbn=25.0[ppm] Boron concentration in normal cells
#define RBEthrm 2.9   //RBE=2.9 Relative biologic effectiveness for thermal neutron
#define RBEfast 2.4   //RBE=2.4 Relative biologic effectiveness for fast neutron
#define CBEn 1.34     //CBEn=1.35 Compound biological effectiveness (CBE) of other normal cells
#define CBEs 2.5     //CBEs=1.35 Compound biological effectiveness (CBE) of skin cells
#define CBEm 4.9     //CBEm=1.35 Compound biological effectiveness (CBE) of mucosa cells
#define CBEt 4.0      //CBEt=3.8  Compound biological effectiveness (CBE) of tumor cells
#define TNR 3.5       //Tumor/normal ratio of boron concentration

#define BYTE 256		// 1 byte

typedef struct _Seng_tab
{	
	char buf[BYTE];
} Seng_tab;

int main(int argc, char *argv[])
{
	FILE *fp, *fpo;

	Seng_tab *eng_tab;
	
	char buf[BYTE];
	char filename[BYTE], outfile[BYTE];
	char dum[BYTE];

	int flg;
	int i,j,k,l,n;
	int num_tab, num_nker, num_gker, num_spe;
	
	float val, val1, val2, val3, val4;
	float *eng_nker, *eng_gker;
	float *hydkerma, *nitkerma, *othkerma, *nkerma, *gkerma;
	float el, eu, nf, gf, err;
	float **hydd, **nitd, **othd, **nd, **gd, **Nbd, **Sbd, **Mbd, **Tbd;

// counting number of data
	sprintf(filename,"eng-list.txt");
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("%s can't be opened.\n",filename);
		return -100;
	}
	i = 0;
	fgets(dum,256,fp);
	while(fgets(dum,256,fp)!=NULL)
	{
		i++;
	}
	num_tab = i; //number of table data
	fclose(fp);
	
	sprintf(filename,"nkerma_AI.txt");
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("%s can't be opened.\n",filename);
		return -100;
	}
	i = 0;
	fgets(dum,256,fp);
	while(fgets(dum,256,fp)!=NULL)
	{
		i++;
	}
	num_nker = i; //number of neutron kerma
	fclose(fp);

	sprintf(filename,"gkerma_AI.txt");
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("%s can't be opened.\n",filename);
		return -100;
	}
	i = 0;
	fgets(dum,256,fp);
	while(fgets(dum,256,fp)!=NULL)
	{
		i++;
	}
	num_gker = i; //number of gamma-ray kerma
	fclose(fp);

//allocating data array
	hydkerma = (float *)malloc(sizeof(float)*(float)num_nker);
	nitkerma = (float *)malloc(sizeof(float)*(float)num_nker);
	othkerma = (float *)malloc(sizeof(float)*(float)num_nker);
	nkerma = (float *)malloc(sizeof(float)*(float)num_nker);
	gkerma = (float *)malloc(sizeof(float)*(float)num_gker);
	eng_tab = (Seng_tab *)malloc(sizeof(Seng_tab)*num_tab);
	eng_nker = (float *)malloc(sizeof(float)*num_nker);
	eng_gker = (float *)malloc(sizeof(float)*num_gker);

//reading data
	sprintf(filename,"eng-list.txt");
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("%s can't be opened.\n",filename);
		return -100;
	}
	i=0;
	fgets(dum,256,fp);
	while(fgets(dum,256,fp)!=NULL)
	{
		if(i>num_tab)
		{
			printf("Reading process reached EOF.\n");
			return -999;
		}
		sscanf(dum,"%s", eng_tab[i].buf);
		i++;
	}
	fclose(fp);

	sprintf(filename,"nkerma_AI.txt");
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("%s can't be opened.\n",filename);
		return -100;
	}
	i=0;
	fgets(dum,256,fp);
	while(fgets(dum,256,fp)!=NULL)
	{
		if(i>num_nker)
		{
			return -999;
		}

		sscanf(dum,"%f %f %f %f %f", &val, &val1, &val2, &val3, &val4);
		eng_nker[i]   = val;
		hydkerma[i]   = val1;
		nitkerma[i]   = val2;
		othkerma[i]   = val3;
		nkerma[i]     = val4;

		i++;
	}
	fclose(fp);

	sprintf(filename,"gkerma_AI.txt");
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("%s can't be opened.\n",filename);
		return -100;
	}
	i=0;
	fgets(dum,256,fp);
	while(fgets(dum,256,fp)!=NULL)
	{
		if(i>num_gker)
		{
			return -999;
		}
		
		sscanf(dum,"%f %f %f %f %f", &val, &val1, &val2, &val3, &val4);
		eng_gker[i]   = val;
		gkerma[i]     = val4;

		i++;
	}
	fclose(fp);

//loading table data
	int ne_tab=200; //number of energy bins in table data
	int num_tally=200; //number of tallies in table data
	//hydrogen dose, nitrogon dose, other neutron dose, neutron dose, boron dose
	float hydd_d, nitd_d, othd_d, nd_d, Nbd_d, Sbd_d, Mbd_d, Tbd_d;
	hydd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		hydd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) hydd[i][j]=0.0;		//hydd[input energy][depth]
	}
	nitd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		nitd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) nitd[i][j]=0.0;		//nitd[input energy][depth]
	}
	othd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		othd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) othd[i][j]=0.0;		//othd[input energy][depth]
	}
	nd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		nd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) nd[i][j]=0.0;		//nd[input energy][depth]
	}
	Nbd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		Nbd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) Nbd[i][j]=0.0;		//Nbd[input energy][depth]
	}
	Sbd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		Sbd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) Sbd[i][j]=0.0;		//Sbd[input energy][depth]
	}
	Mbd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		Mbd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) Mbd[i][j]=0.0;		//Mbd[input energy][depth]
	}
	Tbd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		Tbd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) Tbd[i][j]=0.0;		//Tbd[input energy][depth]
	}

	char path2tab[]="../results";
	for (i=0;i<num_tab;i++){	//input energy loop
		l=0;
		sprintf(filename,"%s/%sMeV-nspe.dat",path2tab,eng_tab[i].buf);
		if((fp=fopen(filename,"r"))==NULL)
		{
			printf("%s can't be opened.\n",filename);
			return -100;
		}
		printf("%s is successfully opened.\n",filename);
		flg=0;
		while(l<num_tally)
		{
			while(flg==0)
			{
				sprintf(buf,"\n");
				fgets(dum,256,fp);
				sscanf(dum,"%s %s",&buf,&buf);
				if(strncmp(buf,"e-lower",7)==0)
				{
					flg=1;
				}
			}

			int kl, ku;
			float ave;
			float hydkerma_e, nitkerma_e, othkerma_e, nkerma_e, bkerma_e;
			hydd_d=0.0;
			nitd_d=0.0;
			othd_d=0.0;
			nd_d=0.0;
			Nbd_d=0.0;
			Sbd_d=0.0;
			Mbd_d=0.0;
			Tbd_d=0.0;
			for(j=0;j<ne_tab;j++)
			{
				
				float nd_de;
				fgets(dum,256,fp);
				sscanf(dum,"%f %f %f %f",&el,&eu,&nf,&err);

				kl = 0;
				ku = 1;
				ave = exp( (log(el)+log(eu))*0.5 );
				for(k=0;k<num_nker;k++)
				{
					if(ave>eng_nker[k])
					{
						kl = k;
						ku = k+1;
					}
				}
				if(ku>=num_nker)
				{
					ku = num_nker-1;
					kl = num_nker-2;
				}
				hydkerma_e = exp( log(hydkerma[kl])+(log(hydkerma[ku])-log(hydkerma[kl]))/(log(eng_nker[ku])-log(eng_nker[kl]))*(log(ave)-log(eng_nker[kl])) );
				nitkerma_e = exp( log(nitkerma[kl])+(log(nitkerma[ku])-log(nitkerma[kl]))/(log(eng_nker[ku])-log(eng_nker[kl]))*(log(ave)-log(eng_nker[kl])) );
				othkerma_e = exp( log(othkerma[kl])+(log(othkerma[ku])-log(othkerma[kl]))/(log(eng_nker[ku])-log(eng_nker[kl]))*(log(ave)-log(eng_nker[kl])) );
				nkerma_e = exp( log(nkerma[kl])+(log(nkerma[ku])-log(nkerma[kl]))/(log(eng_nker[ku])-log(eng_nker[kl]))*(log(ave)-log(eng_nker[kl])) );
				bkerma_e = pow(10.0, -0.496*(6.0+log10(ave))-1.839);	//Boron kerma coeff. [pGy cm**2/ugB-10]
				hydd_d+=RBEfast*hydkerma_e*nf;                         //Hydrogen dose (1H(n,n')p)
				nitd_d+=RBEthrm*nitkerma_e*nf;                         //Nitrogen dose (14N(n,p)14C)
				othd_d+=RBEfast*othkerma_e*nf;                         //Other neutron RBE dose (12C(n,alpha)9Be and 16O(n,alpha)13C)
				nd_d+=(RBEfast*hydkerma_e+RBEthrm*nitkerma_e+RBEfast*othkerma_e)*nf;   //Neutron RBE dose 
				Nbd_d+=Cbn*CBEn*bkerma_e*nf;	               //Boron-dose in other normal cells
				Sbd_d+=Cbn*CBEs*bkerma_e*nf;	               //Boron-dose in skin cells
				Mbd_d+=Cbn*CBEm*bkerma_e*nf;	               //Boron-dose in mucosa cells
				Tbd_d+=Cbn*TNR*CBEt*bkerma_e*nf;               //Boron-dose in tumor cells
				flg=0;
			}
			hydd[i][l]=hydd_d;	//i:mono-eng, l:depth
			nitd[i][l]=nitd_d;
			othd[i][l]=othd_d;
			nd[i][l]=nd_d;
			Nbd[i][l]=Nbd_d;
			Sbd[i][l]=Sbd_d;
			Mbd[i][l]=Mbd_d;
			Tbd[i][l]=Tbd_d;
			l++;
		}
		fclose(fp);
	}
	//secondary gamma-ray dose
	float gd_d;
	gd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		gd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) gd[i][j]=0.0;		//gd[input energy][depth]
	}
	for (i=0;i<num_tab;i++){	//input energy loop
		l=0;
		sprintf(filename,"%s/%sMeV-gspe.dat",path2tab,eng_tab[i].buf);
		if((fp=fopen(filename,"r"))==NULL)
		{
			printf("%s can't be opened.\n",filename);
			return -100;
		}
		printf("%s is successfully opened.\n",filename);
		flg=0;
		while(l<num_tally)
		{
			while(flg==0)
			{
				sprintf(buf,"\n");
				fgets(dum,256,fp);
				sscanf(dum,"%s %s",&buf,&buf);
				if(strncmp(buf,"e-lower",7)==0)
				{
					flg=1;
				}
			}

			int kl, ku;
			float ave;
			float gkerma_e;
			gd_d=0.0;
			for(j=0;j<ne_tab;j++)
			{
				
				float gd_de;
				fgets(dum,256,fp);
				sscanf(dum,"%f %f %f %f",&el,&eu,&gf,&err);

				kl = 0;
				ku = 1;
				ave = exp( (log(el)+log(eu))*0.5 );
				for(k=0;k<num_gker;k++)
				{
					if(ave>eng_gker[k])
					{
						kl = k;
						ku = k+1;
					}
				}
				if(ku>=num_gker)
				{
					ku = num_gker-1;
					kl = num_gker-2;
				}
				gkerma_e = exp( log(gkerma[kl])+(log(gkerma[ku])-log(gkerma[kl]))/(log(eng_gker[ku])-log(eng_gker[kl]))*(log(ave)-log(eng_gker[kl])) );
				gd_d+=gkerma_e*gf;     //Gamma-ray dose by neutron (1H(n,gamma)2H)
				flg=0;
			}
			gd[i][l]=gd_d;	//i:mono-eng, l:depth
			l++;
		}
		fclose(fp);
	}
//This is the end of common discription
//write file
	char path2ndose[]="./ndose25ppm";
	for(i=0;i<num_tab;i++){
		sprintf(filename,"%s/%03d-ndose.dat",path2ndose,i);
		fp=fopen(filename,"w");
		fprintf(fp,"#Eng=%sMeV\n#column1:depth[mm]\n#column2:hydrogen-dose[pGy-eq/source]\n#column3:nitrogen-dose[pGy-eq/source]\n#column4:other-N-dose-eq[pGy-eq/source]\n#column5:N-dose[pGy-eq/source]\n#column6:G-dose[pGy/source]\n#column7:B-dose in other normal cells[pGy-eq/source]\n#column8:B-dose in skin cells[pGy-eq/source]\n#column9:B-dose in mucosa cells[pGy-eq/source]\n#column10:B-dose in tumor cells[pGy-eq/source]\n",eng_tab[i].buf);
		fprintf(fp,"#constants:Cb(N)=%.1fppm, RBE for thermal neutron=%.1f, RBE for fast neutron=%.1f, CBE(N)=%.2f, CBE(S)=%.1f, CBE(M)=%.1f, CBE(T)=%.1f, T/N=%.1f\n", Cbn, RBEthrm, RBEfast, CBEn, CBEs, CBEm, CBEt, TNR);
		for(j=0;j<num_tally;j++){
			fprintf(fp,"%i\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\n", j, hydd[i][j], nitd[i][j], othd[i][j], nd[i][j], gd[i][j], Nbd[i][j], Sbd[i][j], Mbd[i][j], Tbd[i][j]);
		}
		fclose(fp);
	}


// release memory
	free(hydkerma);
	free(nitkerma);
	free(othkerma);
	free(nkerma);
	free(gkerma);
	free(eng_tab);
	free(eng_nker);
	free(eng_gker);
	free(hydd);
	for(i=0;i<num_tab;i++){
		free(hydd[i]);
	}
	free(nitd);
	for(i=0;i<num_tab;i++){
		free(nitd[i]);
	}
	free(othd);
	for(i=0;i<num_tab;i++){
		free(othd[i]);
	}
	free(nd);
	for(i=0;i<num_tab;i++){
		free(nd[i]);
	}
	free(Nbd);
	for(i=0;i<num_tab;i++){
		free(Nbd[i]);
	}
	free(Sbd);
	for(i=0;i<num_tab;i++){
		free(Sbd[i]);
	}
	free(Mbd);
	for(i=0;i<num_tab;i++){
		free(Mbd[i]);
	}
	free(Tbd);
	for(i=0;i<num_tab;i++){
		free(Tbd[i]);
	}
	free(gd);
	for(i=0;i<num_tab;i++){
		free(gd[i]);
	}
	
//	fclose(fpo);
	return 1;
}
