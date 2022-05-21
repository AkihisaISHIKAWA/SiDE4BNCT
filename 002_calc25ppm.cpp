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
	float ene;
} Seng_tab;

int main(int argc, char *argv[])
{
	FILE *fp, *outfp;

	Seng_tab *eng_tab;
	
	char buf[BYTE];
	char filename[BYTE], outfile[BYTE];
	char dum[BYTE];

	int flg;
	int i,j,k,l,n,id;
	int num_tab, num_spe;
	
	float val, val1, val2, val3, val4, val5, val6, val7, val8;
	float el, eu, nf, err;
	float spe_el, spe_eu, spe_nf, spe_err;
	float hydd_e, nitd_e, othd_e, nd_e, gd_e, Nbd_e, Sbd_e, Mbd_e, Tbd_e;
	float **hydd, **nitd, **othd, **nd, **gd, **Nbd, **Sbd, **Mbd, **Tbd;
	float *hydGy, *nitGy, *othGy, *nGy, *gGy, *NbGy, *SbGy, *MbGy, *TbGy, *NGy, *SGy, *MGy, *TGy;

//	flg=0;
//	if(argc>1)
//	{
//		sscanf(argv[1],"%s",&flnm0);
//	}

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
	num_tab = i;
	fclose(fp);
	
//allocating data array
	eng_tab = (Seng_tab *)malloc(sizeof(Seng_tab)*num_tab);

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
		sscanf(dum,"%f", &eng_tab[i].ene);
		i++;
	}
	fclose(fp);


//loading table data
	int num_tally=200;
	float hydd_d;
	hydd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		hydd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) hydd[i][j]=0.0;		//hydd[input energy][depth]
	}
	float nitd_d;
	nitd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		nitd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) nitd[i][j]=0.0;		//nitd[input energy][depth]
	}
	float othd_d;
	othd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		othd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) othd[i][j]=0.0;		//othd[input energy][depth]
	}
	float nd_d;
	nd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		nd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) nd[i][j]=0.0;		//nd[input energy][depth]
	}
	float Nbd_d;
	Nbd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		Nbd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) Nbd[i][j]=0.0;		//Nbd[input energy][depth]
	}
	float Sbd_d;
	Sbd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		Sbd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) Sbd[i][j]=0.0;		//Sbd[input energy][depth]
	}
	float Mbd_d;
	Mbd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		Mbd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) Mbd[i][j]=0.0;		//Mbd[input energy][depth]
	}
	float Tbd_d;
	Tbd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		Tbd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) Tbd[i][j]=0.0;		//Tbd[input energy][depth]
	}
	gd = (float **)malloc(sizeof(float)*num_tab);
	for(i=0;i<num_tab;i++){
		gd[i]=(float *)malloc(sizeof(float)*num_tally);
		for(j=0;j<num_tally;j++) gd[i][j]=0.0;		//gd[input energy][depth]
	}

//read file
	char path2ndose[]="./ndose25ppm";
	for(i=0;i<num_tab;i++){
		sprintf(filename,"%s/%03d-ndose.dat",path2ndose,i);
		fp=fopen(filename,"r");
		fgets(dum,256,fp);	//skip header part #Eng=XXXMeV
		fgets(dum,256,fp);	//skip header part #column1:depth[mm]
		fgets(dum,256,fp);	//skip header part #column2:hydrogen-dose[pGy-eq/source]
		fgets(dum,256,fp);	//skip header part #column3:nitrogen-dose[pGy-eq/source]
		fgets(dum,256,fp);	//skip header part #column4:other-N-dose-eq[pGy-eq/source]
		fgets(dum,256,fp);	//skip header part #column5:N-dose[pGy-eq/source]
		fgets(dum,256,fp);	//skip header part #column6:G-dose[pGy/source]
		fgets(dum,256,fp);	//skip header part #column7:B-dose in other normal cells[pGy-eq/source]
		fgets(dum,256,fp);	//skip header part #column8:B-dose in skin cells[pGy-eq/source]
		fgets(dum,256,fp);	//skip header part #column9:B-dose in mucosa cells[pGy-eq/source]
		fgets(dum,256,fp);	//skip header part #column10:B-dose in tumor cells[pGy-eq/source]
		fgets(dum,256,fp);	//skip header part #constants:Cb(N)=25.0ppm, RBE for thermal neutron=2.9, RBE for fast neutron=2.4, CBE(N)=1.34, CBE(S)=2.50, CBE(M)=4.90, CBE(T)=4.0, T/N=3.5
		for(j=0;j<num_tally;j++)
		{
			fgets(dum,256,fp);
			sscanf(dum, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", &id, &val, &val1, &val2, &val3, &val4, &val5, &val6, &val7, &val8);
			hydd[i][j]  = val;
			nitd[i][j]  = val1;
			othd[i][j]  = val2;
			nd[i][j]    = val3;
			gd[i][j]    = val4;
			Nbd[i][j]   = val5;
			Sbd[i][j]   = val6;
			Mbd[i][j]   = val7;
			Tbd[i][j]   = val8;
		}
		fclose(fp);
	}
printf("Finished reading process!\n");

//allocation dose-depth memory
	hydGy=(float *)malloc(sizeof(float)*num_tally);
	nitGy=(float *)malloc(sizeof(float)*num_tally);
	othGy=(float *)malloc(sizeof(float)*num_tally);
	nGy=(float *)malloc(sizeof(float)*num_tally);
	gGy=(float *)malloc(sizeof(float)*num_tally);
	NbGy=(float *)malloc(sizeof(float)*num_tally);
	SbGy=(float *)malloc(sizeof(float)*num_tally);
	MbGy=(float *)malloc(sizeof(float)*num_tally);
	TbGy=(float *)malloc(sizeof(float)*num_tally);
	NGy=(float *)malloc(sizeof(float)*num_tally);
	SGy=(float *)malloc(sizeof(float)*num_tally);
	MGy=(float *)malloc(sizeof(float)*num_tally);
	TGy=(float *)malloc(sizeof(float)*num_tally);
	for(i=0;i<num_tally;i++){   //i:depth
		hydGy[i]=0.0;
		nitGy[i]=0.0;
		othGy[i]=0.0;
		nGy[i]=0.0;
		gGy[i]=0.0;
		NbGy[i]=0.0;
		SbGy[i]=0.0;
		MbGy[i]=0.0;
		TbGy[i]=0.0;
		NGy[i]=0.0;
		SGy[i]=0.0;
		MGy[i]=0.0;
		TGy[i]=0.0;
	}

//load input energy spectrum
	float ave_eng;
	int jl, ju;
//	char path2spe[]="./spectra";
	char path2spe[]="C:/Lab/001_Study/003_Prog/phantom-dose-new-ver4/001_main/spectra";
	sprintf(filename,"%s/NUANS-LiF(nat)-nspe-fc-2.2e5_tc_9e11.dat",path2spe);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("%s can't be opened.\n",filename);
		return -100;
	}
	i = 0;
	flg = 0;
	fgets(dum,256,fp);
	while(fgets(dum,256,fp)!=NULL && flg == 0)
	{
		if( strlen(dum)>5)
		{
			sscanf(dum,"%f %f %f %f",&spe_el, &spe_eu, &spe_nf, &spe_err);
			ave_eng = exp((log(spe_el)+log(spe_eu))*0.5);
			
			//search energy
			jl = 0;		//lower side energy bin for interpolating
			ju = 0;		//upper side energy bin for interpolating
			for(j=0;j<num_tab;j++)
			{
				if(eng_tab[j].ene < ave_eng)
				{
					jl = j;
					ju = j+1;
				}
			}
			if(ju+1 >= num_tab)
			{
				ju = num_tab - 1;
				jl = num_tab - 2;
			}
			
//interpolating
			for(k=0;k<num_tally;k++)	//loop for depth
			{
				hydd_e = exp( log(hydd[jl][k])+(log(hydd[ju][k])-log(hydd[jl][k]))/(log(eng_tab[ju].ene)-log(eng_tab[jl].ene))*(log(ave_eng)-log(eng_tab[jl].ene)) )*spe_nf;
				nitd_e = exp( log(nitd[jl][k])+(log(nitd[ju][k])-log(nitd[jl][k]))/(log(eng_tab[ju].ene)-log(eng_tab[jl].ene))*(log(ave_eng)-log(eng_tab[jl].ene)) )*spe_nf;
				othd_e = exp( log(othd[jl][k])+(log(othd[ju][k])-log(othd[jl][k]))/(log(eng_tab[ju].ene)-log(eng_tab[jl].ene))*(log(ave_eng)-log(eng_tab[jl].ene)) )*spe_nf;
				nd_e = exp( log(nd[jl][k])+(log(nd[ju][k])-log(nd[jl][k]))/(log(eng_tab[ju].ene)-log(eng_tab[jl].ene))*(log(ave_eng)-log(eng_tab[jl].ene)) )*spe_nf;
				gd_e = exp( log(gd[jl][k])+(log(gd[ju][k])-log(gd[jl][k]))/(log(eng_tab[ju].ene)-log(eng_tab[jl].ene))*(log(ave_eng)-log(eng_tab[jl].ene)) )*spe_nf;
				Nbd_e = exp( log(Nbd[jl][k])+(log(Nbd[ju][k])-log(Nbd[jl][k]))/(log(eng_tab[ju].ene)-log(eng_tab[jl].ene))*(log(ave_eng)-log(eng_tab[jl].ene)) )*spe_nf;
				Sbd_e = exp( log(Sbd[jl][k])+(log(Sbd[ju][k])-log(Sbd[jl][k]))/(log(eng_tab[ju].ene)-log(eng_tab[jl].ene))*(log(ave_eng)-log(eng_tab[jl].ene)) )*spe_nf;
				Mbd_e = exp( log(Mbd[jl][k])+(log(Mbd[ju][k])-log(Mbd[jl][k]))/(log(eng_tab[ju].ene)-log(eng_tab[jl].ene))*(log(ave_eng)-log(eng_tab[jl].ene)) )*spe_nf;
				Tbd_e = exp( log(Tbd[jl][k])+(log(Tbd[ju][k])-log(Tbd[jl][k]))/(log(eng_tab[ju].ene)-log(eng_tab[jl].ene))*(log(ave_eng)-log(eng_tab[jl].ene)) )*spe_nf;

				hydGy[k] += hydd_e;
				nitGy[k] += nitd_e;
				othGy[k] += othd_e;
				nGy[k] += nd_e;
				gGy[k] += gd_e;
				NbGy[k] += Nbd_e;
				SbGy[k] += Sbd_e;
				MbGy[k] += Mbd_e;
				TbGy[k] += Tbd_e;
			}
		}
		else
		{
			flg = 1;
		}

		i++;
	}
	fclose(fp);

//write file
	sprintf(outfile,"./depVSndose25ppm-fc-2.2e5_tc_9e11.dat");
	fp=fopen(outfile,"w");
	fprintf(fp,"#inputfile:%s\n#column1:depth[mm]\n#column2:hydrogen-dose[Gy-eq/h]\n#column3:nitrogen-dose[Gy-eq/h]\n#column4:other-N-dose[Gy-eq/h]\n#column5:N-dose[Gy-eq/h]\n#column6:G-dose[Gy/h]\n#column7:B-dose in other normal cells[Gy-eq/h]\n#column8:B-dose in skin cells[Gy-eq/h]\n#column9:B-dose in mucosa cells[Gy-eq/h]\n#column10:B-dose in tumor cells[Gy-eq/h]\n#column11:Total dose in other normal cells[Gy-eq/h]\n#column12:Total dose in skin cells[Gy-eq/h]\n#column13:Total dose in mucosa cells[Gy-eq/h]\n#column14:Total dose in tumor cells[Gy-eq/h]\n",filename);
	fprintf(fp,"#constants:Cb(N)=%.1fppm, RBE for thermal neutron=%.1f, RBE for fast neutron=%.1f, CBE(N)=%.2f, CBE(S)=%.1f, CBE(M)=%.1f, CBE(T)=%.1f, T/N=%.1f\n", Cbn, RBEthrm, RBEfast, CBEn, CBEs, CBEm, CBEt, TNR);
	for(i=0;i<num_tally;i++){
		hydGy[i]*=3.6E-9*PI*5.0*5.0;
		nitGy[i]*=3.6E-9*PI*5.0*5.0;
		othGy[i]*=3.6E-9*PI*5.0*5.0;
		nGy[i]*=3.6E-9*PI*5.0*5.0;
		gGy[i]*=3.6E-9*PI*5.0*5.0;
		NbGy[i]*=3.6E-9*PI*5.0*5.0;
		SbGy[i]*=3.6E-9*PI*5.0*5.0;
		MbGy[i]*=3.6E-9*PI*5.0*5.0;
		TbGy[i]*=3.6E-9*PI*5.0*5.0;
		NGy[i]=nGy[i]+gGy[i]+NbGy[i];
		SGy[i]=nGy[i]+gGy[i]+SbGy[i];
		MGy[i]=nGy[i]+gGy[i]+MbGy[i];
		TGy[i]=nGy[i]+gGy[i]+TbGy[i];
		fprintf(fp,"%d\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\n",i,hydGy[i],nitGy[i],othGy[i],nGy[i],gGy[i],NbGy[i],SbGy[i],MbGy[i],TbGy[i],NGy[i],SGy[i],MGy[i],TGy[i]);
	}
	fclose(fp);	

// release memory
	free(eng_tab);

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
	free(gd);
	for(i=0;i<num_tab;i++){
		free(gd[i]);
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

	free(hydGy);
	free(nitGy);
	free(othGy);
	free(nGy);
	free(gGy);
	free(NbGy);
	free(SbGy);
	free(MbGy);
	free(TbGy);

	return 1;
}
