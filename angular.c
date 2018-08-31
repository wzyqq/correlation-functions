#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "mpi.h"
#include "allvars.h"
#define dim 3
#define dim_coordinate 2
#define theta_max 16
#define lens (theta_max)
#define root  0
#define rr0 -10000
#define theta_low 0.00079 //10**(-3.1)
#define theta_high 1.26 //10**(0.1)
#define pi 3.141592654
#define ang2rad pi/180.0
#define rad2ang 180.0/pi

void pair_count(int n, int tree_n, double *position, double *tree_position, struct kdtree *tree, long long *nn ,double thetamax2)
{
	double loc_temp;
	int i, j, di,loc;
	int nfound;
	int nalloc;
	struct kdtree_result *nnresults=NULL;
	double find[3];
	double distance;	
	double theta;


	for(i=0; i<n; i++)
	{
		for(di=0; di<dim; di++)
		{
			find[di] = position[i*dim+di];
//			printf("find=%lf\n", find[di]);
		}
//		printf("haha\n");
		nalloc = tree_n;
//		printf("nalloc=%d\n", nalloc);
		nnresults = (struct kdtree_result*)malloc(nalloc*sizeof(struct kdtree_result));
		kdtree_r_nearest(&tree,find,thetamax2,&nfound,nalloc,nnresults);
		if(nfound>0)
		{
			for(j=0; j<nfound; j++)
			{	
				distance = sqrt(nnresults[j].dis);
				theta = 2*asin(distance*0.5)*rad2ang;
//				printf("distance=%lf\n", distance);
				if((theta<theta_high)&&(theta>=theta_low))
				{						
					loc_temp = (log10(theta)*5+15.5);
					if(loc_temp>=0)
					{
						loc=(int)loc_temp;
						if(loc<theta_max)
						{	
							nn[loc]=nn[loc]+1;
						}
					}
				}
			}
		}
		free(nnresults);
	}

// 		for(j=0; j<theta_max; j++)
// 		{
// 			printf("%d\t%d\t%lld\n", i, j, nn[j]);
// 		}
// 	printf("nnfinish\n");	
}

int main(int argc, char **argv)
{	
	MPI_Init(&argc, &argv);

	int i,j,myid;
	int di;
	int nfound;
	int nalloc;
	int subsample;
	int ngalaxy;
	int rgalaxy;
	int loc;
	double d_divide_r;
	long long *dd=NULL;
	long long *dr=NULL;
	long long *rr=NULL;
	double *dposition=NULL;
	double *rposition=NULL;
	double *alldposition=NULL;
	double *allrposition=NULL;
	double *dcoordinate=NULL;
	double *rcoordinate=NULL;
	double *alldcoordinate=NULL;
	double *allrcoordinate=NULL;	
	double find[3];
	double thetamax2;
	double distance;
	double *correlation=NULL;
	double *w=NULL;
	long long *dd_patch_patch=NULL;
	long long *dr_patch_patch=NULL;
	long long *rr_patch_patch=NULL;	
	long long *dd_patch_all=NULL;
	long long *dr_patch_all=NULL;
	long long *rr_patch_all=NULL;	
	long long *dd_all_all=NULL;
	long long *dr_all_all=NULL;
	long long *rr_all_all=NULL;
	long long *jackdd=NULL;
	long long *jackdr=NULL;
	long long *jackrr=NULL;
	int jackngalaxy;
	int jackrgalaxy;
	struct kdtree *ddtree=NULL;
	struct kdtree *rrtree=NULL;
	struct kdtree *ddpatch_tree=NULL;
	struct kdtree *rrpatch_tree=NULL;
	struct kdtree *ddall_tree=NULL;
	struct kdtree *rrall_tree=NULL;	
	char file_data[200];
	char file_random[200];
	char data_length[200];
	char random_length[200];
	char file_data_all[200];
	char file_random_all[200];
	int data_length_all;
	int random_length_all;
	char filesuffix[20];
	char file_ind[200] = "all_ind";
	double start,end,cost;
	long long *result_jackrr;		
	double *result_correlation;
	double theta;
	FILE *parafile,*f20,*f30,*f40,*f41,*f50,*f51;
	FILE *f_all;
	FILE *f_ind;
	FILE *f_rr;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	start =clock();

//--------------------------------------------------------------------------read data	
	parafile = fopen("infile", "r");

	fscanf(parafile, "%s",  file_data_all);	
	fscanf(parafile, "%s",  file_random_all);	
	fscanf(parafile, "%d",  &data_length_all);	
	fscanf(parafile, "%d",  &random_length_all);
	fscanf(parafile, "%d",  &subsample);
	fscanf(parafile, "%lf", &thetamax2);	
	fscanf(parafile, "%s",  file_data);	
	fscanf(parafile, "%s",  file_random);
	fscanf(parafile, "%s",  data_length);	
	fscanf(parafile, "%s",  random_length);
	fclose(parafile);



	f20 = fopen(file_data_all, "r");
	f30 = fopen(file_random_all, "r");
	if(f20 ==NULL)printf("cao");
	if(f30 ==NULL)printf("cao");
	alldcoordinate = (double*)malloc(dim_coordinate* data_length_all*sizeof(double));
	allrcoordinate = (double*)malloc(dim_coordinate* random_length_all*sizeof(double));	
	alldposition = (double*)malloc(dim* data_length_all*sizeof(double));
	allrposition = (double*)malloc(dim* random_length_all*sizeof(double));


//	printf("good\n");  fflush(stdout);
	for(i=0; i<data_length_all; i++)
	{	
		for(j=0; j<dim_coordinate; j++)
		{
			fscanf(f20, "%lf", &alldcoordinate[i*dim_coordinate+j]);
//			printf("%lf\n", alldcoordinate[i*dim_coordinate+j]);
		}
		alldposition[i*dim+0] = cos(alldcoordinate[i*dim_coordinate+1]*ang2rad)*cos(alldcoordinate[i*dim_coordinate+0]*ang2rad);
		alldposition[i*dim+1] = cos(alldcoordinate[i*dim_coordinate+1]*ang2rad)*sin(alldcoordinate[i*dim_coordinate+0]*ang2rad);
		alldposition[i*dim+2] = sin(alldcoordinate[i*dim_coordinate+1]*ang2rad);		
//		printf("%.8lf\t%.8lf\t%.8lf\n", alldposition[i*dim+0], alldposition[i*dim+1], alldposition[i*dim+2]);
	}
	fclose(f20);
//	printf("good1\n");

	for(i=0; i<random_length_all; i++)
	{	
		for(j=0; j<dim_coordinate; j++)
		{
			fscanf(f30, "%lf", &allrcoordinate[i*dim_coordinate+j]);
//			printf("%lf\n", allrcoordinate[i*dim_coordinate+j]);			
		}		
		allrposition[i*dim+0] = cos(allrcoordinate[i*dim_coordinate+1]*ang2rad)*cos(allrcoordinate[i*dim_coordinate+0]*ang2rad);
		allrposition[i*dim+1] = cos(allrcoordinate[i*dim_coordinate+1]*ang2rad)*sin(allrcoordinate[i*dim_coordinate+0]*ang2rad);
		allrposition[i*dim+2] = sin(allrcoordinate[i*dim_coordinate+1]*ang2rad);
//		printf("%.8lf\t%.8lf\t%.8lf\n", allrposition[i*dim+0], allrposition[i*dim+1], allrposition[i*dim+2]);		
	}
	fclose(f30);
//	printf("good2\n");

//--------------------------------------------------------------------------create kdtree
	ddall_tree    = kdtree_create(alldposition, data_length_all, FALSE, TRUE);
	rrall_tree    = kdtree_create(allrposition, random_length_all, FALSE, TRUE);
//--------------------------------------------------------------------------create finish

	result_jackrr = (long long*)malloc(subsample*theta_max*sizeof(long long));
	result_correlation = (double*)malloc(subsample*theta_max*sizeof(double));
	dd_patch_patch = (long long*)malloc(theta_max*sizeof(long long));
	dr_patch_patch = (long long*)malloc(theta_max*sizeof(long long));
	rr_patch_patch = (long long*)malloc(theta_max*sizeof(long long));	
	dd_patch_all   = (long long*)malloc(theta_max*sizeof(long long));
	dr_patch_all   = (long long*)malloc(theta_max*sizeof(long long));
	rr_patch_all   = (long long*)malloc(theta_max*sizeof(long long));
	dd_all_all     = (long long*)malloc(theta_max*sizeof(long long));
	dr_all_all     = (long long*)malloc(theta_max*sizeof(long long));
	rr_all_all     = (long long*)malloc(theta_max*sizeof(long long));
	jackdd = (long long*)malloc(theta_max*sizeof(long long));
	jackdr = (long long*)malloc(theta_max*sizeof(long long));
	jackrr = (long long*)malloc(theta_max*sizeof(long long));
	correlation = (double*)malloc(theta_max*sizeof(double));

	sprintf(filesuffix, "_%d", myid);
	
//	printf("%s\n", file_data); fflush(stdout);
	f40 = fopen(strcat(file_data,filesuffix) ,"r"); 
	f50 = fopen(strcat(file_random,filesuffix) ,"r"); 
	f41 = fopen(strcat(data_length,filesuffix) ,"r"); 
	f51 = fopen(strcat(random_length,filesuffix) ,"r");

//	printf("%s\n", file_data);
//	printf("ff%s\n", filesuffix);
	if(f40 ==NULL)printf("cao40");
	if(f50 ==NULL)printf("cao50");
	if(f41 ==NULL)printf("cao41");
	if(f51 ==NULL)printf("cao51");
	fscanf(f41, "%d", &ngalaxy);
	fscanf(f51, "%d", &rgalaxy);
	fclose(f41);
	fclose(f51);

	dcoordinate = (double*)malloc(dim_coordinate* ngalaxy*sizeof(double));
	rcoordinate = (double*)malloc(dim_coordinate* rgalaxy*sizeof(double));	
	dposition = (double*)malloc(dim*ngalaxy*sizeof(double));
	rposition = (double*)malloc(dim*rgalaxy*sizeof(double));

	for(i=0; i<ngalaxy; i++)
	{
		for(j=0; j<dim_coordinate; j++)
		{
			fscanf(f40, "%lf", &dcoordinate[i*dim_coordinate+j]);
//			printf("%lf\n", dcoordinate[i*dim_coordinate+j]);			
		}
		dposition[i*dim+0] = cos(dcoordinate[i*dim_coordinate+1]*ang2rad)*cos(dcoordinate[i*dim_coordinate+0]*ang2rad);
		dposition[i*dim+1] = cos(dcoordinate[i*dim_coordinate+1]*ang2rad)*sin(dcoordinate[i*dim_coordinate+0]*ang2rad);
		dposition[i*dim+2] = sin(dcoordinate[i*dim_coordinate+1]*ang2rad);
//		printf("%.8lf\t%.8lf\t%.8lf\n", dposition[i*dim+0], dposition[i*dim+1], dposition[i*dim+2]);		
	}
	fclose(f40);
//	printf("good3\n");

	for(i=0; i<rgalaxy; i++)
	{
		for(j=0; j<dim_coordinate; j++)
		{
			fscanf(f50, "%lf", &rcoordinate[i*dim_coordinate+j]);
//			printf("%lf\n", rcoordinate[i*dim_coordinate+j]);						
		}
		rposition[i*dim+0] = cos(rcoordinate[i*dim_coordinate+1]*ang2rad)*cos(rcoordinate[i*dim_coordinate+0]*ang2rad);
		rposition[i*dim+1] = cos(rcoordinate[i*dim_coordinate+1]*ang2rad)*sin(rcoordinate[i*dim_coordinate+0]*ang2rad);
		rposition[i*dim+2] = sin(rcoordinate[i*dim_coordinate+1]*ang2rad);	
//		printf("%.8lf\t%.8lf\t%.8lf\n", rposition[i*dim+0], rposition[i*dim+1], rposition[i*dim+2]);				
	}
	fclose(f50);
//	printf("good4\n");


	ddpatch_tree  = kdtree_create(dposition, ngalaxy, FALSE, TRUE);
	rrpatch_tree  = kdtree_create(rposition, rgalaxy, FALSE, TRUE);	
//--------------------------------------------------------------------------read finish
	
//--------------------------------------------------------------------------find dd nearest node
//	printf("dd_patch_patch\n");
	pair_count(ngalaxy, ngalaxy,            dposition, dposition,    ddpatch_tree, dd_patch_patch, thetamax2);
//	printf("dd_patch_all\n");
	pair_count(ngalaxy, data_length_all,    dposition, alldposition, ddall_tree,   dd_patch_all,   thetamax2);
//	printf("dr_patch_patch\n");	
	pair_count(ngalaxy, rgalaxy,            dposition, rposition,    rrpatch_tree, dr_patch_patch, thetamax2);
//	printf("dr_patch_all\n");
	pair_count(ngalaxy, random_length_all,  dposition, allrposition, rrall_tree,   dr_patch_all,   thetamax2);
//	printf("rr_patch_patch\n");	
	pair_count(rgalaxy, rgalaxy,            rposition, rposition,    rrpatch_tree, rr_patch_patch, thetamax2);
//	printf("rr_patch_all\n");
	pair_count(rgalaxy, random_length_all,  rposition, allrposition, rrall_tree,   rr_patch_all,   thetamax2);	

//--------------------------------------------------------------------------find rr finish

	free(dposition);
 	free(rposition);

//--------------------------------------------------------------------------destroy kdtree	
	kdtree_destroy(&ddall_tree);
	kdtree_destroy(&rrall_tree);
	kdtree_destroy(&ddpatch_tree);
	kdtree_destroy(&rrpatch_tree);	
//	printf("jie\n");
//--------------------------------------------------------------------------destroy finish		

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(dd_patch_all, dd_all_all, lens, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(dr_patch_all, dr_all_all, lens, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(rr_patch_all, rr_all_all, lens, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if(myid==root)
	{	
		f_all =	fopen("alldddrrr","w");				
		for(j=0; j<theta_max; j++)
		{
			fprintf(f_all,"%d\t%lld\t%lld\t%lld\n", j, dd_all_all[j], dr_all_all[j], rr_all_all[j]);
		}
		fclose(f_all);
	}
	MPI_Barrier(MPI_COMM_WORLD);


	f_ind = fopen(strcat(file_ind,filesuffix) ,"w"); 

	jackngalaxy = data_length_all - ngalaxy;
	jackrgalaxy = random_length_all - rgalaxy;
	d_divide_r = jackngalaxy*1.0/jackrgalaxy;

	for(j=0; j<theta_max; j++)
	{
		jackdd[j] = dd_all_all[j] - 2*dd_patch_all[j] + dd_patch_patch[j];
		jackdr[j] = dr_all_all[j] - 2*dr_patch_all[j] + dr_patch_patch[j];
		jackrr[j] = rr_all_all[j] - 2*rr_patch_all[j] + rr_patch_patch[j];
		if(jackrr[j]!=0)
		{
			correlation[j] = (jackdd[j]-
				2*jackdr[j]*d_divide_r+
				jackrr[j]*d_divide_r*d_divide_r)/
			(jackrr[j]*d_divide_r*d_divide_r);
		}
		else
		{
			correlation[j] = rr0;
		}
		
	fprintf(f_ind, "%d\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n", j, dd_patch_patch[j],
	dr_patch_patch[j],rr_patch_patch[j],dd_patch_all[j],dr_patch_all[j],rr_patch_all[j],jackdd[j],jackdr[j],jackrr[j]);
//	printf("myid=%d,corr[%d][%d]=%.18lf\n",myid, i,j, correlation[j]);
	}
	fclose(f_ind);


	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(jackrr, theta_max, MPI_DOUBLE, result_jackrr, theta_max, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
		
	if(myid==root)
	{	
		f_rr = fopen("allrr","w");
		for(i=0 ;i<subsample*theta_max; i++)
		{
			fprintf(f_rr, "%lld\n", result_jackrr[i]);
		}	
	}


	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Gather(correlation, theta_max, MPI_DOUBLE, result_correlation, theta_max, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==root)
	{	
		for(i=0 ;i<subsample*theta_max; i++)
		{
			printf("%.10lf\n",result_correlation[i]);
		}
		
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	fflush(stdout);

	free(dd_patch_patch);
	free(dr_patch_patch);
	free(rr_patch_patch);
	free(dd_patch_all);
	free(dr_patch_all);
	free(rr_patch_all);
	free(dd_all_all);
	free(dr_all_all);
	free(rr_all_all);
	free(jackdd);
	free(jackdr);
	free(jackrr);
	free(correlation);
	free(alldposition);
	free(allrposition);
	free(result_correlation);

	MPI_Finalize();
	end =clock();
	cost = (end -start)/CLOCKS_PER_SEC ;
	printf("time=%lf\n",cost);

	return 0;
}
