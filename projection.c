#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "mpi.h"
#include "allvars.h"
#define dim 3
#define number 5
#define perpendicular_max 28
#define parallel_max 500
#define lens (perpendicular_max*parallel_max)
#define root  0

void pair_count(int n, int tree_n, double *position, double *tree_position, struct kdtree *tree, long long *nn ,double rmax2)
{

	int i, j, di;
	int nfound;
	int nalloc;
	int para_loc,per_loc,loc;
	double dx,dy,dz,adx,ady,adz,para,per;
	struct kdtree_result *nnresults=NULL;
	double find[3];
	double distance;	

	for(i=0; i<n; i++)
	{
		for(di=0; di<dim; di++)
		{
			find[di] = position[i*dim+di];
//		printf("find=%lf\n", find[di]);
		}
//		printf("haha\n");
		nalloc = tree_n;
//		printf("nalloc=%d\n", nalloc);
		nnresults = (struct kdtree_result*)malloc(nalloc*sizeof(struct kdtree_result));
		kdtree_r_nearest(&tree,find,rmax2,&nfound,nalloc,nnresults);
		if(nfound>0)
		{
			for(j=0; j<nfound; j++)
			{	
				distance = sqrt(nnresults[j].dis);
//				printf("distance=%lf\n", distance);
				loc = nnresults[j].idx;
				dx	= position[i*dim]   - tree_position[loc*dim];
				dy 	= position[i*dim+1] - tree_position[loc*dim+1];
				dz 	= position[i*dim+2] - tree_position[loc*dim+2];
				adx = position[i*dim]   + tree_position[loc*dim];
				ady = position[i*dim+1] + tree_position[loc*dim+1];
				adz = position[i*dim+2] + tree_position[loc*dim+2];
				para= ((dx*adx) + (dy*ady) + (dz*adz))/sqrt(adx*adx+ady*ady+adz*adz);
				per = sqrt((dx*dx + dy*dy + dz*dz) - para*para);
				//dx	= position[i*dim]   - tree_position[loc*dim];
				//dy 	= position[i*dim+1] - tree_dposition[loc*dim+1];
				//dz 	= position[i*dim+2] - tree_dposition[loc*dim+2];
//				para = dx;
//				per= sqrt(dy*dy + dz*dz);
//				printf("%lf\t%lf\n", para,per);
				if((fabs(para)<parallel_max)&&(per<57)&&(per>=0.0891))
				{						
					para_loc=(int)(para/2.0+parallel_max/2.0);
					per_loc=(int)((log10(per)+1.05)/0.1);
//					printf("para_loc%d\n", para_loc);
//					printf("per_loc%d\n", per_loc);
					if((per_loc>=0 )&&(per_loc<perpendicular_max))
					{	
						nn[para_loc*perpendicular_max+per_loc]=nn[para_loc*perpendicular_max+per_loc]+1; //attention!There's some problem in the first point(per_loc),because (int)(-0.5)=0
					}
				}
			}
		}
		free(nnresults);
	}
// 	for(i=0; i<parallel_max; i++)
// 	{
// 		for(j=0; j<perpendicular_max; j++)
// 		{
// 			printf("%d\t%d\t%lld\n", i, j, nn[i*perpendicular_max+j]);
// 		}
// 	}
// 	printf("nnfinish\n");	
}

int main(int argc, char **argv)
{	
	MPI_Init(&argc, &argv);

	int i,j,myid,m,low_limit,high_limit;
	int subsample;
	int ngalaxy=0;
	int rgalaxy=0;
	double rmax2;
	double d_divide_r;
	double *dposition=NULL;
	double *rposition=NULL;
	double *alldposition=NULL;
	double *allrposition=NULL;
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
	double *result_wp;
	FILE *parafile,*f20,*f30,*f40,*f41,*f50,*f51;
	FILE *f_all;
	FILE *f_ind;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	start =clock();

//--------------------------------------------------------------------------read data	
	parafile = fopen("infile", "r");

	fscanf(parafile, "%s",  file_data_all);	
	fscanf(parafile, "%s",  file_random_all);	
	fscanf(parafile, "%d",  &data_length_all);	
	fscanf(parafile, "%d",  &random_length_all);
	fscanf(parafile, "%d",  &subsample);
	fscanf(parafile, "%lf", &rmax2);	
	fscanf(parafile, "%s",  file_data);	
	fscanf(parafile, "%s",  file_random);
	fscanf(parafile, "%s",  data_length);	
	fscanf(parafile, "%s",  random_length);
	fclose(parafile);



	f20 = fopen(file_data_all, "r");
	f30 = fopen(file_random_all, "r");
	if(f20 ==NULL)printf("cao");
	if(f30 ==NULL)printf("cao");
	alldposition = (double*)malloc(dim* data_length_all*sizeof(double));
	allrposition = (double*)malloc(dim* random_length_all*sizeof(double));

//	printf("good\n");  
	fflush(stdout);
	for(i=0; i<data_length_all; i++)
	{	
		for(j=0; j<dim; j++)
		{
			fscanf(f20, "%lf", &alldposition[i*dim+j]);
//			printf("%lf\t",alldposition[i*dim+j]);
		}
//		printf("\n");
	}
	fclose(f20);

	for(i=0; i<random_length_all; i++)
	{	
		for(j=0; j<dim; j++)
		{
			fscanf(f30, "%lf", &allrposition[i*dim+j]);
//			printf("%lf\t",allrposition[i*dim+j]);
		}
//		printf("\n");
	}
	fclose(f30);

//--------------------------------------------------------------------------create kdtree

	ddall_tree    = kdtree_create(alldposition, data_length_all, FALSE, TRUE);
	rrall_tree    = kdtree_create(allrposition, random_length_all, FALSE, TRUE);
//--------------------------------------------------------------------------create finish

	result_wp = (double*)malloc(subsample*perpendicular_max*sizeof(double));
	dd_patch_patch = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	dr_patch_patch = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	rr_patch_patch = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));	
	dd_patch_all   = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	dr_patch_all   = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	rr_patch_all   = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	dd_all_all     = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	dr_all_all     = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	rr_all_all     = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	jackdd = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	jackdr = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	jackrr = (long long*)malloc(parallel_max*perpendicular_max*sizeof(long long));
	correlation = (double*)malloc(parallel_max*perpendicular_max*sizeof(double));
	w = (double*)malloc(perpendicular_max*sizeof(double));

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

	dposition = (double*)malloc(dim*ngalaxy*sizeof(double));
	rposition = (double*)malloc(dim*rgalaxy*sizeof(double));

	for(i=0; i<ngalaxy; i++)
	{
		for(j=0; j<dim; j++)
		{
			fscanf(f40, "%lf", &dposition[i*dim+j]);
//			printf("%lf\t", dposition[i*dim+j]);
		}
//		printf("\n");		
	}
	fclose(f40);

	for(i=0; i<rgalaxy; i++)
	{
		for(j=0; j<dim; j++)
		{
			fscanf(f50, "%lf", &rposition[i*dim+j]);
//			printf("%lf\t", rposition[i*dim+j]);
		}
//		printf("\n");		
	}
	fclose(f50);

	ddpatch_tree  = kdtree_create(dposition, ngalaxy, FALSE, TRUE);
	rrpatch_tree  = kdtree_create(rposition, rgalaxy, FALSE, TRUE);	
//--------------------------------------------------------------------------read finish
	
//--------------------------------------------------------------------------find nearest node
//	printf("dd_patch_patch\n");
	pair_count(ngalaxy, ngalaxy,            dposition, dposition,    ddpatch_tree, dd_patch_patch, rmax2);
//	printf("dd_patch_all\n");
	pair_count(ngalaxy, data_length_all,    dposition, alldposition, ddall_tree,   dd_patch_all,   rmax2);
//	printf("dr_patch_patch\n");	
	pair_count(ngalaxy, rgalaxy,            dposition, rposition,    rrpatch_tree, dr_patch_patch, rmax2);
//	printf("dr_patch_all\n");
	pair_count(ngalaxy, random_length_all,  dposition, allrposition, rrall_tree,   dr_patch_all,   rmax2);
//	printf("rr_patch_patch\n");	
	pair_count(rgalaxy, rgalaxy,            rposition, rposition,    rrpatch_tree, rr_patch_patch, rmax2);
//	printf("rr_patch_all\n");
	pair_count(rgalaxy, random_length_all,  rposition, allrposition, rrall_tree,   rr_patch_all,   rmax2);	

//--------------------------------------------------------------------------find  finish

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
		for(i=0; i<parallel_max; i++)
		{ 
			for(j=0; j<perpendicular_max; j++)
			{
				fprintf(f_all,"%d\t%d\t%lld\t%lld\t%lld\n", i, j, dd_all_all[i*perpendicular_max+j], 
					dr_all_all[i*perpendicular_max+j], rr_all_all[i*perpendicular_max+j]);
			}
		}
		fclose(f_all);
	}
	MPI_Barrier(MPI_COMM_WORLD);


	f_ind = fopen(strcat(file_ind,filesuffix) ,"w"); 

	jackngalaxy = data_length_all - ngalaxy;
	jackrgalaxy = random_length_all - rgalaxy;
	d_divide_r = jackngalaxy*1.0/jackrgalaxy;
	for(i=0; i<parallel_max; i++)
	{
		for(j=0; j<perpendicular_max; j++)
		{
			jackdd[i*perpendicular_max+j] = dd_all_all[i*perpendicular_max+j] - 2*dd_patch_all[i*perpendicular_max+j] + dd_patch_patch[i*perpendicular_max+j];
			jackdr[i*perpendicular_max+j] = dr_all_all[i*perpendicular_max+j] - 2*dr_patch_all[i*perpendicular_max+j] + dr_patch_patch[i*perpendicular_max+j];
			jackrr[i*perpendicular_max+j] = rr_all_all[i*perpendicular_max+j] - 2*rr_patch_all[i*perpendicular_max+j] + rr_patch_patch[i*perpendicular_max+j];
			if(jackrr[i*perpendicular_max+j]!=0)
			{
				correlation[i*perpendicular_max+j] = (jackdd[i*perpendicular_max+j]-
					2*jackdr[i*perpendicular_max+j]*d_divide_r+
					jackrr[i*perpendicular_max+j]*d_divide_r*d_divide_r)/
				(jackrr[i*perpendicular_max+j]*d_divide_r*d_divide_r);
			}
			else
			{
				correlation[i*perpendicular_max+j] = -10000;
			}
			
		fprintf(f_ind, "%d\t%d\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n", i, j, dd_patch_patch[i*perpendicular_max+j],
		dr_patch_patch[i*perpendicular_max+j],rr_patch_patch[i*perpendicular_max+j],dd_patch_all[i*perpendicular_max+j],
		dr_patch_all[i*perpendicular_max+j],rr_patch_all[i*perpendicular_max+j],jackdd[i*perpendicular_max+j],
		jackdr[i*perpendicular_max+j],jackrr[i*perpendicular_max+j]);
//		printf("myid=%d,corr[%d][%d]=%.18lf\n",myid, i,j, correlation[i*perpendicular_max+j]);
		}
	}
	fclose(f_ind);

	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	for(m=1; m<=number; m++)
	{
		low_limit=parallel_max/2-50*m;
		high_limit=parallel_max/2+50*m;
		for(j=0; j<perpendicular_max; j++)	
		{ 
			w[j]=0;
		}
		for(j=0; j<subsample*perpendicular_max; j++)	
		{ 
			result_wp[j]=0;
		}
		for(j=0; j<perpendicular_max; j++)			
		{
			for(i=low_limit; i<high_limit; i++)		
			{ 	
				if(correlation[i*perpendicular_max+j]!=-10000)
				{
					w[j] = w[j]+correlation[i*perpendicular_max+j]*2;
				}
				else
				{
					w[j] = 0;
					break;
				}
			}
//			printf("%lf\n", w[j]);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(w, perpendicular_max, MPI_DOUBLE, result_wp, perpendicular_max, MPI_DOUBLE, root, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if(myid==root)
		{	
			for(i=0 ;i<subsample*perpendicular_max; i++)
			{
				printf("%.10lf\n",result_wp[i]);
			}
			
		}
		//MPI_Barrier(MPI_COMM_WORLD);
		fflush(stdout);
	}
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
	free(w);
	free(alldposition);
	free(allrposition);
	free(result_wp);

	MPI_Finalize();
	end =clock();
	cost = (end -start)/CLOCKS_PER_SEC ;
	printf("time=%lf\n",cost);

	return 0;
}
