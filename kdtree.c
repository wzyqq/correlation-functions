/*************************************************************************
	> File Name: kdtree.c
	> Author: 
	> Mail: 
	> Created Time: 2018年03月03日 星期六 22时17分53秒
 ************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allvars.h"

struct kdtree *kdtree_create(double *input_data, int n, int sort, int rearrange)
{
    int i,j;
    struct kdtree *mr;
    mr = malloc(sizeof(struct kdtree));
    mr->the_data = input_data;
    mr->n = n;
//	printf("good\n");
    build_tree(&mr);

	mr->sort = sort;
	mr->rearrange = rearrange;

	if(mr->rearrange)
	{
		mr->rearranged_data = malloc(dim*(mr->n)*sizeof(double));
		for(i=0; i<n; i++)
        {
        	for(j=0; j<dim; j++)
        	{
        		mr->rearranged_data[i*dim+j] = mr->the_data[mr->ind[i]*dim+j];
//                printf("arr=%lf\n", mr->rearranged_data[i*dim+j]);
			}
		}
	}
    else
    {
		{
			mr->rearranged_data=NULL;
		}
	}
    // for(j=0;j<mr->n;j++)
    // {
    //     printf("%d\n", mr->ind[j]);
    // }
    return mr;
}
void build_tree(struct kdtree **tp)
{
    int j;
    struct tree_node *dummy=NULL;

    (*tp)->ind = malloc((*tp)->n*sizeof(int));
    for(j=0;j<(*tp)->n;j++)
    {
        (*tp)->ind[j]=j;
//		printf("%d\n", (*tp)->ind[j]);
    }
//	printf("good\n");
    (*tp)->root=malloc(sizeof(struct tree_node));
//	printf("%d\n",(*tp)->n-1 );
    (*tp)->root=build_tree_for_range(tp,0,(*tp)->n-1,&dummy);
}
struct tree_node *build_tree_for_range(struct kdtree **tp,int l,int u,struct tree_node **parent)
{
    struct tree_node *res;
    int i,c,m;
    int recompute;
    double average;
    double extent[dim];

    res        = malloc(sizeof(struct tree_node));    
    res->box   = malloc(dim*sizeof(struct interval));
    res->left  = malloc(sizeof(struct tree_node));
    res->right = malloc(sizeof(struct tree_node));

//	printf("good\n");
    if(u<l)
    {
        res = NULL;
        return res;
    }
    if((u-l)<=bucket_size)
    {
        for(i=0; i<dim; i++)
        {
            spread_in_coordinate(tp,i,l,u,&res->box[i]);
//			printf("finalnode%lf\t%lf\n", res->box[i].lower, res->box[i].upper);
        }
        res->cut_dim = 0;
        res->cut_val = 0;
        res->l = l;
        res->u = u;
//        printf("finalnode l=%d\t u=%d\n", res->l, res->u);
        res->left =NULL;
        res->right =NULL;
        return res; 
    }
    else
    {
        for(i=0;i<dim;i++)
        {
            recompute=TRUE;
//			printf("good1\n" );
            if((*parent)!=NULL)
            {
//				printf("biejia\n");
                if(i!=(*parent)->cut_dim)
                {
                    recompute=FALSE;
                }
            }
            if(recompute)
            {
//				printf("good2\n");
                spread_in_coordinate(tp,i,l,u,&res->box[i]);
//				printf("%lf\t%lf\n", res->box[i].lower, res->box[i].upper);
            }
            else
            {
//				printf("wocao\n");
                res->box[i]=(*parent)->box[i];
//				printf("qigua%lf\t%lf\n", res->box[i].lower, res->box[i].upper);
            }
        }
        for(i=0;i<dim;i++)
        {
			extent[i] = res->box[i].upper-res->box[i].lower;
//			printf("extent[%d]=%lf\n", i, extent[i]);
        }
        c = find_maxmum(extent);
//		printf("c = %d\n", c);
		// for(i=0; i<(*tp)->n*dim; i++)
		// {
		// 	printf("(*tp)->the_data[%d]=%lf\n", i, (*tp)->the_data[i]);
		// }
		// for(i=0; i<(*tp)->n; i++)
		// {
		// 	printf("(*tp)->ind[%d]=%d\n", i, (*tp)->ind[i]);
		// }
//		printf("l=%d\tu=%d\n", l ,u);
        average = sum((*tp)->the_data, (*tp)->ind, c, l, u)*1.0/(u-l+1);
//		printf("average=%lf\n", average);
        res->cut_val = average;
        m = select_on_coordinate_value((*tp)->the_data,(*tp)->ind,c,average,l,u);
//		printf("m=%d\n", m);
        res->cut_dim = c;
        res->l = l;
        res->u = u;
        res->left  = build_tree_for_range(tp,l,m,&res);
        res->right = build_tree_for_range(tp,m+1,u,&res);
        if(res->right==NULL)
        {
//        	printf("111\n");
            for(i=0;i<dim;i++)
            {
                res->box[i] = res->left->box[i];
            }
            res->cut_val_left = res->left->box[c].upper;
            res->cut_val = res->cut_val_left;
        }
        else if(res->left==NULL)
        {
//        	printf("222\n");
            for(i=0; i<dim; i++)
            {
                res->box[i] = res->right->box[i];
            }
            res->cut_val_right = res->right->box[c].upper;
            res->cut_val = res->cut_val_right;
        }
        else
        {
//			printf("333\n");
            res->cut_val_right = res->right->box[c].lower;
            res->cut_val_left = res->left->box[c].upper;
            res->cut_val = (res->cut_val_left + res->cut_val_right)/2.0;
            for(i=0; i<dim; i++)
            {
            	res->box[i].upper = max(res->left->box[i].upper,res->right->box[i].upper);
           		res->box[i].lower = min(res->left->box[i].lower,res->right->box[i].lower);
//           		printf("box=%lf\t%lf\n", res->box[i].upper,res->box[i].lower);
        	}
//        	printf("%lf\t%lf\t%lf\n", res->cut_val_right, res->cut_val_left, res->cut_val);
        }
    }
    return res;
}


void spread_in_coordinate(struct kdtree **tp,int c,int l,int u,struct interval *interv)
{
    double last,lmax,lmin,t,smin,smax;
    int i,ulocal;
    double *v;
    int *ind;

    v = (*tp)->the_data;
    ind = (*tp)->ind;
    smin = v[c+ind[l]*dim];
    smax = smin;
//	printf("smin=%lf\n", smin);

    ulocal = u;
//	printf("ulocal=%d\n",ulocal);
    for(i=l+2; i<=ulocal; i=i+2)
    { 
        lmin = v[c+ind[i-1]*dim];
        lmax = v[c+ind[i]*dim];
//		printf("lmin=%lf\tlmax=%lf\n", lmin, lmax);
//		printf("i=%d\n", i);
        if (lmin>lmax)
        {
            t = lmin;
            lmin = lmax;
            lmax = t;
        }
        if (smin>lmin) 
        {
            smin = lmin;
        }
        if (smax<lmax)
        {
            smax = lmax;
        }
        // printf("smin=%lf\n", smin);
        // printf("smax=%lf\n", smax);
    }
//	printf("woi=%d\n", i);
    if(i==(ulocal+1))
    {
        last = v[c+ind[ulocal]*dim];
//		printf("last=%lf\n", last);
        if (smin>last) 
        {
            smin = last;
        } 
        if (smax<last)
        {
            smax = last;
        }
//		printf("smin=%lf\n", smin);
//		printf("smax=%lf\n", smax);
    }
    interv->lower = smin;
    interv->upper = smax;
}


int select_on_coordinate_value(double *v,int *ind,int c,double alpha,int li,int ui)
{
	int i;
    int lb,rb;
    int tmp;
    int res;
    lb = li;
    rb = ui;
    while (lb < rb)
    {
        if ( v[c+ind[lb]*dim] <= alpha )
        {
            lb = lb+1;
        }
        else
        {
            tmp = ind[lb]; 
            ind[lb] = ind[rb]; 
            ind[rb] = tmp;
            rb = rb-1;
        }
    }
    if (v[c+ind[lb]*dim] <= alpha)
    {   
        res = lb;
    }
    else
    {
        res = lb-1;
    }
    // for(i=li; i<=ui; i++)
    // {
    // 	printf("zheshishenmea=%lf\n", v[c+ind[i]*dim]);
    // }
    return res;
}


double max(double a, double b)
{
	if(a<=b)
	{
		return b;
	}
	else
	{
		return a;
	}
}

double min(double a, double b)
{
	if(a<=b)
	{
		return a;
	}
	else
	{
		return b;
	}
}
int find_maxmum(double *array)
{
	int i,l,t_loc,last_loc;
	double lmin,lmax;
	double smin,smax;
    double t,last;
	int lmin_loc,lmax_loc;
	int smin_loc,smax_loc;
    smax = array[0];
    smin = array[0];
    smin_loc = 0;
    smax_loc = 0;

    l = 0;

    for(i=l+2; i<=dim-1; i=i+2)
    { 
        lmin = array[i-1];
        lmax = array[i];
        lmin_loc = i-1;
        lmax_loc = i;
        if (lmin>lmax)
        {
            t = lmin;
            lmin = lmax;
            lmax = t;
            t_loc = lmin_loc;
            lmin_loc = lmax_loc;
            lmax_loc = t_loc;
        }
        if (smin>lmin) 
        {
            smin = lmin;
            smin_loc = lmin_loc;
        }
        if (smax<lmax)
        {
            smax = lmax;
        	smax_loc = lmax_loc;
        }
    }
    if(i==dim)
    {
        last = array[dim-1];
        last_loc = dim-1;
        if (smin>last) 
        {
            smin = last;
            smin_loc = last_loc;
        } 
        if (smax<last)
        {
            smax = last;
            smax_loc = last_loc;
        }
    }
    return smax_loc;
}
double sum(double *array, int *ind, int c, int location1, int location2)
{
	int i;
	double result=0;
//	printf("location1=%d\tlocation2=%d\n", location1,location2);
	for(i=location1; i<=location2; i++)
	{
		result = result + array[ind[i]*dim+c];
	}
	return result;
}



void kdtree_r_nearest(struct kdtree **tp,double *qv,double r2,int *nfound,int nalloc,struct kdtree_result *results) 
{
//	int i;
	sr = malloc(sizeof(struct tree_search_record));
    sr->qv = qv;
    sr->ballsize = r2;
    sr->nfound = 0;
    sr->results = results;
    sr->nalloc = nalloc;
    sr->overflow = FALSE; 
    sr->ind = (*tp)->ind;
    sr->rearrange= (*tp)->rearrange;
//    printf("good\n");
    if ((*tp)->rearrange)
    {
//    	printf("lei\n");
    	sr->data = (*tp)->rearranged_data;
    }
    else
    {
//    	printf("yuan\n");
    	sr->data = (*tp)->the_data;
		// for(i=0; i<(*tp)->n*dim; i++)
		// {
		// 	printf("(*tp)->the_data[%d]=%lf\n", i, (*tp)->the_data[i]);
		// }
    }

    search(&(*tp)->root);
    *nfound = sr->nfound;
    if ((*tp)->sort)
    {
//    	printf("buhuiba\n");
		kdtree_sort_results(*nfound, results);
    }
    if (sr->overflow)
    {
    	printf("KD_TREE_TRANS: warning! return from kdtree_r_nearest found more neighbors");
    }
}


void search(struct tree_node **node)
{
    struct tree_node *ncloser, *nfarther;
    int cut_dim, i;
    double qval, dis;
    double ballsize;
    double *qv;
    struct interval *box; 

    if (((*node)->left==NULL) || ((*node)->right==NULL))
    {
//    	printf("zuihou\n");
		process_terminal_node_fixedball(node);
    }
    else
    {
//    	printf("hao\n");
		qv = sr->qv;
		cut_dim = (*node)->cut_dim;
		qval = qv[cut_dim];
//		printf("qval=%lf\n", qval);
        if (qval < (*node)->cut_val)
        {
            ncloser = (*node)->left;
            nfarther = (*node)->right;
//            dis = pow((*node)->cut_val_right - qval,2.0);
             dis = ((*node)->cut_val_right - qval) * ((*node)->cut_val_right - qval);
        }
        else
        {
            ncloser = (*node)->right;
            nfarther = (*node)->left;
//            dis = pow((*node)->cut_val_left - qval,2.0);
            dis = ((*node)->cut_val_left - qval) *((*node)->cut_val_left - qval);
        }

        if (ncloser!=NULL)  
        {
            search(&ncloser);
        }
        if (nfarther!=NULL)
        {
            ballsize = sr->ballsize;
            if (dis <= ballsize) 
            {
                box = (*node)->box;
                for(i=0; i<dim; i++)
                {
                    if (i != cut_dim)
                    {
                        dis = dis + dis2_from_bnd(qv[i],box[i].lower,box[i].upper);
                        if (dis > ballsize)
                        {   
                            return;
                        }
                    }
                }

            search(&nfarther);
            }
        }
    }
}


double dis2_from_bnd(double x,double amin,double amax)
{
    double res;
    if (x > amax)
    {
//       res = pow(x-amax,2.0);
        res = (x-amax)*(x-amax);
        return res;
    }
    else
    {   
        if (x < amin)
        { 
//            res = pow(amin-x,2.0);
            res = (amin-x)* (amin-x);
            return res;
        }
        else
        {
          res = 0.0;
          return res;
        }
    }
    return res;
}




void process_terminal_node_fixedball(struct tree_node **node)
{   
    double *qv;
    int *ind;
    double *data;
    int nfound;
    int i, indexofi, k;
    double ballsize, sd;
    int rearrange;

    qv = sr->qv;
    ballsize = sr->ballsize; 
//    printf("ballsize=%lf\n", ballsize);
    rearrange = sr->rearrange;
    ind = sr->ind;
    data = sr->data;
    nfound = sr->nfound;

    i=(*node)->l;
    while(i<=(*node)->u)
    {
        if (rearrange)
        {
            sd = 0.0;
            for(k=0; k<dim; k++)
            {  
//                sd = sd + pow((data[k+dim*i] - qv[k]),2.0);
                sd = sd + (data[k+dim*i] - qv[k])*(data[k+dim*i] - qv[k]);
                if (sd>ballsize) 
                {    
                    goto mainloop;
                }
            }
//            printf("sd1=%lf\n", sd);

            indexofi = ind[i]; 
        }
        else
        {
            indexofi = ind[i];
            sd = 0.0;
            for(k=0; k<dim; k++)
            {
//                sd = sd + pow((data[k+dim*indexofi] - qv[k]),2.0);
                sd = sd + (data[k+dim*indexofi] - qv[k])*(data[k+dim*indexofi] - qv[k]);
                if (sd>ballsize)
                { 
                    goto mainloop;
                }
            }
//            printf("sd2=%lf\n", sd);
        }

        nfound = nfound+1;
//        printf("nfound=%d\n", nfound);
        if (nfound > sr->nalloc)
        {
            sr->overflow = TRUE;
        }
        else
        {
            sr->results[nfound-1].dis = sd;
            sr->results[nfound-1].idx = indexofi;
        }
	mainloop:i++;
    }
    sr->nfound = nfound;
}
void kdtree_sort_results(int nfound,struct kdtree_result *results)
{
    if (nfound >= 1) 
    {
    	qsort(results, nfound, sizeof(struct kdtree_result), cmp);
    }
    return;
}

int cmp(const void *a ,const void *b)
{
    return (((struct kdtree_result *)a)->dis > ((struct kdtree_result *)b)->dis) ? 1 : -1;
}


void destroy_node(struct tree_node **np)
{
    if ((*np)->left!=NULL)
    {
        destroy_node(&((*np)->left));
        (*np)->left=NULL;
    }
    if ((*np)->right!=NULL)
    {
        destroy_node(&((*np)->right));
        (*np)->right=NULL;
    }
    if ((*np)->box!=NULL) 
    {
        free((*np)->box);
    }
        free((*np));
}

void kdtree_destroy(struct kdtree **tp)
{
    destroy_node(&((*tp)->root));
    free ((*tp)->ind);
    (*tp)->ind=NULL;

    if ((*tp)->rearrange) 
    {
       free((*tp)->rearranged_data);
       (*tp)->rearranged_data=NULL;
    }
    free((*tp));
}