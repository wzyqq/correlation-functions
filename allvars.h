#define bucket_size 12
//#define bucket_size 12
#define TRUE 1
#define FALSE 0
#define dim 3

struct interval
{
    double lower,upper;
};
struct tree_node
{
    int cut_dim;
    double cut_val;
    double cut_val_left,cut_val_right;
    int l,u;
    struct tree_node *left,*right;
    struct interval *box;
};
struct kdtree
{
    int n;
    double *the_data;
    int *ind;
    int sort;
    int rearrange;
    double *rearranged_data;
    struct tree_node *root;
};
struct kdtree_result
{
    double dis;
    int idx;
};


struct tree_search_record
{       
    int nfound;
    double ballsize;
    int nalloc;
    int rearrange; 
    int overflow;
    double *qv; 
    struct kdtree_result  *results; 
    double *data;
    int *ind;
}*sr;



struct kdtree *kdtree_create(double *input_data, int n, int sort, int rearrange);
void build_tree(struct kdtree **tp);
struct tree_node *build_tree_for_range(struct kdtree **tp,int l,int u,struct tree_node **parent);
void spread_in_coordinate(struct kdtree **tp,int c,int l,int u,struct interval *interv);
int select_on_coordinate_value(double *v,int *ind,int c,double alpha,int li,int ui);
double max(double a, double b);
double min(double a, double b);
int find_maxmum(double *array);
double sum(double *array, int *ind, int c, int location1, int location2);
void kdtree_r_nearest(struct kdtree **tp,double *qv,double r2,int *nfound,int nalloc,struct kdtree_result *results);
void search(struct tree_node **node);
double dis2_from_bnd(double x,double amin,double amax);
void process_terminal_node_fixedball(struct tree_node **node);
void kdtree_sort_results(int nfound,struct kdtree_result *results);
int cmp(const void *a ,const void *b);
void destroy_node(struct tree_node **np);
void kdtree_destroy(struct kdtree **tp);





 

