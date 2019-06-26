
void       read_image(FILE *,DATA_TYPE *,int,int,int,int,double,double);
char       max(int ,int ,int ,DATA_TYPE * ,DATA_TYPE *);
void       get_stamps(DATA_TYPE *, char *map);
void       get_stamp_vectors(DATA_TYPE *,DATA_TYPE *);
double     power(double ,int ), **build_matrix();
void       xy_conv_stamp(DATA_TYPE *,int,int,char);
void       allocate_vectors(), make_bg_vectors(int, int);
void       build_scprod(DATA_TYPE *);
void       lubksb(double **,int ,int *,double *);
void       ludcmp(double **,int ,int *,double *);
DATA_TYPE  *kernel_vector(int,int,int,int,char *);
void       make_kernel(int, int, double *);
void       kernel_convolve(DATA_TYPE *, double *);
double     get_background(int, int, double *);
void       cut_stamp(int , int ,DATA_TYPE *,DATA_TYPE *,double *);
char       *defect_map(DATA_TYPE *,DATA_TYPE *);
void       check_stamps();
void       build_matrix0();
void       build_scprod0();
char       check_again(double *,DATA_TYPE *);
void       make_model(int, double *);
void       write_header(FILE *);
void       write_image(DATA_TYPE *, FILE *,int,int);
void       write_image2(DATA_TYPE *, FILE *,int,int);
char       read_config(char *);
FILE       *read_header(char *,int *,double *,double *);
void       quick_sort(DATA_TYPE *,int *,int n);
int        get_index(char *,char *);

