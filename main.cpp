
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>

#define root 0

// return 1 if in set, 0 otherwise
int inset(double real, double img, int maxiter){
    double z_real = real;
    double z_img = img;

    for(int iters = 0; iters < maxiter; iters++){
        double z2_real = z_real*z_real-z_img*z_img;
        double z2_img = 2.0*z_real*z_img;
        z_real = z2_real + real;
        z_img = z2_img + img;
        if(z_real*z_real + z_img*z_img > 4.0) return 0;
    }
    return 1;
}

// count the number of points in the set, within the region
void mandelbrotSetCount(double real_lower, double real_upper, double img_lower, double img_upper, int num, int maxiter){

    int nnodes, //number of MPI nodes in the computation
            rank, //node number
            chunk; //number of complex handled by each node

    double start_r, end_r; //start,end real region for this node

    int buffer_p[2];
    buffer_p[0]=num;
    buffer_p[1]=maxiter;

    MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Barrier(MPI_COMM_WORLD);

    printf("nproc:%d,rank:%d\n",nnodes,rank);

// Distribute the data

//    MPI_Bcast(&maxiter,1,MPI_INT,root,MPI_COMM_WORLD);
//    MPI_Bcast(&num,1,MPI_INT,root,MPI_COMM_WORLD);

    chunk=num/nnodes;

    double real_len=real_upper-real_lower;
    start_r=real_lower+rank*real_len/nnodes;
    end_r=start_r+(rank+1)*real_len/nnodes;

    int global_count=0,local_count=0;

    double real_step = (real_upper-real_lower)/num;
    double img_step = (img_upper-img_lower)/num;

//    int pNum=omp_get_num_procs();
//
//    printf("Processor number is %d\n",pNum);

//#pragma omp parallel for schedule(dynamic)
    for(int real=0; real<chunk; real++){
        for(int img=0; img<num; img++){
            local_count+=inset(start_r+real*real_step,img_lower+img*img_step,maxiter);
        }
    }

//#pragma omp barrier

    printf("No.%d node counted %d numbers\n", rank,local_count);

    MPI_Reduce(&local_count,&global_count,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

    if(rank==root)
        printf("%d\n", global_count);
}

// main
int main(int argc, char *argv[]){


    double real_lower;
    double real_upper;
    double img_lower;
    double img_upper;
    int num;
    int maxiter;
    // each six input arguments is a region
    int num_regions = (argc-1)/6;

    double stamp;

    MPI_Init(&argc,&argv);

    stamp=MPI_Wtime();

//#pragma omp parallel for
    for(int region=0;region<num_regions;region++){
        // scan the arguments
        sscanf(argv[region*6+1],"%lf",&real_lower);
        sscanf(argv[region*6+2],"%lf",&real_upper);
        sscanf(argv[region*6+3],"%lf",&img_lower);
        sscanf(argv[region*6+4],"%lf",&img_upper);
        sscanf(argv[region*6+5],"%i",&num);
        sscanf(argv[region*6+6],"%i",&maxiter);
        mandelbrotSetCount(real_lower,real_upper,img_lower,img_upper,num,maxiter);
    }

    printf("time_cost:%.16g\n",MPI_Wtime()-stamp);

    MPI_Finalize ();
    
    return EXIT_SUCCESS;
}