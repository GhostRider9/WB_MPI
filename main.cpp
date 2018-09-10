
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>

#define root 0
#define slice 20
#define tag_do_work 1
#define tag_done_work 2
#define tag_shutdown 4

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


void masterDoWork(double real_lower, double real_upper,int num){

    int nnodes;
    MPI_Comm_size(MPI_COMM_WORLD,&nnodes);

    double real_step=(real_upper-real_lower)/num;
    double gap=slice*real_step;

    double buffer[2];
    double temp_lower;
    double temp_upper=real_lower;

    bool allDone= false;

    printf("master distributes the work\n");

    MPI_Status status;
    int count=1;

    while(1){

        MPI_Recv(NULL,0,MPI_INT,MPI_ANY_SOURCE,tag_done_work,MPI_COMM_WORLD,&status);
        //printf("Master received node %d asks for work\n",status.MPI_SOURCE);

        if(allDone){
            //printf("Send no work signal to node %d\n", status.MPI_SOURCE);
            MPI_Ssend(NULL,0,MPI_INT,status.MPI_SOURCE,tag_shutdown,MPI_COMM_WORLD);
            count++;
            if(count==nnodes){
                break;
            }
        }else{
            temp_lower=temp_upper;
            temp_upper+=gap;
            buffer[0]=temp_lower;
            buffer[1]=temp_upper;
            if(temp_upper>=real_upper){
                temp_upper=real_upper;
                allDone= true;
            }
            MPI_Ssend(buffer,2,MPI_DOUBLE,status.MPI_SOURCE,tag_do_work,MPI_COMM_WORLD);
        }
    }

    //printf("Master shutdown\n");
}


int slaveDoWork(double real_lower, double real_upper, double img_lower, double img_upper, int num, int maxiter){
    double buffer[2];
    double start_r, end_r;
    double img_step = (img_upper-img_lower)/num;
    double real_step = (real_upper-real_lower)/num;
    int local_count=0;
    bool noWork= false;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Status status;
    MPI_Ssend(NULL,0,MPI_INT,root,tag_done_work,MPI_COMM_WORLD);

    while(!noWork){
        MPI_Probe (root, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if(status.MPI_TAG==tag_do_work){
            MPI_Recv(buffer,2,MPI_DOUBLE,root,tag_do_work,MPI_COMM_WORLD,&status);

            //printf("Node %d receive data from root, data %lf,%lf\n", rank,buffer[0],buffer[1]);

            start_r=buffer[0];
            end_r=buffer[1];

            double real_temp=start_r;

            for(int real=1; real_temp<end_r; real++){
//                printf("%lf",real_temp);
                for(int img=0; img<num; img++){
                    local_count+=inset(real_temp,img_lower+img*img_step,maxiter);
                    real_temp=start_r+real*real_step;
                }
            }
            //printf("No.%d node counted %d numbers\n", rank,local_count);

            MPI_Ssend(NULL,0,MPI_INT,root,tag_done_work,MPI_COMM_WORLD);
            //printf("Node %d asks for work\n",rank);

        } else if(status.MPI_TAG==tag_shutdown){
            MPI_Recv(NULL,0,MPI_INT,root,tag_shutdown,MPI_COMM_WORLD,&status);
            //printf("No work need to do\n");
            noWork=true;
        }
    }
    return local_count;
}

// count the number of points in the set, within the region
void mandelbrotSetCount(double real_lower, double real_upper, double img_lower, double img_upper, int num, int maxiter){

    int nnodes, //number of MPI nodes in the computation
            rank; //node number

    MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Barrier(MPI_COMM_WORLD);

    printf("nnodes:%d,rank:%d\n",nnodes,rank);

    int global_count=0,local_count=0;

    if(rank==root){
        masterDoWork(real_lower,real_upper,num);
    }
    else{
        //printf("slave node %d do work\n",rank);
        local_count=slaveDoWork(real_lower,real_upper,img_lower,img_upper,num,maxiter);
        //printf("No.%d node counted %d numbers\n", rank,local_count);
    }

//    MPI_Barrier(MPI_COMM_WORLD);
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
    int provided;

    int rank;

    // enable multi-thread support
    MPI_Init_thread (& argc , &argv , MPI_THREAD_FUNNELED ,& provided);

//    MPI_Init(&argc,&argv);

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

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(rank==root){
        printf("time_cost:%.16g\n",MPI_Wtime()-stamp);
    }

    MPI_Finalize ();

    return EXIT_SUCCESS;
}