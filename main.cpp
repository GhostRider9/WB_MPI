
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>
#include <math.h>
#include "edge_tracking.cpp"

#define root 0
#define slice 200
#define tag_do_work 1
#define tag_done_work 2
#define tag_shutdown 3


#define PERIODICITY_CHECKING_ENABLED 0
#define CARDIOID_BULB_CHECKING 0
#define IMG_AXIS_SYMMETRY 1

#define EPSINON 0.00001

static int state=0;

int inset(double real, int maxiter){
    double z_real = real;
    double z_img = 0;
    for(int iters = 0; iters < maxiter; iters++){
        double z2_real = z_real*z_real-z_img*z_img;
        double z2_img = 2.0*z_real*z_img;
        z_real = z2_real + real;
        z_img = z2_img;
        if(z_real*z_real + z_img*z_img > 4.0) return 0;
    }
    return 1;
}

int realAxisCount(float real_lower,float real_upper,int num,int maxiter){
    int count=0;
    double real_step = (real_upper-real_lower)/num;
    for(int real=0; real<num; real++){
        count+=inset(real_lower+real*real_step,maxiter);
    }
    return count;
}


void masterDoWork(double real_lower, double real_upper,int num){

//#pragma omp parallel
//    {
//        printf("Master threads: %d\n",omp_get_num_threads());
//    }

    int nprocess;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocess);

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
        //printf("Master received process %d asks for work\n",status.MPI_SOURCE);

        if(allDone){
            //printf("Send no work signal to process %d\n", status.MPI_SOURCE);
            MPI_Ssend(NULL,0,MPI_INT,status.MPI_SOURCE,tag_shutdown,MPI_COMM_WORLD);
            count++;
            if(count==nprocess){
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


int slaveDoWork(double real_lower, double real_upper, double img_lower, double img_upper, int img_num, int maxiter){

//#pragma omp parallel
//    {
//        printf("Slave threads: %d\n",omp_get_num_threads());
//    }


    double buffer[2];
    double start_r, end_r;

    int local_count=0;
    bool noWork= false;


    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Status status;
    MPI_Ssend(NULL,0,MPI_INT,root,tag_done_work,MPI_COMM_WORLD);


    EdgeTracker* et = new EdgeTracker();

    while(!noWork){
        MPI_Probe (root, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if(status.MPI_TAG==tag_do_work){
            MPI_Recv(buffer,2,MPI_DOUBLE,root,tag_do_work,MPI_COMM_WORLD,&status);

//            printf("Process %d receive data from root, data %lf,%lf\n", rank,buffer[0],buffer[1]);

            start_r=buffer[0];
            end_r=buffer[1];

//            double real_temp=start_r;

//            for(int real=1; real_temp<end_r; real++){
////                printf("%lf",real_temp);
//                for(int img=0; img<num; img++){
//                    local_count+=inset(real_temp,img_lower+img*img_step,maxiter);
//                    real_temp=start_r+real*real_step;
//                }
//            }
//            printf("%lf ",start_r);
            local_count += et->pointsInRegion(start_r,end_r,img_lower,img_upper,maxiter,slice,img_num);
//            printf("No.%d process counted %d numbers\n", rank,local_count);

            MPI_Ssend(NULL,0,MPI_INT,root,tag_done_work,MPI_COMM_WORLD);
            //printf("Process %d asks for work\n",rank);

        } else if(status.MPI_TAG==tag_shutdown){
            MPI_Recv(NULL,0,MPI_INT,root,tag_shutdown,MPI_COMM_WORLD,&status);
            //printf("No work need to do\n");
            noWork=true;
        }
    }
    return local_count;
}



// count the number of points in the set, within the region
int mandelbrotSetCount(double real_lower, double real_upper, double img_lower, double img_upper, int real_num, int img_num, int maxiter){

    int nprocess, //number of MPI process in the computation
            rank; //process number

    MPI_Comm_size(MPI_COMM_WORLD,&nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Barrier(MPI_COMM_WORLD);

    printf("nprocess:%d,rank:%d\n",nprocess,rank);

    int global_count=0,local_count=0;

    if(rank==root){
        masterDoWork(real_lower,real_upper,real_num);
    }
    else{
        //printf("slave process %d do work\n",rank);
        local_count=slaveDoWork(real_lower,real_upper,img_lower,img_upper,img_num,maxiter);
        //printf("No.%d process counted %d numbers\n", rank,local_count);
    }

//    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&local_count,&global_count,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

    return global_count;
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



//#pragma omp parallel for
    for(int region=0;region<num_regions;region++){

        int rank;

        // enable multi-thread support
        MPI_Init_thread (& argc , &argv , MPI_THREAD_FUNNELED ,& provided);

//    MPI_Init(&argc,&argv);

        stamp=MPI_Wtime();

        int global_count=0;
        // scan the arguments
        sscanf(argv[region*6+1],"%lf",&real_lower);
        sscanf(argv[region*6+2],"%lf",&real_upper);
        sscanf(argv[region*6+3],"%lf",&img_lower);
        sscanf(argv[region*6+4],"%lf",&img_upper);
        sscanf(argv[region*6+5],"%i",&num);
        sscanf(argv[region*6+6],"%i",&maxiter);

        double img_length = img_upper-img_lower;
        double real_length = real_upper-real_lower;
        double img_step=img_length/num;
        int temp_num1, temp_num2;

        if(IMG_AXIS_SYMMETRY && fmod(num,img_length)<EPSINON){
            if(img_upper>0 && img_lower<0){
                if(abs(img_lower+img_upper)<EPSINON){
                    state=1;
                }
                else {
                    if(abs(img_upper)<abs(img_lower)){
                        double temp=img_lower;
                        img_lower=-img_upper;
                        img_upper=-temp;
                    }
                    state=2;
                }
            }
        }


//        printf("State: %d\n",state);
        switch(state){
            case 0:
                global_count=mandelbrotSetCount(real_lower,real_upper,img_lower,img_upper,num,num,maxiter);
                break;
            case 1:
                temp_num1=num/2;
                global_count= -realAxisCount(real_lower,real_upper,num,maxiter)
                                  +2*mandelbrotSetCount(real_lower,real_upper,0,img_upper,num,temp_num1,maxiter);
                break;
            case 2:
                temp_num1=(int)(-img_lower*num/img_length);
                temp_num2=(int)((img_upper+img_lower)*num/img_length);

                printf("%d,%d\n",temp_num1,temp_num2);
                printf("%lf,%lf\n",(-img_lower*num/img_length),((img_upper+img_lower)*num/img_length));

                global_count= -realAxisCount(real_lower,real_upper,num,maxiter)
                                  +2*mandelbrotSetCount(real_lower,real_upper,0,-img_lower,num,temp_num1,maxiter)
                                  +mandelbrotSetCount(real_lower,real_upper,-img_lower,img_upper,num,temp_num2,maxiter);
                break;
        }

        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

        if(rank==root){
            printf("%d\n",global_count);
            printf("time_cost:%.16g\n",MPI_Wtime()-stamp);
        }

        MPI_Finalize ();

    }



    return EXIT_SUCCESS;
}