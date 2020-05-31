#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include "srcdef.cc"
#include <iostream>

int main (int argc, char **argv){

        int i, j, pix, rank, size,num;

        double tx_run,ty_run,detpol_run,alpha_run,beta_run,epeak_run,G4new_px_run,G4new_py_run,G4new_pz_run;
                double G4new_vx_run,G4new_vy_run,G4new_vz_run;


 char builddir[100]="/mnt/home/project/ccztipoc/astrosat_mass_model_pol/CZTIMassModel-build";
 char outdir[200]="/mnt/home/project/ccztipoc/astrosat_mass_model_pol/Results/Latest_Energy_100_400/pol_alpha/";





 int N_ANG = 1;
 int N_GRB = 100;
 int ngrb_run;
  double thetax[1]={167.29};


    double thetay[1]={159.22};


    double detpol[1]={0};


    double alpha[1]={-0.97};
    double beta[1]={-2.9257};
    double epeak[1]={800.};


        double G4new_px[1]={0};

        double G4new_py[1]={0};

        double G4new_pz[1]={0};



    double G4new_vx[1]={-0.347146};

        double G4new_vy[1]={0.914831};

        double G4new_vz[1]={-0.206334};

        int satIllumFlag=1;

        int polFlag=0;

   int verbose=0;
        long nevt=1E6;


   long NRUNS=N_ANG*N_GRB;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);


  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


        for(i=0;i<NRUNS;i++)
        {
                if((i%size)== 0)
                {
                        num = i + rank;


                        tx_run=thetax[(num % N_ANG)];
                        ty_run=thetay[(num % N_ANG)];
                        detpol_run=detpol[(num % N_ANG)];

                                                ngrb_run=(int)(num/N_ANG);
                                                cout<<"HAHAH "<<ngrb_run<<"\n";

                                                alpha_run=alpha[(num % N_ANG)];
                                                beta_run=beta[(num % N_ANG)];
                                                epeak_run=epeak[(num % N_ANG)];

                                                G4new_px_run=G4new_px[(num % N_ANG)];
                                                G4new_py_run=G4new_py[(num % N_ANG)];
                                                G4new_pz_run=G4new_pz[(num % N_ANG)];

                        G4new_vx_run=G4new_vx[(num % N_ANG)];
                        G4new_vy_run=G4new_vy[(num % N_ANG)];
                        G4new_vz_run=G4new_vz[(num % N_ANG)];

                        char RUN_STR[40];
                        char cp_command[1000];
                                        char cd_command[1000];
                        char rm_command[1000];

                                                char gps_fnam[1000];
                        char verbose_fnam[1000];
                        char run_command[2000];

                                                                                                                                                                                         
                        printf("%d\t%d\t%d\n",num,ngrb_run,num%N_ANG);



                        sprintf(RUN_STR,"N%03d_P%06.1lf_TX%06.1lf_TY%06.1lf",ngrb_run,detpol_run,tx_run,ty_run);

                        cout<<outdir<<RUN_STR<<"\n";


                                    sprintf(rm_command,"rm -r %s/%s",outdir,RUN_STR);
                        sprintf(cp_command,"cp -r %s %s/%s",builddir,outdir,RUN_STR);

                                sprintf(cd_command,"cd  %s/%s",outdir,RUN_STR);

                        system(rm_command);
                        system(cp_command);


                        sprintf(gps_fnam,"%s/%s/gps_%s.mac",outdir,RUN_STR,RUN_STR);
                        sprintf(verbose_fnam,"%s/%s/verbose.out",outdir,RUN_STR);

                                                cout<<gps_fnam<<"\n";

                        srcDefn gpsDef;
                        gpsDef.computeGeantSource(tx_run,ty_run, satIllumFlag, polFlag,detpol_run);
                        gpsDef.writeMonoEneGpsFile(gps_fnam, verbose, nevt,alpha_run,beta_run,epeak_run,G4new_px_run, G4new_py_run, G4new_pz_run,G4new_vx_run,G4new_vy_run,G4new_vz_run);

                        sprintf(run_command,"%s; %s/%s/CZTIMM %s %d > %s",cd_command,outdir,RUN_STR,gps_fnam,rank,verbose_fnam);

                        cout<<run_command<<"\n";


                                                system(run_command);

                }

        }

        MPI_Finalize();
        return 0;
}


