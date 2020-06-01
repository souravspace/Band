#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include "CZTIMM_SRCDEF_POL_Band_new.cc"
#include <iostream>

int main (int argc, char **argv){

        int i, j, pix, rank, size,num;
        double tx_run,ty_run,detpol_run,alpha_run,beta_run,epeak_run;

 char builddir[100]="/mnt/home/project/ccztipoc/CZTIMM/code/mass_model/CZTIMassModel-build_pol_band_test";
 char outdir[200]="/mnt/home/project/ccztipoc/gitdown/Test/Results";      


 int N_ANG = 1;                          
 int P_ANG = 1;                          
 int N_GRB = 25;                         
 int ngrb_run;

	double thetax[1]={167.29204};    
 	double thetay[1]={159.22561};      
	double detpol[1]={0};            

	double alpha[1]={-0.98};        
        double beta[1]={-2.67};
	double epeak[1]={1014.9};


	int satIllumFlag=1;      
        int polFlag=0;                 
        int verbose=0;

        int e_min = 1;
        int e_max = 10000;
        long nevt=4.0E6;                


long NRUNS=N_ANG*N_GRB*P_ANG;

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
                        detpol_run=detpol[(num % P_ANG)];

						ngrb_run=(int)(num/1);


						alpha_run=alpha[(num % N_ANG)];
						beta_run=beta[(num % N_ANG)];
						epeak_run=epeak[(num % N_ANG)];
						

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
                        gpsDef.computeGeantSource(tx_run,ty_run,alpha_run,beta_run,epeak_run,e_min,e_max,nevt,satIllumFlag,polFlag,detpol_run);
                        gpsDef.writeMonoEneGpsFile(gps_fnam, verbose, nevt);
                     
                        sprintf(run_command,"%s; %s/%s/CZTIMM %s %d > %s",cd_command,outdir,RUN_STR,gps_fnam,rank,verbose_fnam);

                        cout<<run_command<<"\n";


//	                system(run_command);
						
                }

        }

        MPI_Finalize();
        return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////


