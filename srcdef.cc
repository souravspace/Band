#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;



class srcDefn
{
        public:

                double thetax,thetay;
                int pol;
                int satIllum;
                double detpolang;


                double vx,vy,vz;
                double theta,phi;
                double polx,poly,polz;

                double G_posx,G_posy,G_posz;
                double G_vx,G_vy,G_vz;
                double G_polx,G_poly,G_polz;


                double G_Rot1[3];
                double G_Rot2[3];

                double G_halfx;
                double G_halfy;

                void computeGeantSource(double tx,double ty, int satIllumFlag, int polFlag, double detpol);

                void writeMonoEneGpsFile(char *filename, int verbose, long nevt, double alpha_run, double beta_run, double epeak_run, double G4new_px_run, double G4new_py_run, double G4new_pz_run,double G4new_vx_run,double G4new_vy_run,double G4new_vz_run);
};

void srcDefn::computeGeantSource(double tx,double ty, int satIllumFlag, int polFlag, double detpol_run)
{

        thetax=tx;
        thetay=ty;
        satIllum= satIllumFlag;
        pol=polFlag;
        detpolang=detpol_run;
        double PI=3.1415;

        double R=100.0;
        double th_x=thetax*(PI/180.);
        double th_y=thetay*(PI/180.);

    vz=1./sqrt(1.+tan(th_x)*tan(th_x)+tan(th_y)*tan(th_y));
    vx=vz*tan(th_x);
    vy=vz*tan(th_y);

    vx=fabs(vx);
    vy=fabs(vy);

    if (thetax < 0.0) vx=fabs(vx)*(-1.0);

    if (thetay < 0.0) vy=fabs(vy)*(-1.0);

    vz=vx/tan(th_x);
    vz=fabs(vz);

    if (fabs(thetax) > 90.0 || fabs(thetay) > 90.0) vz=fabs(vz)*(-1.0);


    G_vx=-vy;
    G_vy=-vz;
    G_vz=-vx;



    G_Rot1[0]=vx/sqrt(1-vz*vz);
    G_Rot1[1]=0;
    G_Rot1[2]=-vy/sqrt(1-vz*vz);

    G_Rot2[0]=vz*vy/sqrt(1-vz*vz);
    G_Rot2[1]=-sqrt(vx*vx+vy*vy);
    G_Rot2[2]=vz*vx/sqrt(1-vz*vz);

        if(satIllum==0)
        {

        G_posx=-G_vx*R;
        G_posy=-G_vy*R;
        G_posz=-G_vz*R;
                G_halfx=50.0;
                G_halfy=50.0;
        }
        else
        {

        G_posx=-G_vx*R;
        G_posy=-G_vy*R;
        G_posz=-G_vz*R;
        G_halfx=100.0;
        G_halfy=100.0;
        }




        if(pol==1 && detpol_run!=0)
        {

    double chi=detpol_run*(PI/180.);


    double rb=1.0/(vx*cos(chi)+vy*sin(chi));

    double xb=rb*cos(chi);
    double yb=rb*sin(chi);
    double zb=0.0;




    polx=xb-vx;
    poly=yb-vy;
    polz=zb-vz;

    double p=sqrt(polx*polx+poly*poly+polz*polz);



    polx/=p;
    poly/=p;
    polz/=p;


    G_polx=poly;
    G_poly=polz;
    G_polz=polx; 
    
    G_polz=-(G_polx*G_vx+G_poly*G_vy)/G_vz;

        }
        else
        {
                G_polx=0;
                G_poly=0;
                G_polz=0;
                G_Rot1[0]=0;
                G_Rot1[1]=0;
                G_Rot1[2]=1;
                G_Rot2[0]=1;
                G_Rot2[1]=0;
                G_Rot2[2]=0;


        }

        return;

}

void srcDefn::writeMonoEneGpsFile(char *filename, int verbose, long nevt, double alpha_run,double beta_run, double epeak_run, double G4new_px_run, double G4new_py_run, double G4new_pz_run,double G4new_vx_run,double G4new_vy_run,double G4new_vz_run)
{

        FILE *fgps=fopen(filename,"w");

fprintf(fgps,"/gps/particle gamma \n\
/gps/spectrum/type Band\n\
/gps/spectrum/emin 0.100 MeV\n\
/gps/spectrum/emax 0.400 MeV\n\
/gps/spectrum/Norm 1\n\
/gps/spectrum/alpha_indx %lf \n\
/gps/spectrum/beta_indx %lf \n\
/gps/spectrum/Epeak %lf \n\
/gps/pos/type Plane \n\
/gps/pos/shape Square \n\
/gps/pos/centre %lf %lf %lf cm\n\
/gps/direction %lf %lf %lf \n\
/gps/pos/rot1 %lf %lf %lf \n\
/gps/pos/rot2 %lf %lf %lf \n\
/gps/polarization %lf %lf %lf \n\
/gps/pos/halfx %lf cm\n\
/gps/pos/halfy %lf cm\n\
/tracking/verbose %d\n\
/gps/verbose %d\n\
/run/beamOn %ld\n",
alpha_run,\
beta_run,\
epeak_run,\
G_posx,G_posy,G_posz,\
G4new_vx_run,G4new_vy_run,G4new_vz_run,\
G_Rot1[0],G_Rot1[1],G_Rot1[2],\
G_Rot2[0],G_Rot2[1],G_Rot2[2],\
G4new_px_run,G4new_py_run,G4new_pz_run,\
G_halfx,G_halfy,\
verbose,verbose,nevt);

        fclose(fgps);

        return;
}



