#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

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

                double spec[37];
                double ener[37]= {1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};
                double norm;
                double photons ;
                double phot_tot;  

                void computeGeantSource(double tx,double ty,double alpha_run, double beta_run, double epeak_run, int e_min_run, int e_max_run, long nevt_run, int satIllumFlag, int polFlag, double detpol);
                void writeMonoEneGpsFile(char *filename, int verbose, long nevt);


};



void srcDefn::computeGeantSource(double tx,double ty,double alpha_run, double beta_run, double epeak_run,int e_min_run, int e_max_run,long nevt_run, int satIllumFlag, int polFlag, double detpol_run)
{
        thetax=tx;
        thetay=ty;
        satIllum= satIllumFlag;
        pol=polFlag;
        detpolang=detpol_run;
        double PI=3.1415;
        double enval;

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

/*
////To check band normalization and input spectrum
std::ofstream ofiles;
ofiles.open("input_band_spec.dat", std::ios::app);
phot_tot = 0;
for (int ie = e_min_run ; ie <=e_max_run ; ie++)
{
enval = ie*1.;
if ((alpha_run - beta_run)*epeak_run >= enval)
{
photons = pow ((enval/100.), alpha_run)*exp (-enval/epeak_run) ;
}
else
{
photons = pow (((alpha_run - beta_run) * epeak_run/100.), (alpha_run - beta_run))* exp (beta_run - alpha_run) * pow ((enval/100.), beta_run); 
}
ofiles << ie << "\t\t"<< photons<<endl;
phot_tot = phot_tot + photons ;
}
ofiles.close();
norm = (nevt_run*1.)/phot_tot;
cout <<"Total photons is " <<phot_tot<<endl;
*/


double enn,phot;

for (int i = 0 ; i <37 ; i++)
{
enn = ener[i];
if ((alpha_run - beta_run)*epeak_run >= enn)
{
phot = pow ((enn/100.), alpha_run)*exp (-enn/epeak_run) ;
}
else
{
phot = pow (((alpha_run - beta_run) * epeak_run/100.), (alpha_run - beta_run))* exp (beta_run - alpha_run) * pow ((enn/100.), beta_run); 
}
spec[i] = phot ;
}

return;
}

void srcDefn::writeMonoEneGpsFile(char *filename, int verbose, long nevt)
{
FILE *fgps=fopen(filename,"w");
fprintf(fgps,"/gps/particle gamma \n\
/gps/pos/type Plane \n\
/gps/pos/shape Square \n\
/gps/pos/centre %lf %lf %lf cm\n\
/gps/direction %lf %lf %lf \n\
/gps/pos/rot1 %lf %lf %lf \n\
/gps/pos/rot2 %lf %lf %lf \n\
/gps/polarization %lf %lf %lf \n\
/gps/pos/halfx %lf cm\n\
/gps/pos/halfy %lf cm\n\
/gps/ene/type Arb \n\
/gps/hist/type arb \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/point %lf %lf \n\
/gps/hist/inter Lin \n\
/tracking/verbose %d\n\
/gps/verbose %d\n\
/run/printProgress 10000 \n\
/run/beamOn %ld\n",
G_posx,G_posy,G_posz,\
G_vx,G_vy,G_vz,\
G_Rot1[0],G_Rot1[1],G_Rot1[2],\
G_Rot2[0],G_Rot2[1],G_Rot2[2],\
G_polx,G_poly,G_polz,\
G_halfx,G_halfy,\
1./1000.,spec[0],\
2./1000.,spec[1],\
3./1000.,spec[2],\
4./1000.,spec[3],\
5./1000.,spec[4],\
6./1000.,spec[5],\
7./1000.,spec[6],\
8./1000.,spec[7],\
9./1000.,spec[8],\
10./1000.,spec[9],\
20./1000.,spec[10],\
30./1000.,spec[11],\
40./1000.,spec[12],\
50./1000.,spec[13],\
60./1000.,spec[14],\
70./1000.,spec[15],\
80./1000.,spec[16],\
90./1000.,spec[17],\
100./1000.,spec[18],\
200./1000.,spec[19],\
300./1000.,spec[20],\
400./1000.,spec[21],\
500./1000.,spec[22],\
600./1000.,spec[23],\
700./1000.,spec[24],\
800./1000.,spec[25],\
900./1000.,spec[26],\
1000./1000.,spec[27],\
2000./1000.,spec[28],\
3000./1000.,spec[29],\
4000./1000.,spec[30],\
5000./1000.,spec[31],\
6000./1000.,spec[32],\
7000./1000.,spec[33],\
8000./1000.,spec[34],\
9000./1000.,spec[35],\
10000./1000.,spec[36],\
verbose,verbose,nevt);

fclose(fgps);

return;
}



