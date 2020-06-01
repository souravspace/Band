#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

/* Function to computes source definition for Geant4 simulations 
 * from thetax and thetay. 
 *
 * Function to produces GPS macro file from this definition 
 *
 *
 * Usage: Create instance of srcDefn class in your code and 
 * use the methods of class to create macro file. For details 
 * see the class definition below. (Sample main code at the end)
 *
 * Based on grbpospol.pro  
 *
 * Mithun NPS (22/12/2016)
 *
 *
 * Edited on 06/09/18 by Aarthy to input new px, py, pz geant given by Mithun from Dipankar's 
 * method of conversion of det to sky angle. Since vectors are directly available the G_polxyz 
 * calculation is not used here
 *
 *
 * Edited on 03/04/19 by Aarthy. Correct calculation of Px Py Pz, no longer given as input directly
 * Formalism in grbpospol by Mithun implemented here
 *
 * */


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
		void writeMonoEneGpsFile(char *filename, int verbose, long nevt, double alpha_run, double beta_run, double epeak_run);
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
	
    //CHANGE FROM CZTI xyz to Geant4 xyz 
    G_vx=-vy;
    G_vy=-vz;
    G_vz=-vx;

	//printf("HERE %lf\t%lf\t%lf\n",vx,vy,vz);

	// Plane rotations (As given by Sujay Mate)
 
	G_Rot1[0]=vx/sqrt(1-vz*vz);  //cos(phi) 
    G_Rot1[1]=0;                 // 0 
    G_Rot1[2]=-vy/sqrt(1-vz*vz); // -sin(phi)
	
    G_Rot2[0]=vz*vy/sqrt(1-vz*vz); // cos(theta)*sin(phi)
    G_Rot2[1]=-sqrt(vx*vx+vy*vy);  // -sin(theta)
    G_Rot2[2]=vz*vx/sqrt(1-vz*vz); // cos(theta)*cos(phi)

	if(satIllum==0) 
	{
		//Postion and Length for only CZTI Illumination
    	G_posx=-G_vx*R;
    	G_posy=-G_vy*R;
    	G_posz=-G_vz*R;
		G_halfx=50.0;
		G_halfy=50.0;
	}
	else
	{
		//Postion and Length for Full Satellite Illumination
        G_posx=-G_vx*R;
        G_posy=-G_vy*R;
        G_posz=-G_vz*R;	
        G_halfx=100.0;
        G_halfy=100.0;

		/* *************THIS IS INCORRECT (KEPT AS PLACE HOLDER) : NEEDS TO BE DONE ******************** */		
	}
		

    //Computing polarization vector direction (if required)

	if(pol==1 && detpol_run!=0)
	{
    //Point of intersection of projection of polarization vector in CZTI plane (given by cos(chi) i + sin(chi) j)
    //and the source plane Point (xb,yb,zb)

    double chi=detpol_run*(PI/180.);


    double rb=1.0/(vx*cos(chi)+vy*sin(chi));

    double xb=rb*cos(chi);
    double yb=rb*sin(chi);
    double zb=0.0;

    //Polarization vector - direction numbers (Line passing through P(xb,yb,zb) and (vx,vy,vz))

//    polx=xb-vx;
//   poly=yb-vy;
//    polz=zb-vz;
	
	polx=xb;
	poly=yb;
	polz=zb-(1./vz);

    double p=sqrt(polx*polx+poly*poly+polz*polz);

    // Normalize to get direction cosines

    polx/=p;
    poly/=p;
    polz/=p;

    //Change to Geant4 coordinate system

    G_polx=poly;
    G_poly=polz;
    G_polz=polx;

    // To take care of round of errors
//    G_polz=-(G_polx*G_vx+G_poly*G_vy)/G_vz;

	double pz_Geanta=-(G_polx*G_vx+G_poly*G_vy)/G_vz;
    if(abs(G_polz-pz_Geanta) > 0.0001)
	{
		printf("SOMETHING IS WRONG **** %f\n",G_polz-pz_Geanta);
		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",detpol_run,vx,vy,vz,polx,poly,polz);
    }
    else
	    printf("ALL GOOD %f\n", G_polz-pz_Geanta);
	    G_polz=pz_Geanta;
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

void srcDefn::writeMonoEneGpsFile(char *filename, int verbose, long nevt, double alpha_run, double beta_run, double epeak_run)
{

	FILE *fgps=fopen(filename,"w");

fprintf(fgps,"/gps/particle gamma \n\
/gps/spectrum/type Band\n\
/gps/spectrum/emin 0.100 MeV\n\
/gps/spectrum/emax 2.00 MeV\n\
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
G_vx,G_vy,G_vz,\
G_Rot1[0],G_Rot1[1],G_Rot1[2],\
G_Rot2[0],G_Rot2[1],G_Rot2[2],\
G_polx,G_poly,G_polz,\
G_halfx,G_halfy,\
verbose,verbose,nevt);

	fclose(fgps);

	return;
}
