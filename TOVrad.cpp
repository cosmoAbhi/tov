#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <math.h>
#include <sstream>
using namespace std;

double EOS(double, int);//returns P as a function of rho
double EOSinv(double, int);//returns rho as a function of P
double EOSder(double, double, int); //returns the value of dp/d(rho)
double mDer(double, double);
double PDer(double, double, double, double);
double betaDer(double, double, double, double, double, double, int);
double rstart=0.001;
double rstop=14.000;
double rstep=0.001;
double r=0.001;
double K=101.4496;//7.29935;
double tau=2.0;//5.0/3.0;
double a0=1.0;
double expLambda,g;
//all units are in km
//arrays for storing values of dependent variables
double m[14000];
double rho[14000];
double P[14000];
double y[14000];
double H[14000];
double beta[14000];
double Ye_BE[14000];
//arrays for storing intermediate values required to compute dependents by RK4
double mRK[4];
double PRK[4];
double betaRK[4];
double mass, press, dens, Beta, h;
//arrays to store tabulated EOS entries
int ser[1000];
double nB[1000];
double rhoEOS[1000];
double pEOS[1000];
int EOScount=0;//variable to keep a count on the no. of tabulated EOS entries
int choice;//variable to specify EOS choice
double frac;//variable for carrying out interpolation

int main()
{
    cout<<"Specify EOS choice: 1 for polytrope, 2 for EOS.txt: ";
    cin>>choice;
    if(choice==2)
    {
        ifstream infile;
        infile.open("EOS.txt");

        while (!infile.eof())
        {
            infile>> ser[EOScount];
            infile>> nB[EOScount];
            infile>> rhoEOS[EOScount];
            infile>> pEOS[EOScount];
            rhoEOS[EOScount]*=7.4237e-19;
            pEOS[EOScount++]*=8.2601e-40;
        }
        infile.close();
    }
    ofstream write;
    write.open("TOVrad_output.txt");

    rho[0]=0.0012355;//0.000741111;//km^(-2) central density
    P[0]=EOS(rho[0], choice);
    m[0]=4/3*M_PI*pow(rstart, 3.0)*rho[0];
    H[0]=a0*pow(rstart, 2.0);//a0 is arbitrarily chosen
    beta[0]=2*a0*rstart;
    y[0]=r*beta[0]/H[0];
    g=-2.0*PDer(rstart, m[0], P[0], rho[0])/(P[0]+rho[0]);
    Ye_BE[0]=0.1*(1-exp(-rho[0]/0.0006))+0.1*exp(-rho[0]/0.00002);
    write<<"    r       m           P           rho         H                   y           Ye_BE       g"<<endl;
    write<<r<<" "<< m[0]<<" "<<P[0]<<" "<<rho[0]<<" "<<H[0]<<" "<<y[0]<<" "<<Ye_BE[0]<<" "<<g<<endl;
    for(int i=1;i<14000;i++)
    {
        if(r==rstop)
            break;
        //RK4 step 1: calculating 'k1's and m, P, rho, beta, H values to be used in step 2 (refer to http://mathworld.wolfram.com/Runge-KuttaMethod.html)
        mRK[0]=rstep*mDer(r, rho[i-1]);
        PRK[0]=rstep*PDer(r, m[i-1], P[i-1], rho[i-1]);
        betaRK[0]=rstep*betaDer(r, m[i-1], P[i-1], rho[i-1], beta[i-1], H[i-1], choice);
        mass=m[i-1]+mRK[0]/2.0;
        press=P[i-1]+PRK[0]/2.0;
        Beta=beta[i-1]+betaRK[0]/2.0;
        dens=EOSinv(press, choice);
        h=H[i-1]+Beta*rstep;
        //RK4 step 2: calculating 'k2's and values for step 3
        mRK[1]=rstep*mDer(r+rstep/2.0, dens);
        PRK[1]=rstep*PDer(r+rstep/2.0, mass, press, dens);
        betaRK[1]=rstep*betaDer(r+rstep/2.0, mass, press, dens, Beta, h, choice);
        mass=m[i-1]+mRK[1]/2.0;
        press=P[i-1]+PRK[1]/2.0;
        Beta=beta[i-1]+betaRK[1]/2.0;
        dens=EOSinv(press, choice);
        h=H[i-1]+Beta*rstep;
        //RK4 step 3: calculating 'k3's and values for step 4
        mRK[2]=rstep*mDer(r+rstep/2.0, dens);
        PRK[2]=rstep*PDer(r+rstep/2.0, mass, press, dens);
        betaRK[2]=rstep*betaDer(r+rstep/2.0, mass, press, dens, Beta, h, choice);
        mass=m[i-1]+mRK[2];
        press=P[i-1]+PRK[2];
        Beta=beta[i-1]+betaRK[2];
        dens=EOSinv(press, choice);
        h=H[i-1]+Beta*rstep;
        //RK4 step 4: calculating 'k4's
        mRK[3]=rstep*mDer(r+rstep, dens);
        PRK[3]=rstep*PDer(r+rstep, mass, press, dens);
        betaRK[3]=rstep*betaDer(r+rstep, mass, press, dens, Beta, h, choice);
        //calculating RK4 final results for this step
        m[i]=m[i-1]+mRK[0]/6.0+mRK[1]/3.0+mRK[2]/3.0+mRK[3]/6.0;
        P[i]=P[i-1]+PRK[0]/6.0+PRK[1]/3.0+PRK[2]/3.0+PRK[3]/6.0;
        if(P[i]<=0)
            break;
        beta[i]=beta[i-1]+betaRK[0]/6.0+betaRK[1]/3.0+betaRK[2]/3.0+betaRK[3]/6.0;
        rho[i]=EOSinv(P[i], choice);
        H[i]=H[i-1]+beta[i]*rstep;
        r=r+rstep;
        y[i]=r*beta[i]/H[i];
        g=-2.0*PDer(r,m[i],P[i],rho[i])/(P[i]+rho[i]);
        if(rho[i]>=1e-4)
            Ye_BE[i]=0.1*(1-exp(-rho[i]/0.0006));
        else if(rho[i]>=1e-6&&rho[i]<1e-4)
            Ye_BE[i]=0.1*(1-exp(-rho[i]/0.0006))+0.1*exp(-rho[i]/0.00002);

        if(isnan(m[i])==false)
            write<<r<<" "<< m[i]<<" "<<P[i]<<" "<<rho[i]<<" "<<H[i]<<" "<<y[i]<<" "<<Ye_BE[i]<<" "<<g<<endl;
    }
    write.close();
}
double EOS(double dens, int choice)
{
    switch(choice)
    {
        case 1: return K*pow(dens, tau);
        case 2: for(int i=1;i<EOScount;i++)
                {
                    if(rhoEOS[i-1]==dens)
                        return pEOS[i-1];
                    if(dens>rhoEOS[i-1]&&dens<rhoEOS[i])
                    {
                        frac=(dens-rhoEOS[i-1])/(rhoEOS[i]-rhoEOS[i-1]);
                        return (frac*(pEOS[i]-pEOS[i-1])+pEOS[i-1]);
                    }

                }
                if(rhoEOS[EOScount-1]==dens)
                    return pEOS[EOScount-1];
    }
    return 0.0;
}
double EOSinv(double pres, int choice)
{
    switch(choice)
    {
        case 1: return pow(pres/K, 1/tau);
        case 2: for(int i=1;i<EOScount;i++)
                {
                    if(pEOS[i-1]==pres)
                        return rhoEOS[i-1];
                    if(pres>pEOS[i-1]&&pres<pEOS[i])
                    {
                        frac=(pres-pEOS[i-1])/(pEOS[i]-pEOS[i-1]);
                        //cout<<frac<<endl;
                        //cout<<(frac*(rhoEOS[i]-rhoEOS[i-1])+rhoEOS[i-1])<<endl;
                        return (frac*(rhoEOS[i]-rhoEOS[i-1])+rhoEOS[i-1]);
                    }

                }
                if(pEOS[EOScount-1]==pres)
                    return rhoEOS[EOScount-1];
    }
    return 0.0;
}
double EOSder(double P, double rho, int choice)
{
    switch(choice)
    {
        case 1: return tau*P/rho;
        case 2: for(int i=1;i<EOScount;i++)
                {
                    if(pEOS[i-1]==P)
                        return (pEOS[i]-pEOS[i-1])/(rhoEOS[i]-rhoEOS[i-1]);
                    if(P>pEOS[i-1]&&P<pEOS[i])
                    {
                        frac=(P-pEOS[i-1])/(pEOS[i]-pEOS[i-1]);
                        if(frac<0.5)
                            return (P-pEOS[i-1])/(rho-rhoEOS[i-1]);
                        else
                            return (pEOS[i]-P)/(rhoEOS[i]-rho);
                    }

                }
                if(pEOS[EOScount-1]==P)
                    return (pEOS[EOScount-1]-pEOS[EOScount-2])/(rhoEOS[EOScount-1]-rhoEOS[EOScount-2]);;
    }
    return 0.0;
}
double mDer(double r, double rho)
{
    return 4.0*M_PI*pow(r, 2.0)*rho;
}
double PDer(double r, double m, double P, double rho)
{
    return -1.0*(rho+P)*(m+4.0*M_PI*pow(r, 3.0)*P)/(pow(r, 2.0)-2.0*m*r);
}
double betaDer(double r, double m, double P, double rho, double beta, double H, int choice)
{
    expLambda=r/(r-2*m);
    return 2*expLambda*H*(-2*M_PI*(5*rho+9*P+(rho+P)/EOSder(P, rho, choice))+3/pow(r, 2.0)+2*expLambda*pow((m/pow(r, 2.0)+4*M_PI*r*P),2.0))+2*beta/r*expLambda*(-1+m/r+2*M_PI*pow(r, 2.0)*(rho-P));
}
