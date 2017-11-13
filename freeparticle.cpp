//  freeparticle.cpp
//  You must link this code with the fftw3 library. On Unix systems, link with -lfftw3 -lm.

#include "freeparticle.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <fftw3.h>
#include <complex>
#include <stdlib.h>
#include <time.h>

using namespace std;

#define numslices ((tf-t0)/(dt*div))
#define N 2048  // Finite difference resolution *number of points
#define dt 0.001  // Time step in program seconds
#define L 20.  //  Box length

#define RE 0  //Placeholder for real number designator
#define IM 1 //Placeholder for imaginary number designator

#define t0 0. // initial program time
#define tf 50. // program end time

#define div 100  // down-scale factor

#define AMPLITUDE 1.


// Physical dimensionality variables
double x;
double x0;
double p;
double k;

double t; //this variable stores the time

double psi[2][2][N]; //this stores the wavefunction
double chi[2][2][N]; //this stores the wavefunction in fourier space

double h = L/N;
double dx = h;

double sigma = 0.5;
int slicenum = 0; //intitalizes output field iterator
int num = 0;
double A;

double realsum;
double complexsum;
double probabilitysum;

const int numberslices = tf/(dt*div);
double norms[numberslices];

double alpha = 0.6; //This will be our measurement parameter

int guess[numberslices]; //this value will hold which of our functions we are using to make a measurement (measurement not always made every step but we randomly choose everytime to calculate our probability)

double probofmeasure[numberslices]; //stores probability of each measurement at weach step
double roll[numberslices];//this stores the value to be compared with our probability

int measurement[numberslices];

double leftfunction(double pos)
{
    if (pos < 0.0) {

        return alpha;
    }

    if (pos >= 0.0) {

        return 0.0;
    }
}

double rightfunction(double pos)
{
    if (pos < 0.0) {

        return 0.0;
    }

    if (pos >= 0.0) {

        return alpha;
    }
}

void updatewfrm()
{
    printf("\n commencing updating waveform\n");
    double updatedpsi[2][N];

    if (guess[num] == 0) {
        for (int index = 0; index < N; index++) {

            x = (index*dx) - (L/2);
            updatedpsi[RE][index] = psi[0][RE][index]*leftfunction(x);
            updatedpsi[IM][index] = psi[0][IM][index]*leftfunction(x);
        }
    }
    else
    {
        for (int index = 0; index < N; index++) {

            x = (index*dx) - (L/2);
            updatedpsi[RE][index] = psi[0][RE][index]*rightfunction(x);
            updatedpsi[IM][index] = psi[0][IM][index]*rightfunction(x);
        }

    }


    double psisum = 0.;
    double updatedpsisum = 0.;

    for (int index = 0; index < N; index++) {
        psisum += psi[0][RE][index]*psi[0][RE][index] + psi[0][IM][index]*psi[0][IM][index];
        updatedpsisum += updatedpsi[RE][index]*updatedpsi[RE][index] + updatedpsi[IM][index]*updatedpsi[IM][index];


    }
    printf("\n Psisum = %lf   Updatedsum = %lf\n", psisum, updatedpsisum);
    printf("\nDifference in sums is  %lf\n", (psisum-updatedpsisum));
    double normsum = 0.0;

    for (int index=0; index < N; index++) {
        normsum += updatedpsi[RE][index]*updatedpsi[RE][index] + updatedpsi[IM][index]*updatedpsi[IM][index];
    }

    printf("\nnormalization sum =   %lf\n", normsum);

    normsum = sqrt(normsum);

    for (int index=0; index < N; index++) {
        psi[0][RE][index] = updatedpsi[RE][index] /normsum;
        psi[0][IM][index] = updatedpsi[IM][index] /normsum;
    }
    normsum = 0.0;

    return;
}

void normalize()
{

    double normsum = 0.0;

    for (int index=0; index < N; index++) {
        normsum += psi[0][RE][index]*psi[0][RE][index] + psi[0][IM][index]*psi[0][IM][index];
    }

    normsum = sqrt(normsum);

    for (int index=0; index < N; index++) {
        psi[0][RE][index] = psi[0][RE][index]/normsum;
        psi[0][IM][index] = psi[0][IM][index]/normsum;
    }
    normsum = 0.0;

    return;
}


void renormalize()
{
    printf("\n commencing renormilaztion\n");

    return;
}

double calcprobabilityofmeasure(int whichfunc){

    double sum = 0.0;

    if (whichfunc == 0) {
        for (int index = 0; index < N; index++) {
            x = index * dx - (L/2);
            sum += abs((psi[0][RE][index]-psi[0][IM][index])*(psi[0][RE][index]+psi[0][IM][index]) * leftfunction(x)*leftfunction(x));
        }
    }

    if (whichfunc == 1) {
        for (int index = 0; index < N; index++) {
            x = index * dx - (L/2);
            sum += abs((psi[0][RE][index]-psi[0][IM][index])*(psi[0][RE][index]+psi[0][IM][index]) * rightfunction(x)*rightfunction(x));
        }
    }
    return sum;
}

void roller(int first){

    roll[first]=(double)rand()/RAND_MAX;//generate random real between 0 and 1

    return;
}

void calcnorm(int first)
{
    double temp = 0.0;
    for (int index = 0; index < N; index++)
    {
        temp += psi[0][RE][index]*psi[0][RE][index]+psi[0][IM][index]*psi[0][IM][index];
    }
    norms[first] = temp;
}

void outputnorm()
{
    static FILE *slicenorm;
    static char name[500];

    sprintf(name,"./slices/norm/slices_norm.dat");
    slicenorm=fopen(name,"w");

    for (int index = 0; index < numberslices; index++)
    {

        fprintf(slicenorm,"%d  %lf", index, norms[index]);
        fprintf(slicenorm,"\n");
    }

    fclose(slicenorm);

}


void outputfield(int first)//outputs the field values
{
    static FILE *slicefield;
    static char name[500];

    sprintf(name,"./slices/slices_fields_%d.dat", first);
    slicefield=fopen(name,"w");

    double psiprob[N];

    for (int index = 0; index < N; index++)
    {
        psiprob[index]= psi[0][RE][index]*psi[0][RE][index] + psi[0][IM][index]*psi[0][IM][index];
    }

    for (int index = 0; index < N; index++)
    {

        fprintf(slicefield,"%d  %lf  %lf  %lf", index, psi[0][RE][index], psi[0][IM][index], psiprob[index]);
        fprintf(slicefield,"\n");
    }

    fclose(slicefield);
}

void outputmeasurement()
{
    static FILE *slicemeas;
    static char name[500];

    sprintf(name,"./slices/norm/slices_meas.dat");
    slicemeas=fopen(name,"w");

    fprintf(slicemeas,"Slice index; Measurement taken? (1 yes, 0 no); Which side?; Prob of measurement; Rolled Value;\n");
    for (int index = 0; index < numberslices; index++)
    {

        fprintf(slicemeas,"%d  %d  %d  %lf  %lf", index,measurement[index], guess[index], probofmeasure[index],roll[index]);
        fprintf(slicemeas,"\n");
    }

    fclose(slicemeas);

}


double potential(double x)
{
    double U = 85.0;
    return U * pow((( 1.0 - cos(6.0*3.141592*(0.65*x-L/2.)/L))/2.), 2.0);
}

int main ()
{
    t=0.0;

    srand(time(NULL)); //Sets time seed

    fftw_complex *in, *out, *in2,*out2;
    fftw_plan plan, plan2;//plan will be foward and plan2 will be backward

    printf("\n***********************************\n\nSetting gaussian wavefunction... \n");
    for (int index = 0; index < N; index++)
    {
        x = (index*dx) - (L/2);
        k = 80*3.141592625/L;
        A = AMPLITUDE;
        x0 = -L/4;


        psi[0][RE][index] = A*cos(k*(x-x0)) * exp(-((x-x0)*(x-x0))/(4*sigma*sigma));//Sets gaussian real part of psi
        psi[0][IM][index] = A*sin(k*(x-x0)) * exp(-((x-x0)*(x-x0))/(4*sigma*sigma));//Sets gaussian for the imaginary parts of psi
        normalize();

    }
    printf("\nInitial conditions set \n");



    double realsum = 0.0;
    double complexsum = 0.0;
    double probabilitysum = 0.0;

    for (int index = 0; index < N; index++)
    {
        realsum += psi[0][RE][index];
        complexsum += psi[0][IM][index];
        probabilitysum += psi[0][RE][index] * psi[0][RE][index] + psi[0][IM][index] * psi[0][IM][index];
    }

    printf("\n\n\nRealsum = %lf\nComplexsum = %lf\n\nProbabilitysum = %lf\n\n",realsum,complexsum,probabilitysum);



    //create both fftw plans to be used to time evolve our wavefunction
    //Forward DFT
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


    //Backward DFT
    in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    plan2 = fftw_plan_dft_1d(N, in2, out2, FFTW_BACKWARD, FFTW_ESTIMATE);

    printf("FFTw plans set!\n\n\n Ready to begin\n");

    printf("Beginning time-evolution\n");



    //time evolution begins
    while (t<tf)
    {

        if (num % div == 0)
        {
            //outputfield(slicenum);
            slicenum++;

        }

        for (int index=0; index < N; index++)//update phase of position space wavefunction
        {
            x = index * dx - (L/2);
            psi[1][RE][index] = psi[0][RE][index] * cos(potential(x) * dt) + psi[0][IM][index]*sin(potential(x) * dt);
            psi[1][IM][index] = psi[0][IM][index] * cos(potential(x) * dt) - psi[0][RE][index] * sin(potential(x) * dt);
        }

        for (int index = 0; index < N; index++)//moves wavefunction to in[] for FFTw foward plan "plan"
        {
            in[index][0] = psi[1][RE][index];
            in[index][1] = psi[1][IM][index];
        }

        fftw_execute(plan);//transform now stored in out in FFTw format

        for (int index = 0; index < N; index++) {
            chi[0][RE][index] = out[index][0];
            chi[0][IM][index] = out[index][1];
        }

        for (int index = 0; index < N; index++)//here we update the phases in momentum space
        {
            p = ((2*3.1415926535)/L)*(((index + (N/2)) % N) - N/2);
            chi[1][RE][index] = chi[0][IM][index]*sin((dt*p*p)/2) + chi[0][RE][index]*cos((dt*p*p)/2);
            chi[1][IM][index] = chi[0][IM][index]*cos((dt*p*p)/2) - chi[0][RE][index]*sin((dt*p*p)/2);
        }

        for (int index = 0; index < N; index++)//moves wavefunction to in2[] for FFTw backwards plan "plan2"
        {
            in2[index][0] = chi[1][RE][index];
            in2[index][1] = chi[1][IM][index];
        }

        fftw_execute(plan2);

        for (int index = 0; index < N; index++) {
            psi[0][RE][index] = out2[index][0];
            psi[0][IM][index] = out2[index][1];
        }

        for (int index = 0; index < N; index++) //this loop accounts for unnormalized DFT after fwd and bkwd trnsfms
        {
            psi[0][RE][index] = psi[0][RE][index]/N;
            psi[0][IM][index] = psi[0][IM][index]/N;
        }





        calcnorm(num);
        guess[num] = rand()%2;
        probofmeasure[num]=calcprobabilityofmeasure(guess[num]);
        roller(num); //commence roll

        if (probofmeasure[num] > roll[num])
        {
            measurement[num] = 1;
            updatewfrm();
            //renormalize();


        }
        else{
            measurement[num] = 0;
        }

        num++;
        t += dt;


    }
    outputnorm();
    outputmeasurement();
    printf("*** time-evolution has been terminated ***\n");
    fftw_destroy_plan(plan);
    fftw_destroy_plan(plan2);

    fftw_free(in); fftw_free(in2); fftw_free(out); fftw_free(out2);



    printf("\n Number of Slices Outputted: %d \n", numberslices);
    return 0;



}
