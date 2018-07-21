// Copyright 2013 by Carsten Richter
// Contact: carsten.richter@desy.de and
//          carsten.richter@physik.tu-freiberg.de
//
// This file is part of pyxrr.
//
// pyxrr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pyxrr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pyxrr inside the 'copying.txt' file.
// If not, see <http://www.gnu.org/licenses/>.



#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>

#define pi 3.141592653589793 //2384626433832795

struct MyMatrix { double complex m1;
                  double complex m2;
                  double complex m3;
                  double complex m4;
                 };






void interface_mltply(struct MyMatrix *Matrix,
                      const unsigned short this, 
                      const unsigned short next, 
                      double complex n[],
                      double complex q[],
                      double complex sintheta[],
                      double sigma,
                      double dnext,
                      const unsigned char pol) {
    /*
        combined calculation of interface reflectivity and multiplication onto
        given matrix.
    */
    // Fresnelformel mit Rauigkeits Debye-Waller-like Faktoren
    double complex r_tb, exp_iphi, M1, M2, M3, M4;
    if(pol==0) r_tb=( q[this]-q[next])
                      *cexp(-q[this]*q[next]*sigma*sigma/2)
                      /(q[this]+q[next]); //senkrecht polarisiert
    else r_tb=(n[this]*sintheta[next]-n[next]*sintheta[this])
              *cexp(-q[this]*q[next]*sigma*sigma/2)
              /(n[this]*sintheta[next]+n[next]*sintheta[this]); //parallel polarisiert
    // Phase des Folgelayers
    exp_iphi=cexp(I * dnext * q[next] / 2);
    // Matrix fuer Interface
    //printf("%i -> %i : dnext=%f, rough=%f, phase=(%f,%f)\n", this, next,dnext, sigma, creal(exp_iphi), cimag(exp_iphi));
    M1=((*Matrix).m1 + (*Matrix).m3 * r_tb) / exp_iphi;
    M2=((*Matrix).m2 + (*Matrix).m4 * r_tb) / exp_iphi;
    M3=((*Matrix).m1 * r_tb + (*Matrix).m3) * exp_iphi;
    M4=((*Matrix).m2 * r_tb + (*Matrix).m4) * exp_iphi;
    /*
    printf("M = [[%f+i*%f, %f+i*%f], [%f+i*%f, %f+i*%f]]\n", 
            creal(M1), cimag(M1),
            creal(M1), cimag(M1),
            creal(M1), cimag(M1),
            creal(M1), cimag(M1));
    */
    (*Matrix).m1=M1;
    (*Matrix).m2=M2;
    (*Matrix).m3=M3;
    (*Matrix).m4=M4;
}

struct MyMatrix mmult2dim(struct MyMatrix *matrix1, struct MyMatrix *matrix2) {
    /*
        simple matrix multiplication of 2x2 Matrices (MyMatrix)
    */
    struct MyMatrix swap;
    swap.m1=(*matrix1).m1*(*matrix2).m1 + (*matrix1).m3*(*matrix2).m2;
    swap.m2=(*matrix1).m2*(*matrix2).m1 + (*matrix1).m4*(*matrix2).m2;
    swap.m3=(*matrix1).m1*(*matrix2).m3 + (*matrix1).m3*(*matrix2).m4;
    swap.m4=(*matrix1).m2*(*matrix2).m3 + (*matrix1).m4*(*matrix2).m4;
    return swap;
}


double complex amplitude(double theta,
                 double lambda,
                 unsigned int periods[],
                 unsigned int groupsize[],
                 double *thickness, 
                 short periodic[],
                 double complex n[],
                 double roughness[],
                 unsigned int ngroups, 
                 unsigned int nlayers,
                 unsigned int ninterfaces,
                 const unsigned char pol) {
    /* 
       h  : group idx
       i  : layer idx
       j  : layer idx in group
       l  : interface idx
       jL : layer idx in multilayer group (could be just j)
       i_anf : idx of first layer of group
       l_anf : idx of first interface of group
       next : idx of next layer
       jN : index of period in multilayer
    */
    unsigned short h, i, j, l, jL, i_anf, l_anf, next;
    unsigned int jN;
    double complex q[nlayers], sintheta[nlayers], dnext;
    struct MyMatrix Matrix, MMultilayer;
    //double sinalpha=sin((90.0-theta)*pi/180);
    double sinalpha = cos(theta*pi/180);

    //printf("%f\n", theta);
    for(i=0; i<nlayers; i++){
        sintheta[i] = csqrt(1 - (n[0]/n[i])*(n[0]/n[i]) * sinalpha*sinalpha);
        q[i] = 4*pi*n[i]*sintheta[i]/lambda; // momentum transfer
    }

    // init Matrix
    Matrix.m1=1.0;
    Matrix.m2=0.0;
    Matrix.m3=0.0;
    Matrix.m4=1.0;
    i=0; // layer index (without substrate)
    l=1; // interface index
    // compute product
    for ( h=0; h<ngroups; h++)  {
        i_anf=i; // first layer of the group
        l_anf=l;

        if (periodic[h] && periods[h]>1) {
            // prepare multilayer
            MMultilayer.m1=1.0;
            MMultilayer.m2=0.0;
            MMultilayer.m3=0.0;
            MMultilayer.m4=1.0;
            next = i_anf;
            dnext = thickness[next];
            // first the last interface of the group:
            interface_mltply(&MMultilayer,
                             i_anf+groupsize[h]-1,
                             i_anf,
                             n,
                             q,
                             sintheta,
                             roughness[l_anf],//[l+groupsize[h]-1], 
                             thickness[i_anf],
                             pol);
            l++;
        }
        // generate multilayer, example: 3 layers 1,2,3:
        for (j=0; j<groupsize[h]-1; j++) { //all except last layer of group
            next = i+1;
            interface_mltply(&Matrix,
                             i,
                             next,
                             n,
                             q,
                             sintheta, 
                             roughness[l],
                             thickness[next],
                             pol); //makes the multilayer part (M1*M2)
            if (periodic[h] && periods[h]>1) interface_mltply(&MMultilayer,
                                                               i,
                                                               next,
                                                               n,
                                                               q,
                                                               sintheta,
                                                               roughness[l],
                                                               thickness[next],
                                                               pol); //the actually periodic part Multilayer=M3*M1*M2
            i++; // consume a layer
            l++; // consume an interface
        }



        // finally: last Layer of Multilayer has the first as successor
        if(periods[h]>1) {
            //iterate through periods
            if (periodic[h]) {
                //currently the multilayer is M3*M1*M2
                //decomposition into binary format to speed up
                jN = (periods[h]-1);
                while(jN > 0) {
                    if ((jN%2)==1) Matrix=mmult2dim(&Matrix, &MMultilayer);
                    jN = jN/2;
                    MMultilayer = mmult2dim(&MMultilayer, &MMultilayer);
                }
            }
            else {
                for ( jN=1; jN<(periods[h]); jN++)  { // The non-periodic case -- later
                    // generate multilayer, example: 3 layers 1,2,3:
                    for ( jL=0; jL<groupsize[h]; jL++)  {
                        if(jL==0) next = i_anf;
                        else next = i+1;
                        //printf("%i %f\n", jN, theta);

                        //if(dlen[next]>1) dnext = d[next][jN]; //later
                        //else dnext = d[next];
                        dnext = thickness[next];

                        interface_mltply(&Matrix,
                                         i,
                                         next,
                                         n,
                                         q,
                                         sintheta,
                                         roughness[l],
                                         dnext,
                                         pol); //hier ist dann Matrix=M1*M2*M3
                        if(jL==0) {
                            i=i_anf;
                            l=l_anf;
                        }
                        else {
                            i++;
                            l++;
                        }
                    }
                }
            }
        //l++; //gleiche Schicht (M3) aber naechste Grenzflaeche
        //printf("%i %i \n", i, nlayers);
        }

        if(i<(nlayers-1)) {
            // last layer of periodic group (interface to next group)
            // but only if we are not yet in the substrate
            // since there is no following layer i+1
            next = i+1;
            if(i<(nlayers-2)) dnext = thickness[next];
            else dnext = 0;
            //printf("last block, h=%i, size=%i, ", h, groupsize[h]);
            interface_mltply(&Matrix,
                             i,
                             next,
                             n,
                             q,
                             sintheta, 
                             roughness[l],
                             dnext,
                             pol);
            i++;
            l++;
        }
        //printf("%i %i %i\n", i, groupsize[h], j);
    }
    //printf("%i %i \n", i, l);
    /*
    printf("M1=%f+i%f, M2=%f+i%f, M=%f+i%f, ",
         creal(Matrix.m1),
         cimag(Matrix.m1),
         creal(Matrix.m2),
         cimag(Matrix.m2),
         creal(Matrix.m2/Matrix.m1),
         cimag(Matrix.m2/Matrix.m1)
        );
    */
    return Matrix.m2/Matrix.m1;
}


void reflectivity(double complex *output,
                  double *theta,
                  double *thickness,
                  double *roughness,
                  double complex *n,
                  double lambda,
                  double polarization,
                  unsigned int *periods,
                  unsigned int *groupsize,
                  short *periodic,
                  unsigned int ntheta,
                  unsigned int ngroups,
                  unsigned int nlayers,
                  unsigned int ninterfaces,
                  unsigned int numthreads) {
/*
 *    Calculates the Reflectivity of a Multilayer.
 *
 *    Inputs:\n \
 *     - output: output array for the result
 *     - theta: array of glancing angles
 *     - thickness: array of layer thicknesses in angstrom
 *     - roughness: array of interface roughnesses in angstrom
 *     - n: refractive index for each layer
 *     - lambda: x-ray wavelength in angstrom
 *     - polarization: x-ray wavelength in angstrom
 *     - periods: contains number of periods for each group
 *     - groupsize: array containing number of unique Layers per group
 *     - ntheta: lengths of theta array
 *     - ngroups: total number of groups
 *     - nlayers: total number of layers
 *     - ninterfaces: total number of interfaces
 *     - numthreads: number of CPUs to use for calculation
 *
 */
 
    int i;

    if (numthreads<=0) numthreads = omp_get_num_procs();
    else if (numthreads>omp_get_num_procs()) numthreads = omp_get_num_procs();
    omp_set_num_threads(numthreads);
    
    /* Do the calculation. */
    // case sigma amplitude (polarization == -1)
    if ( fabs(polarization + 1) < 0.001 ) { 
        #pragma omp parallel for
        for ( i=0; i<ntheta; i++) {
            output[i] = amplitude(theta[i],
                                  lambda,
                                  periods,
                                  groupsize,
                                  thickness,
                                  periodic,
                                  n,
                                  roughness,
                                  ngroups,
                                  nlayers,
                                  ninterfaces,
                                  0);
        }
    }
    // case pi amplitude (polarization == -2)
    else if ( fabs(polarization + 2) < 0.001 ) { 
        #pragma omp parallel for
        for ( i=0; i<ntheta; i++) {
            output[i] = amplitude(theta[i],
                                  lambda,
                                  periods,
                                  groupsize,
                                  thickness,
                                  periodic,
                                  n,
                                  roughness,
                                  ngroups,
                                  nlayers,
                                  ninterfaces,
                                  1);
        }
    }
    // case pi intensity (polarization == 1)
    else if ( polarization > .999 ) { 
        #pragma omp parallel for
        for ( i=0; i<ntheta; i++) {
            output[i] = cabs(amplitude(theta[i],
                                  lambda,
                                  periods,
                                  groupsize,
                                  thickness,
                                  periodic,
                                  n,
                                  roughness,
                                  ngroups,
                                  nlayers,
                                  ninterfaces,
                                  1));
            output[i] = output[i]*output[i];
        }
    }
    // case sigma intensity (polarization == 0)
    else if ( fabs(polarization) < .001 ) { 
        #pragma omp parallel for
        for ( i=0; i<ntheta; i++) {
            output[i] = cabs(amplitude(theta[i],
                                  lambda,
                                  periods,
                                  groupsize,
                                  thickness,
                                  periodic,
                                  n,
                                  roughness,
                                  ngroups,
                                  nlayers,
                                  ninterfaces,
                                  0));
            output[i] = output[i]*output[i];
        }
    }
    // case total intensity (polarization == 0)
    else {
        double r_s, r_p;
        #pragma omp parallel for private(r_s, r_p)
        for ( i=0; i<ntheta; i++) {
            r_s = cabs(amplitude(theta[i],
                                  lambda,
                                  periods,
                                  groupsize,
                                  thickness,
                                  periodic,
                                  n,
                                  roughness,
                                  ngroups,
                                  nlayers,
                                  ninterfaces,
                                  0));
            r_p = cabs(amplitude(theta[i],
                                  lambda,
                                  periods,
                                  groupsize,
                                  thickness,
                                  periodic,
                                  n,
                                  roughness,
                                  ngroups,
                                  nlayers,
                                  ninterfaces,
                                  1));
            output[i] = (r_s*r_s)*(1-polarization) + (r_p*r_p)*polarization;
        }
    }


}


