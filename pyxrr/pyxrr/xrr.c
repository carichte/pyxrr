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
#include "Python.h"
#include <numpy/arrayobject.h>
#include <math.h>
#include <complex.h>


#define pi 3.141592653589793 //2384626433832795

struct MyMatrix { double complex m1;
                  double complex m2;
                  double complex m3;
                  double complex m4;
                 };

/* ==== Set up the methods table ====================== */
int  not_doublevector(PyArrayObject *vec)  {
    if (vec->descr->type_num != NPY_DOUBLE || vec->nd != 1)  {
        PyErr_SetString(PyExc_ValueError, 
        "In not_doublevector: array must be of type Float and 1 dimensional (n).");
        return 1;
    }
    return 0;
}

int  not_intvector(PyArrayObject *vec)  {
    if (vec->descr->type_num != NPY_LONG || vec->nd != 1)  {
        PyErr_SetString(PyExc_ValueError, 
        "In not_intvector: array must be of type Integer and 1 dimensional (n).");
        return 1;
    }
    return 0;
}

double *py_floatvector_to_Carrayptrs(PyArrayObject *arrayin)  {
    //int n;
    //n=arrayin->dimensions[0];
    return (double *) arrayin->data;  /* pointer to arrayin data as double */
}

int *py_intvector_to_Carrayptrs(PyArrayObject *arrayin)  {
    int n;
    n=arrayin->dimensions[0];
    int *swap = PyArray_DATA(arrayin);
    //swap = arrayin->data;
    int *result = PyArray_DATA(arrayin);
    int i=0;
    int j=0;
    while(j<n) {
        //printf("%i\n", swap[i]);
        if (swap[i] != 0)  {
            result[j]=swap[i];
            j++;
        }
        i++;
    }
    return (int *) result;
}


PyArrayObject *pyvector(PyObject *objin) {
     return (PyArrayObject *) PyArray_ContiguousFromObject(objin,
                              NPY_DOUBLE, 1,1);
}


void interface_mltply(struct MyMatrix *Matrix, const short int this, 
                      const short int next, 
                      double complex n[], double complex q[], 
                      double complex sintheta[], double sigma, double dnext,
                      const unsigned char pol) {
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
    M1=((*Matrix).m1 + (*Matrix).m3 * r_tb) / exp_iphi;
    M2=((*Matrix).m2 + (*Matrix).m4 * r_tb) / exp_iphi;
    M3=((*Matrix).m1 * r_tb + (*Matrix).m3) * exp_iphi;
    M4=((*Matrix).m2 * r_tb + (*Matrix).m4) * exp_iphi;
    (*Matrix).m1=M1;
    (*Matrix).m2=M2;
    (*Matrix).m3=M3;
    (*Matrix).m4=M4;
}

struct MyMatrix mmult2dim(struct MyMatrix *matrix1, struct MyMatrix *matrix2) {
    struct MyMatrix swap;
    swap.m1=(*matrix1).m1*(*matrix2).m1 + (*matrix1).m3*(*matrix2).m2;
    swap.m2=(*matrix1).m2*(*matrix2).m1 + (*matrix1).m4*(*matrix2).m2;
    swap.m3=(*matrix1).m1*(*matrix2).m3 + (*matrix1).m3*(*matrix2).m4;
    swap.m4=(*matrix1).m2*(*matrix2).m3 + (*matrix1).m4*(*matrix2).m4;
    return swap;
}

struct MyMatrix interface(int this, int next, int this_sigma, double complex n[], 
                          double complex q[], double complex sintheta[],
                          double sigma[], double d[], unsigned char pol) {
    double complex r_tb, exp_iphi;
    struct MyMatrix result;
    // Fresnelformel mit Rauigkeits Debye-Waller-like Faktoren
    if(pol==0) r_tb=(q[this]-q[next])
                    *cexp(-q[this]*q[next]*sigma[this_sigma]*sigma[this_sigma]/2)
                    /(q[this]+q[next]); //senkrecht polarisiert
    else r_tb=(n[this]*sintheta[next]-n[next]*sintheta[this])
                    *cexp(-q[this]*q[next]*sigma[this_sigma]*sigma[this_sigma]/2)
                    /(n[this]*sintheta[next]+n[next]*sintheta[this]); //parallel polarisiert
    // Phase des Folgelayers
    exp_iphi=cexp(I*d[next]*q[next]/2);
    // Matrix fuer Interface
    result.m1=1/exp_iphi;
    result.m2=r_tb/exp_iphi;
    result.m3=r_tb*exp_iphi;
    result.m4=exp_iphi;
    return result;
}



double amplitude(double theta, double lambda, int N[], int layercount[], double** d, 
                 int dlen[], double n[], double k[], double sigma[], int groups, 
                 int layers, int sigma_dim, const unsigned char pol) {
    short int h, i, j, l, m, jN, jL, i_anf, l_anf, next; //diverse indizes
    const double complex n_0=1.0-n[0]+I*k[0];
    double complex q[layers+1], n_compl[layers+1], sintheta[layers+1], dnext;
    struct MyMatrix Matrix;
    //double sinalpha=sin((90.0-theta)*pi/180);
    double sinalpha=cos(theta*pi/180);
    for(m=0; m<=layers; m++){
        n_compl[m]=1.0-n[m]+I*k[m]; //komplexe Brechzahlen
        sintheta[m]=csqrt(1-(n_0/n_compl[m])*(n_0/n_compl[m])*sinalpha*sinalpha);
        q[m]=4*pi*n_compl[m]*sintheta[m]/lambda; // Impulsuebertrag
    }
    // Matrix initialisieren
    Matrix.m1=1.0;
    Matrix.m2=0.0;
    Matrix.m3=0.0;
    Matrix.m4=1.0;
    i=0; // layer index (without substrate)
    l=0; // interface index
    //Produkt ausrechnen
    for ( h=0; h<groups; h++)  {
        i_anf=i;
        l_anf=l;
        
        // generiert einen Multilayer, Beispiel: 3 Schichten 1,2,3:
        for (j=0; j<layercount[h]; j++) { 
            if(j<(layercount[h]-1)) {
                //all except last layer of group
                next = i+1;
                dnext = d[next][0];
                interface_mltply(&Matrix, i, next, n_compl, q, sintheta, 
                                 sigma[l], dnext, pol);
                i++;
                l++;
            }
            
            // letzter Layer der Multischicht hat ersten Layer als Folgeschicht:
            else if(N[h]>1) {
                //iterating through periods
                for ( jN=1; jN<(N[h]); jN++)  {
                    
                    // generiert einen Multilayer, Beispiel: 3 Schichten 1,2,3
                    for ( jL=0; jL<layercount[h]; jL++)  {
                        if(jL==0) next = i_anf;
                        else next = i+1;
                        
                        if(dlen[next]>1) {
                            if (!(dlen[next] == N[h])) {
                                printf("Length mismach: %i %i", dlen[next], N[h]);
                                PyErr_SetString(PyExc_ValueError, 
                                "Assertion: len(d)+1 = len(delta) = len(beta).");
                            return 0;
                            }
                            dnext = d[next][jN];
                        }
                        else dnext = d[next][0];
                        
                        interface_mltply(&Matrix, i, next, n_compl, q,
                                         sintheta, sigma[l], dnext, pol);
                                        //hier ist dann Matrix=M1*M2*M3
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
                l++; //gleiche Schicht (M3) aber naechste Grenzflaeche
                //printf("%i %i \n", i, layers);
            }
        }
        //Die Letzte Schicht des Multilayers (Interface zu naechster Gruppe)
        if(i<layers) { // nur falls wir noch nicht im Substrat sind
            // dort gibts naemlich keinen folgelayer i+1
            next = i+1;
            if(i<(layers-1)) dnext = d[next][0];
            else dnext=0;
            interface_mltply(&Matrix, i, next, n_compl, q, sintheta, 
                             sigma[l], dnext, pol);
            //Matrix=mmult2dim(&Matrix, &M);
            i++;
            l++;
        }
        //printf("%i %i %i\n", i, layercount[h], j);
    }
    //printf("%i %i \n", i, l);
    return cabs(Matrix.m2/Matrix.m1);
}


double amplitude_fast(double theta, double lambda, int N[], int layercount[], double d[], 
                      double n[], double k[], double sigma[], int groups, int layers, 
                      int sigma_dim, unsigned char pol) {
    int i, l, h,j,m, jN; //diverse indizes
    const double complex n_0=1.0-n[0]+I*k[0];
    double complex q[layers+1], n_compl[layers+1], sintheta[layers+1];
    struct MyMatrix M, Matrix, M_multilayer;
    //double sinalpha=sin((90.0-theta)*pi/180);
    double sinalpha=cos(theta*pi/180);
    for(m=0; m<=layers; m++){
        n_compl[m]=1.0-n[m]+I*k[m]; //komplexe Brechzahlen
        sintheta[m]=csqrt(1-(n_0/n_compl[m])*(n_0/n_compl[m])*sinalpha*sinalpha);
        q[m]=4*pi*n_compl[m]*sintheta[m]/lambda; // Impulsuebertrag
    }
    // Matrix initialisieren
    Matrix.m1=1.0;
    Matrix.m2=0.0;
    Matrix.m3=0.0;
    Matrix.m4=1.0;
    i=0;
    l=0;
    for ( h=0; h<groups; h++)  {
        M_multilayer.m1=1.0;
        M_multilayer.m2=0.0;
        M_multilayer.m3=0.0;
        M_multilayer.m4=1.0;
        // generiert einen Multilayer, Beispiel: 3 Schichten 1,2,3:
        for ( j=0; j<layercount[h]; j++)  {
            if(j<(layercount[h]-1)) {
                M=interface(i, i+1, l, n_compl, q, sintheta, sigma, d, pol);
                // Wichtig fuer die Reihenfolge:
                //  zuletzt muss Matrix fuer letztes Interface des ML stehen
                M_multilayer=mmult2dim(&M_multilayer, &M); //macht Multilayer=M1*M2
                Matrix=mmult2dim(&Matrix, &M); //gesamtschicht: Matrix=M1*M2
                i++;
                l++;
            }
            // letzter Layer der Multischicht hat ersten Layer als Folgeschicht:
            else if(N[h]>1) {
                M=interface(i, i+1-layercount[h], l, n_compl, q, sintheta, sigma, d, pol);
                M_multilayer = mmult2dim(&M, &M_multilayer);
                //hier ist dann Multilayer=M3*M1*M2
                l++;
                //jetzt zerlegung in binaeres Format um Rechenzeit zu sparen
                jN = (N[h]-1);
                while(jN > 0) {
                    if ((jN%2)==1) Matrix=mmult2dim(&Matrix, &M_multilayer);
                    jN = jN/2;
                    M_multilayer = mmult2dim(&M_multilayer, &M_multilayer);
                }
                /*for ( jN=0; jN<(N[h]-1); jN++)  {
                    Matrix=mmult2dim(&Matrix, &M_multilayer);
                    //hier: Matrix=M1*M2*(M3*M1*M2)**(N-1)
                }*/
            }
        }
        //Die Letzte Schicht des Multilayers (Interface zu naechster Gruppe)
        if(i<layers) { // nur falls wir noch nicht im Substrat sind
            // dort gibts naemlich keinen folgelayer i+1
            M=interface(i, i+1, l, n_compl, q, sintheta, sigma, d, pol);
            Matrix=mmult2dim(&Matrix, &M);
            i++;
            l++;
        }
        //printf("%i %i %i\n", i, layercount[h], j);
    }
    //printf("%i %i \n", i, l);
    return cabs(Matrix.m2/Matrix.m1);
}



static PyObject *reflectivity(PyObject *self, PyObject *args)  {
    PyArrayObject *theta_range, *r_values, *N, *layercount, *n, *k, *sigma;
    PyObject* dseq;
    double *cin, *cout, r_s, r_p, lambda, *n_in, *k_in, *sigma_in;
    double **d;
    
    float polarization;
    // The C vectors to be created to point to the 
    //   python vectors, cin and cout point to the row
    //   of vecin and vecout, respectively
    int i, dim, N_dim, d_dim, n_dim, k_dim, sigma_dim, dims[2], 
        *N_ptr, *layercount_ptr, *dlen, lc_dim, cpus, fast = 1;
    /* Parse tuples separately since args will differ between C fcns */
    if (!PyArg_ParseTuple(args, "O!O!O!OO!O!O!df", 
                            &PyArray_Type, &theta_range, 
                            &PyArray_Type, &N, 
                            &PyArray_Type, &layercount, 
                                           &dseq, // list of arrays!
                            &PyArray_Type, &n, 
                            &PyArray_Type, &k, 
                            &PyArray_Type, &sigma, 
                            &lambda, 
                            &polarization)) return NULL;
    if (NULL == theta_range)  return NULL;
    if (NULL == N)  return NULL;
    if (NULL == layercount)  return NULL;
    if (NULL == dseq)  return NULL;
    if (NULL == n)  return NULL;
    if (NULL == k)  return NULL;
    if (NULL == sigma)  return NULL;
    /* Check that object input is 'double' type and a vector
    Not needed if python wrapper function checks before call to this routine */
    if (not_doublevector(theta_range)) return NULL;
    if (not_intvector(N)) return NULL;
    if (not_intvector(layercount)) return NULL;
    if(!dseq) return 0;
    if (not_doublevector(n)) return NULL;
    if (not_doublevector(k)) return NULL;
    if (not_doublevector(sigma)) return NULL;
    
    /* Get the dimension of the input */
    dim=dims[0]=theta_range->dimensions[0];
    N_dim=N->dimensions[0];
    lc_dim=layercount->dimensions[0];
    n_dim=n->dimensions[0];
    k_dim=k->dimensions[0];
    sigma_dim=sigma->dimensions[0];
    
    dseq = PySequence_Fast(dseq, "argument must be iterable");
    d_dim = PySequence_Fast_GET_SIZE(dseq); // get length of given list of d arrays
    d = malloc(sizeof(double *) * d_dim); // allocate memory for array of pointers
    dlen = malloc(sizeof(int) * d_dim); // allocate memory for array of int
    
    if((!d) || (!dlen)) {
        Py_DECREF(dseq);
        return PyErr_NoMemory();
    }
    
    
    /* check dimensions of the input */
    if (!((d_dim+1 == n_dim) && (d_dim+1 == k_dim))) {
        PyErr_SetString(PyExc_ValueError, 
                        "Assertion: len(d)+1 = len(delta) = len(beta).");
        return NULL;
    }
    if (!(N_dim == lc_dim)) {
        PyErr_SetString(PyExc_ValueError,
                        "N and LayerCount have to be of equal length.");
        return NULL;
    }
    
    
    for(i=0; i < d_dim; i++) {
        //PyObject *item = PyList_GET_ITEM(dseq, i);
        PyArrayObject* aitem = pyvector(PyList_GET_ITEM(dseq, i));
        if(!aitem) {
            Py_DECREF(dseq);
            free(d);
            return 0;
        }
        if (not_doublevector(aitem)) {
            Py_DECREF(dseq);
            free(d);
            PyErr_SetString(PyExc_TypeError, 
                            "all items must be arrays of float64");
            return 0;
        }
        dlen[i] = aitem->dimensions[0];
        d[i] = py_floatvector_to_Carrayptrs(aitem);
        Py_DECREF(aitem);
        //printf("%i %i \n", i, dlen[i]);
        
    }
    Py_DECREF(dseq);
    
    /* Make a new double vector of same dimension */
    //r_values=(PyArrayObject *) PyArray_SimpleNew(1,(npy_intp*)(dims),NPY_DOUBLE);
    r_values=(PyArrayObject *) PyArray_FromDims(1,dims,NPY_DOUBLE);
    
    /* Change contiguous arrays into C *arrays   */
    cin = py_floatvector_to_Carrayptrs(theta_range);
    //cin = theta_range->data;  /* would this also work? */
    //int *N_ptr = PyArray_DATA(N);
    N_ptr= py_intvector_to_Carrayptrs(N);
    layercount_ptr = py_intvector_to_Carrayptrs(layercount);
    //int *layercount_ptr=layercount->data;
    n_in=py_floatvector_to_Carrayptrs(n);
    k_in=py_floatvector_to_Carrayptrs(k);
    sigma_in=py_floatvector_to_Carrayptrs(sigma);
    cout=py_floatvector_to_Carrayptrs(r_values);

    // Is the Multilayer strictly periodic?
    for (i=0; i<d_dim; i++) fast = fast && (dlen[i] == 1);
    
    
    /* Do the calculation. */
    
    
    cpus = omp_get_num_procs();
    omp_set_num_threads(cpus);
    
    if(fast) {
        double *d_in;
        d_in = malloc(sizeof(double) * d_dim);
        if(!d_in) {
            Py_DECREF(dseq);
            return PyErr_NoMemory();
        }
        for (i=0; i<d_dim; i++) d_in[i] = d[i][0];
        
        if ( polarization == 0 ) { 
            #pragma omp parallel for private(r_s)
            for ( i=0; i<dim; i++) {
                r_s=amplitude_fast(cin[i], lambda, N_ptr, layercount_ptr, d_in,
                              n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 0);
                cout[i]=r_s*r_s;
            }
        }
        else if ( polarization == 1 ) { 
            #pragma omp parallel for private(r_p)
            for ( i=0; i<dim; i++) {
                r_p=amplitude_fast(cin[i], lambda, N_ptr, layercount_ptr, d_in,
                              n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 1);
                cout[i]=r_p*r_p;
            }
        }
        else {
            #pragma omp parallel for private(r_s, r_p)
            for ( i=0; i<dim; i++) {
                r_s=amplitude_fast(cin[i], lambda, N_ptr, layercount_ptr, d_in,
                              n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 0);
                r_p=amplitude_fast(cin[i], lambda, N_ptr, layercount_ptr, d_in,
                              n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 1);
                cout[i]=r_s*r_s*(1-polarization) +  r_p*r_p*polarization;
            }
        }
    }
    else {
        if ( polarization == 0 ) { 
            #pragma omp parallel for private(r_s)
            for ( i=0; i<dim; i++) {
                r_s=amplitude(cin[i], lambda, N_ptr, layercount_ptr, d, dlen,
                              n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 0);
                cout[i]=r_s*r_s;
            }
        }
        else if ( polarization == 1 ) { 
            #pragma omp parallel for private(r_p)
            for ( i=0; i<dim; i++) {
                r_p=amplitude(cin[i], lambda, N_ptr, layercount_ptr, d, dlen,
                              n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 1);
                cout[i]=r_p*r_p;
            }
        }
        else {
            #pragma omp parallel for private(r_s, r_p)
            for ( i=0; i<dim; i++) {
                r_s=amplitude(cin[i], lambda, N_ptr, layercount_ptr, d, dlen,
                              n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 0);
                r_p=amplitude(cin[i], lambda, N_ptr, layercount_ptr, d, dlen,
                              n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 1);
                cout[i]=(r_s*r_s)*(1-polarization) + (r_p*r_p)*polarization;
            }
        }
    }
    free(d);
    return PyArray_Return(r_values);
}

static PyMethodDef methods[] = {
    {"reflectivity", reflectivity, METH_VARARGS, 
    "reflectivity(theta_range, N, LayerCount, d, delta, beta, sigma, lambda, pol) \n \
    \n \
    Calculates the Reflectivity of a Multilayer.\n \
    \n \
    \n \
    Inputs:\n \
     - theta_range: array of glancing angles theta\n \
     - N array: with number of periods for each group\n \
     - LayerCount: array containing number of unique Layers per group\n \
     - d: array containing the d-spacings (Angstrom) of all unique layers except substrate\n \
     - delta & beta: arrays containing the optical constants of all unique layers\n \
     - sigma: array containing the RMS roughness parameters (Angstrom) of all unique interfaces\n \
     - lambda: Wavelength (Angstrom, float)\n \
     - pol: part of Radiation that is parallel-polarized (float)\n \
            (0 => perpendicular, 1=> parallel, 0.5=> unpolarized)\n\n \
      len(N) = len(LayerCount)\n \
      len(d)+1 = len(delta) = len(beta)\n \
      len(sigma) > len(d)\n\
      \n \
     All values ordered from top to bottom!\n\
      \n\
      \n\
    Outputs:\n \
     - array of reflectivity values for each theta"},
    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};

void initxrr() {
    (void) Py_InitModule3("xrr", methods, 
    "Contains only one function, which is for calculating X-Ray-Reflectivity (xrr).\n\
    It is supposed to be used by the python xrr wrapper (pyxrr)\n\n\
    See xrr.reflectivity for more details.");
    import_array();
}