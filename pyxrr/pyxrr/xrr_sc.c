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



//#include <omp.h>
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

int not_intvector(PyArrayObject *vec)  {
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
    int n, i=0;
    n=arrayin->dimensions[0];
    int *result = PyArray_DATA(arrayin);
    
    //PyArray_ISSIGNED(arrayin) &&
    if (PyArray_ITEMSIZE(arrayin)==sizeof(int)) ;
    else if (PyArray_ITEMSIZE(arrayin)==sizeof(long)) {
        long *swap = PyArray_DATA(arrayin);
        for (i=0; i<n; i++) { 
            result[i] = (int) swap[i];
        }
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
    /*
        combined calculation of interface reflectivity and multiplication onto
        given matrix.
    */
    // Fresnelformel mit Rauigkeits Debye-Waller-like Faktoren
    double complex r_tb, exp_iphi, M1, M2, M3, M4;
    //printf("%i -> %i\n", this, next);
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



double amplitude(double theta, double lambda, int N[], int layercount[], double** d, 
                 int dlen[], double n[], double k[], double sigma[], int groups, 
                 int layers, int sigma_dim, const unsigned char pol) {
    short int h, i, j, l, m, jN, jL, i_anf, l_anf, next, periodic=1; //diverse indizes
    const double complex n_0=1.0-n[0]+I*k[0];
    double complex q[layers+1], n_compl[layers+1], sintheta[layers+1], dnext;
    struct MyMatrix Matrix, MMultilayer;
    //double sinalpha=sin((90.0-theta)*pi/180);
    double sinalpha=cos(theta*pi/180);
    char myerror[80];
    
    
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
        
        //check dimensions:
        periodic = 1;
        //printf("%i %i\n",h, N[h]);
        for (j=0; j<layercount[h]; j++) { 
            if((j+i_anf)==layers) continue;
            else if (!(dlen[j+i_anf] == N[h] || dlen[j+i_anf]<=1)) {
                sprintf(myerror, 
                        "Length mismach: len(d)=%i!=%i=N for group %i.",
                        dlen[j+i_anf], N[h], h);
                puts(myerror);
                PyErr_SetString(PyExc_ValueError, myerror);
                return 0;
            }
            else if(dlen[j+i_anf]>1) periodic=0;
        }
        //printf("%i %i %i %f\n", h, layercount[h], periodic, theta);
        
        if (periodic && N[h]>1) {
            MMultilayer.m1=1.0;
            MMultilayer.m2=0.0;
            MMultilayer.m3=0.0;
            MMultilayer.m4=1.0;
            next = i_anf;
            dnext = d[next][0];
            // first the last interface of the group:
            interface_mltply(&MMultilayer, i_anf+layercount[h]-1, i_anf, 
                             n_compl, q, sintheta, sigma[l+layercount[h]-1], 
                             dnext, pol);
            //M = interface(i, i+1-layercount[h], l, n_compl, q, sintheta, sigma, d, pol);
        }
        // generiert einen Multilayer, Beispiel: 3 Schichten 1,2,3:
        for (j=0; j<layercount[h]-1; j++) { 
            //all except last layer of group
            next = i+1;
            dnext = d[next][0];
            interface_mltply(&Matrix, i, next, n_compl, q, sintheta, 
                             sigma[l], dnext, pol); //macht Multilayer=M1*M2
            if (periodic) interface_mltply(&MMultilayer, i, next, n_compl, q,
                    sintheta, sigma[l], dnext, pol); //macht Multilayer=M3*M1*M2
            i++;
            l++;
        }
        // letzter Layer der Multischicht hat ersten Layer als Folgeschicht:
        if(N[h]>1) {
            //iterating through periods
            if (periodic) {
                //hier ist dann Multilayer=M3*M1*M2
                //jetzt zerlegung in binaeres Format um Rechenzeit zu sparen
                jN = (N[h]-1);
                while(jN > 0) {
                    if ((jN%2)==1) Matrix=mmult2dim(&Matrix, &MMultilayer);
                    jN = jN/2;
                    MMultilayer = mmult2dim(&MMultilayer, &MMultilayer);
                }
            }
            else {
                for ( jN=1; jN<(N[h]); jN++)  {
                    // generiert einen Multilayer, Beispiel: 3 Schichten 1,2,3
                    for ( jL=0; jL<layercount[h]; jL++)  {
                        if(jL==0) next = i_anf;
                        else next = i+1;
                        
                        if(dlen[next]>1) dnext = d[next][jN];
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
            }
        l++; //gleiche Schicht (M3) aber naechste Grenzflaeche
        //printf("%i %i \n", i, layers);
        }
        
        if(i<layers) {
            // last layer of periodic stack (interface to next group)
            // but only if we are not yet in the substrate
            // since there is no following layer i+1
            next = i+1;
            if(i<(layers-1)) dnext = d[next][0];
            else dnext=0;
            interface_mltply(&Matrix, i, next, n_compl, q, sintheta, 
                             sigma[l], dnext, pol);
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
        *dlen, lc_dim;
    int *N_ptr, *layercount_ptr;
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
            free(dlen);
            return 0;
        }
        if (not_doublevector(aitem)) {
            Py_DECREF(dseq);
            free(d);
            free(dlen);
            PyErr_SetString(PyExc_TypeError, 
                            "all items of ``d'' must be arrays of float64");
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
    //cin = (double *) theta_range->data;  /* would this also work? */
    N_ptr= py_intvector_to_Carrayptrs(N);
    layercount_ptr = py_intvector_to_Carrayptrs(layercount);
    n_in=py_floatvector_to_Carrayptrs(n);
    k_in=py_floatvector_to_Carrayptrs(k);
    sigma_in=py_floatvector_to_Carrayptrs(sigma);
    cout=py_floatvector_to_Carrayptrs(r_values);

    
    /* Do the calculation. */
    //cpus = omp_get_num_procs();
    //omp_set_num_threads(cpus);
    //omp_set_num_threads(1);
    
    if ( polarization == 0 ) { 
        //#pragma omp parallel for private(r_s)
        for ( i=0; i<dim; i++) {
            r_s=amplitude(cin[i], lambda, N_ptr, layercount_ptr, d, dlen,
                          n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 0);
            cout[i]=r_s*r_s;
        }
    }
    else if ( polarization == 1 ) { 
        //#pragma omp parallel for private(r_p)
        for ( i=0; i<dim; i++) {
            r_p=amplitude(cin[i], lambda, N_ptr, layercount_ptr, d, dlen,
                          n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 1);
            cout[i]=r_p*r_p;
        }
    }
    else {
        //#pragma omp parallel for private(r_s, r_p)
        for ( i=0; i<dim; i++) {
            r_s=amplitude(cin[i], lambda, N_ptr, layercount_ptr, d, dlen,
                          n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 0);
            r_p=amplitude(cin[i], lambda, N_ptr, layercount_ptr, d, dlen,
                          n_in, k_in, sigma_in, N_dim, d_dim, sigma_dim, 1);
            cout[i]=(r_s*r_s)*(1-polarization) + (r_p*r_p)*polarization;
        }
    }
    free(d);
    free(dlen);
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
     - d: list of arrays containing the d-spacings (Angstrom) of all unique layers except substrate.\n \
          The array shape must be (1,) or (N,) where N is the number ob periods of the corresponding group.\
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