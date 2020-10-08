#include <stdlib.h>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <time.h>

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "myHeader.cpp"
#include <Eigen/Core>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include <unistd.h>

#include <thread>

//#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;

//const double PI  =3.141592653589793238462;
/*
Chuan Gao C++
*/

/* 
usage
./BicMix --nf 100 --y ./AOAS/data.txt --out result --sep space
*/

int main(int argc, char *argv[])
{

    clock_t start, end;
    double cpu_time_used;

    start = clock();

    // declare variables
    int nf = 100, s_n = 0, d_y = 0, i_seed = 0, n_itr = 5001, write_itr = 50, interval = 500;
    double a = 0.5, b = 0.5, c = 0.5, d = 0.5, g = 1, h = 1, alpha = 1, beta = 1;
    double tol = 1e-3;

    string file_y, dir_out, sep = "tab";
    string lam_method_in = "matrix", x_method_in = "dense";
    stringstream ss;

    // read in argument
    string s_nf = "--nf", s_y = "--y", s_out = "--dir_out", s_sep = "--sep", s_a = "--a", s_b = "--b", s_c = "--c", s_d = "--d", s_e = "--e", s_f = "--f", s_seed = "--seed", s_itr = "--itr", s_interval = "--interval", s_write_itr = "--write_itr", s_tol = "--tol", s_lam_method = "--lam_method", s_x_method = "x_method";

    for (int i = 0; i < argc; i++)
    {
        if (s_nf.compare(argv[i]) == 0){nf = atoi(argv[i + 1]);}
        if (s_y.compare(argv[i]) == 0){file_y = argv[i + 1];}
        if (s_out.compare(argv[i]) == 0){dir_out = argv[i + 1];}
        if (s_sep.compare(argv[i]) == 0){sep = argv[i + 1];}
        if (s_a.compare(argv[i]) == 0){a = atof(argv[i + 1]);}
        if (s_b.compare(argv[i]) == 0){b = atof(argv[i + 1]);}

        if (s_c.compare(argv[i]) == 0){c = atof(argv[i + 1]);}
        if (s_d.compare(argv[i]) == 0){d = atof(argv[i + 1]);}

        if (s_e.compare(argv[i]) == 0){g = atof(argv[i + 1]);}
        if (s_f.compare(argv[i]) == 0){h = atof(argv[i + 1]);}

        if (s_interval.compare(argv[i]) == 0){interval = atoi(argv[i + 1]);}
        if (s_write_itr.compare(argv[i]) == 0){write_itr = atoi(argv[i + 1]);}
        if (s_seed.compare(argv[i]) == 0){i_seed = atoi(argv[i + 1]);}
        if (s_itr.compare(argv[i]) == 0){n_itr = atoi(argv[i + 1]);}
        if (s_tol.compare(argv[i]) == 0){tol = atof(argv[i + 1]);}
        if (s_lam_method.compare(argv[i]) == 0){lam_method_in = argv[i + 1];}
        if (s_x_method.compare(argv[i]) == 0){x_method_in = argv[i + 1];}
    }

    //dir_out.replace(dir_out.end(),dir_out.end(),"/","");
    // convert directory name to char_array
    int n_char = dir_out.length();
    // declaring character array
    char char_array[n_char + 1];
    // copying the contents of the
    // string to char array
    strcpy(char_array, dir_out.c_str());
    DIR *dir = opendir(char_array);
    if (dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else
    {
        /* Directory does not exist. */
        printf("Can't open results directory, stopping. \n");
        exit(0);
    }

    // calculate the sample size and the gene numbers
    cal_y_dimensions(file_y, sep, s_n, d_y);

    // write command into file for later references, also write the dimension of the gene expression matrix
    ss.str("");
    ss.clear();
    ss << dir_out << "/command.txt";
    ofstream f_com(ss.str().c_str());
    if (f_com.is_open())
    {
        for (int i = 0; i < argc; i++)
        {
            f_com << argv[i] << " ";
        }
        f_com << endl;
    }
    f_com << "Y matrix has dimension of " << s_n << " by " << d_y << endl
          << endl;
    f_com << "Starting analysis using factor number of " << nf << endl;
    //f_com << "Command is written in command.txt" << endl;
    f_com << "The total number of runs is set to " << n_itr << endl;
    f_com << "Results will be written for every " << write_itr << " iterations" << endl;
    f_com << "Convergence is reached if the total number of sparse elements remain unchanged for " << interval << " iterations" << endl;
    f_com.close();

    cout << "Starting analysis using factor number of " << nf << endl;
    cout << "Details of the run parameters can be found in command.txt" << endl;
    cout << "The total number of runs is set to " << n_itr << endl;
    cout << "Results will be written for every " << write_itr << " iterations" << endl;
    cout << "Convergence is reached if the total number of sparse elements remain unchanged for " << interval << " iterations" << endl;

    // read in the Y matrix
    MatrixXd Y = MatrixXd::Constant(s_n, d_y, 0);
    
    read_y(file_y, sep, Y);

    //cout << Y.block(0,0,10,10);
    cout << "a " << a << endl;
    cout << "b " << b << endl;
    cout << "c " << c << endl;
    cout << "d " << d << endl;
    cout << "g " << g << endl;
    cout << "h " << h << endl;

    cout << "lam_method " << lam_method_in << endl;
    cout << "x_method " << x_method_in << endl;

    VectorXd mu = VectorXd::Constant(s_n, 0);

    //cal_mu(Y, dir_out, s_n, d_y);


    Eigen::initParallel();

    double tol_in = (double)tol;

    long seed;
    //seed = time (NULL) * getpid();
    //seed = 1000;

    seed = (long)i_seed;

    ss.str("");
    ss.clear();
    ss << dir_out << "/seed";
    ofstream f_seed (ss.str().c_str());
    f_seed << seed << endl;
    f_seed.close();
    

    VectorXd PSI=VectorXd::Constant(s_n,1);
    VectorXd PSI_INV=VectorXd::Constant(s_n,1);
   
    int nt=nf;
    
    MatrixXd LAM=MatrixXd::Constant(s_n,nt,0);
    MatrixXd THETA=MatrixXd::Constant(s_n,nf,1);
    MatrixXd DELTA=MatrixXd::Constant(s_n,nf,1);
    VectorXd PHI = VectorXd::Constant(nf,1);
    VectorXd TAU = VectorXd::Constant(nf,1);
    
    double nu = 1;
    double ETA = 1;
    double GAMMA = 1;
    
    VectorXd count_lam = VectorXd::Constant(nf,0);
    VectorXd index = VectorXd::Constant(nf,0);
    
    double nmix=2;
    double zi = double(1)/nmix;
    MatrixXd Z = MatrixXd::Constant(nmix,nf,zi);
    MatrixXd logZ = MatrixXd::Constant(nmix,nf,log(zi));
    MatrixXd LOGV = MatrixXd::Constant(nmix,1,log(zi));
    
    
    // fill in the lambda matrix
    
    init_lam(LAM, s_n, nf, seed);
    
    //VectorXd lam_count_v = VectorXd::Constant(n_itr,0);
    
    //declare and initialize parameters related to X
    double XI=1,VARSIG=1,OMEGA=1;
    
    MatrixXd EX=MatrixXd::Constant(nt,d_y,0);
    //MatrixXd TEX=MatrixXd::Constant(d_y,nt,0);
    MatrixXd EXX=MatrixXd::Constant(nt,nt,0);
    
    VectorXd KAPPA = VectorXd::Constant(nf,1);
    VectorXd LAMX = VectorXd::Constant(nf,1);
    MatrixXd RHO=MatrixXd::Constant(nf,d_y,1);
    MatrixXd SIGMA=MatrixXd::Constant(nf,d_y,1);
    
    VectorXd count_x = VectorXd::Constant(nf,0);
    
    MatrixXd O = MatrixXd::Constant(nmix,nf,zi);
    MatrixXd logO = MatrixXd::Constant(nmix,nf,log(zi));
    MatrixXd LOGVO = MatrixXd::Constant(nmix,1,log(zi));
    
    MatrixXd LPL = MatrixXd::Constant(nf,nf,0);
    //MatrixXd vLXL = MatrixXd::Constant(s_n,s_n,0);
    //MatrixXd partR = MatrixXd::Constant(nf,d_y,0);
    //MatrixXd partL = MatrixXd::Constant(s_n,nf,0);
    
    // fill in the EX matrix
    init_ex(EX, d_y, nf, seed);
    
    EXX=EX*EX.transpose();
    
    //VectorXd x_count_v = VectorXd::Constant(n_itr,0);
    //MatrixXd LAM_T=LAM.transpose();
    
    LPL.setZero();
    for(int i=0;i<s_n;i++){
        LPL += LAM.transpose().col(i)*PSI_INV(i)*LAM.row(i);
    }
    
    VectorXd det_psi = VectorXd::Constant(n_itr,0);
    
    //cout << "You passed method " << method << endl;
    
    //string lam_method_in = lam_method;
    
    for(int itr=0;itr<=n_itr;itr++){
        
        cal_lam_all(LAM, Y,EX,PSI_INV,EXX,Z,LPL,THETA,PHI, s_n, d_y, nf,
                    a, b, c, d, g, h, GAMMA, ETA, nu, TAU, DELTA, alpha, beta, lam_method_in);
        
        cal_z(logZ, LOGV, LAM,  THETA,  DELTA, PHI,  Z, nf, s_n, a, b, alpha, beta,true);
        
        //if(x_method.compare("bicmix") == 0){
        //cout << "You passed bicmix" << endl;
        cal_ex_all( LAM,  Y, EX, PSI_INV, EXX, O, LPL, SIGMA, LAMX,  s_n,  d_y,  nf,
                   a,  b,  c,  d,  g,  h,  VARSIG,  OMEGA,  XI,  KAPPA,  RHO,  logO,  LOGVO,  alpha,  beta, x_method_in);
        //}
        //if(method.compare("sfamix") == 0){
        //cal_ex_simple(EX, EXX, LAM, PSI_INV, Y, s_n, nf, d_y);
        //}
        
        
        
        // count the number of non-zero values in each row of the x matrix
        int count_nonzero = 0;
        count_nonzero = count_nonzero_lam_ex(count_x, count_lam, index, LAM, EX, PHI, LAMX, nf, s_n, d_y);
        
        // remove factors, loadings that are exclusively zeros, and assign to new matrix
        if(count_nonzero != nf){
            nf=count_nonzero;
            nt=nf;
            //red_dim(EX, EXX, LAM, THETA, DELTA, PHI, TAU, Z, logZ, count_lam, index, KAPPA, LAMX,
            //        RHO, SIGMA, count_x, O, logO, LPL, partR, partL, nf, nt, s_n, d_y, nmix, zi);
            red_dim(EX, EXX, LAM, THETA, DELTA, PHI, TAU, Z, logZ, count_lam, index, KAPPA, LAMX,
                    RHO, SIGMA, count_x, O, logO, LPL, nf, nt, s_n, d_y, nmix, zi,true);
            //red_dim(EX_merge, EXX_merge, LAM_merge, THETA_merge, DELTA_merge, PHI_merge, TAU_merge, Z_merge, logZ_merge, count_lam, index, KAPPA_merge, LAMX_merge,
            //       RHO_merge, SIGMA_merge, count_x, O_merge, logO_merge, LPL_merge, nfcov, nt, s_n, d_y, nmix, zi,false);
            
        }
        
        MatrixXd LX=LAM*EX;
        cal_psi(PSI, LX, Y, LAM, EXX, s_n, d_y);
        inv_psi_vec(PSI,PSI_INV,s_n);
        
        for(int i=0;i<s_n;i++){
            det_psi(itr) = det_psi(itr)+log(PSI(i));
        }
        
        if(itr%10==0){
            cout << "itr " << itr << endl;
            
            cout << "number of factors " << nf << endl;
            cout << "count_lam" << endl << count_lam.transpose() << endl;
            //cout << "number of betas " << ncov << endl;
            //cout << "count_beta" << endl << count_lam_cov.transpose() << endl;
            cout << "count_x" << endl << count_x.transpose() << endl;
        }
        if(itr % write_itr == 0){
            //write_final_hidden(dir_out, LAM, Z, EX, EXX, O, PSI, lam_count_v, LAMX, PHI, itr, seed);

            // write LAM
            string sin;
            sin = dir_out + "/LAM_" + std::to_string(itr);
            write_file <MatrixXd> (LAM,sin);
            
            sin = dir_out + "/Z_" + std::to_string(itr);
            write_file <MatrixXd> (Z,sin);
            
            sin = dir_out + "/EX_" + std::to_string(itr);
            write_file <MatrixXd> (EX,sin);
            
            sin = dir_out + "/EXX_" + std::to_string(itr);
            write_file <MatrixXd> (EXX,sin);
            
            sin = dir_out + "/PSI_" + std::to_string(itr);
            write_file <VectorXd> (PSI,sin);

		if(x_method_in.compare("sparse") == 0){
			sin = dir_out + "/O_" + std::to_string(itr);
            		write_file <MatrixXd> (O,sin);
            
		}
        }
        if(itr>10){
            if(abs(det_psi(itr) - det_psi(itr-1)) < tol_in){
                //itr_at = itr;
                //write_final_beta(dir_out, LAM_cov, PSI, itr, seed);
                //write_final_hidden(dir_out, LAM, Z, EX, EXX, O, PSI, lam_count_v, LAMX, PHI, itr, seed);
                string sin = dir_out + "/itr";
                write_file <int> (itr,sin);
                
                // write LAM
                sin = dir_out + "/LAM";
                write_file <MatrixXd> (LAM,sin);
                
                sin = dir_out + "/Z";
                write_file <MatrixXd> (Z,sin);
                
                sin = dir_out + "/EX";
                write_file <MatrixXd> (EX,sin);
                
                sin = dir_out + "/EXX";
                write_file <MatrixXd> (EXX,sin);
                
                sin = dir_out + "/PSI";
                write_file <VectorXd> (PSI,sin);

		if(x_method_in.compare("sparse") == 0){
			sin = dir_out + "/O";
            		write_file <MatrixXd> (O,sin);
            
		}
                
                break;
            }
        }
        
        
    }
    
    
}

