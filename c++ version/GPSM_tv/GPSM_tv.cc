#include <iostream>
#include <opencv2/opencv.hpp>
#include <string>
#include <fstream>
#include <cmath>
#include "function.h"
using namespace cv;
using namespace std;

int main()
{
    int i, j;

    // define parameters
    int phase_n = 4; // number of mask phase
    float s_phi = 2; // step of amplitude mask optimization
    float s_theta = 0.01; // step of phase mask optimization
    float a = 80; // steepness of sigmoid function
    float t_r = 0.5; // threshold of sigmoid function
    float t_m = 0.5; // threshold of discretization of mask
    float gamma_r_D = 0.01; // weight of discrete penalty function of amplitude mask
    float gamma_a_D = 0.001; // weight of discrete penalty function of phase mask
    float gamma_r_TV = 0.1; // weight of wavelet penalty function of amplitude mask
    float gamma_a_TV = 0.001; // weight of wavelet penalty function of phase mask

    Mat pz; // desired output pattern
    string filepath_pz="../Matrices/txt_version/pz_u.txt";
    pz = loadMat(N,N,filepath_pz);

    Mat ra; // initial phase pattern
    string filepath_ra="../Matrices/txt_version/ra_u.txt";
    ra = loadMat(N,N,filepath_ra);
    
    // show pz
    Mat pz_norm;
    pz.convertTo(pz_norm,CV_32F,256);
    namedWindow("pz", 0);
    resizeWindow("pz", window_size,window_size);
    imshow("pz", pz_norm);
    //waitKey();
    
    // Gradient matrix of cost function of \phi(amplitude)
    Mat dr = Mat::zeros(Size(N,N), CV_32F);
    
    // Gradient matrix of cost function of \theta(phase)
    Mat da = Mat::zeros(Size(N,N), CV_32F);

    // Gradient matrix of discretization penalty of \phi(amplitude)
    Mat dr_D = Mat::zeros(Size(N,N), CV_32F);
    
    // Gradient matrix of discretization penalty of \theta(phase)
    Mat da_D = Mat::zeros(Size(N,N), CV_32F);
    
    // Gradient matrix of wavelet penalty of \phi(amplitude)
    Mat dr_TV = Mat::zeros(Size(N,N), CV_32F);

    // Gradient matrix of wavelet penalty of \theta(phase)
    Mat da_TV = Mat::zeros(Size(N,N), CV_32F);
   
    // output pattern error in each iteration
    float convergence[maxloop];
    
    // Index of iteration number
    int count = 0;
    
    // Output pattern error of optimized pole-level mask
    float sum6 = 100;

    // Output pattern error of optimized complex-valued mask
    float sum8 = 100;

    // amplitude impulse response of coherent imaging system
    Mat h;
    h = Gaussian_kernel(11, 14);
    
    // rotate matrix h clockwise 180 degrees to get matrix g
    // because h is rotationally symmertirc, g is the same as h
    Mat g = h;

    // initialize amplitude matrix, r = \phi
    Mat rr = Mat::ones(Size(N,N), CV_32FC1);
    for(i=0;i!=N;i++){
        for(j=0;j!=N;j++){
            rr.at<float>(i,j) = pz.at<float>(i,j) == 0? M_PI*4.0/5.0 : M_PI/5.0;
        }
    }
    
    // maske pattern
    Mat m = Mat::zeros(Size(N,N), CV_32FC1);

    // initial errors
    float cun = 1000;
    
    // dr,da
    Mat delta_r = Mat::zeros(Size(N,N), CV_32FC1); // dr * s_phi
    Mat delta_a = Mat::zeros(Size(N,N), CV_32FC1); // da * s_theta
    
    // cos(rr),cos(ra),sin(ra)
    Mat cos_rr = Mat::zeros(Size(N,N), CV_32FC1); // cos(rr)
    Mat cos_ra = Mat::zeros(Size(N,N), CV_32FC1); // cos(ra)
    Mat sin_ra = Mat::zeros(Size(N,N), CV_32FC1); // sin(ra)
    
    // real part of m, imaginary part of m
    Mat mr = Mat::zeros(Size(N,N), CV_32FC1); // real part of m
    Mat mi = Mat::zeros(Size(N,N), CV_32FC1); // imaginary part of m
    
    // abs(m)
    Mat mmo = Mat::zeros(Size(N,N), CV_32FC1); // Amplitude pattern of m
    Mat mr_squared = Mat::zeros(Size(N,N), CV_32FC1); // mr^2
    Mat mi_squared = Mat::zeros(Size(N,N), CV_32FC1); // mi^2
    
    // output pattern of optimized pole-level mask
    Mat viccbin = Mat::zeros(Size(N,N), CV_32FC1);
    
    // mainloop
    while ((sum6 > min_error) && (count < maxloop)){
        dr.convertTo(delta_r,CV_32FC1,s_phi); // equals to dr*s_phi
        da.convertTo(delta_a,CV_32FC1,s_theta); // equals to da * s_theta
        rr = rr - delta_r;
        ra = ra - delta_a;
        
        cos_rr = cos_mat(rr);
        cos_ra = cos_mat(ra);
        sin_ra = sin_mat(ra);
        
        cos_rr += 1;
        // m = mr + i * mi = 0.5*(1+cos(rr))*exp(i*ra)
        mr = cos_rr.mul(cos_ra);
        mi = cos_rr.mul(sin_ra);
        mr.convertTo(mr,CV_32FC1,0.5); 
        mi.convertTo(mi,CV_32FC1,0.5);
        pow(mr,2,mr_squared);
        pow(mi,2,mi_squared);
        sqrt(mr_squared + mi_squared, mmo);
        
        // four-phase PSM
        if (phase_n == 4){
            // pole-level mask, discretize the value to 0 or 1
            Mat viccone = NotSmaller(mmo, t_m);        
            // phase of pi/4
            Mat vicctwo = AndMat(NotSmaller(mr, 0), NotSmaller(mi, 0));
            Mat vicctwo_real, vicctwo_imag;
            // real part of vicctwo * exp(i*pi/4)
            vicctwo.convertTo(vicctwo_real, CV_32FC1,cos(M_PI/4.0));
            // imaginary part of vicctwo * exp(i*pi/4)
            vicctwo.convertTo(vicctwo_imag, CV_32FC1, sin(M_PI/4.0));
            
            // phase of pi*3/4
            Mat viccthree = AndMat(NotLarger(mr, 0), NotSmaller(mi, 0));
            Mat viccthree_real, viccthree_imag;
            // real part of viccthree * exp(i*pi*3/4)
            viccthree.convertTo(viccthree_real, CV_32FC1,cos(M_PI*3.0/4.0));
            // imaginary part of viccthree * exp(i*pi*3/4)
            viccthree.convertTo(viccthree_imag, CV_32FC1, sin(M_PI*3.0/4.0));
        
            // phase of pi*5/4
            Mat viccfour = AndMat(NotLarger(mr, 0), NotLarger(mi, 0));
            Mat viccfour_real, viccfour_imag;
            // real part of viccfour * exp(i*pi*5/4)
            viccfour.convertTo(viccfour_real, CV_32FC1,cos(M_PI*5.0/4.0));
            // imaginary part of viccfour * exp(i*pi*5/4)
            viccfour.convertTo(viccfour_imag, CV_32FC1, sin(M_PI*5.0/4.0));

            // phase of pi*7/4
            Mat viccfive = AndMat(NotSmaller(mr, 0), NotLarger(mi, 0));
            Mat viccfive_real, viccfive_imag;
            // real part of viccfive * exp(i*pi*7/4)
            viccfive.convertTo(viccfive_real, CV_32FC1,cos(M_PI*7.0/4.0));
            // imaginary part of viccfive * exp(i*pi*7/4)
            viccfive.convertTo(viccfive_imag, CV_32FC1, sin(M_PI*7.0/4.0));

            // phase pattern of mask pattern
            Mat viccsix_real, viccsix_imag;
            viccsix_real = vicctwo_real + viccthree_real + 
                           viccfour_real + viccfive_real;
            viccsix_imag = vicctwo_imag + viccthree_imag +
                            viccfour_imag + viccfive_imag;
            
            // pole-level mask pattern
            Mat viccin_real, viccin_imag;
            viccin_real = viccone.mul(viccsix_real);
            viccin_imag = viccone.mul(viccsix_imag);
            
            // caculate the convolution between viccin and h;
            Mat viccout_real, viccout_imag;
            filter2D(viccin_real,viccout_real,viccin_real.depth(),h);
            filter2D(viccin_imag,viccout_imag,viccin_imag.depth(),h);
            
            Mat viccout_real_squared, viccout_imag_squared;
            pow(viccout_real, 2.0, viccout_real_squared);
            pow(viccout_imag, 2.0, viccout_imag_squared);
            
            // amplitude of viccout
            Mat viccout_abs;
            sqrt(viccout_real_squared + viccout_imag_squared,viccout_abs);

            // output pattern of pole-level mask
            viccbin = NotSmaller(viccout_abs, t_r);
            
            // caculate the error between output pattern and desired pattern
            Mat diff_abs; // absolute value of the difference between 2 matrix
            absdiff(pz, viccbin, diff_abs);

            // caculate the sum of all elements of diff_abs
            sum6 = sum(diff_abs)[0]; // sum() return type is cv::Scalar

            convergence[count] = sum6;
        }
        // two-phase PSM
        else if(phase_n == 2){
            // pole-level mask, discretize the value to 0 or 1
            Mat viccone = NotSmaller(mmo, t_m);
            // phase of 0
            Mat vicctwo = NotSmaller(mr, 0);
            // phase of pi
            Mat viccthree = NotLarger(mr, 0);
            
            viccthree.convertTo(viccthree, CV_32FC1,-1.0);
            // Phase pattern of mask pattern
            Mat viccfour = vicctwo + viccthree;
            Mat viccin = viccone.mul(viccfour);
            
            // caculate the convolution between viccin and h;
            Mat viccout = Mat::zeros(Size(N,N), CV_32FC1);
            filter2D(viccin,viccout,viccin.depth(),h);
            Mat viccout_abs = abs(viccout);

            // output pattern of pole-level mask
            viccbin = NotSmaller(viccout_abs, t_r);
            Mat diff_abs = Mat::zeros(Size(N,N), CV_32FC1);
            absdiff(pz, viccbin, diff_abs);

            sum6 = sum(diff_abs)[0];

            convergence[count] = sum6;
        } 
        if (cun > sum6){
            cun = sum6;
        }
        
        count += 1; // increase iterations
        
        cout << "iterations: " << count 
             << ", errors: " << cun << endl;
        // convolution between continuous mask and filter h
        Mat mid1_real = Mat::zeros(Size(N,N), CV_32FC1);
        Mat mid1_imag = Mat::zeros(Size(N,N), CV_32FC1);
        filter2D(mr, mid1_real, mr.depth(),h);
        filter2D(mi, mid1_imag, mi.depth(),h);

        Mat mid1_real_squared = Mat::zeros(Size(N,N), CV_32FC1);
        Mat mid1_imag_squared = Mat::zeros(Size(N,N), CV_32FC1);
        pow(mid1_real,2.0,mid1_real_squared);
        pow(mid1_imag,2.0,mid1_imag_squared);
        
        // amplitude
        Mat mid1mo = Mat::zeros(Size(N,N), CV_32FC1);
        sqrt(mid1_real_squared + mid1_imag_squared, mid1mo);
        
        Mat z = Mat::zeros(Size(N,N), CV_32FC1);
        Mat matrix_ones = Mat::ones(Size(N,N),CV_32FC1);
        Mat mid1mo_mul_a = Mat::zeros(Size(N,N), CV_32FC1);
        Mat exp_mid_value = Mat::zeros(Size(N,N), CV_32FC1);
        // -1 * a * mid1mo
        mid1mo.convertTo(mid1mo_mul_a, CV_32FC1, -a); 
        // exp(-1*a*mid1mo + a*t_r)
        exp(mid1mo_mul_a + a*t_r, exp_mid_value);
        exp_mid_value = exp_mid_value + 1.0;
        
        z = matrix_ones / exp_mid_value;

        Mat mid3 = Mat::zeros(Size(N,N), CV_32FC1);
        // (pz-z) .* z .* (1-z) .* mid1_real .*(1./mid1mo)
        mid3 = ((((pz - z).mul(z)).mul(1-z)).mul(mid1_real)).mul(matrix_ones/mid1mo);
        
        Mat mid5 = Mat::zeros(Size(N,N), CV_32FC1);
        filter2D(mid3, mid5, mid3.depth(),g); 

        Mat mid7 = Mat::zeros(Size(N,N), CV_32FC1);
        // (pz-z) .* z .* (1-z) .* mid1_imag .*(1./mid1mo)
        mid7 = ((((pz - z).mul(z)).mul(1-z)).mul(mid1_imag)).mul(matrix_ones/mid1mo);
        
        Mat mid9 = Mat::zeros(Size(N,N), CV_32FC1);
        filter2D(mid7, mid9, mid7.depth(),g);

        // gradient of the discretization penalty corresponding to amplitude (phi)
        Mat sin_rr = Mat::zeros(Size(N,N),CV_32FC1);
        sin_rr = sin_mat(rr);
        dr_D = sin_rr.mul(1.0+cos_rr);
        dr_D.convertTo(dr_D,CV_32FC1,-0.5);
        
        // gradient of the discretization penalty corresponding to phase
        if(phase_n == 4){
            Mat four_mul_ra = Mat::zeros(Size(N,N),CV_32FC1);
            ra.convertTo(four_mul_ra, CV_32FC1,4.0);
            
            Mat cos_ra_pi = Mat::zeros(Size(N,N),CV_32FC1);
            cos_ra_pi = cos_mat(four_mul_ra - M_PI*1.5);
            
            Mat sin_ra_pi = Mat::zeros(Size(N,N),CV_32FC1);
            sin_ra_pi = sin_mat(four_mul_ra - M_PI*0.5);
            
            da_D = (sin_ra_pi + 1.0).mul(cos_ra_pi);
            da_D.convertTo(da_D,CV_32FC1,8.0);
        }
        else if (phase_n == 2){
            Mat two_mul_ra = Mat::zeros(Size(N,N),CV_32FC1);
            ra.convertTo(two_mul_ra, CV_32FC1,2.0);
            
            Mat cos_ra_pi = Mat::zeros(Size(N,N),CV_32FC1);
            cos_ra_pi = cos_mat(two_mul_ra - M_PI*0.5);
            
            Mat sin_ra_pi = Mat::zeros(Size(N,N),CV_32FC1);
            sin_ra_pi = sin_mat(two_mul_ra - M_PI*0.5);
            
            da_D = (sin_ra_pi + 1.0).mul(cos_ra_pi);
            da_D.convertTo(da_D,CV_32FC1,4.0);
        }
        
        // gradient of total variation penalty
        Mat f = Mat::zeros(Size(N,N),CV_32FC1);
        Mat m_pz_real = Mat::zeros(Size(N,N),CV_32FC1);
        absdiff(mr,pz,m_pz_real);

        pow(m_pz_real,2.0,m_pz_real);
        pow(mi,2.0,mi_squared);
        sqrt(m_pz_real + mi_squared, f); // f=abs(m-pz), m is a complex matrix
        // right shift of f
        Mat f_right = Mat::zeros(Size(N,N),CV_32FC1);
        for(i=0;i!=N;i++){
            for(j=1;j!=N;j++){
                f_right.at<float>(i,j) = f.at<float>(i,j-1);
            }
        }

        // up shift of f
        Mat f_up = Mat::zeros(Size(N,N),CV_32FC1);
        for(i=0;i!=N-1;i++){
            for(j=0;j!=N;j++){
                f_right.at<float>(i,j) = f.at<float>(i+1,j);
            }
        }
        
        Mat f1 = Mat::zeros(Size(N,N),CV_32FC1);
        f1 = sign_mat(f - f_right);
        for(i=0;i!=N;i++){
            f1.at<float>(i,0) = 0;
        }

        Mat f2 = Mat::zeros(Size(N,N),CV_32FC1);
        f2 = sign_mat(f - f_up);
        for(j=0;j!=N;j++){
            f2.at<float>(N-1,j) = 0;
        }
        // left shift of f1
        Mat f1_left = Mat::zeros(Size(N,N),CV_32FC1);
        for(i=0;i!=N;i++){
            for(j=0;j!=N-1;j++){
                f1_left.at<float>(i,j) = f1.at<float>(i, j+1);
            }
        }
        // down shift of f2
        Mat f2_down = Mat::zeros(Size(N,N),CV_32FC1);
        for(i=1;i!=N;i++){
            for(j=0;j!=N;j++){
                f2_down.at<float>(i,j) = f2.at<float>(i-1,j);
            }
        }
        
        Mat f11 = f1 - f1_left;
        for(i=0;i!=N;i++){
            f11.at<float>(i,N-1) = 0;
        }

        Mat f22 = f2 - f2_down;
        for(j=0;j!=N;j++){
            f22.at<float>(0,j) = 0;
        }

        Mat f3 = f11 + f22;

        // gradient of total variation penalty corresponding to amplitude
        Mat real_m_z_exp = Mat::zeros(Size(N,N),CV_32FC1);
        real_m_z_exp = (mr - z).mul(cos_ra) + mi.mul(sin_ra);

        Mat abs_m_z = Mat::zeros(Size(N,N),CV_32FC1); // abs(m-z)
        Mat m_z_real = Mat::zeros(Size(N,N),CV_32FC1);
        absdiff(mr,z,m_z_real);
        pow(m_z_real,2.0,m_z_real);
        pow(mi,2.0,mi_squared);
        sqrt(m_z_real + mi_squared, abs_m_z); 

        dr_TV = ((f3.mul(sin_rr)).mul(real_m_z_exp)) / abs_m_z;
        dr_TV.convertTo(dr_TV,CV_32FC1,-0.5);

        // gradient of total variation penalty corresponding to theta
        real_m_z_exp = mi.mul(cos_ra) - (mr - z).mul(sin_ra);
        da_TV = ((f3.mul(1.0 + cos_rr)).mul(real_m_z_exp)) / abs_m_z;
        da_TV.convertTo(da_TV,CV_32FC1,0.5);
        
        // gradient of overall cost function
        Mat dr_item1 = Mat::zeros(Size(N,N),CV_32FC1);
        Mat dr_item2 = Mat::zeros(Size(N,N),CV_32FC1);
        Mat dr_item3 = Mat::zeros(Size(N,N),CV_32FC1);
        Mat dr_item4 = Mat::zeros(Size(N,N),CV_32FC1);
        sin_rr.convertTo(sin_rr, CV_32FC1,a); // sin_rr = a*sin(rr)
        dr_item1 = (sin_rr.mul(cos_ra)).mul(mid5); 
        dr_item2 = (sin_rr.mul(sin_ra)).mul(mid9);
        dr_D.convertTo(dr_item3, CV_32FC1,gamma_r_D);
        dr_TV.convertTo(dr_item4,CV_32FC1,gamma_r_TV);

        dr = dr_item1 + dr_item2 + dr_item3 + dr_item4;

        Mat da_item1 = Mat::zeros(Size(N,N),CV_32FC1);
        Mat da_item2 = Mat::zeros(Size(N,N),CV_32FC1);
        Mat da_item3 = Mat::zeros(Size(N,N),CV_32FC1);
        Mat da_item4 = Mat::zeros(Size(N,N),CV_32FC1);
        
        cos_rr = cos_rr + 1.0;
        cos_rr.convertTo(cos_rr,CV_32FC1,a);
        da_item1 = (cos_rr.mul(sin_ra)).mul(mid5);
        da_item2 = (cos_rr.mul(cos_ra)).mul(mid9);
        da_D.convertTo(da_item3,CV_32FC1,gamma_a_D);
        da_TV.convertTo(da_item4,CV_32FC1,gamma_a_TV);

        da = da_item1 + da_item2 + da_item3 + da_item4;
    }

    // display part
    // output pattern of initial mask pz
    Mat output_pz = Mat::zeros(Size(N,N),CV_32FC1);
    filter2D(pz,output_pz,pz.depth(),h);
    output_pz = NotSmaller(output_pz, t_r);
    // scaled to 0-255
    output_pz.convertTo(output_pz,CV_32F,256);
    namedWindow("output pattern of mask pz", 0);
    resizeWindow("output pattern of mask pz", window_size,window_size);
    imshow("output pattern of mask pz", output_pz);
    
    // Magnitude of optimized comolex-valued mask
    Mat magnitude_mask = Mat::zeros(Size(N,N),CV_32FC1);
    magnitude_mask = cos_mat(rr) + 1.0;
    magnitude_mask.convertTo(magnitude_mask,CV_32FC1,0.5);
    namedWindow("Magnitude of optimized complex-valued mask",0);    
    resizeWindow("Magnitude of optimized complex-valued mask",window_size,window_size);
    imshow("Magnitude of optimized complex-valued mask",magnitude_mask);
    
    // Magnitude of optimized pole-level mask
    cos_rr = cos_mat(rr);
    cos_ra = cos_mat(ra);
    sin_ra = sin_mat(ra);
    cos_rr += 1;
    // m = mr + i * mi = 0.5*(1+cos(rr))*exp(i*ra)
    mr = cos_rr.mul(cos_ra);
    mi = cos_rr.mul(sin_ra);
    mr.convertTo(mr,CV_32FC1,0.5); 
    mi.convertTo(mi,CV_32FC1,0.5);
    pow(mr,2,mr_squared);
    pow(mi,2,mi_squared);
    sqrt(mr_squared + mi_squared, mmo);
    
    Mat pole_mask = NotSmaller(mmo, t_m);
    pole_mask.convertTo(pole_mask,CV_32FC1,256);
    namedWindow("optimized pole-level mask",0);
    resizeWindow("optimized pole-level mask",window_size,window_size);
    imshow("optimized pole-level mask",pole_mask);

    // output pattern of optimized pole-level mask
    Mat output_m = Mat::zeros(Size(N,N),CV_32FC1);
    viccbin.convertTo(output_m,CV_32FC1,256);
    namedWindow("output pattern of optimized pole-level mask",0);
    resizeWindow("output pattern of optimized pole-level mask",window_size,window_size);
    imshow("output pattern of optimized pole-level mask",output_m);
    
    // convergence plot
    string title = "convergence plot";
    plot_lines(title, convergence, count,6,20);
    waitKey(0);

    return 0;
}






















