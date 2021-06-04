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
    float step = 0.5; // step of gradient descent method
    float a = 80; // steepness of sigmoid function
    float t_r = 0.5; // threshold of sigmoid function
    float t_m = 0.5; // threshold of discretization of mask
    float gamma_D = 0.01; // weight of discrete penalty 
    float gamma_TV = 0.025; // weight of wavelet penalty

    Mat pz; // desired output pattern
    string filepath_pz="../Matrices/txt_version/pz_f.txt";
    pz = loadMat(N,N,filepath_pz);
    
    // show pz
    Mat pz_norm;
    pz.convertTo(pz_norm,CV_32F,256);
    namedWindow("pz", 0);
    resizeWindow("pz", window_size,window_size);
    imshow("pz", pz_norm);
    //waitKey();
    
    // Gradient matrix of cost function of
    Mat d = Mat::zeros(Size(N,N), CV_32F);    

    // Gradient matrix of discretization penalty
    Mat d_D = Mat::zeros(Size(N,N), CV_32F);
    
    // Gradient matrix of total variation penalty
    Mat d_TV = Mat::zeros(Size(N,N), CV_32F);

   
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
    Mat r = Mat::ones(Size(N,N), CV_32FC1);
    for(i=0;i!=N;i++){
        for(j=0;j!=N;j++){
            r.at<float>(i,j) = pz.at<float>(i,j) == 0? M_PI*4.0/5.0 : M_PI/5.0;
        }
    }
    
    // maske pattern
    Mat m = Mat::zeros(Size(N,N), CV_32FC1);

    Mat delta_d = Mat::zeros(Size(N,N), CV_32FC1); // dr * s_phi
    
    // cos(r),sin(r)
    Mat cos_r = Mat::zeros(Size(N,N), CV_32FC1); // cos(rr)
    Mat sin_r = Mat::zeros(Size(N,N), CV_32FC1); // sin(ra)   
    
    // output pattern of optimized pole-level mask
    Mat viccbin = Mat::zeros(Size(N,N), CV_32FC1);
    
    // mainloop
    while ((sum6 > min_error) && (count < maxloop)){
        d.convertTo(delta_d,CV_32FC1,step); // equals to dr*s_phi
        r = r - delta_d;
        cos_r = cos_mat(r);
        sin_r = sin_mat(r);
        
        cos_r += 1;
        // m = (1 + cos(r)) / 2
        cos_r.convertTo(m,CV_32FC1,0.5); 
        
        // binary mask
        Mat viccin = NotSmaller(m, t_m);
        Mat viccout = Mat::zeros(Size(N,N),CV_32FC1);
        filter2D(viccin,viccout,viccin.depth(),h);
        
        viccbin = NotSmaller(viccout, t_r);
        
        Mat pz_sub_viccbin = Mat::zeros(Size(N,N),CV_32FC1);
        absdiff(pz, viccbin, pz_sub_viccbin);
        sum6 = sum(pz_sub_viccbin)[0];
        convergence[count] = sum6;
        count += 1; // increase iterations
        cout << "iterations: " << count 
             << ", errors: " << sum6 << endl;
        
        // caculate z
        Mat mid1 = Mat::zeros(Size(N,N),CV_32FC1);
        // convolution between m and h
        filter2D(m,mid1,m.depth(),h);

        Mat z = Mat::zeros(Size(N,N), CV_32FC1);
        Mat matrix_ones = Mat::ones(Size(N,N),CV_32FC1);
        Mat mid1_mul_a = Mat::zeros(Size(N,N), CV_32FC1);
        Mat exp_mid_value = Mat::zeros(Size(N,N), CV_32FC1);
        // -1 * a * mid1mo
        mid1.convertTo(mid1_mul_a, CV_32FC1, -a); 
        // exp(-1*a*mid1mo + a*t_r)
        exp(mid1_mul_a + a*t_r, exp_mid_value);
        exp_mid_value = exp_mid_value + 1.0;
        z = matrix_ones / exp_mid_value;

        Mat mid3 = Mat::zeros(Size(N,N), CV_32FC1);
        // (pz-z) .* z .* (1-z) .* mid1_real .*(1./mid1mo)
        mid3 = ((pz - z).mul(z)).mul(1-z);
        
        Mat mid5 = Mat::zeros(Size(N,N), CV_32FC1);
        filter2D(mid3, mid5, mid3.depth(),g); 

        Mat mid7 = Mat::zeros(Size(N,N), CV_32FC1);
        

        // gradient of the discretization penalty corresponding to amplitude (phi)
        sin_r = sin_mat(r);
        m.convertTo(d_D,CV_32FC1,-8.0);
        d_D = d_D + 4.0;
        d_D.convertTo(d_D,CV_32FC1,-0.5);
        d_D = d_D.mul(sin_r);
     
        // gradient of total variation penalty
        Mat f = Mat::zeros(Size(N,N),CV_32FC1);
        absdiff(m,pz,f);

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

        d_TV = (f3.mul(sin_r)).mul(sign_mat(m-pz));
        d_TV.convertTo(d_TV,CV_32FC1,-0.5);
        
        // gradient of overall cost function
        Mat d_item1 = Mat::zeros(Size(N,N),CV_32FC1);
        Mat d_item2 = Mat::zeros(Size(N,N),CV_32FC1);
        Mat d_item3 = Mat::zeros(Size(N,N),CV_32FC1);
        sin_r.convertTo(sin_r, CV_32FC1,2*a); // sin_rr = 2*a*sin(rr)
        d_item1 = sin_r.mul(mid5); 
        d_D.convertTo(d_item2, CV_32FC1,gamma_D);
        d_TV.convertTo(d_item3,CV_32FC1,gamma_TV);

        d = d_item1 + d_item2 + d_item3;

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
    
    
    // Magnitude of optimized pole-level mask  
    Mat pole_mask = NotSmaller(m, t_m);
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






















