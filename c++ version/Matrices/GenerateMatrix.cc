#include <iostream>
#include <opencv2/opencv.hpp>
#include <math.h>
using namespace cv;
using namespace std;
int main()
{
    // generate matrix pz_f
    Mat pz_f = Mat::zeros(Size(80, 80), CV_8UC1);
    int i, j;
    for(j=20;j!=60;j++){
        for(i=12;i!=22;i++){
            pz_f.at<uchar>(i,j) = 1;
        }
        for(i=27;i!=37;i++){
            pz_f.at<uchar>(i,j) = 1;
        }
        for(i=42;i!=52;i++){
            pz_f.at<uchar>(i,j) = 1;
        }
        for(i=57;i!=67;i++){
            pz_f.at<uchar>(i,j) = 1;
        }

    }
    
    // show pz_f
    Mat pz_f_converted;
    //convertTo(dst_mat, type, scale, shift)
    pz_f.convertTo(pz_f_converted, CV_8UC1,256);
    namedWindow("pz_f", 0);
    resizeWindow("pz_f", 640,640);
    imshow("pz_f", pz_f_converted);
    
    // preserce pz_f as pz_f.xml
    FileStorage fs("pz_f.xml",FileStorage::WRITE);
    fs << "pz_f" <<pz_f;
    fs.release();
    // waitKey();
    
    // generate ra_f
    Mat ra_f = Mat::zeros(Size(80,80), CV_32F);
    for(j=0;j!=80;j++){
        for(i=0;i!=24;i++){
            ra_f.at<float>(i,j) = M_PI * 6.0 / 5.0;
        }
        for(i=24;i!=39;i++){
            ra_f.at<float>(i,j) = M_PI / 5.0;
        }
        for(i=39;i!=54;i++){
            ra_f.at<float>(i,j) = M_PI * 6.0 / 5.0;
        }
        for(i=54;i!=80;i++){
            ra_f.at<float>(i,j) = M_PI / 5.0;
        }
    }

    // show ra_f
    Mat ra_f_converted;
    ra_f.convertTo(ra_f_converted, CV_32F, 1.0);
    namedWindow("ra_f", 0);
    resizeWindow("ra_f",640,640);
    imshow("ra_f", ra_f_converted);
    
    // preserve ra_f
    FileStorage fs_ra_f("ra_f.xml", FileStorage::WRITE);
    fs_ra_f << "ra_f" << ra_f;
    fs_ra_f.release();
    //waitKey();

    // generate matrix pz_u
    Mat pz_u = Mat::zeros(Size(80, 80), CV_8UC1);
    for(j=20;j!=60;j++){
        for(i=20;i!=35;i++){
            pz_u.at<uchar>(i,j) = 1;
        }
    }
    for(i=35;i!=60;i++){
        for(j=20;j!=36;j++){
            pz_u.at<uchar>(i,j) = 1;
        }
        for(j=45;j!=60;j++){
            pz_u.at<uchar>(i,j) = 1;
        }
    }
    // show pz_u
    Mat pz_u_converted;
    //convertTo(dst_mat, type, scale, shift)
    pz_u.convertTo(pz_u_converted, CV_8UC1,256);
    namedWindow("pz_u", 0);
    resizeWindow("pz_u", 640,640);
    imshow("pz_u", pz_u_converted);
    
    // preserce pz_f as pz_f.xml
    FileStorage fs_pz_u("pz_u.xml",FileStorage::WRITE);
    fs_pz_u << "pz_u" << pz_u;
    fs_pz_u.release();
    //waitKey();
    
    // generate ra_u
    Mat ra_u = Mat::ones(Size(80,80), CV_32F);
    ra_u  = ra_u * M_PI * 7.0 / 4.0;
    for(j=20;j!=60;j++){
        for(i=20;i!=35;i++){
            ra_u.at<float>(i,j) = M_PI / 4.0;
        }
    }
    for(j=0;j!=20;j++){
        for(i=35;i!=80;i++){
            ra_u.at<float>(i,j) = M_PI * 5.0 / 4.0;
        }
    }
    for(j=20;j!=35;j++){
        for(i=35;i!=60;i++){
            ra_u.at<float>(i,j) = M_PI * 3.0 / 4.0;
        }

        for(i=60;i!=80;i++){
            ra_u.at<float>(i,j) = M_PI * 5.0 / 4.0;
        }

    }
    for(j=35;j!=40;j++){
        for(i=35;i!=80;i++){
            ra_u.at<float>(i,j) = M_PI * 5.0 / 4.0;
        }
    }
    for(j=40;j!=45;j++){
        for(i=35;i!=80;i++){
            ra_u.at<float>(i,j) = M_PI / 4.0;
        }
    }
    for(j=45;j!=60;j++){
        for(i=60;i!=80;i++){
            ra_u.at<float>(i,j) = M_PI / 4.0;
        }

    }
    for(j=60;j!=80;j++){
        for(i=35;i!=80;i++){
            ra_u.at<float>(i,j) = M_PI / 4.0;
        }
    }
    // show ra_u
    Mat ra_u_converted;
    ra_u.convertTo(ra_u_converted, CV_32F, 1.0);
    namedWindow("ra_u", 0);
    resizeWindow("ra_u",640,640);
    imshow("ra_u", ra_u_converted);
    
    // preserve ra_u
    FileStorage fs_ra_u("ra_u.xml", FileStorage::WRITE);
    fs_ra_u << "ra_u" << ra_u;
    fs_ra_u.release();
    //waitKey();
    
    // generate scale2
    Mat scale2 = Mat::ones(Size(80,80), CV_32F);
    scale2 = scale2 * 0.7;
    for(j=35;j!=45;j++){
        for(i=35;i!=60;i++){
            scale2.at<float>(i,j) = 1.6;
        }
    }

    // show scale2
    Mat scale2_converted;
    scale2.convertTo(scale2_converted, CV_32F, 1.0);
    namedWindow("scale2", 0);
    resizeWindow("scale2",640,640);
    imshow("scale2", scale2_converted);
    
    // preserve ra_f
    FileStorage fs_scale2("scale2.xml", FileStorage::WRITE);
    fs_scale2 << "scale2" << scale2;
    fs_scale2.release();
    waitKey();
    /* read matrix
    FileStorage fs("pz_f.xml", FileStorage::READ);
    Mat readed_mat;
    fs["pz_f"] >> readed_mat;
    fs.release();
`   */
    return 0;
}

