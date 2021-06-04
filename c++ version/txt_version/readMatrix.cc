#include <iostream>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <fstream>
using namespace cv;
using namespace std;
int main()
{
    int i ,j;
    Mat pz_t = Mat::zeros(Size(50,50),CV_32F);
    ifstream infile;
    infile.open("pz_t.txt");
    for(i=0;i!=50;i++){
        for(j=0;j!=50;j++){
            infile >> pz_t.at<float>(i,j);
        }
    }
    //cout << pz_t << endl;
    infile.close();
    // show pz_t
    Mat pz_t_converted;
    pz_t.convertTo(pz_t_converted,CV_32F,256);
    namedWindow("pz_t", 0);
    resizeWindow("pz_t", 640, 640);
    imshow("pz_t", pz_t_converted);
    waitKey();
    
    return 0;
}

