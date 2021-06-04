#ifndef FUNCTION_H
#define FUNCTION_H
#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <string>
#include <cmath>

#define N 80
#define maxloop 27
#define min_error 44 
#define window_size 640
// read files to matrix
cv::Mat loadMat(int rows, int cols, std::string& filepath);

// generate Gaussian kernel
cv::Mat Gaussian_kernel(int size, int sigma);

// compare cv:Mat and a float number
// if the element in corresponding position is 
// not smaller than the number
// then the element is assigned to 1;
// else, it's assigned to 0;
// return a matrix with 0 and 1;
cv::Mat NotSmaller(cv::Mat m, float num);

// m[i,j] <= num, return 1, else return 0
cv::Mat NotLarger(cv::Mat m, float num);

// operator & between two Mat
cv::Mat AndMat(cv::Mat m, cv::Mat n);

// cos() of matrix
cv::Mat cos_mat(cv::Mat m);

// sin() of matrix
cv::Mat sin_mat(cv::Mat m);

// plot line graphs
void plot_lines(std::string title, float *arr, int arr_len, 
                int interval_x, int interval_y);

// find the maximum vaule of an array
template <typename T>
T max_value(const T *arr, int len);

// find the minimum vaule of an array
template <typename T>
T min_value(const T *arr, int len);

// sign() return 1 is x > 0;
// return 0 if x == 0;
// return -1 if x < 0;
cv::Mat sign_mat(cv::Mat m);
#endif
