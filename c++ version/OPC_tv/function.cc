#include "function.h"

cv::Mat loadMat(int rows, int cols, std::string& filepath){
    int i, j;
    cv::Mat matrix = cv::Mat::zeros(cv::Size(rows, cols), CV_32F);
    std::ifstream infile;
    // the type of parameter of open() is const char*
    // using .c_str() to transfer string to const char*
    infile.open(filepath.c_str(),std::ios::in);
    if(!infile.is_open())
    {
        std::cout << "no such file!" << std::endl;
        return matrix;
    }
    for(i=0;i!=rows;i++){
        for(j=0;j!=cols;j++){
            infile >> matrix.at<float>(i,j);
        }
    }
    infile.close();
    return matrix;
}

// generate Gaussian_kernel, the sama as fspecial('gaussian',size,sigma) in matlab
cv::Mat Gaussian_kernel(int size, int sigma){
    int m = size / 2;
    cv::Mat kernel(size, size, CV_32FC1);
    float s = 2 * sigma * sigma;
    float sum = 0.0;
    int i, j;
    for(i = 0; i != size; i++){
        for(j = 0; j !=size; j++){
            int x = i - m;
            int y = j - m;
            kernel.ptr<float>(i)[j] = exp(-(x*x + y*y) / s);
            sum += kernel.ptr<float>(i)[j];
        }
    }
    for(i = 0; i != size; i++){
        for(j = 0; j !=size; j++){
            kernel.ptr<float>(i)[j] /= sum ;
        }
    }

    return kernel;
}

// compare cv:Mat and a float number
// if the element in corresponding position is 
// not smaller than the number
// then the element is assigned to 1;
// else, it's assigned to 0;
// return a matrix with 0 and 1;
cv::Mat NotSmaller(cv::Mat m, float num){
    cv::Mat result = cv::Mat::zeros(cv::Size(m.rows,m.cols),CV_32FC1);
    int i, j;
    for(i=0;i!=m.rows;i++){
        for(j=0;j!=m.cols;j++){
            result.at<float>(i,j) = (m.at<float>(i,j) >= num) ? 1 :0;
        }
    }
    return result;
}

// m[i,j] <= num, return 1, else return 0
cv::Mat NotLarger(cv::Mat m, float num){
    cv::Mat result = cv::Mat::zeros(cv::Size(m.rows,m.cols),CV_32FC1);
    int i, j;
    for(i=0;i!=m.rows;i++){
        for(j=0;j!=m.cols;j++){
            result.at<float>(i,j) = (m.at<float>(i,j) <= num) ? 1 :0;
        }
    }
    return result;
}

// return 1 only if m(i,j) = n(i,j) = 1, else return 0
cv::Mat AndMat(cv::Mat m, cv::Mat n){
    if((m.rows!= n.rows) || (m.cols != n.cols)){
        std::cout << "the two matrix has different sizes!" << std::endl;
        exit(0); // exit procedure
    }
    else{
        cv::Mat result = cv::Mat::zeros(cv::Size(m.rows,m.cols), CV_32FC1);
        int i, j;
        for(i=0; i!=m.rows;i++){
            for(j=0;j!=m.cols;j++){
                result.at<float>(i,j) = ((m.at<float>(i,j) == 1) && 
                                         (n.at<float>(i,j) == 1))? 1 : 0;
            }
        }
        return result;
    }
}
// cos() of matrix
cv::Mat cos_mat(cv::Mat m){
    int i, j;
    cv::Mat result = cv::Mat::zeros(cv::Size(m.rows,m.cols),CV_32FC1);
    for(i=0;i!=m.rows;i++){
        for(j=0;j!=m.cols;j++){
            result.at<float>(i,j) = cos(m.at<float>(i,j));
        }
    }
    return result;
}

// sin() of matrix
cv::Mat sin_mat(cv::Mat m){
    int i, j;
    cv::Mat result = cv::Mat::zeros(cv::Size(m.rows,m.cols),CV_32FC1);
    for(i=0;i!=m.rows;i++){
        for(j=0;j!=m.cols;j++){
            result.at<float>(i,j) = sin(m.at<float>(i,j));
        }
    }
    return result;
}

// plot line graphs
void plot_lines(std::string title,float *arr, int arr_len,
                int interval_x, int interval_y){
    cv::Mat plot = cv::Mat::zeros(cv::Size(640,640), CV_32FC1);
    plot.setTo(255); // set background color as white
    
    int i , point_space;
    float min_num, max_num;
    min_num = min_value(arr, arr_len);
    max_num = max_value(arr, arr_len);
    float *arr_norm = new float[arr_len];
    if(min_num == max_num){
        for(i=0;i!=arr_len;i++){
            *(arr_norm + i) = 320;
        }
    }
    else{
        for(i=0;i!=arr_len;i++){
            // normalized to [30, 600]
            *(arr_norm + i) = *(arr + i) * 570.0 / max_num + 30.0; 
                              
        }
    }

    point_space = 560 / (arr_len-1); // [50,600]

    std::vector<cv::Point> points;
    for(i=0;i!=arr_len;i++){
        // cv::Point(i, *(arr+i));
        // flip up and down
        cv::Point p1(50 + i*point_space, 640 - *(arr_norm+i));
        points.push_back(p1);        
    }
    cv::polylines(plot,points,false,cv::Scalar(0,255,0),2);
    cv::namedWindow(title,0);
    // draw coordinate axis
    std::vector<cv::Point> vertex;
    vertex.push_back(cv::Point(50,30));
    vertex.push_back(cv::Point(50,610));
    vertex.push_back(cv::Point(610,610));
    vertex.push_back(cv::Point(610,30));
    cv::polylines(plot,vertex,true,cv::Scalar(0,255,0),1);
    
    std::vector<cv::Point> axis_x; // horizontal coordinate position
    std::vector<cv::Point> axis_y; // vertical coordinate position

    int x, y;
    int x_num = arr_len / interval_x + 2;
    int y_num = max_num / interval_y + 2;
    // plot horizontal axis
    for(x=0;x!=x_num;x++){
        //char value[10];
        std::string value_x = std::to_string(x*interval_x);
        axis_x.push_back(cv::Point(x*interval_x*point_space + 40,630));
        cv::putText(plot, value_x,axis_x[x],cv::FONT_ITALIC,
                    0.5,cv::Scalar(0,0,255),1,8);

        cv::Point p0(x*interval_x*point_space+50,610);
        cv::Point p2(x*interval_x*point_space+50, 615);
        cv::line(plot,p0, p2, cv::Scalar(0,255,255),1);
    }

    // plot vertical axis
    for(y=1;y!=y_num;y++){
        std::string value_y = std::to_string(y*interval_y);
        
        float Py_y = 640 - y*interval_y*570.0/max_num - 30;
        
        // value
        cv::Point Py(10,Py_y + 5);
        cv::putText(plot,value_y,Py,cv::FONT_ITALIC,
                    0.5,cv::Scalar(0,0,255),1,8);
        // short lines
        cv::Point Py_1(45,Py_y);
        cv::Point Py_2(50,Py_y);
        cv::line(plot,Py_1,Py_2,cv::Scalar(0,0,255),1);
    }
    cv::imshow(title,plot);
    cv::waitKey(0);
    delete [] arr_norm;
}

// find the maximum vaule of an array
template <typename T>
T max_value(const T *arr, int len){
    T result = 0;
    for(int i =0; i!=len;i++){
        if (*(arr+i)>result){
            result = *(arr+i);
        }
    }
    return result;
}

// find the minimum vaule of an array
template <typename T>
T min_value(const T *arr, int len){
    T result = 0;
    for(int i =0; i!=len;i++){
        if (*(arr+i)<result){
            result = *(arr+i);
        }
    }
    return result;
}

// sign() return 1 is x > 0;
// return 0 if x == 0;
// return -1 if x < 0;
cv::Mat sign_mat(cv::Mat m){
    cv::Mat result = cv::Mat::zeros(cv::Size(m.rows,m.cols),CV_32FC1);
    for(int i=0;i!=m.rows;i++){
        for(int j=0;j!=m.cols;j++){
            if(m.at<float>(i,j) > 0){
                result.at<float>(i,j) = 1;
            }
            else{
                if(m.at<float>(i,j) == 0){
                    result.at<float>(i,j) = 0;
                }
                else{
                    result.at<float>(i,j) = -1;
                }

            }
        }
    }
    return result;
}
