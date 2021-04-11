# 基于离散罚函数和小波罚函数的相移掩膜优化算法 （GPSM_wa）
**Generalized Phase-Shift Mask Optimazation Algorithm based on discretization penalty and wavelet penalty**

GPSM_wa算法在相干照明系统中对N×N的目标版图做广义梯度的相移掩膜优化。

该算法产生优化后的四相或二相相移掩膜，并包含了离散罚函数及局部小波罚函数。

掩膜板上的不同区域可以被赋予不同的小波罚函数权重。

若所有区域的权重均为1，则为全局小波罚函数。

## 1. Matlab版代码及说明

> *function [] = GPSM_wa(N, pz, ra, phase_n, s_phi, s_theta, a, t_r, t_m, gamma_r_D, gamma_a_D, gamma_r_WA, gamma_a_WA, scale, epsilon, maxloop);*

函数一共有16个参数，每个参数的含义如下：

*N*： Mask的维度；

*pz*：目标输出版图；

*ra*：Mask的初始相位版图；

*phase_n*：掩膜的不同相位的数量，在这里为2或4；

*s_phi*：Mask幅度版图优化的步长；

*s_theta*：Mask相位版图优化的步长；

*a*：sigmoid函数的陡度（steepness）；

*t_r*：sigmoid函数的处理阈值；

*t_m*：Mask的全局阈值；

*gamma_r_D*：Mask幅度版图的相应离散罚函数权重；

*gamma_a_D*: Mask相位版图的相应离散罚函数权重；

*gamma_r_WA*：Mask幅度版图的相应小波罚函数权重；

*gamma_a_WA*：Mask相位版图的相应小波罚函数权重；

*scale*：局部小波罚函数的区域性权重；

*epsilon*：输出版图的最大可容忍误差；

*maxloop*：最大迭代次数；
> *dr=zeros(N,N);   %Gradient of the cost function corresponding to \phi*

> *da=zeros(N,N);   %Gradient of the cost function corresponding to \theta*

> *dr_D=zeros(N,N);   %Gradient of the discretization penalty corresponding to \phi*

> *da_D=zeros(N,N);   %Gradient of the discretization penalty corresponding to \theta*

> *dr_WA=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \phi*

> *da_WA=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \theta*

> *convergence=zeros(maxloop,1);   %Output pattern error in each iteration*

> *count=0;   %Index of iteration number*

> *sum6=100;   %Output pattern error corresponding to the optimized pole-level mask*

> *sum8=100;   %Output pattern error corresponding to the optimized complex-valued mask*

以上为基本变量的定义，包括梯度，每次迭代后的误差，迭代次数和输出版图误差

> h = fspecial('gaussian',11, 14);

h: 振幅冲激响应函数，fspecial('gaussian', hsize, sigma)返回大小为hsize且标准偏差为sigma（正）的旋转对称高斯低通滤波器。

hsize可以是一个向量，指定h中的行数和列数，也可以是一个标量，在这种情况下，h是一个方阵。

因此h为11×11的标准差为14的高斯型低通滤波器矩阵。

h的物理意义为光学系统对于输入的Mask的响应，通过Mask与h做卷积，即可得到光学系统的输出。

> for ii = 1:11
>     for j = 1:11
>         h1((ii-1)*11+j) = h(ii, j);
>     end
> end
> for ii = 1:11
>     for j=1:11
>         g(ii,j)=h1((11-ii)*11+(12-j));
>     end
> end
