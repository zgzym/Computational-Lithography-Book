# 基于离散罚函数和小波罚函数的相移掩膜优化算法 （GPSM_wa）
**Generalized Phase-Shift Mask Optimazation Algorithm based on discretization penalty and wavelet penalty**

### 本文的matlab代码来源于《Computational Lithography》一书；

GPSM_wa算法在相干照明系统中对N×N的目标版图做广义梯度的相移掩膜优化。

该算法产生优化后的四相或二相相移掩膜，并包含了离散罚函数及局部小波罚函数。

掩膜板上的不同区域可以被赋予不同的小波罚函数权重。

若所有区域的权重均为1，则为全局小波罚函数。

## 1. Matlab版代码及说明

> function [] = GPSM_wa(N, pz, ra, phase_n, s_phi, s_theta, a, t_r, t_m, gamma_r_D, gamma_a_D, gamma_r_WA, gamma_a_WA, scale, epsilon, maxloop);

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
> dr=zeros(N,N);   %Gradient of the cost function corresponding to \phi

> da=zeros(N,N);   %Gradient of the cost function corresponding to \theta

> dr_D=zeros(N,N);   %Gradient of the discretization penalty corresponding to \phi

> da_D=zeros(N,N);   %Gradient of the discretization penalty corresponding to \theta

> dr_WA=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \phi

> da_WA=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \theta

> convergence=zeros(maxloop,1);   %Output pattern error in each iteration

> count=0;   %Index of iteration number

> sum6=100;   %Output pattern error corresponding to the optimized pole-level mask

> sum8=100;   %Output pattern error corresponding to the optimized complex-valued mask

以上为基本变量的定义，包括梯度，每次迭代后的误差，迭代次数和输出版图误差

> h = fspecial('gaussian',11, 14);

h: 振幅冲激响应函数，fspecial('gaussian', hsize, sigma)返回大小为hsize且标准偏差为sigma（正）的旋转对称高斯低通滤波器。

hsize可以是一个向量，指定h中的行数和列数，也可以是一个标量，在这种情况下，h是一个方阵。

因此h为11×11的标准差为14的高斯型低通滤波器矩阵。

h的物理意义为光学系统对于输入的Mask的响应，通过Mask与h做卷积，即可得到光学系统的输出。

> for ii = 1:11
> 
>     for j = 1:11
>     
>         h1((ii-1)*11+j) = h(ii, j);
>         
>     end
> end
> 
> for ii = 1:11
> 
>     for j=1:11
>     
>         g(ii,j)=h1((11-ii)*11+(12-j));
>         
>     end
>     
> end

*h1: h的一维向量展开；
*g: 将h旋转180度

> rr = pi*4/5*(pz==0) + pi/5*(pz==1);

*rr : 初始化Mask的振幅，将目标版图作为初始Mask版图；

> m=zeros(N,N);   %Mask pattern
> cun=1000; % 初始误差

下面的while循环是算法主体部分：

> while (sum6>epsilon) & (count<maxloop)
> 
>    count=count+1; 
>    
>    rr=rr-s_phi*dr;   %Update
>    
>    ra=ra-s_theta*da;   %Update
>    
>    m=0.5.*(1+cos(rr)).*exp(i.*ra);   %Calculate continuous mask pattern
>    
>    mr=real(m);   %Real part of continuous mask pattern
>    
>    mi=imag(m);   %Imaginary part of continuous mask pattern
>    
>    mmo=abs(m);   %Amplitude pattern of continuous mask pattern
>    
>    %%%%%%Quantize the complex-valued mask to pole-level mask%%%%%
>    
>    if (phase_n==4)   %Four-phase PSM
>    
>        viccone=mmo>t_m;   %Transparent area on the mask
>        
>        vicctwo=(mr>=0)&(mi>=0);   %Area with phase of pi/4
>        
>        vicctwo=exp(i*pi/4)*vicctwo;
>        
>        viccthree=(mr<=0)&(mi>=0);   %Area with phase of pi*3/4
>        
>        viccthree=exp(i*pi*3/4)*viccthree;
>        
>        viccfour=(mr<=0)&(mi<=0);   %Area with phase of pi*5/4
>        
>        viccfour=exp(i*pi*5/4)*viccfour;
>        
>        viccfive=(mr>=0)&(mi<=0);   %Area with phase of pi*7/4
>        
>        viccfive=exp(i*pi*7/4)*viccfive;
>        
>        viccsix=vicctwo+viccthree+viccfour+viccfive;   %Phase pattern of mask pattern
>        
>        viccin=viccone.*viccsix;   %Pole-level mask pattern
>        
>    elseif (phase_n==2)   %Two-phase PSM
>    
>        viccone=mmo>t_m;   %Transparent area on the mask
>        
>        vicctwo=mr>0;   %Area with phase of 0
>        
>        viccthree=mr<=0;   %Area with phase of pi
>        
>        viccthree=-1*viccthree;
>        
>        viccfour=vicctwo+viccthree;   %Phase pattern of mask pattern
>        
>       viccin=viccone.*viccfour;   %Pole-level mask pattern
>       
>    end
>    
>     viccout=imfilter(viccin,h);
>     
>    viccbin=abs(viccout)>t_r;   %Output pattern of pole-level mask
>
>    sum6=sum(sum(abs(pz-viccbin)));  
>      
>    convergence(count,1)=sum6;
>    
>    if cun>sum6
>       
>         cun=sum6;
>    
>    end
>    
>    disp(cun);
>   
>    mid1=imfilter(m,h);   %Convolution between continuous mask and low-pass filter
>    
>    mid1mo=abs(mid1);   %Convolution between continuous mask amplitude and low-pass filter
>    
>    mid1r=imfilter(mr,h);   %Convolution between real part of continuous mask amplitude and low-pass filter 
>    
>    mid1i=imfilter(mi,h);   %Convolution between imaginary part of continuous mask amplitude and low-pass filter
>    
>    z=1./ (  1+exp(-1*a*(mid1mo)+a*t_r)  ); 
>    
>    mid3=( pz-z ).*z.*(1-z).*mid1r.*(1./mid1mo);   
>    
>    mid5=imfilter(mid3,g);   
>    
>    mid7=( pz-z ).*z.*(1-z).*mid1i.*(1./mid1mo);   
>    
>    mid9=imfilter(mid7,g);
>    
>    %%%%%%Gradient of the discretization penalty corresponding to \phi%%%%%%  
>    
>    dr_D=(-0.5)*sin(rr).*(1+cos(rr));
>    
>    %%%%%%Gradient of the discretization penalty corresponding to \theta%%%%%% 
>    
>    if (phase_n==4)   %Four-phase PSM
>    
>        da_D=8*( sin(4*ra-pi*3/2) + 1 ).*cos(4*ra-pi*3/2);
>    
>    elseif (phase_n==2)   %Two-phase PSM
>    
>        da_D=4.*( sin(2.*ra-pi/2)+1 ).*cos(2.*ra-pi/2);
>    
>    end
>
>    %%%%%%Gradient of wavelet penaly corresponding to \phi%%%%%%
>    
>    for ii=0:N/2-1
>    
>        for jj=0:N/2-1
>        
>            dr_WA(ii*2+1,jj*2+1)= scale(ii*2+1,jj*2+1) * (-1)*sin(rr(ii*2+1,jj*2+1))*real( exp((-i)*ra(ii*2+1,jj*2+1)) *( 3*m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - 
>                                 
>                                   m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) );
>            
>            dr_WA(ii*2+1,jj*2+2)= scale(ii*2+1,jj*2+2) * (-1)*sin(rr(ii*2+1,jj*2+2))*real( exp((-i)*ra(ii*2+1,jj*2+2)) *( 3*m(ii*2+1,jj*2+2) - m(ii*2+1,jj*2+1) -               >                                  
>                                  m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) );
>            
>            dr_WA(ii*2+2,jj*2+1)= scale(ii*2+2,jj*2+1) * (-1)*sin(rr(ii*2+2,jj*2+1))*real( exp((-i)*ra(ii*2+2,jj*2+1)) *( 3*m(ii*2+2,jj*2+1) - m(ii*2+1,jj*2+1) -               >                                  
>                                  m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+2) ) );
>            
>            dr_WA(ii*2+2,jj*2+2)= scale(ii*2+2,jj*2+2) * (-1)*sin(rr(ii*2+2,jj*2+2))*real( exp((-i)*ra(ii*2+2,jj*2+2)) *( 3*m(ii*2+2,jj*2+2) - m(ii*2+1,jj*2+1) -               >                                   
>                                  m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) ) );
>       
>         end
>   
>     end
>    
>    %%%%%%Gradient of wavelet penaly corresponding to \theta%%%%%%
>    
>    for ii=0:N/2-1
>       
>         for jj=0:N/2-1
>            
>            da_WA(ii*2+1,jj*2+1)= scale(ii*2+1,jj*2+1) *  (1+cos(rr(ii*2+1,jj*2+1)))*real( (-i)*exp((-i)*ra(ii*2+1,jj*2+1)) *( 3*m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - >
>            
>            m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) );
>            
>            da_WA(ii*2+1,jj*2+2)= scale(ii*2+1,jj*2+2) *  (1+cos(rr(ii*2+1,jj*2+2)))*real( (-i)*exp((-i)*ra(ii*2+1,jj*2+2)) *( 3*m(ii*2+1,jj*2+2) - m(ii*2+1,jj*2+1) - >
>            
>            m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) );
>            
>            da_WA(ii*2+2,jj*2+1)= scale(ii*2+2,jj*2+1) *  (1+cos(rr(ii*2+2,jj*2+1)))*real( (-i)*exp((-i)*ra(ii*2+2,jj*2+1)) *( 3*m(ii*2+2,jj*2+1) - m(ii*2+1,jj*2+1) - 
>            
>            m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+2) ) );
>            
>            da_WA(ii*2+2,jj*2+2)= scale(ii*2+2,jj*2+2) *  (1+cos(rr(ii*2+2,jj*2+2)))*real( (-i)*exp((-i)*ra(ii*2+2,jj*2+2)) *( 3*m(ii*2+2,jj*2+2) - m(ii*2+1,jj*2+1) -           >                                  
>            m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) ) );
>            
>        end
>        
>    end
>
>    %%%%%%Gradient of overall cost function%%%%%% 
>    
>    dr=a*sin(rr).*cos(ra).*mid5 + a*sin(rr).*sin(ra).*mid9 + gamma_r_D*dr_D + gamma_r_WA*dr_WA;
>    
>    da=2*a*0.5*(1+cos(rr)).*sin(ra).*mid5 - 2*a*0.5*(1+cos(rr)).*cos(ra).*mid9 + gamma_a_D*da_D + gamma_a_WA*da_WA;
>    
> end
