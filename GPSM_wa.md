# 基于离散罚函数和小波罚函数的相移掩膜优化算法 （GPSM_wa）
**Generalized Phase-Shift Mask Optimazation Algorithm based on discretization penalty and wavelet penalty**

### 本文的matlab代码来源于《Computational Lithography》一书；

GPSM_wa算法在相干照明系统中对N×N的目标版图做广义梯度的相移掩膜优化。

该算法产生优化后的四相或二相相移掩膜，并包含了离散罚函数及局部小波罚函数。

掩膜板上的不同区域可以被赋予不同的小波罚函数权重。

若所有区域的权重均为1，则为全局小波罚函数。

离散罚函数与小波罚函数的作用在于降低mask版图的复杂性，

但一般是以增加曝光版图误差为代价。

本文采用的小波变换为1级哈尔小波变换，

采用更深一层的哈尔小波变换可以进一步降低mask版图的复杂性，但可能会增加误差。

## 1. Matlab版代码及说明
<table><tr><td> function [] = GPSM_wa(N, pz, ra, phase_n, s_phi, s_theta, a, t_r, t_m, gamma_r_D, gamma_a_D, gamma_r_WA, gamma_a_WA, scale, epsilon, maxloop);</td></tr></table>

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
<table><tr><td>
dr=zeros(N,N);   %Gradient of the cost function corresponding to \phi

da=zeros(N,N);   %Gradient of the cost function corresponding to \theta

dr_D=zeros(N,N);   %Gradient of the discretization penalty corresponding to \phi

da_D=zeros(N,N);   %Gradient of the discretization penalty corresponding to \theta

dr_WA=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \phi

da_WA=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \theta

convergence=zeros(maxloop,1);   %Output pattern error in each iteration

count=0;   %Index of iteration number

sum6=100;   %Output pattern error corresponding to the optimized pole-level mask

sum8=100;   %Output pattern error corresponding to the optimized complex-valued mask
</td></tr></table>

以上为基本变量的定义，包括梯度，每次迭代后的误差，迭代次数和输出版图误差

<table><tr><td>
h = fspecial('gaussian',11, 14);
</td></tr></table>

h: 振幅冲激响应函数，fspecial('gaussian', hsize, sigma)返回大小为hsize且标准偏差为sigma（正）的旋转对称高斯低通滤波器。

hsize可以是一个向量，指定h中的行数和列数，也可以是一个标量，在这种情况下，h是一个方阵。

因此h为11×11的标准差为14的高斯型低通滤波器矩阵。

h的物理意义为光学系统对于输入的Mask的响应，通过Mask与h做卷积，即可得到光学系统的输出。

<table><tr><td>
for ii = 1:11
 
     for j = 1:11
     
         h1((ii-1)*11+j) = h(ii, j);
         
     end
 end
 
 for ii = 1:11
 
     for j=1:11
     
         g(ii,j)=h1((11-ii)*11+(12-j));
         
     end
     
 end
</td></tr></table>

*h1: h的一维向量展开；
*g: 将h旋转180度，不是转置！

<table><tr><td>
rr = pi*4/5*(pz==0) + pi/5*(pz==1);
</td></tr></table>

*rr : 初始化Mask的振幅，将目标版图作为初始Mask版图；

<table><tr><td>
m=zeros(N,N);   %Mask pattern
cun=1000; % 初始误差
</td></tr></table>

下面的while循环是算法主体部分：

<table><tr><td>
while (sum6>epsilon) & (count<maxloop)  // 当误差大于所设定的误差epsilon并且循环次数小于最大循环测试maxloop时继续循环
                                                                        \
   count=count+1;  // 循环次数+1
    
   rr=rr-s_phi*dr;   // 更新振幅rr，s_phi为上次循环后计算的振幅的梯度，dr为步长
    
   ra=ra-s_theta*da;   // 更新相位ra，s_theta为上次循环后计算的相位的梯度，da为步长
    
   m=0.5.*(1+cos(rr)).*exp(i.*ra);   // 将mask的振幅的限制条件（0-1）用cos函数替代
</td></tr></table>
 
 振幅表达式如下：
 
 ![2fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/2.png)
 
计算mask pattern的实部，虚部和振幅
 
 <table><tr><td>  
   mr=real(m);   %Real part of continuous mask pattern
   
   mi=imag(m);   %Imaginary part of continuous mask pattern
    
   mmo=abs(m);   %Amplitude pattern of continuous mask pattern
 </td></tr></table>
 
 根据mask pattern的实部，虚部的正负判断相位，根据振幅是否超过阈值t_m判断为透明区域或非透明区域
 
 <table><tr><td>  
   %%%%%%Quantize the complex-valued mask to pole-level mask%%%%%
    
   if (phase_n==4)   %Four-phase PSM
    
        viccone=mmo>t_m;   %Transparent area on the mask，幅度大于t_m为1
        
        vicctwo=(mr>=0)&(mi>=0);   %Area with phase of pi/4
        
        vicctwo=exp(i*pi/4)*vicctwo;
        
        viccthree=(mr<=0)&(mi>=0);   %Area with phase of pi*3/4
        
        viccthree=exp(i*pi*3/4)*viccthree;
        
        viccfour=(mr<=0)&(mi<=0);   %Area with phase of pi*5/4
        
        viccfour=exp(i*pi*5/4)*viccfour;
        
        viccfive=(mr>=0)&(mi<=0);   %Area with phase of pi*7/4
        
        viccfive=exp(i*pi*7/4)*viccfive;
        
        viccsix=vicctwo+viccthree+viccfour+viccfive;   % MASK的相位版图
        
        viccin=viccone.*viccsix;   %Pole-level mask pattern
        
    elseif (phase_n==2)   %Two-phase PSM
    
        viccone=mmo>t_m;   %Transparent area on the mask
        
        vicctwo=mr>0;   %Area with phase of 0
        
        viccthree=mr<=0;   %Area with phase of pi
        
        viccthree=-1*viccthree;
        
        viccfour=vicctwo+viccthree;   %Phase pattern of mask pattern
        
       viccin=viccone.*viccfour;   %Pole-level mask pattern
       
    end
  </td></tr></table>
  
将pole-level mask pattern 与冲激响应函数h做卷积
  
<table><tr><td> 
    viccout=imfilter(viccin,h); 
</td></tr></table>
  
根据输出版图的振幅（绝对值）与sigmoid函数阈值t_r（显影阈值）比较，将其转为输出版图

<table><tr><td> 
    viccbin=abs(viccout)>t_r;   %Output pattern of pole-level mask
</td></tr></table>

计算输出版图viccbin与目标版图pz之间的误差sum6

<table><tr><td>
    sum6=sum(sum(abs(pz-viccbin)));  
      
    convergence(count,1)=sum6;
    
    if cun>sum6
       
         cun=sum6;
    
    end
    
    disp(cun);
</td></tr></table>   

将连续型mask m与冲激响应函数h做卷积：

<table><tr><td>
mid1=imfilter(m,h);   %Convolution between continuous mask and low-pass filter
    
mid1mo=abs(mid1);   %卷积后的绝对值，振幅
    
mid1r=imfilter(mr,h);   %Convolution between real part of continuous mask amplitude and low-pass filter 
     
mid1i=imfilter(mi,h);   %Convolution between imaginary part of continuous mask amplitude and low-pass filter
    
z=1./ (  1+exp(-1*a*(mid1mo)+a*t_r)  ); % sigmoid函数
</td></tr></table>

sigmoid函数表达式如下：
 ![fig3](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/3.png)
 
接下来计算cost function的梯度，首先计算几个中间值：

<table><tr><td>
    mid3=( pz-z ).*z.*(1-z).*mid1r.*(1./mid1mo);   
</td></tr></table>

mid3的表达式为：
 ![fig4](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/4.png)

其中，mid1r为HmR，即m的实部与h的卷积，1./mid1mo即为T(m)，T(m)表达式为：
 ![fig5](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/5.png)

<table><tr><td>
    mid5=imfilter(mid3,g);   
</td></tr></table>

mid5为mid3与g的卷积，g为h的转置，因此mid5的表达式为：
![fig6](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/6.png)

<table><tr><td>
    mid7=( pz-z ).*z.*(1-z).*mid1i.*(1./mid1mo);   
</td></tr></table>

mid7的表达式为：
![fig7](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/7.png)

<table><tr><td>
    mid9=imfilter(mid7,g);
</td></tr></table>

mid9为mid7与g的卷积，其表达式为：
![fig8](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/8.png)

离散罚函数用于降低mask版图的复杂性，将振幅和相位都限制在
预先设定的几个离散值
<table><tr><td>
    %%%%%%Gradient of the discretization penalty corresponding to \phi%%%%%%  
    
    dr_D=(-0.5)*sin(rr).*(1+cos(rr));
    
    %%%%%%Gradient of the discretization penalty corresponding to \theta%%%%%% 
    
    if (phase_n==4)   %Four-phase PSM
    
        da_D=8*( sin(4*ra-pi*3/2) + 1 ).*cos(4*ra-pi*3/2);
    
    elseif (phase_n==2)   %Two-phase PSM
    
        da_D=4.*( sin(2.*ra-pi/2)+1 ).*cos(2.*ra-pi/2);
    
    end
</td></tr></table>

离散罚函数振幅项的梯度为：

![fig9](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/9.png)

对于四相的相移掩膜，其离散罚函数的相位项的梯度为：

![fig10](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/10.png)

对于二相的相移掩膜，其离散罚函数的相位项的梯度为：

![fig11](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/11.png)

计算局域小波罚函数的振幅项的梯度：

<table><tr><td>
%%%%%%Gradient of wavelet penaly corresponding to \phi%%%%%%
    
for ii=0:N/2-1
    
    for jj=0:N/2-1
        
     dr_WA(ii*2+1,jj*2+1)= scale(ii*2+1,jj*2+1) * (-1)*sin(rr(ii*2+1,jj*2+1))*real( exp((-i)*ra(ii*2+1,jj*2+1)) *( 3*m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - 
                                 
                           m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) );
            
     dr_WA(ii*2+1,jj*2+2)= scale(ii*2+1,jj*2+2) * (-1)*sin(rr(ii*2+1,jj*2+2))*real( exp((-i)*ra(ii*2+1,jj*2+2)) *( 3*m(ii*2+1,jj*2+2) - m(ii*2+1,jj*2+1) -                                            
                                  m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) );
            
     dr_WA(ii*2+2,jj*2+1)= scale(ii*2+2,jj*2+1) * (-1)*sin(rr(ii*2+2,jj*2+1))*real( exp((-i)*ra(ii*2+2,jj*2+1)) *( 3*m(ii*2+2,jj*2+1) - m(ii*2+1,jj*2+1) -                                  
                                  m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+2) ) );
            
     dr_WA(ii*2+2,jj*2+2)= scale(ii*2+2,jj*2+2) * (-1)*sin(rr(ii*2+2,jj*2+2))*real( exp((-i)*ra(ii*2+2,jj*2+2)) *( 3*m(ii*2+2,jj*2+2) - m(ii*2+1,jj*2+1) -                                      
                                  m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) ) );
       
     end
   
end
</td></tr></table>

dr_WA的表达式为：
 
![fig12](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/12.png) 

计算局域小波罚函数的相位项的梯度：

<table><tr><td>
%%%%%%Gradient of wavelet penaly corresponding to \theta%%%%%%

for ii=0:N/2-1
   
	 for jj=0:N/2-1
		
		da_WA(ii*2+1,jj*2+1)= scale(ii*2+1,jj*2+1) *  (1+cos(rr(ii*2+1,jj*2+1)))*real( (-i)*exp((-i)*ra(ii*2+1,jj*2+1)) *( 3*m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - >
		
		m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) );
		
		da_WA(ii*2+1,jj*2+2)= scale(ii*2+1,jj*2+2) *  (1+cos(rr(ii*2+1,jj*2+2)))*real( (-i)*exp((-i)*ra(ii*2+1,jj*2+2)) *( 3*m(ii*2+1,jj*2+2) - m(ii*2+1,jj*2+1) - >
		
		m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) );
		
		da_WA(ii*2+2,jj*2+1)= scale(ii*2+2,jj*2+1) *  (1+cos(rr(ii*2+2,jj*2+1)))*real( (-i)*exp((-i)*ra(ii*2+2,jj*2+1)) *( 3*m(ii*2+2,jj*2+1) - m(ii*2+1,jj*2+1) - 
		
		m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+2) ) );
		
		da_WA(ii*2+2,jj*2+2)= scale(ii*2+2,jj*2+2) *  (1+cos(rr(ii*2+2,jj*2+2)))*real( (-i)*exp((-i)*ra(ii*2+2,jj*2+2)) *( 3*m(ii*2+2,jj*2+2) - m(ii*2+1,jj*2+1) -                                             
		m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) ) );
		
	end
	
end
</td></tr></table>

da_WA的表达式为：

![fig13](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/13.png) 

其中，scale(i,j)为局域罚函数的权重矩阵。

最后，计算整个cost function的梯度：

<table><tr><td>
%%%%%%Gradient of overall cost function%%%%%% 

dr=a*sin(rr).*cos(ra).*mid5 + a*sin(rr).*sin(ra).*mid9 + gamma_r_D*dr_D + gamma_r_WA*dr_WA;

da=2*a*0.5*(1+cos(rr)).*sin(ra).*mid5 - 2*a*0.5*(1+cos(rr)).*cos(ra).*mid9 + gamma_a_D*da_D + gamma_a_WA*da_WA;

end
</td></tr></table>

dr的前两项表达式为：

![fig14](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/14.png)

da的前两项表达式为：

![fig15](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/15.png)



