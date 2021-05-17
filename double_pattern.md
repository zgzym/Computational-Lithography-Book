# 基于两次曝光的PSM优化算法 （double_pattern）

*Double patterning optimization using two generalized PSMs in coherent imaging system.*

### 本文的matlab代码来源于《Computational Lithography》一书；

PSM_svd算法用于相干成像系统中采用2个广义PSM进行2次曝光来进行梯度优化，适用于相干成像系统。

### 本算法产生优化后的二相相移掩膜并包含了离散罚函数和小波罚函数。

## Matlab算法

```python
	function [] = double_pattern(s_phi_one, s_theta_one, s_phi_two, s_theta_two, a, t_r, t_m, gamma_r_D_one, gamma_a_D_one, 
	
				  gamma_r_WA_one, gamma_a_WA_one, gamma_r_D_two, gamma_a_D_two, gamma_r_WA_two, gamma_a_WA_two, epsilon, maxloop);
```
s_phi_one: 第1个mask的强度（φ）的步长

s_theta_one: 第1个mask的相位（θ）的步长

s_phi_two: 第2个mask的强度（φ）的步长

s_theta_two: 第2个mask的相位（θ）的步长

a: sigmoid函数的陡度

t_r: 处理阈值

t_m: mask的全局阈值，用于将连续型mask(cos(m))离散化为-1,0,1

gamma_r_D_one: 第1个mask强度矩阵的离散罚函数权重

gamma_a_D_one: 第1个mask相位矩阵的离散罚函数权重

gamma_r_WA_one: 第1个mask强度矩阵的小波罚函数权重

gamma_a_WA_one: 第1个mask相位矩阵的小波罚函数权重

gamma_r_D_two: 第2个mask强度矩阵的离散罚函数权重

gamma_a_D_two: 第2个mask相位矩阵的离散罚函数权重

gamma_r_WA_two: 第2个mask强度矩阵的小波罚函数权重

gamma_a_WA_two: 第2个mask相位矩阵的小波罚函数权重

epsilon：输出版图的可容忍误差

maxloop：最大迭代次数

初始化参数：
```python
	N=80;   %Mask dimension
	
	dr_one=zeros(N,N);   %Gradient of the cost function corresponding to \phi for the first mask
	
	da_one=zeros(N,N);   %Gradient of the cost function corresponding to \theta for the first mask
	
	dr_D_one=zeros(N,N);   %Gradient of the discretization penalty corresponding to \phi for the first mask
	
	da_D_one=zeros(N,N);   %Gradient of the discretization penalty corresponding to \theta for the first mask
	
	dr_WA_one=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \phi for the first mask
	
	da_WA_one=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \theta for the first mask
	
	dr_two=zeros(N,N);   %Gradient of the cost function corresponding to \phi for the second mask
	
	da_two=zeros(N,N);   %Gradient of the cost function corresponding to \theta for the second mask
	
	dr_D_two=zeros(N,N);   %Gradient of the discretization penalty corresponding to \phi for the second mask
	
	da_D_two=zeros(N,N);   %Gradient of the discretization penalty corresponding to \theta for the second mask
	
	dr_WA_two=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \phi for the second mask
	
	da_WA_two=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \theta for the second mask
	
	% epsilon=10;   %Tolerable output pattern error
	
	% maxloop=100;   %Maximum iteration number
	
	convergence=zeros(maxloop,1);   %Output pattern error in each iteration
	
	count=0;   %Index of iteration number
	
	sum6=100;   %Output pattern error corresponding to the optimized pole-level mask
	
	cun=100;   %The lowest output pattern error that has been obtained

```
相干成像系统的幅度冲激响应函数：
```python
	h=fspecial('gaussian',11,14);
```

h的转置(g):
```python
	for ii=1:11
		
		for j=1:11
			
			h1((ii-1)*11+j)=h(ii,j);
		
		end
	
	end
	
	for ii=1:11
		
		for j=1:11
			
			g(ii,j)=h1((11-ii)*11+(12-j));
		
		end
	
	end
```

定义desired output pattern：

```python
	pz=zeros(N,N);
	
	for ii=21:60
		
		for j=21:60
			
			pz(ii,j)=1;
		
		end
	
	end
	
	for ii=36:60
		
		for j=36:45
			
			pz(ii,j)=0;
		
		end
	
	end
```

初始化mask1的强度版图和相位版图：
```python
	%%%%%%The initialization of \phi, where r=\phi for the first mask%%%%%%
	
	rr_one=ones(N,N)*pi*4/5;
	
	for ii=21:35
		
		for j=21:60
			
			rr_one(ii,j)=pi/5;
		
		end
	
	end

	%%%%%%The initialization of \theta, where r=\theta for the first mask%%%%%%
	
	ra_one=ones(N,N)*pi/5;
```

初始化mask2的强度版图和相位版图：
```python
	%%%%%%The initialization of \phi, where r=\phi for the second mask%%%%%%
	
	rr_two=ones(N,N)*pi*4/5;
	
	for ii=36:60
		
		for j=21:35
			
			rr_two(ii,j)=pi/5;
		
		end
	
	end
	
	for ii=36:60
		
		for j=46:60
			
			rr_two(ii,j)=pi/5;
		
		end
	
	end
	
	%%%%%%The initialization of \theta, where ra=\theta for the second mask%%%%%%
	
	for ii=1:80
		
		for j=1:40
			
			ra_two(ii,j)=6*pi/5;
		
		end
	
	end
	
	for ii=1:80
		
		for j=41:80
			
			ra_two(ii,j)=pi/5;
		
		end
	
	end
```

下面是主循环：
```python
	%%%%%%Double-patterning for generalized PSM optimization in coherent imaging system%%%%%
	
	m_one=zeros(N,N);   %The first mask pattern
	
	m_two=zeros(N,N);   %The second mask pattern
	
	while (sum6>epsilon) & (count<maxloop) % 误差大于epsilon或迭代次数小于maxloop时继续迭代
		
		count=count+1; 
		
		rr_one=rr_one-s_phi_one*dr_one;   %Update
		
		ra_one=ra_one-s_theta_one*da_one;   %Update
		
		rr_two=rr_two-s_phi_two*dr_two;   %Update
		
		ra_two=ra_two-s_theta_two*da_two;   %Update
```

计算复数域的mask pattern：
```python
	m_one=0.5.*(1+cos(rr_one)).*exp(i.*ra_one);   %Calculate the first complex-valued mask pattern
    
	mr_one=real(m_one);   %Real part of the first complex-valued mask pattern
    
	mi_one=imag(m_one);   %Imaginary part of the first complex-valued mask pattern
    
	mmo_one=abs(m_one);   %Amplitude pattern of the first complex-valued mask pattern
    
	m_two=0.5.*(1+cos(rr_two)).*exp(i.*ra_two);   %Calculate the second complex-valued mask pattern
    
	mr_two=real(m_two);   %Real part of the second complex-valued mask pattern
    
	mi_two=imag(m_two);   %Imaginary part of the second complex-valued mask pattern
    
	mmo_two=abs(m_two);   %Amplitude pattern of the second complex-valued mask pattern
```

根据t_m（全局阈值）将连续型mask（mmo_one和mmo_two）量化为Pole-level mask：
```python

    %%%%%%Quantize the first complex-valued mask to pole-level mask%%%%%%  
    
	viccone_one=mmo_one>t_m;
    
	vicctwo_one=mr_one>0;
    
	viccthree_one=mr_one<=0;
    
	viccthree_one=-1*viccthree_one;
    
	viccfour_one=vicctwo_one+viccthree_one;
    
	viccin_one=viccone_one.*viccfour_one;   %The first pole-level mask
    
    viccout_one=imfilter(viccin_one,h);
    
	viccbin_one=abs(viccout_one)>t_r;   %Output pattern of the first pole-level mask
    
    %%%%%%Quantize the second complex-valued mask to pole-level mask%%%%%%
    
	viccone_two=mmo_two>t_m;
    
	vicctwo_two=mr_two>0;
    
	viccthree_two=mr_two<=0;
    
	viccthree_two=-1*viccthree_two;
    
	viccfour_two=vicctwo_two+viccthree_two;
    
	viccin_two=viccone_two.*viccfour_two;   %The second pole-level mask

    viccout_two=imfilter(viccin_two,h);
    
	viccbin_two=abs(viccout_two)>t_r;   %Output pattern of the second pole-level mask
```

计算误差：
```python
	sum6=sum(sum(abs(pz-((viccbin_one+viccbin_two)>=1))));   %Output pattenr error
    
	convergence(count,1)=sum6;
    
	if cun>sum6   %The lowest output pattern error that has been obtained
       
	   cun=sum6;
    
	end
```

计算mask与冲激响应函数h的矩阵，并计算相应的输出版图
```python
    mid1_one=imfilter(m_one,h);   %Convolution between the first complex-valued mask and low-pass filter
    
	mid1mo_one=abs(mid1_one);   %Convolution between the first complex-valued mask amplitude and low-pass filter
    
	mid1r_one=imfilter(mr_one,h);   %Convolution between real part of the first complex-valued mask amplitude and low-pass filter
    
	mid1i_one=imfilter(mi_one,h);   %Convolution between imaginary part of the first complex-valued mask amplitude and low-pass filter
    
	mid1_two=imfilter(m_two,h);   %Convolution between the second complex-valued mask and low-pass filter
    
	mid1mo_two=abs(mid1_two);   %Convolution between the second complex-valued mask amplitude and low-pass filter
    
	mid1r_two=imfilter(mr_two,h);   %Convolution between real part of the second complex-valued mask amplitude and low-pass filter
    
	mid1i_two=imfilter(mi_two,h);   %Convolution between imaginary part of the second complex-valued mask amplitude and low-pass filter
   
    z_one=1./ (  1+exp(-1*a*(mid1mo_one)+a*t_r)  ); % mask1的输出版图
    
	z_two=1./ (  1+exp(-1*a*(mid1mo_two)+a*t_r)  ); % mask2的输出版图
```
其中， z_one和z_two的表达式分别：

![1fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/double_pattern_1.png)

两个版图叠加后的版图：
```python
	F=0.5*(tanh(z_one+z_two-1)+1);   %Overall photoresist effect of the two exposures
```

其中，tanh是用于近似替代阶跃函数U(·)，因为阶跃函数求导后会引入Dirac冲激项，不方便进行进一步的分析，两个函数如下图所示：

![2fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/double_pattern_2.png)

下面计算mask1和mask2的cost function的梯度：
```python
	%%%%%%Calculations for the gradient of the first mask%%%%%%
    
	mid3_one=(pz-F).*(sech(z_one+z_two-1)).^2.*z_one.*(1-z_one).*mid1r_one.*(1./mid1mo_one);   
    
	mid5_one=imfilter(mid3_one,g);  
    
	mid7_one=(pz-F).*(sech(z_one+z_two-1)).^2.*z_one.*(1-z_one).*mid1i_one.*(1./mid1mo_one);   
    
	mid9_one=imfilter(mid7_one,g);
    
	%%%%%%Calculations for the gradient of the second mask%%%%%%
    
	mid3_two=(pz-F).*(sech(z_one+z_two-1)).^2.*z_two.*(1-z_two).*mid1r_two.*(1./mid1mo_two);   
    
	mid5_two=imfilter(mid3_two,g);   
    
	mid7_two=(pz-F).*(sech(z_one+z_two-1)).^2.*z_two.*(1-z_two).*mid1i_two.*(1./mid1mo_two);   
    
	mid9_two=imfilter(mid7_two,g);
```

计算离散罚函数及小波罚函数的梯度：
```python
    %%%%%%Gradient of discretization penaly of the first mask%%%%%%
    
	dr_D_one=(-0.5)*sin(rr_one).*(1+cos(rr_one));
    
	da_D_one=8*( sin(4*ra_one-pi*3/2)+1 ).*cos(4*ra_one-pi*3/2);
 
    %%%%%%Gradient of discretization penaly of the second mask%%%%%%
    
	dr_D_two=(-0.5)*sin(rr_two).*(1+cos(rr_two));
    
	da_D_two=8*( sin(4*ra_two-pi*3/2)+1 ).*cos(4*ra_two-pi*3/2);
    
    %%%%%%Gradient of wavelet penalty of the first mask%%%%%%
    
	for ii=0:N/2-1
        
		for jj=0:N/2-1
            
			dr_WA_one(ii*2+1,jj*2+1)=  -1*sin(rr_one(ii*2+1,jj*2+1))*real( exp((-i)*ra_one(ii*2+1,jj*2+1)) *( 3*m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+1) - m_one(ii*2+2,jj*2+2) ) );
            
			dr_WA_one(ii*2+1,jj*2+2)=  -1*sin(rr_one(ii*2+1,jj*2+2))*real( exp((-i)*ra_one(ii*2+1,jj*2+2)) *( 3*m_one(ii*2+1,jj*2+2) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+2,jj*2+1) - m_one(ii*2+2,jj*2+2) ) );
            
			dr_WA_one(ii*2+2,jj*2+1)=  -1*sin(rr_one(ii*2+2,jj*2+1))*real( exp((-i)*ra_one(ii*2+2,jj*2+1)) *( 3*m_one(ii*2+2,jj*2+1) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+2) ) );
            
			dr_WA_one(ii*2+2,jj*2+2)=  -1*sin(rr_one(ii*2+2,jj*2+2))*real( exp((-i)*ra_one(ii*2+2,jj*2+2)) *( 3*m_one(ii*2+2,jj*2+2) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+1) ) );
        
		end
    
	end
    
	for ii=0:N/2-1
        
		for jj=0:N/2-1
            
			da_WA_one(ii*2+1,jj*2+1)=  (1+cos(rr_one(ii*2+1,jj*2+1)))*real( (-i)*exp((-i)*ra_one(ii*2+1,jj*2+1)) *( 3*m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+1) - m_one(ii*2+2,jj*2+2) ) );
            
			da_WA_one(ii*2+1,jj*2+2)=  (1+cos(rr_one(ii*2+1,jj*2+2)))*real( (-i)*exp((-i)*ra_one(ii*2+1,jj*2+2)) *( 3*m_one(ii*2+1,jj*2+2) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+2,jj*2+1) - m_one(ii*2+2,jj*2+2) ) );
            
			da_WA_one(ii*2+2,jj*2+1)=  (1+cos(rr_one(ii*2+2,jj*2+1)))*real( (-i)*exp((-i)*ra_one(ii*2+2,jj*2+1)) *( 3*m_one(ii*2+2,jj*2+1) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+2) ) );
            
			da_WA_one(ii*2+2,jj*2+2)=  (1+cos(rr_one(ii*2+2,jj*2+2)))*real( (-i)*exp((-i)*ra_one(ii*2+2,jj*2+2)) *( 3*m_one(ii*2+2,jj*2+2) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+1) ) );
        
		end
    
	end
    
    %%%%%%Gradient of wavelet penalty of the second mask%%%%%%
    
	for ii=0:N/2-1
        
		for jj=0:N/2-1
            
			dr_WA_two(ii*2+1,jj*2+1)=  -1*sin(rr_two(ii*2+1,jj*2+1))*real( exp((-i)*ra_two(ii*2+1,jj*2+1)) *( 3*m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+1) - m_two(ii*2+2,jj*2+2) ) );
            
			dr_WA_two(ii*2+1,jj*2+2)=  -1*sin(rr_two(ii*2+1,jj*2+2))*real( exp((-i)*ra_two(ii*2+1,jj*2+2)) *( 3*m_two(ii*2+1,jj*2+2) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+2,jj*2+1) - m_two(ii*2+2,jj*2+2) ) );
            
			dr_WA_two(ii*2+2,jj*2+1)=  -1*sin(rr_two(ii*2+2,jj*2+1))*real( exp((-i)*ra_two(ii*2+2,jj*2+1)) *( 3*m_two(ii*2+2,jj*2+1) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+2) ) );
            
			dr_WA_two(ii*2+2,jj*2+2)=  -1*sin(rr_two(ii*2+2,jj*2+2))*real( exp((-i)*ra_two(ii*2+2,jj*2+2)) *( 3*m_two(ii*2+2,jj*2+2) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+1) ) );
        
		end
    
	end
    
	for ii=0:N/2-1
        
		for jj=0:N/2-1
            
			da_WA_two(ii*2+1,jj*2+1)=  (1+cos(rr_two(ii*2+1,jj*2+1)))*real( (-i)*exp((-i)*ra_two(ii*2+1,jj*2+1)) *( 3*m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+1) - m_two(ii*2+2,jj*2+2) ) );
            
			da_WA_two(ii*2+1,jj*2+2)=  (1+cos(rr_two(ii*2+1,jj*2+2)))*real( (-i)*exp((-i)*ra_two(ii*2+1,jj*2+2)) *( 3*m_two(ii*2+1,jj*2+2) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+2,jj*2+1) - m_two(ii*2+2,jj*2+2) ) );
            
			da_WA_two(ii*2+2,jj*2+1)=  (1+cos(rr_two(ii*2+2,jj*2+1)))*real( (-i)*exp((-i)*ra_two(ii*2+2,jj*2+1)) *( 3*m_two(ii*2+2,jj*2+1) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+2) ) );
            
			da_WA_two(ii*2+2,jj*2+2)=  (1+cos(rr_two(ii*2+2,jj*2+2)))*real( (-i)*exp((-i)*ra_two(ii*2+2,jj*2+2)) *( 3*m_two(ii*2+2,jj*2+2) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+1) ) );
        
		end
    
	end
```

最后计算总的cost function的梯度：
```python
	%%%%%%Gradient of the overall cost function for the first mask%%%%%%
    
	dr_one=0.5*( a*sin(rr_one).*cos(ra_one).*mid5_one + a*sin(rr_one).*sin(ra_one).*mid9_one ) + gamma_r_D_one*dr_D_one + gamma_r_WA_one*dr_WA_one;
    
	da_one=a*0.5*(1+cos(rr_one)).*sin(ra_one).*mid5_one - a*0.5*(1+cos(rr_one)).*cos(ra_one).*mid9_one + gamma_a_D_one*da_D_one + gamma_a_WA_one*da_WA_one;      

    %%%%%%Gradient of the overall cost function for the second mask%%%%%%
    
	dr_two=0.5*( a*sin(rr_two).*cos(ra_two).*mid5_two + a*sin(rr_two).*sin(ra_two).*mid9_two ) + gamma_r_D_two*dr_D_two + gamma_r_WA_two*dr_WA_two;
    
	da_two=a*0.5*(1+cos(rr_two)).*sin(ra_two).*mid5_two - a*0.5*(1+cos(rr_two)).*cos(ra_two).*mid9_two + gamma_a_D_two*da_D_two + gamma_a_WA_two*da_WA_two; 

end
```
其中，mask1与mask2的强度矩阵的梯度为：

![3fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/double_pattern_3.png)

其中，mask1与mask2的相位矩阵的梯度为：

![4fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/double_pattern_4.png)

![5fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/double_pattern_5.png)
