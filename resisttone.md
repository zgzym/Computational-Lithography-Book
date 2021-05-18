# 基于光刻胶极性反转的两相位相移掩膜优化算法 （resisttone）

### 本文的matlab代码来源于《Computational Lithography》一书；

proc_dct算法用于部分相干成像系统中的相移掩膜优化，采用SVD模型和光刻胶极性反转法结合梯度下降法。

本算法包含了离散罚函数和小波罚函数，同时可以选择是否进行二维离散余弦变换处理来降低mask的复杂度。

## Matlab算法

```python
	function [] = resisttone(N_mask, desire_pattern, distribution, pixel, k, NA, lamda, sigma, order, 
							 threshold, flag, step, a, t_r, t_m, gamma_D, gamma_WA, epsilon, maxloop);
```
N_mask: mask的维度

desire_pattern: Desired output pattern

distribution: 光刻胶的分布

pixel: 像素尺寸(nm)

k: 工艺常数

NA: 数值孔径

lamda: 波长(nm)

sigma: 部分相干因子

order：Bessel函数的阶数

threshold: 离散余弦变换谱的高频分量的截止阈值，保留的低频分量个数为：(threshold - 1) * (threshold - 2) / 2

flag：flag=1时进行2D DCT，flag=0不进行2D DCT

step：步长

a：sigmoid函数的陡度

t_r: 一阶相干近似的处理阈值

t_m: mask的全局阈值（用于离散化mask)

gamma_D：离散罚函数的权重

gamma_WA: 小波罚函数的权重

epsilon：输出版图的最小误差

maxloop：最大循环次数


调用SOC()函数计算TCC：
```python	                
	midway=(N_mask+1)/2;   %Middle point of mask
	
	TCC=zeros(N_mask^2,N_mask^2);   %Transmission cross coefficients
	
	m_trinary_new=zeros(N_mask,N_mask);   %Trinary mask after the optimal post-processing

	%%%%%%Calculate the transmission cross coefficients%%%%%%
	[TCC] = SOCS(N_mask, pixel, k, NA, lamda, midway, sigma, order);
```

对TCC做SVD分解：
```python
	%%%%%%Singular value decomposition of the partially coherent imaging system%%%%%%
	[U,S,V]=svd(TCC);   %Singular value decomposition
	
	h_1_fre=reshape(U(1:N_mask^2,1:1),N_mask,N_mask);
	
	h_1=(fftshift(ifft2(ifftshift((h_1_fre)))));   %The impulse response of the first order approximation 
	
	sum_eigenvalue=(sum(sum(S)));   %The summation of the eigenvalues 
```

初始化参数：
```python
	t_r_real=0;   %Global threshold of photoresist effect for all of the eigen value
	
	d=zeros(N_mask,N_mask);   %Gradient of the cost function
	
	d_D=zeros(N_mask,N_mask);   %Gradient of the discretization penalty
	
	d_WA=zeros(N_mask,N_mask);   %Gradient of the wavelet penalty

	convergence=zeros(maxloop,1);   %Output pattern error in each iteration
	
	count=0;   %Index of iteration number
	
	sum6=100;   %Output pattern error corresponding to the optimized trinary mask
	
	sum8=100;   %Output pattern error corresponding to the optimized real-valued mask
```

计算冲激响应函数h_1的转置：
```python
	%%%%%%the amplitude impulse response of the partially coherent imaging system%%%%%%
	h_1_1=zeros(N_mask,N_mask);
	
	h_1_1=h_1;

	for ii=1:N_mask
		
		for jj=1:N_mask
			
			h_1_vector((ii-1)*N_mask+jj)=h_1_1(ii,jj);
		
		end
	
	end
	
	for ii=1:N_mask
		
		for jj=1:N_mask
			
			g_1(ii,jj)=h_1_vector((N_mask-ii)*N_mask+(N_mask+1-jj)); %inverse vector
		
		end
	
	end
```

下面为主循环：
```python
	%%%%%%The desired pattern is two-bars with photoresist reversion%%%%%%
	pz=desire_pattern-distribution+1; % 本行代码引入光刻胶极性反转的影响

	%%%%%%The initialization of \theta, where r=\theta%%%%%%
	r=pi/2*(pz==0) + pi/5*(pz==1);

	%%%%%%PSM optimization in partially coherent imaging system%%%%%%
	m=zeros(N_mask,N_mask);   %Mask pattern
	
	while (sum6>epsilon) & (count<maxloop)
		
		count=count+1;
		
		%%%%%%Calculate pattern error%%%%%%
		m=cos(r);   %Gray mask
		
		m_trinary_p=m>t_m;
		
		m_trinary_n=-1*(m<(-1*t_m));
		
		m_trinary=m_trinary_p+m_trinary_n;   %Trinary mask
		
		aerial=zeros(N_mask,N_mask);   %Aerial image 
		
		aerial=(  abs(imfilter(double(m_trinary),h_1_1)).^2   );
		
		z_trinary=aerial>t_r;   %Binary output pattern
		
		sum6=sum(sum(abs(abs(pz)-z_trinary)));   %Output pattern error of trinary mask 
		
		convergence(count,1)=sum6; 
	  
		%%%%%%Gradient of cost function%%%%%%
		mid1=abs(imfilter(double(m),h_1_1)).^2;
		
		z=1./(  1+exp(-a*mid1+a*t_r)  ); 
		
		mid3=(pz-z).*z.*(1-z);   
		
		mid4=mid3.*imfilter(double(m),h_1_1);
		
		mid4_4=mid3.*imfilter(double(m),conj(h_1_1));
		
		mid5=real(imfilter(double(mid4),conj(g_1))+imfilter(double(mid4_4),g_1));
		
		%%%%%%Gradient of discretization penaly%%%%%%  
		d_D=( (-18)*m.^3+2*m ).*((-1)*sin(r));
		
		%%%%%%Gradient of wavelet penalty%%%%%%
		for ii=0:(N_mask/2-1)
			
			for jj=0:(N_mask/2-1)
				
				d_WA(ii*2+1,jj*2+1)= ( 3*m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+1,jj*2+1));
				
				d_WA(ii*2+1,jj*2+2)= ( 3*m(ii*2+1,jj*2+2) - m(ii*2+1,jj*2+1) - m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+1,jj*2+2));
				
				d_WA(ii*2+2,jj*2+1)= ( 3*m(ii*2+2,jj*2+1) - m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+2,jj*2+1));
				
				d_WA(ii*2+2,jj*2+2)= ( 3*m(ii*2+2,jj*2+2) - m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) ) * (-0.5)*sin(r(ii*2+2,jj*2+2));
			
			end
		
		end
	   
		%%%%%%%Calculate whole revision vector%%%%%%%%%%
		d=2*a*mid5.*sin(r) + gamma_D*d_D + gamma_WA*d_WA;
		
		r=r-step*d;   %Update

	end
```
上述代码的相关公式可参考PSM_svd.md

设置光刻胶显影阈值：
```python
	t_r_real=t_r*sum_eigenvalue;   %Global threshold of photoresist effect for all of the eigen value
```

当flag=1时，做二维离散余弦变换：
```python
	[m_trinary_new] = proc_dct(N_mask, pz, m, t_r, t_r_real, t_m, TCC, threshold); 
```
