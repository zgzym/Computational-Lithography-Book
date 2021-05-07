# 基于傅里叶级数展开模型的OPC优化算法 （OPC_acaa）

*Fourier series expansion model*

### 本文的matlab代码来源于《Computational Lithography》一书；

OPC_fse算法用于在部分相干成像系统中采用傅里叶级数展开模型对二元掩膜进行梯度优化，

之前所介绍的GPSM_wa、GPSM_tv、OPC_tv和PSM_tv等都是用于相干成像系统。

傅里叶级数展开模型将部分相干光源分解为一系列相干子光源的叠加。

## 本算法产生优化后的二元掩膜并包含了离散罚函数和全局小波罚函数。

## Matlab算法

```python
	function [] = OPC_fse(N, pz, N_filter, pixel, k, NA, lamda, order, sigma_large_inner, sigma_large_outer, step, a, t_r, t_m, gamma_D, gamma_WA, epsilon, maxloop);
```
N: mask的维度

pz: Desired output pattern

N_filter: 幅度冲激响应矩阵维度

pixel: 像素大小 (nm)

k：工艺常数

NA: 数值孔径

lamda: 波长 (nm)

order: Bessel函数的阶数

sigma_large_inner: 内部分相干因子

sigma_large_outer：外部分相干因子

step：步长

a：sigmoid函数的陡度（steepness)

t_r: sigmoid函数的处理阈值

t_m: mask的全局阈值（显影阈值）

gamma_D：离散罚函数的权重

gamma_WA: 全局小波罚函数的权重

epsilon：输出版图的可容忍误差

maxloop：最大迭代次数

```python
	D=pixel*N; % 版图的尺寸（边长）
	
	D_C_1=lamda/2/sigma_large_outer/NA;   %环形光源外环的相干长度
	
	D_C_2=lamda/2/sigma_large_inner/NA;   %环形光源内环的相干长度
	
	omega_0=pi/D; %ω0
	
	N_coherence=floor(2*D/(2*D_C_1)+1)+2;   %环形照明光源的尺寸
	
	midway_coherence=(N_coherence+1)/2;   %环形照明光源的中心点
```

上述代码的计算公式位于文献第26页。

接下来初始化迭代的矩阵的参数：

```python
	d=zeros(N,N);   %Gradient of the cost function
	
	d_D=zeros(N,N);   %Gradient of the discretization penalty
	
	d_WA=zeros(N,N);   %Gradient of the wavelet penalty
	
	convergence=zeros(maxloop,1);   %Output pattern error in each iteration
	
	count=0;   %Index of iteration number
	
	sum6=10000;   %Output pattern error corresponding to the Fourier series expansion model
	
	sum8=10000;   %Output pattern error corresponding to the average coherent approximation model
```

生成幅度冲激响应函数：

```python
	h=zeros(N_filter,N_filter); %幅度冲激响应矩阵 
	
	radius=0;
	
	midway=(N_filter+1)/2; %幅度冲激响应矩阵的中心点坐标
	
	for row=1:N_filter
		
		for column=1:N_filter
			
			radius=pixel*sqrt( (row-midway)^2 + (column-midway)^2 ); % 距中心点的距离
			
			if (radius<=(midway)*pixel)
				
				argument=2*pi*radius*NA/lamda;
				
				if (radius==0)
					
					h(row,column)=h(row-1,column);
				
				else
					
					h(row,column)=besselj(order,argument)/argument; 
				
				end
			
			end
		
		end	
	
	end
	
	h=h/sum(sum(h));   %Normalization
	
	figure 
	
	surf(h);
	
	figure
	
	surf(abs(fftshift(fft2(h)))); % 将h做傅里叶变换后将零频移至矩阵中心

	for ii=1:N_filter
		
		for jj=1:N_filter
			
			h_vector((ii-1)*N_filter+jj)=h(ii,jj); % h的向量化
		
		end
	
	end
	
	for ii=1:N_filter
		
		for jj=1:N_filter
			
			g(ii,jj)=h_vector((N_filter-ii)*N_filter+(N_filter+1-jj)); % g为h的转置
		
		end
	
	end
	
	figure
	
	surf(g);
```

初始化theta:

```python
%%%%%%The initialization of \theta, where r=\theta%%%%%%
r=pi*4/5*(pz==0) + pi/5*(pz==1);
```

生成环形照明光源的像素化矩阵表示：

```python
	radius_1=D/(2*D_C_1);   %环形照明外环的半径（矩阵索引值）
	
	radius_2=D/(2*D_C_2);   %环形照明内环的半径（矩阵索引值）
	
	yita=zeros(N_coherence,N_coherence);   %环形照明光源矩阵
	
	% 判断像素点距中心的距离
	for row=1:N_coherence
	
		for column=1:N_coherence
		
			radius=pixel*sqrt( (row-midway_coherence)^2 + (column-midway_coherence)^2 );
			
			if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
			
				yita(row,column)=1;
			
			end
			
		end
	
	end
	
	yita=yita/sum(sum(yita)); % 归一化
	
	figure
	
	imshow(yita>0);
```

下面为主循环：
```python
	%%%%%%OPC optimization in partially coherent imaging system%%%%%%
	
	m=zeros(N,N);   %Mask pattern
	
	while (sum6>epsilon) & (count<maxloop)   
	   
	   count=count+1; 
	   
	   r=r-step*d;   %Update
	   
	   %%%%%%Calculate pattern error%%%%%%
	   
	   m=(1+cos(r))/2;   %Grey mask
	   
	   m_binary=m>t_m;   %Binary mask
	   
	   aerial=zeros(N,N);   %Aerial image，光刻胶表面光强分布矩阵
	   
	   for p=1:N_coherence
		   
		   for q=1:N_coherence
			   
			   radius=pixel*sqrt( (p-midway_coherence)^2 + (q-midway_coherence)^2 );    
			   
			   if (radius<=radius_1*pixel) & (radius>=radius_2*pixel) % 位于环形照明光源内部时
				   
				   exponential=zeros(N_filter,N_filter); % has the same dimension as the filter h
				   
				   for row=1:N_filter
					   
					   for column=1:N_filter
						   
						   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
						   
						   exponential(row,column)=exp(i*omega_0*argument);
					   
					   end
				   
				   end
				   
				   aerial=aerial+yita(p,q)* abs(  imfilter(double(m_binary),h.*exponential)  ).^2;
			   
			   end
		   
		   end
	   
	   end
	   
	   z_binary=aerial>t_r;   %Binary output pattern，根据全局阈值生成显影后的版图
	   
	   sum6=sum(sum(abs(pz-double(z_binary)))); % 计算输出版图的误差
	   
	   convergence(count,1)=sum6;
```

上述代码的计算公式位于文献第26页。

其中，aerial的表达式为：

![1fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_fse_1.png)

yita(p, q)即为； 

![2fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_fse_2.png)

double(m_binary)为M(r),exponential表达式为：

![3fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_fse_3.png)

接着计算cost function的梯度：
```python
	%%%%%%Gradient of cost function%%%%%%
	   
	   mid1=zeros(N,N);
	   
	   for p=1:N_coherence
		   
		   for q=1:N_coherence
			   
			   radius=pixel*sqrt( (p-midway_coherence)^2 + (q-midway_coherence)^2 );    
			   
			   if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
				   
				   exponential=zeros(N_filter,N_filter);
				   
				   for row=1:N_filter
					   
					   for column=1:N_filter
						   
						   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
						   
						   exponential(row,column)=exp(i*omega_0*argument);
					   
					   end
				   
				   end
				   
				   mid1=mid1+yita(p,q)* abs(imfilter(double(m),h.*exponential)).^2;
			   
			   end
		   
		   end
	   
	   end
	   
	   z=1./ (  1+exp(-a*mid1+a*t_r)  ); 
	   
	   mid3=( pz-z ).*z.*(1-z);   
	   
	   mid5=zeros(N,N);
	   
	   for p=1:N_coherence
		   
		   for q=1:N_coherence
			   
			   radius=pixel*sqrt( (p-midway_coherence)^2 + (q-midway_coherence)^2 );    
			   
			   if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
				   
				   exponential=zeros(N_filter,N_filter);
				   
				   for row=1:N_filter
					   
					   for column=1:N_filter
						   
						   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
						   
						   exponential(row,column)=exp(i*omega_0*argument);
					   
					   end
				   
				   end
				   
				   mid5=mid5+yita(p,q)* real(imfilter( mid3.*imfilter(double(m),h.*exponential) , conj(g.*exponential))+
				        
						imfilter( mid3.*imfilter(double(m),conj(h.*exponential)), g.*exponential));
			   
			   end
		   
		   end
	   
	   end
```

其中，mid1的表达式为：

![4fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_fse_4.png)

z的表达式为：

![5fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_fse_5.png)

mid3的表达式为：

![6fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_fse_6.png)

mid5的表达式为：

![7fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_fse_7.png)

## mid5的代码与文献的计算公式有一点出入，在加法的第2项中，H^m没有复共轭，但是代码里用conj()函数取了复共轭，暂时没有验算是否需要取复共轭。

最后计算离散罚函数与全局小波罚函数的梯度：
```python
	 %%%%%%Gradient of discretization penaly%%%%%%  
	   
	   d_D=( (-8)*m+4 )*(-0.5).*sin(r);
	   
	   %%%%%%Gradient of wavelet penaly%%%%%%
	   
	   for ii=0:(N/2-1)
		   
		   for jj=0:(N/2-1)
			   
			   d_WA(ii*2+1,jj*2+1)= ( 3*m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+1,jj*2+1));
			   
			   d_WA(ii*2+1,jj*2+2)= ( 3*m(ii*2+1,jj*2+2) - m(ii*2+1,jj*2+1) - m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+1,jj*2+2));
			   
			   d_WA(ii*2+2,jj*2+1)= ( 3*m(ii*2+2,jj*2+1) - m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+2,jj*2+1));
			   
			   d_WA(ii*2+2,jj*2+2)= ( 3*m(ii*2+2,jj*2+2) - m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) ) * (-0.5)*sin(r(ii*2+2,jj*2+2));
		   
		   end
	   
	   end
	 
	   %%%%%%Gradient of overall cost function%%%%%%
	   
	   d=a*mid5.*(0.5*sin(r))+gamma_D*d_D +gamma_WA*d_WA;

	   disp(count);
	   
	   disp(sum6);
	
	end
```python
