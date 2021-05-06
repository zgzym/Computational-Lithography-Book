# 基于平均相干近似模型的OPC优化算法 （OPC_acaa）

*average coherence approximation algorithm*

### 本文的matlab代码来源于《Computational Lithography》一书；

OPC_acaa算法用于在部分相干成像系统中采用平均相干近似模型对二元掩膜进行梯度优化，

之前所介绍的GPSM_wa、GPSM_tv、OPC_tv和PSM_tv等都是用于相干成像系统。

平均相干近似模型将部分相干光源分解为一个相干光源和一个完全不相干光源的叠加。

## 本算法产生优化后的二元掩膜并包含了离散罚函数和全局小波罚函数。

## Matlab算法

```python
	function [] = OPC_acaa(N, pz, N_filter, pixel, k, NA, lamda, order, sigma_large_inner, sigma_large_outer, step, a, t_r, tr_approx, t_m, gamma_D, gamma_WA, epsilon, maxloop);
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

tr_approx：平均相干近似模型的处理阈值

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

接着冲激响应函数的幅度矩阵：

```python
	a_half_d=2*D_C_1;
	
	d=lamda/(2*NA);
	
	h_C=zeros(N_filter,N_filter);
	
	h_I=zeros(N_filter,N_filter);
	
	h_total=zeros(N_filter,N_filter);
	
	radius=0;
	
	midway=(N_filter+1)/2;%低通滤波器h的中心点索引值(midway, midway)
```

根据坐标判断距中心点的距离(radius)，其中besselj(order, argument)返回自变量为argument的order阶一类bessel函数值
```python
	for row=1:N_filter
		
		for column=1:N_filter
			
			radius=pixel*sqrt( (row-midway)^2 + (column-midway)^2 );
			
			if (radius<=(midway)*pixel)
				
				argument=2*pi*radius*NA/lamda;
				
				if (radius==0)
					
					h_total(row,column)=h_total(row-1,column);
				
				else
					
					h_total(row,column)=besselj(order,argument)/argument;
				
				end
			
			end
		
		end
	
	end
```

h_total的计算公式为：

![1fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_1.png)

```python
	h_total=h_total/sum(sum(h_total)); %归一化
	
	integral=sum(sum(h_total.^2))*pixel^2; % 计算Norm^2，在SVD中，该值表示kernel对于Model的重要程度

	vic=zeros(N+1,N+1);
	
	% vic的中心区域(N_filter)为h_total，其余区域为0
	vic((N+2)/2-midway+1:(N+2)/2+midway-1,(N+2)/2-midway+1:(N+2)/2+midway-1)=h_total;
	
	%fft2()为2D傅里叶变换，fftshift()将零频分量移至矩阵中心位置
	f_coe_fft=fftshift(fft2(vic.^2));
```

下图分别为fftshift()之前和之后的频率分布图：

![2fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_2.png)

![3fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_3.png)

将位于环形照明光源之外的点的fft系数设为0：
```python
	for row=1:N
	
		for column=1:N
		
			radius=pixel*sqrt( (row-(N+2)/2)^2 + (column-(N+2)/2)^2 );
			
			if (radius>radius_1*pixel) || (radius<radius_2*pixel)
			
				f_coe_fft(row,column)=0;
				
			end
			
		end
	
	end	
```
```python
	f_coe=real(ifft2(ifftshift(f_coe_fft)));
	
	f_coe=a_half_d^2*f_coe/(integral);
```

f_coe的表达式为：

![4fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_4.png)

其中，integral即为：

![5fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_5.png)

a_half_d^2为：

![6fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_6.png)

下面计算相干光源部分的占比：
```python
	% 
	f_coe=f_coe((N+2)/2-midway+1:(N+2)/2+midway-1,(N+2)/2-midway+1:(N+2)/2+midway-1);
```
```python
	h_C=(sqrt(f_coe).*h_total);   %Equivalent coherent component
	
	h_I=(sqrt(1-f_coe)).*h_total;   %Equivalent incoherent component
```
其中，h_C为等价幅度冲激响应函数的相干分量，h_I为非相干分量，其表达式为：

![7fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_7.png)

将幅度冲激响应矩阵向量化：
```python
	for ii=1:N_filter
		
		for jj=1:N_filter
			
			h_C_vector((ii-1)*N_filter+jj)=h_C(ii,jj);
			
			h_I_vector((ii-1)*N_filter+jj)=h_I(ii,jj);
		
		end
	
	end
	
	for ii=1:N_filter   %h的转置，g_I为h_I的转置的平方
		
		for jj=1:N_filter
			
			g_C(ii,jj)=h_C_vector((N_filter-ii)*N_filter+(N_filter+1-jj)); %inverse vector
			
			g_I(ii,jj)=abs(h_I_vector((N_filter-ii)*N_filter+(N_filter+1-jj)))^2; %inverse vector  have squre!!!!!!!!!
		
		end
	
	end
```

初始化相位和mask：

```python
	%%%%%%The initialization of \theta, where r=\theta%%%%%%
	r=pi*4/5*(pz==0) + pi/5*(pz==1);
	 
	%%%%%%OPC optimization in partially coherent imaging system%%%%%%
	m=zeros(N,N);   %Mask pattern
```

下面是主循环：
```python
	while (sum8>epsilon) & (count<maxloop)
	
	   count=count+1; 
	   
	   r=r-step*d;   %Update
	   
	   %%%%%%Calculate pattern error%%%%%%
	   
	   m=(1+cos(r))/2;   %Grey mask
	   
	   m_binary=m>t_m;   %Binary mask，用于计算离散化的输出版图
	   
	   aerial=zeros(N,N);   %Aerial image 
	   
	   aerial=(abs(imfilter(double(m_binary),h_C)).^2 + imfilter(double(m_binary).^2,abs(h_I).^2));
```
aerial的表达式为：

![8fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_8.png)

```python	   
	   z_binary=aerial>tr_approx;   %Binary output pattern
	   
	   sum8=sum(sum(abs(abs(pz)-z_binary)));
	   
	   convergence(count,1)=sum8;

	   %%%%%%Gradient of cost function%%%%%%
	   
	   % double()用于转换数据精度 
	   
	   mid1=(abs(imfilter(double(m),h_C)).^2 + imfilter(double(m).^2,abs(h_I).^2));
	   
	   z=1./ (1+exp(-a*mid1+a*tr_approx));  
```
z的表达式为：

![9fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_9.png)

```python	   
	   mid3=(pz-z).*z.*(1-z);   
	   
	   mid4=mid3.*imfilter(double(m),h_C);
```

mid4的表达式为：

![10fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_10.png)

```python	   
	   mid4_5=mid3.*imfilter(double(m),conj(h_C));
	   
	   mid5=0.5*(imfilter(double(mid4),conj(g_C))+imfilter(double(mid4_5),g_C));
```
mi5的表达式为：

![11fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_11.png)

```python	   
	   mid7=mid3.*double(m);
	   
	   mid8=imfilter(double(mid7),g_I);
```

mid8的表达式为：

![12fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/OPC_acaa_12.png)

计算离散罚函数和小波罚函数的梯度：

```python	   
	   %%%%%%Gradient of discretization penaly%%%%%%  
	   
	   d_D=((-8)*m+4)*(-0.5).*sin(r);
	   
	   %%%%%%Gradient of wavelet penaly%%%%%%
	   
	   for ii=0:(N/2-1)
		
		for jj=0:(N/2-1)

			   d_WA(ii*2+1,jj*2+1)= ( 3*m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+1,jj*2+1));

			   d_WA(ii*2+1,jj*2+2)= ( 3*m(ii*2+1,jj*2+2) - m(ii*2+1,jj*2+1) - m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+1,jj*2+2));

			   d_WA(ii*2+2,jj*2+1)= ( 3*m(ii*2+2,jj*2+1) - m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+2,jj*2+1));

			   d_WA(ii*2+2,jj*2+2)= ( 3*m(ii*2+2,jj*2+2) - m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) ) * (-0.5)*sin(r(ii*2+2,jj*2+2));

		   end

	   end
```

最后计算总的cost function的梯度:
```python	   
	   %%%%%%Gradient of overall cost function%%%%%%

	   d=2*a*mid5.*sin(r)+2*a*mid8.*sin(r)+gamma_D*d_D +gamma_WA*d_WA;

	end
```
