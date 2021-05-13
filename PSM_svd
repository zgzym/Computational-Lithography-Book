# 基于奇异值分解模型的PSM优化算法 （PSM_svd）

*Two-phase PSM optimization using the singular value decomposition model in partially coherent imaging system*

### 本文的matlab代码来源于《Computational Lithography》一书；

PSM_svd算法用于在部分相干成像系统中采用奇异值分解模型对2种相位的相移掩膜进行梯度优化，适用于部分相干成像系统。


## 本算法产生优化后的二相相移掩膜并包含了离散罚函数。

## Matlab算法

```python
	function [] = PSM_svd(N_mask, pz, r, pixel, k, NA, lamda, sigma, order, step, a, t_r, t_m, gamma_D, epsilon, maxloop);
```
N: mask的维度

pz: Desired output pattern

r: mask的初始相位

pixel: 像素大小 (nm)

k：工艺常数

NA: 数值孔径

lamda: 波长 (nm)

sigma: 部分相干因子

order: Bessel函数的阶数

step：步长

a：sigmoid函数的陡度（steepness)

t_r: 一级相干近似的处理阈值

t_m: mask的全局阈值（显影阈值）

gamma_D：离散罚函数的权重

epsilon：输出版图的可容忍误差

maxloop：最大迭代次数

```python
	midway=(N_mask+1)/2;   % mask的中心点
	
	%% 计算交叉透过系数TCC
	[TCC] = SOCS(N_mask, pixel, k, NA, lamda, midway, sigma, order);
	
```

其中，SOCS()为计算交叉透过系数TCC的函数，代码如下：
```python
	function [TCC] = SOCS(N, pixel, k, NA, lamda, midway, sigma, order);
	
	J=zeros(N,N);   %Effective source
	
	P=zeros(N,N);   %Pupil function = fourier transform of lens shape
	
	TCC=zeros(N^2,N^2);   %Transmission cross-coefficient
```

矩阵P为光瞳函数，即透镜形状的傅里叶变换：
```python
	radius_filter=NA/lamda*pixel*N;
	
	for row=1:N
		
		for column=1:N
			
			if (radius_filter>=sqrt( (row-midway)^2 + (column-midway)^2 ) );
				
				P(row,column)=1/pi/radius_filter^2;
			
			end
		
		end
	
	end
```

其计算公式为：

![1fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_1.png)

幅度冲激响应矩阵：
```python
	h=zeros(N,N);
	
	h=fftshift(ifft2(ifftshift(P)));
```

h即为将P的零频移至矩阵中心。

有效光源：
```python
	radius_illumination=radius_filter*sigma;
	
	for row=1:N
		
		for column=1:N
			
			if (radius_illumination>=sqrt( (row-midway)^2 + (column-midway)^2 ) );
				
				J(row,column)=1/pi/radius_illumination^2;
			
			end
		
		end
	
	end
```
有效光源频域的计算公式为：

![2fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_2.png)

最后计算TCC:
```python
	P_1=P;
	
	P_2=P;
	
	for x=1:N^2
		
		for y=1:N^2
			
			index_1=mod(x-1,N)+1-midway;
			
			index_2=floor((x-1)/N)+1-midway;
			
			index_3=mod(y-1,N)+1-midway;
			
			index_4=floor((y-1)/N)+1-midway;
			
			if (index_1>0)
				
				P_1=[P(index_1+1:N,1:N);zeros(index_1,N)];
			
			end
			
			if (index_1<0)
				
				P_1=[zeros(abs(index_1),N);P(1:N+index_1,1:N)];
			
			end
			
			if (index_2>0)
				
				P_1=[P_1(1:N,index_2+1:N),zeros(N,index_2)];
			
			end
			
			if (index_2<0)
				
				P_1=[zeros(N,abs(index_2)),P_1(1:N,1:N+index_2)];
			
			end

			if (index_3>0)
				
				P_2=[P(index_3+1:N,1:N);zeros(index_3,N)];
			
			end
			
			if (index_3<0)
				
				P_2=[zeros(abs(index_3),N);P(1:N+index_3,1:N)];
			
			end
			
			if (index_4>0)
				
				P_2=[P_2(1:N,index_4+1:N),zeros(N,index_4)];
			
			end
			
			if (index_4<0)
				
				P_2=[zeros(N,abs(index_4)),P_2(1:N,1:N+index_4)];
			
			end
			
			TCC(x,y)=sum(sum(J.*P_1.*P_2));
		
		end
		
		disp(x);
	
	end
```

TCC的计算公式为：

![3fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_3.png)

接下来做奇异值分解
```python
	%%%%%%Singular value decomposition of the partially coherent imaging system%%%%%%
	[U,S,V]=svd(TCC);   %Singular value decomposition， TCC = U * S * V，S为对角矩阵，对角线上的值为本征值
	
	h_1_fre=reshape(U(1:N_mask^2,1:1),N_mask,N_mask); % TCC维度为N_mask^2，h_1_fre为将U的第一列变换为维度为N_mask的方阵
	
	h_1=(fftshift(ifft2(ifftshift((h_1_fre)))));   %The impulse response of the first order approximation 
	
	sum_eigenvalue=(sum(sum(S)));   %The summation of the eigenvalues 
```


svd(TCC)的结果为：（文献中第29页）

![4fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_4.png)

h_1的表达式为：

![5fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_5.png)

初始化参数：
```python
	t_r_real=0;   %Global threshold of photoresist effect for all of the eigen value
	
	d=zeros(N_mask,N_mask);   %Gradient of the cost function
	
	d_D=zeros(N_mask,N_mask);   %Gradient of the discretization penalty
	
	convergence=zeros(maxloop,1);   %Output pattern error in each iteration
	
	count=0;   %Index of iteration number
	
	sum6=100;   %Output pattern error corresponding to the optimized trinary mask
	
	sum8=100;   %Output pattern error corresponding to the optimized real-valued mask
```

计算幅度冲激响应矩阵的转置矩阵g_1:
```python
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
	m=zeros(N_mask,N_mask);   %Mask pattern
	
	while (sum6>epsilon) & (count<maxloop)
		
		count=count+1;
		
		%%%%%%%%%%%%%%%%%%%calculate pattern error%%%%%%%%%%%%%%%%%%
		
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
	   
		%%%%%%%Calculate whole revision vector%%%%%%%%%%
		
		d=2*a*mid5.*sin(r) + gamma_D*d_D;
		
		r=r-step*d;   %Update

	end
```

mid1的表达式为：

![6fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_6.png)

z的表达式为：

![7fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_7.png)

mid3的表达式为：

![8fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_8.png)

mid4的表达式为：

![9fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_9.png)

mid4_4的表达式为：

![10fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_10.png)

2*a*mid5.*sin(r)的表达式为：

![11fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/PSM_svd_11.png)

显影阈值为：
```python
	t_r_real=t_r*sum_eigenvalue;   %Global threshold of photoresist effect for all eigenvalues
```

