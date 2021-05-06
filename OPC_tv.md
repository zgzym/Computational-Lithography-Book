# 基于离散罚函数和全变分罚函数的OPC优化算法 （OPC_tv）

**Generalized OPC Optimazation Algorithm based on discretization penalty and total variation penalty**

### 本文的matlab代码来源于《Computational Lithography》一书；

OPC_tv算法用于在相干成像系统中对于N×N的二元掩膜进行基于梯度的优化，

本算法产生优化后的二元掩膜并包含了离散罚函数和全变分罚函数。

## Matlab算法

```python
	function [] = OPC_tv(N, N_filter, k, pz, s, a, t_r, t_m, gamma_D, gamma_TV, epsilon, maxloop);
```
N: mask的维度

N_filter: 幅度冲激响应矩阵维度

k：工艺常数

pz: 目标版图（desired pattern)

s：步长

a：sigmoid函数的陡度（steepness)

t_r: sigmoid函数的处理阈值

t_m: mask的全局阈值（显影阈值）

gamma_D：离散罚函数的权重

gamma_TV: 全变分罚函数的权重

epsilon：输出版图的可容忍误差

maxloop：最大迭代次数

```python
	d=zeros(N,N);   %cost function的梯度矩阵

	d_D=zeros(N,N);   %离散罚函数的梯度矩阵

	d_TV=zeros(N,N);   %全变分罚函数的梯度矩阵

	convergence=zeros(maxloop,1);   %每一步迭代的输出版图的误差

	count=0;   %迭代次数

	sum6=100;   %连续型掩膜版图的误差

	sum8=100;   %根据显影阈值离散化后的掩膜版图的误差
```

```python
	%%%%%%the amplitude impulse response of the coherent imaging system%%%%%%
	
	h=fspecial('gaussian',N_filter,k); % 高斯型幅度冲激响应函数
	
	%%%%%%The rotation of the amplitude impulse response%%%%%%
	
	for i=1:N_filter
		
		for j=1:N_filter
			
			h1((i-1)*N_filter+j)=h(i,j);
		
		end
	
	end
	
	for i=1:N_filter
		
		for j=1:N_filter
			
			g(i,j)=h1((N_filter-i)*N_filter+(N_filter+1-j));
		
		end
	
	end
```

其中g为h的转置。


```python
	r=pi*4/5*(pz==0) + pi/5*(pz==1); % 初始化mask相位

	m=zeros(N,N);    % Mask pattern
```

接下来的代码为迭代主循环：
```python
	while (sum6>epsilon) & (count<maxloop) % 迭代至误差小于epsilon或者达到最大迭代次数
	
		count=count+1; 
		
		r=r-s*d;   %Update
		
		m=(1+ cos(r))/2;   %Update，连续型mask，这里不存在相位变化
		
		viccin=m>t_m;   %Binary mask，根据全局阈值t_m离散化的mask
		
		viccout=imfilter(viccin,h); % 与冲激响应函数h做卷积，这一步即为成像计算
		
		viccbin=viccout>t_r;   %Output pattern of binary mask，根据sigmoid函数的阈值离散化mask
		
		sum6=sum(sum(abs(pz-viccbin)));  % 计算误差，这里为一阶范数
		
		convergence(count,1)=sum6;  % 记录当前循环的误差

		mid1=imfilter(m,h);
		
		z=1./ (  1+exp(-1*a*mid1+a*t_r)  ); 
		
		mid3=(pz-z).*z.*(1-z);   
		
		mid5=imfilter(mid3,g);

		%%%%%%Gradient of discretization penaly%%%%%%  
		
		d_D=( (-8)*m+4 )*(-0.5).*sin(r);

		%%%%%%Gradient of total variation penaly%%%%%%
		
		% total variation penalty的介绍在文献中93页
		
		f=abs(m-pz); % active pattern，即mask与目标版图之间的总差异，表征mask的复杂度
		
		f_right=zeros(N,N);   %Right shift of f
		
		f_right(1:N,2:N)=f(1:N,1:N-1);
		
		f_up=zeros(N,N);   %Up shift of f
		
		f_up(1:N-1,1:N)=f(2:N,1:N);
		
		f1=sign( f-f_right );
		
		f1(:,1)=0;
		
		f2=sign( f-f_up );
		
		f2(N,:)=0;

		f1_left=zeros(N,N);   %Left shift of f1
		
		f1_left(1:N,1:N-1)=f1(1:N,2:N);
		
		f2_down=zeros(N,N);   %Down shift of f2
		
		f2_down(2:N,N:N)=f2(1:N-1,N:N);
		
		f11=f1-f1_left;
		
		f11(:,N)=0;
		
		f22=f2-f2_down;
		
		f22(1,:)=0;

		d_TV=(f11+f22).*sign(m-pz)*(-0.5).*sin(r);

		%%%%%%Gradient of overall cost function%%%%%%
		
		d=2*a*mid5.*sin(r) +gamma_D*d_D +gamma_TV*d_TV;

	end
```



mid1的计算公式为：

![23fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/23.png)

z的表达式为如下，去掉exp(jθ)项：

![24fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/24.png)

mid3的表达式为：

![25fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/25.png)

mid5的表达式为：

![26fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/26.png)

离散罚函数的梯度公式为：

![27fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/27.png)

全变分罚函数的梯度公式为：

![28fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/28.png)

其中,Qx和Qy分别表示向右和向上的求导
