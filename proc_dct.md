# 基于二维离散余弦变换的后处理优化算法 （proc_dct）

### 本文的matlab代码来源于《Computational Lithography》一书；

proc_dct算法用于同时降低mask版图的复杂性和输出版图的误差。

二维离散余弦变换的介绍见：https://www.cnblogs.com/latencytime/p/10228938.html

### 本算法产生优化后的离散化的mask版图。

## Matlab算法

```python
	function [m_trinary_new] = proc_dct(N_mask, pz, m, t_r, t_r_real, t_m, TCC, threshold);

	error=zeros(101,1);  %The output pattern error versus the number of reserved DCT low frequency components
```
m_trinary_new：简化后的mask版图

N_mask: mask的维度

pz: Desired output pattern

m: 原始mask版图

t_r: 一阶相干近似的处理阈值

t_r_real：部分相干成像系统的处理阈值

t_m: mask的全局阈值（用于离散化mask)

TCC: 交叉透过系数

threshold: 离散余弦变换谱的高频分量的截止阈值，保留的低频分量个数为：(threshold - 1) * (threshold - 2) / 2

调用dct2()函数做二维离散余弦变换：
```python
	%%%%%%2D-DCT of gray optimized mask%%%%%%
	B=dct2(m);   %2D-DCT of gray optimized mask
```

B = dct2(A) 返回矩阵A的二维离散余弦变换. B具有与A相同的维度，B的元素为离散余弦变换系数B(k1,k2)。

保留低频分量：
```python
	B_new=B;   %2D-DCT transform after the cutting off the low frequency components
	
	for ii=102:-1:2 % ii = 102,101,...,2，频率截止阈值
		
		for p=1:N_mask
			
			for q=1:N_mask
				
				if (p+q==ii)
					
					B_new(p,q)=0;
				
				end
			
			end
		
		end
		
		% 做逆二维离散余弦变换，得到去除高频分量后的mask版图
		m_new= idct2(B_new);   %Gray mask after the cutting off the low frequency components
		
		% 根据阈值t_m离散化mask为三元版图-1，0，1
		m_trinary_p=m_new>t_m;
		
		m_trinary_n=-1*(m_new<(-1*t_m));
		
		m_trinary_new=m_trinary_p+m_trinary_n;   %Trinary mask after the cutting off the low frequency components
		
		%%%%Calculate the output pattern error corresponding to the trinary%%%%
		
		%%%%mask after the cutting off the low frequency components%%%%%%%%%%%%
		
		% 计算简化后的mask的输出版图的误差
		aerial=zeros(N_mask,N_mask);
		
		aerial_fre=zeros(N_mask,N_mask);
		
		m_trinary_new_fre=(fftshift(fft2(m_trinary_new)));
		
		for x=1:N_mask^2
			
			for y=1:N_mask^2
				
				index_1=mod(x-1,N_mask)+1;
				
				index_2=floor((x-1)/N_mask)+1;
				
				index_3=mod(y-1,N_mask)+1;
				
				index_4=floor((y-1)/N_mask)+1;
				%采用TCC计算输出版图
				aerial_fre(mod(index_1-index_3,N_mask)+1,mod(index_2-index_4,N_mask)+1)=aerial_fre(mod(index_1-index_3,N_mask)+1,mod(index_2-index_4,N_mask)+1)+TCC(x,y)*(m_trinary_new_fre(index_1,index_2))*conj(m_trinary_new_fre(index_3,index_4));
			
			end
		
		end
		
		aerial=abs(ifft2(aerial_fre))/((N_mask)^2);
		
		z_trinary=aerial>t_r_real;   %Binary output pattern after cutting off the low frequency components
		
		error(103-ii,1)=sum(sum(abs(abs(pz)-z_trinary)));   %Output pattern error after cutting off the low frequency components
		
		disp(error(103-ii,1));
	
	end

```
