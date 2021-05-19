# 光源掩膜协同优化算法

### 本文的matlab代码来源于《Computational Lithography》一书；

SMO_OPC算法用于对原始的环形照明光源进行光源和mask的协同优化。

## Matlab算法

```python
	function [] = smo_OPC(N, pz, N_filter, pixel, k, NA, lamda, order, sigma_large_inner, sigma_large_outer);
```
N: mask的维度

pz: Desired output pattern

N_filter: 冲激响应矩阵的维度

pixel: 像素尺寸(nm)

k: 工艺常数

NA: 数值孔径

lamda: 波长(nm)

order：Bessel函数的阶数

sigma_large_inner：内部部分相干因子

sigma_large_outer: 外部部分相干因子


定义冲激响应矩阵h：
```python	                
	midway=(N_filter+1)/2;   %middle of low pass filter

	%%%%%%the amplitude impulse response of the partially coherent imaging system%%%%%%
	h=zeros(N_filter,N_filter);
	
	radius=0;
	
	for row=1:N_filter
		
		for column=1:N_filter
			
			radius=pixel*sqrt( (row-midway)^2 + (column-midway)^2 );
			
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

```

计算h的转置矩阵g：
```python
	for ii=1:N_filter
		
		for jj=1:N_filter
			
			h_vector((ii-1)*N_filter+jj)=h(ii,jj);
		
		end
	
	end
	
	for ii=1:N_filter
		
		for jj=1:N_filter
			
			g(ii,jj)=h_vector((N_filter-ii)*N_filter+(N_filter+1-jj)); %inverse vector
		
		end
	
	end
```

初始化光源：
```python
	%%%%%%%Initial illumination pattern%%%%%%%
	D=pixel*N;
	
	D_C_1=lamda/2/sigma_large_outer/NA;   %Coherence length
	
	D_C_2=lamda/2/sigma_large_inner/NA;   %Coherence length
	
	omega_0=pi/D;
	
	N_coherence=18;   %Source dimension
	
	midway_coherence=(N_coherence+1)/2;   %Middle point of illumination
	
	radius_1=D/(2*D_C_1);   %Inner radius of annular illumination
	
	radius_2=D/(2*D_C_2);   %Outer radius of annular illumination 
	
	yita=zeros(N_coherence,N_coherence);   %Illumination pattern
	
	for row=1:N_coherence
		
		for column=1:N_coherence
			
			radius=pixel*sqrt( (row-midway_coherence)^2 + (column-midway_coherence)^2 );
			
			if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
				
				yita(row,column)=1;
			
			end
		
		end
	
	end
	
	normalize=sum(sum(yita));   %Normalization factor

	yita_initial=zeros(N_coherence,N_coherence);   %Initial illumination pattern
	
	yita_initial=yita;
```
环形光源初始化的注释可参考OPC_acc或OPC_fse.

迭代参数初始化：
```python
	%%%%%%Initialization of optimization%%%%%%
	count=0;   %Iteration number
	
	m=pz;   %Initial mask
	
	direction_mask=zeros(N,N);   %Cost sensitivity function of mask
	
	direction_source=zeros(N_coherence,N_coherence);   %Cost sensitivity function of source
	
	flag_mask_pre=ones(N,N);   %Locations of the changable pixels on mask in the previous iteration
	
	flag_source_pre=ones(N_coherence,N_coherence);   %Locations of the changable pixels on source in the previous iteration
	
	flag_mask=zeros(N,N);   %Locations of the changable pixels on mask in the current iteration
	
	flag_source=zeros(N_coherence,N_coherence);   %Locations of the changable pixels on source in the current iteration
	
	error=10;   %Output pattern error in the current iteration
```

计算初始化mask（即desired pattern）的输出版图的误差：
```python
	%%%%%%Calculate the output pattern error in the previous iteration%%%%%%
	aerial=zeros(N,N);
	
	normalize=sum(sum(yita));
	
	for p=1:N_coherence
		
		for q=1:N_coherence    
			
			if (yita(p,q)>0)
				
				exponential=zeros(N_filter,N_filter); % has the same dimension as the filter h
				
				for row=1:N_filter
					
					for column=1:N_filter
						
						argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
						
						exponential(row,column)=exp(i*omega_0*argument);
					
					end
				
				end
				
				aerial=aerial+yita(p,q)/ normalize* abs(  imfilter(double(m),h.*exponential)  ).^2;
			
			end
		
		end
	
	end
	
	error_pre=sum(sum((pz-aerial).^2));   %Output pattern error in the previous iteration
```

下面为主循环：
```python
	%%%%%%Source and mask optimization%%%%%%
	while (error > 0)   
	   
	   count=count+1;   %Update the iteration number
	   
	   disp(strcat('Iteration=',num2str(count)));
	   
	   direction_mask=zeros(N,N);   %Reset the cost sensitivity function of mask
	   
	   direction_source=zeros(N_coherence,N_coherence);   %Reset the cost sensitivity function of source
	   
	   flag_mask=zeros(N,N);   %Reset the locations of the changable pixels on mask in the current iteration
	   
	   flag_source=zeros(N_coherence,N_coherence);   %Reset the locations of the changable pixels on source in the current iteration
```

计算aerial image:
```python
	   %%%%%%Calculate the aerial image%%%%%%
	   aerial=zeros(N,N);
	   
	   normalize=sum(sum(yita));
	   
	   for p=1:N_coherence
		   
		   for q=1:N_coherence    
			   
			   if (yita(p,q)>0)
				   
				   exponential=zeros(N_filter,N_filter); 
				   
				   for row=1:N_filter
					   
					   for column=1:N_filter
						   
						   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
						   
						   exponential(row,column)=exp(i*omega_0*argument);
					   
					   end
				   
				   end
				   
				   aerial=aerial+yita(p,q)/ normalize* abs(  imfilter(double(m),h.*exponential)  ).^2;
			   
			   end
		   
		   end
	   
	   end
```

计算mask的成本敏感函数：
```python
	   %%%%%%Calculate the cost sensitivity function of mask%%%%%%
	   for p=1:N_coherence
		   
		   for q=1:N_coherence   
			   
			   if (yita(p,q)>0)
				   
				   exponential=zeros(N_filter,N_filter);
				   
				   for row=1:N_filter
					   
					   for column=1:N_filter
						   
						   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
						   
						   exponential(row,column)=exp(i*omega_0*argument);
					   
					   end
				   
				   end
				   
				   direction_mask=direction_mask+(2)*yita(p,q)/normalize* real(  imfilter( (pz-aerial).*imfilter(double(m),h.*exponential) , conj(g.*exponential) )  +   imfilter( 
				   (pz-aerial).*imfilter(double(m),conj(h.*exponential)) , g.*exponential )    );
			   
			   end
		   
		   end
	   
	   end
```

其中，direction_mask的表达式为：

![1fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/SMO_OPC_1.png)

计算source的cost sensitivity function:
```python
	   %%%%%%Calculate the cost sensitivity function of source%%%%%%
	   normalize=sum(sum(yita));
	   
	   for p=1:N_coherence
		   
		   for q=1:N_coherence   
			   
			   exponential=zeros(N_filter,N_filter);
			   
			   for row=1:N_filter
				   
				   for column=1:N_filter
					   
					   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
					   
					   exponential(row,column)=exp(i*omega_0*argument);
				   
				   end
			   
			   end 
			   
			   direction_source(p,q)=-1*sum(sum( (pz-aerial).* abs(imfilter(double(m),h.*exponential)).^2 )) *(-2)/normalize;   
		   
		   end
	   
	   end
```
其中，direction_source的表达式为：

![2fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/SMO_OPC_2.png)

计算mask版图中可变像素点的位置：
```python
	   %%%%%%Calculate the locations of the changable pixels on mask%%%%%%
	   for ii=2:N-1
		   
		   for jj=2:N-1
			   
			   if (m(ii,jj-1)+m(ii,jj+1)+m(ii-1,jj)+m(ii+1,jj)==0)|(m(ii,jj-1)+m(ii,jj+1)+m(ii-1,jj)+m(ii+1,jj)==4)
				   
				   flag_mask(ii,jj)=0;
			   
			   else
				   
				   flag_mask(ii,jj)=1;
			   
			   end
		   
		   end
	   
	   end
```
计算source版图中可变像素点的位置：
```python
	   %%%%%%Calculate the locations of the changable pixels on source%%%%%%
	   for ii=2:N_coherence-1
		   
		   for jj=2:N_coherence-1
			   
			   if (yita(ii,jj-1)+yita(ii,jj+1)+yita(ii-1,jj)+yita(ii+1,jj)==0)|(yita(ii,jj-1)+yita(ii,jj+1)+yita(ii-1,jj)+yita(ii+1,jj)==4)
				   
				   flag_source(ii,jj)=0;
			   
			   else
				   
				   flag_source(ii,jj)=1;
			   
			   end
		   
		   end
	   
	   end
```

可变像素点为4-boundary pixel，将其翻转不会引入孤立像素点；

4-boundary pixel是指其至少有1个相邻像素点与其值相异。

孤立像素点是指与其4-neighbors的值均不同的像素点；否则，更新可变像素点位置的矩阵。

假如不存在可变像素点，则迭代终止：
```python
	   %%%%%%If no pixel on mask or souce changed, then the program is terminated%%%%%%
	   if(sum(sum(flag_mask~=flag_mask_pre))+sum(sum(flag_source~=flag_source_pre)) == 0)
		   
		   break;
	   
	   else
		   
		   flag_mask_pre=flag_mask;   %Update the locations of the changable pixels on mask
		   
		   flag_source_pre=flag_source;   %Update the locations of the changable pixels on source
	   
	   end
```

寻找梯度值最大的可变像素点来减少输出版图的误差：
```python
		while (sum(sum(flag_mask))+sum(sum(flag_source))>0)   %若mask或source版图中存在可变像素点，则继续循环
		   
		   test_mask=abs(direction_mask.*flag_mask);   %mask版图中可变像素点的梯度大小
		   
		   % find(X) 返回一个包含X中非零元素线性索引值的矢量（线性索引值表示从上往下，再从左往右的索引）.
		   
		   max_index_mask=find(  test_mask >= max(max(test_mask))   );   %找到mask版图中梯度值最大的可变像素点
		   
		   max_index_mask=max_index_mask(1,1);   %选择第一个具有最大梯度的可变像素
		   
		   max_row_mask=mod(max_index_mask,N);   %确定可变像素在mask版图中的行数
		   
		   if (max_row_mask==0)
			   
			   max_row_mask=N;
		   
		   end
		   
		   max_column_mask=floor(max_index_mask/N)+1;   %确定可变像素在mask版图中的列数

		   test_source=abs(direction_source.*flag_source);   %Magnitude of the gradient for the changable pixels on source
		   
		   max_index_source=find(   test_source >= max(max(test_source))   );   %Find the changable pixel with maximum gradient magnitude on source
		   
		   max_index_source=max_index_source(1,1);   %Choose a changable pixel with maximum gradient magnitude on source
		   
		   max_row_source=mod(max_index_source,N_coherence);   %Determine the row of the changable pixel on source
		   
		   if (max_row_source==0)
			   
			   max_row_source=N_coherence;
		   
		   end
		   
		   max_column_source=floor(max_index_source/N_coherence)+1;   %Determine the column of the changable pixel on source
```

首先翻转mask版图上的可变像素：
```python
		   %Change the pixel on mask first%
		   % 若该点的mask版图的梯度大于source版图的梯度，则先翻转mask版图上的可变像素点
		   if (test_mask(max_row_mask,max_column_mask) >= (test_source(max_row_source,max_column_source)))
			   
			   flag_mask(max_row_mask,max_column_mask)=0;   %在后面的迭代过程中不再翻转该点
			   
			   % 若该点在mask版图上值为0并且direction_mask的值大于0，或者在mask版图上值为1并且direction_mask的值小于0，则将该点翻转
			   if ((m(max_row_mask,max_column_mask)==0)&(direction_mask(max_row_mask,max_column_mask)>0)) | 
			       ((m(max_row_mask,max_column_mask)==1)&(direction_mask(max_row_mask,max_column_mask)<0))
				   
				   m(max_row_mask,max_column_mask)=double(~m(max_row_mask,max_column_mask));   %Flip the pixel according to the cost sensitivity and current pixel value
				   
				   %%%%%%Check topological constraint on mask%%%%%%
				   %Flag of the topological constraint on mask
				   %If maxu=0, then the topological constraint on mask is satisfied
				   %If maxu=1, then the topological constraint on mask is not satisfied
				   maxu=0; 
				   % 检查在(max_row_mask,max_column_mask)周围是否存在singular pixel
				   for ii=max(2,max_row_mask-1):min(N-1,max_row_mask+1)
					   if (maxu==1)
						   break;
					   end
					   for jj=max(2,max_column_mask-1):min(N-1,max_column_mask+1)
						   % 若(ii,jj)为singular pixel，则违反了拓扑规则，maxu=1
						   if (m(ii,jj)~=m(ii,jj-1)) & (m(ii,jj)~=m(ii,jj+1)) & (m(ii,jj)~=m(ii-1,jj)) & (m(ii,jj)~=m(ii+1,jj))
							   maxu=1;
							   break;
						   end
					   end
				   end
				   if (maxu==1)   %若违反了拓扑规则，则将该点撤销翻转
					   m(max_row_mask,max_column_mask)=double(~m(max_row_mask,max_column_mask));
					   continue;
				   end
```

计算翻转像素后输出版图的误差是否减小：
```python
	               %%%%%%Check the change of output pattern error%%%%%%
				   normalize=sum(sum(yita));
				   
				   aerial=zeros(N,N);
				   
				   for p=1:N_coherence
					   
					   for q=1:N_coherence    
						   
						   if (yita(p,q)>0)
							   
							   exponential=zeros(N_filter,N_filter); 
							   
							   for row=1:N_filter
								   
								   for column=1:N_filter
									   
									   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
									   
									   exponential(row,column)=exp(i*omega_0*argument);
								   
								   end
							   
							   end
							   
							   aerial=aerial+yita(p,q)/ normalize* abs(  imfilter(double(m),h.*exponential)  ).^2;
						   
						   end
					   
					   end
				   
				   end
				   
				   error=sum(sum((pz-aerial).^2));   %Output pattern error in the current iteration
		
				   if(error>=error_pre)   %If the error is increased, then restore the original pixel values
					   
					   m(max_row_mask,max_column_mask)=double(~m(max_row_mask,max_column_mask));
				   
				   else   %If the error is reduced, then keep the pixel update, and update the output patter error in the previous iteration
					   
					   error_pre=error;
					   
				   end
			  
			  end
```
% 若该点的mask版图的梯度不大于source版图的梯度，则先翻转source版图上的可变像素点：
```python
	       %%%%%%Change the pixel on source first%%%%%%		   
		   else
			   
			   flag_source(max_row_source,max_column_source)=0;   %Do not consider this pixel in the following iterations
			   
			   if ((yita(max_row_source,max_column_source)==0) & (direction_source(max_row_source,max_column_source)>0)) | ((yita(max_row_source,max_column_source)==1) & (direction_source(max_row_source,max_column_source)<0))
				   
				   yita(max_row_source,max_column_source)=double(~yita(max_row_source,max_column_source));   %Flip the pixel according to the cost sensitivity and current pixel value
				   
				   %%%%%%Check topological constraint on source%%%%%%
				   %Flag of the topological constraint on source
				   %If maxu=0, then the topological constraint on source is satisfied
				   %If maxu=1, then the topological constraint on mask is not satisfied
				   maxu=0;   
				   
				   for ii=max(2,max_row_source-1):min(N_coherence-1,max_row_source+1)
					   
					   if (maxu==1)
						   
						   break;
					   
					   end
					   
					   for jj=max(2,max_column_source-1):min(N_coherence-1,max_column_source+1)
						   
						   if (yita(ii,jj)~=yita(ii,jj-1)) & (yita(ii,jj)~=yita(ii,jj+1)) & (yita(ii,jj)~=yita(ii-1,jj)) & (yita(ii,jj)~=yita(ii+1,jj))
							   
							   maxu=1;
							   
							   break;
						   
						   end
					   
					   end
				   
				   end
				   
				   if (maxu==1)   %If the topological constraint on mask is not satisfied, then restore the original pixel values
					   
					   yita(max_row_source,max_column_source)=double(~yita(max_row_source,max_column_source));
					   
					   continue;
				   
				   end
				   
				   %%%%%%Check the change of output pattern error%%%%%%
				   normalize=sum(sum(yita));
				   
				   aerial=zeros(N,N);
				   
				   for p=1:N_coherence
					   
					   for q=1:N_coherence    
						   
						   if (yita(p,q)>0)
							   
							   exponential=zeros(N_filter,N_filter); 
							   
							   for row=1:N_filter
								   
								   for column=1:N_filter
									   
									   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
									   
									   exponential(row,column)=exp(i*omega_0*argument);
								   
								   end
							   
							   end
							   
							   aerial=aerial+yita(p,q)/ normalize* abs(  imfilter(double(m),h.*exponential)  ).^2;
						   
						   end
					   
					   end
				   
				   end
				   
				   error=sum(sum((pz-aerial).^2));   %Output pattern error in the current iteration
				   
				   if(error>=error_pre)   %If the error is increased, then restore the original pixel values 
					   
					   yita(max_row_source,max_column_source)=double(~yita(max_row_source,max_column_source));
				   
				   else   %If the error is reduced, then keep the pixel update, and update the output patter error in the previous iteration
					   
					   error_pre=error;
					   
				   end
			   
			   end
		   
		   end
	   
	   end
	
	end
```
