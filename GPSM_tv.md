# 基于离散罚函数和全变分罚函数的相移掩膜优化算法 （GPSM_tv）

** Generalized Phase-Shift Mask Optimazation Algorithm based on discretization penalty and total variation penalty**

### 本文的matlab代码来源于《Computational Lithography》一书；

基于像素的OPC和PSM优化允许mask版图具有极大的自由度，

因此会导致优化出的mask版图过于复杂，提升了加工难度

因此研究人员引入了全变分罚函数来去除mask版图中过于

细节化的部分。

GPSM_tv与GPSM_wa算法的区别在于将小波罚函数替换成了

全变分罚函数，算法其余部分基本相同，因此在此只介绍

全变分罚函数部分。

## Matlab算法

```matlab
f=abs(m-pz);
```
f 为激活版图（activation pattern)，

f 表示mask版图的复杂度，其计算公式为：

![16fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/16.png)

f 的每个像素表示mask的强度与desired pattern的L1范式，

L1范数是指向量中各个元素绝对值之和。

下面计算全变分罚函数*Rtv

```matlab
f_right=zeros(N,N);   %Right shift of f

f_right(1:N,2:N)=f(1:N,1:N-1); % 将f整体向右平移一列，第一列元素均为0

f_up=zeros(N,N);   %Up shift of f

f_up(1:N-1,1:N)=f(2:N,1:N); % 将f整体向上平移一行，第N行元素均为0
```

全变分罚函数*Rtv的计算公式为;

![17fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/17.png)

||·||代表L1范式，Qx代表水平方向（向右）一阶导数，

Qy代表竖直（向上）一阶导数

Qxf即为f - f_right，Qyf即为f - f_up.

接下来求*Rtv的梯度：

```matlab
f1=sign( f-f_right );

f1(:,1)=0;

f2=sign( f-f_up );

f2(N,:)=0;

f1_left=zeros(N,N);   %Left shift of f1

f1_left(1:N,1:N-1)=f1(1:N,2:N);

f2_down=zeros(N,N);   %Down shift of f2

f2_down(2:N,N:N)=f2(1:N-1,N:N);

f11=f1-f1_left; % 对f1做向左的一阶微分

f11(:,N)=0;

f22=f2-f2_down; % 对f2做向下的一阶微分

f22(1,:)=0;    

f3=(f11+f22);
```

f3的表达式为：

![18fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/18.png)

对于复数域，幅度和相位的全变分罚函数的梯度要分别计算：

```matlab
%%%%%%Gradient of total variation penaly corresponding to \phi%%%%%%

dr_TV=f3*(-1).*sin(rr).*real( (m-z).*exp((-1)*i*ra) ) /2 ./abs(m-z);

%%%%%%Gradient of total variation penaly corresponding to \theta%%%%%%

da_TV=f3.*(1+cos(rr)).*real( (m-z).*exp((-1)*i*ra)*(-1)*i ) /2 ./abs(m-z);
```

其中，rr为幅度，ra为相位。

dr_TV为幅度的梯度，其表达式为：

![19fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/19.png)

W1(m)的表达式为：

![20fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/20.png)

da_TV为相位的梯度，其表达式为：

![21fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/21.png)

W2(m)的表达式为：

![22fd](https://github.com/zgzym/Computational-Lithography-Book/blob/main/images/22.png)
