02 透视投影模型

## 2.1 针孔相机模型

假设某点的世界坐标系：${\mathbf{X}} = {\left[ {X,Y,Z} \right]^T}$，$f$ 为焦距，成像平面上像点为：$[x,y]$，那么显然有：
$$
\frac{X}{Z} = \frac{x}{f},\frac{Y}{Z} = \frac{y}{f}
$$
于是在成像面上的投影点：
$$
x = f\frac{X}{Z},y = f\frac{Y}{Z} \label{2}
$$
或者写为向量形式：
$$
{\bf{x}} = \left[ {\begin{array}{*{20}{c}}
  x \\ 
  y 
\end{array}} \right] = \frac{f}{Z}\left[ {\begin{array}{*{20}{c}}
  X \\ 
  Y 
\end{array}} \right] \label{3}
$$
世界坐标系点：${\bf{X}}_i$ 在相机光学中心 $C=[0,0,0]^T$ 和相应的图像点 ${\bf{x}}_i$，即（$\lambda > 1$）：
$$
{X_i} = \lambda \left[ {{x_i} - C} \right] = \lambda {x_i}
$$

## 2.2 投影矩阵

方程 $\ref{2}$ 和 $\ref{3}$ 描述了笛卡尔坐标域中的非线性变换。使用齐次坐标系，投影转换可以写作线性矩阵方程：
$$
\left[ {\begin{array}{*{20}{c}}
  x \\ 
  y 
\end{array}} \right] = \frac{f}{Z}\left[ {\begin{array}{*{20}{c}}
  X \\ 
  Y 
\end{array}} \right] \equiv \left[ {\begin{array}{*{20}{c}}
  {f \cdot X/Z} \\ 
  {f \cdot Y/Z} \\ 
  1 
\end{array}} \right] \equiv \left[ {\begin{array}{*{20}{c}}
  {fX} \\ 
  {fY} \\ 
  Z 
\end{array}} \right] \equiv \underbrace {\left[ {\begin{array}{*{20}{c}}
  f&0&0&0 \\ 
  0&f&0&0 \\ 
  0&0&1&0 
\end{array}} \right]}_{{M_p}}\left[ {\begin{array}{*{20}{c}}
  X \\ 
  Y \\ 
  Z \\ 
  1 
\end{array}} \right]
$$
或者写地更紧凑一些：
$$
\mathbf{x} = {\hom ^{ - 1}}\left( {{\mathbf{M}_p} \cdot \hom \left( {\mathbf{X}} \right)} \right)
$$
投影矩阵 $M_p$ 可以分解为两个矩阵：$M_f$ 和 $M_0$ 在形式上：
$$
{M_P} = \left[ {\begin{array}{*{20}{c}}
  f&0&0&0 \\ 
  0&f&0&0 \\ 
  0&0&1&0 
\end{array}} \right] = \underbrace {\left[ {\begin{array}{*{20}{c}}
  f&0&0 \\ 
  0&f&0 \\ 
  0&0&1 
\end{array}} \right]}_{{M_f}}\underbrace {\left[ {\begin{array}{*{20}{c}}
  1&0&0&0 \\ 
  0&1&0&0 \\ 
  0&0&1&0 
\end{array}} \right]}_{{M_0}} = {M_f} \cdot {M_0}
$$
其中：$M_f$ 代表（理想）针孔相机的内部与焦距 $f$，$M_0$ 代表相机坐标系和世界坐标系之间的转换。特别地，$M_0$ 经常被称为标准投影矩阵，当如我们假设的一样，光轴与Z轴对齐。

## 2.3 刚体运动

如果照相机有它自己的（非典型的）坐标系，它就会观察到受刚体运动影响的三维点，如A.2所叙述。于是投影矩阵 $M_P$（公式5）被应用来修正（旋转和平移），而非原始点 $\underset{\raise0.3em\hbox{$\smash{\scriptscriptstyle-}$}}{X}  = \hom \left( X \right)$：
$$
x = {\hom ^{ - 1}}\left[ {{M_p}\underset{\raise0.3em\hbox{$\smash{\scriptscriptstyle-}$}}{X} '} \right] = {\hom ^{ - 1}}\left[ {{M_p}{M_{rb}}\hom \left( X \right)} \right]
$$
其中：$M_{rb}$ 代表3D刚体运动，因此，焦距为f的理想针孔相机在刚性运动下的完整透视成像变换可以写为：
$$
x = {\hom ^{ - 1}}\left[ {{M_f} \cdot {M_0} \cdot {M_{rb}} \cdot \hom \left( X \right)} \right]
$$
或者合并$M_0$ 和 $M_{rb}$：
$$
\begin{array}{l}
  \left[ {\begin{array}{*{20}{c}}
  x \\ 
  y 
\end{array}} \right] 

&= {\hom ^{ - 1}}\left( {\left[ {\begin{array}{*{20}{c}}
  f&0&0 \\ 
  0&f&0 \\ 
  0&0&1 
\end{array}} \right]\left[ {\begin{array}{*{20}{c}}
  1&0&0&0 \\ 
  0&1&0&0 \\ 
  0&0&1&0 
\end{array}} \right]\left[ {\begin{array}{*{20}{c}}
  {{r_{11}}}&{{r_{12}}}&{{r_{13}}}&{{t_x}} \\ 
  {{r_{21}}}&{{r_{22}}}&{{r_{23}}}&{{t_y}} \\ 
  {{r_{31}}}&{{r_{32}}}&{{r_{33}}}&{{t_z}} \\ 
  0&0&0&1 
\end{array}} \right]\left[ {\begin{array}{*{20}{c}}
  X \\ 
  Y \\ 
  Z \\ 
  1 
\end{array}} \right]} \right) \hfill \\
   
   
   &= {\hom ^{ - 1}}\left( {\underbrace {\left[ {\begin{array}{*{20}{c}}
  f&0&0 \\ 
  0&f&0 \\ 
  0&0&1 
\end{array}} \right]}_{{M_f}}\underbrace {\left[ {\left. {\begin{array}{*{20}{c}}
  {{r_{11}}}&{{r_{12}}}&{{r_{13}}} \\ 
  {{r_{21}}}&{{r_{22}}}&{{r_{23}}} \\ 
  {{r_{31}}}&{{r_{32}}}&{{r_{33}}} 
\end{array}} \right|\begin{array}{*{20}{c}}
  {{t_x}} \\ 
  {{t_y}} \\ 
  {{t_z}} 
\end{array}} \right]}_{Rt}.\left[ {\begin{array}{*{20}{c}}
  X \\ 
  Y \\ 
  Z \\ 
  1 
\end{array}} \right]} \right) \hfill \\
   
\end{array}
$$

$$
= {\hom ^{ - 1}}\left( {{M_f}\left[ {R|t} \right]\hom \left( X \right)} \right) 
$$

如果 $f=1$，那么 $M_f$ 成为单位矩阵，可以被表示为：
$$
{\bf{x}} = {\hom ^{ - 1}}\left( {\left[ {R|t} \right]\hom \left( X \right)} \right)
$$
在后续，它被称为：标准投影（ “normalized projection”）。

## 2.4 相机内参

在Eqn（11）中进行透视成像变换作为真实相机的模型。特别是，我们需要定义如何通过考虑来将图像平面上的连续 $x/y$ 坐标映射到实际的像素坐标：

- 传感器尺度（可能不同）：$s_x,s_y$在 $x,y$ 方向
- 图像中心：$u_c=(u_c,v_c)$
- 倾斜畸变：$s_{\theta}$ 图像平面（通常是可以忽略不计或零）

最终传感器坐标：$u=[u,v]^T$可以从归一化的图像坐标系 ${\bf{x}}=(x,y)^T$（公式12）计算得到：
$$
\left[ {\begin{array}{*{20}{c}}
  u \\ 
  v 
\end{array}} \right] = {\hom ^{ - 1}}\left( {\left[ {\begin{array}{*{20}{c}}
  {{s_x}}&{{s_\theta }}&{{u_0}} \\ 
  0&{{s_y}}&{{v_0}} \\ 
  0&0&1 
\end{array}} \right]\left[ {\begin{array}{*{20}{c}}
  f&0&0 \\ 
  0&f&0 \\ 
  0&0&1 
\end{array}} \right] \cdot \left[ {\begin{array}{*{20}{c}}
  x \\ 
  y \\ 
  1 
\end{array}} \right]} \right)
$$

$$
\left[ {\begin{array}{*{20}{c}}
  u \\ 
  v 
\end{array}} \right] = {\hom ^{ - 1}}\left( {A \cdot \hom \left( {\bf{x}} \right)} \right)
$$

其中：
$$
A = \left[ {\begin{array}{*{20}{c}}
  {f{s_x}}&{{s_\theta }}&{{u_0}} \\ 
  0&{f{s_y}}&{{v_0}} \\ 
  0&0&1 
\end{array}} \right] = \left[ {\begin{array}{*{20}{c}}
  \alpha &\gamma &{{u_0}} \\ 
  0&\beta &{{v_0}} \\ 
  0&0&1 
\end{array}} \right]
$$
称为相机内参矩阵。将额外的相机内参考虑进去，完整的透视图像变换可以写作：
$$
\left[ {\begin{array}{*{20}{c}}
  u \\ 
  v 
\end{array}} \right] = {\hom ^{ - 1}}\left( {\underbrace {\left[ {\begin{array}{*{20}{c}}
  \alpha &\gamma &{{u_0}} \\ 
  0&\beta &{{v_0}} \\ 
  0&0&1 
\end{array}} \right]}_A\underbrace {\left[ {\begin{array}{*{20}{c}}
  {{r_{11}}}&{{r_{12}}}&{{r_{13}}}&{{t_x}} \\ 
  {{r_{21}}}&{{r_{22}}}&{{r_{23}}}&{{t_y}} \\ 
  {{r_{31}}}&{{r_{32}}}&{{r_{33}}}&{{t_z}} 
\end{array}} \right]}_{W = \left( {R|t} \right)}\left[ {\begin{array}{*{20}{c}}
  {{X_i}} \\ 
  {{Y_i}} \\ 
  {{Z_i}} \\ 
  1 
\end{array}} \right]} \right)
$$

$$
= {\hom ^{ - 1}}\left( {A \cdot W \cdot \hom \left( X \right)} \right)
$$

其中 $A$ 捕捉相机的内在属性（“内在”），$W=(R|t)$ 是投影变换的外部参数（“外部”）。现在我们计算公式17通过两步：

- **Step 1.**  计算标准投影 ${\bf{x}} = (x,y)^T$（如公式12）
  $$
  \left[ {\begin{array}{*{20}{c}}
    x \\ 
    y 
  \end{array}} \right] = {\hom ^{ - 1}}\left( {W \cdot \hom \left( X \right)} \right)
  $$

  $$
   = {\hom ^{ - 1}}\left( {\left[ {\begin{array}{*{20}{c}}
    {{r_{11}}}&{{r_{12}}}&{{r_{13}}}&{{t_x}} \\ 
    {{r_{21}}}&{{r_{22}}}&{{r_{23}}}&{{t_y}} \\ 
    {{r_{31}}}&{{r_{32}}}&{{r_{33}}}&{{t_z}} 
  \end{array}} \right]\left[ {\begin{array}{*{20}{c}}
    {{X_i}} \\ 
    {{Y_i}} \\ 
    {{Z_i}} \\ 
    1 
  \end{array}} \right]} \right)
  $$

  $$
  = \hat P\left( {W,X} \right)
  $$

- **Step 2.** 从规范化的坐标 $\bf{x}$ 到传感器坐标：${\bf{u}}=[u,v]^T$ 通过公式14：
  $$
  \left[ {\begin{array}{*{20}{c}}
    u \\ 
    v 
  \end{array}} \right] = {\hom ^{ - 1}}\left( {A \cdot \hom \left( x \right)} \right) = A' \cdot \hom \left( x \right)
  $$

  $$
   = {\hom ^{ - 1}}\left( {\left[ {\begin{array}{*{20}{c}}
    \alpha &\gamma &{{u_0}} \\ 
    0&\beta &{{v_0}} \\ 
    0&0&1 
  \end{array}} \right]\left[ {\begin{array}{*{20}{c}}
    x \\ 
    y \\ 
    1 
  \end{array}} \right]} \right) = \underbrace {\left[ {\begin{array}{*{20}{c}}
    \alpha &\gamma &{{u_0}} \\ 
    0&\beta &{{v_0}} 
  \end{array}} \right]}_{A'}\left[ {\begin{array}{*{20}{c}}
    x \\ 
    y \\ 
    1 
  \end{array}} \right]
  $$

其中：$A'$ 是 $A$ 上部分 $2 \times 3$ 的子矩阵。通过使用 $A'$，不需要进行笛卡尔坐标系转换。$A$ 和 $A'$ 是2D仿射变换。通过合并这两个步骤，整个从 $3D$ 到 $2D$ 的投影过程总结如下（从世界坐标系 $\bf{X}$ 到传感器坐标：$\bf{u}$）：
$$
u = A'\hom \left( x \right) = A'\hom \left( {{{\hom }^{ - 1}}\left( {W \cdot \hom \left( X \right)} \right)} \right) = P\left( {A,W,X} \right)
$$
我们称 $P(A,W,X)$ 为投影矩阵，将3D点 $X=(X,Y,Z)^T$（世界坐标系）转换到2D传感器点：${\bf{u}} = (u,v)^T$，使用内参$A$ 和外参 $W$。这个函数可以以以下形式分为两个组件函数：
$$
P\left( {A,W,X} \right) = \left[ {\begin{array}{*{20}{c}}
  {{P_x}\left( {A,W,X} \right)} \\ 
  {{P_y}\left( {A,W,X} \right)} 
\end{array}} \right] = \left[ {\begin{array}{*{20}{c}}
  u \\ 
  v 
\end{array}} \right] = {\bf{u}}
$$
我们将在下面的步骤中建立这个符号。

## 2.5 镜头畸变

到目前为止，我们已经依赖于朴素的针孔相机模型，它除了上述的射影变换之外没有显示出任何畸变。真正的相机是用镜头而不是针孔建造的，这些都引入了额外的几何畸变，包括两个主要部分[8，p342]：

- **切向畸变**：由透镜中心从光轴上的位移引起的（这主要是由Eq 13中的变量偏移量$(u_c，v_c)$来处理的。(13)）
- **径向畸变**：由光折射变化引起的径向畸变，这在广角镜头中很明显（“桶畸变”）。

虽然透镜畸真一般是一种复杂的物理现象，但它通常以足够的精度建模为距离透镜中心径向距离$r$的单变量多项式函数$D(r)$ [8，p343](见Sec。2.5.2, Eqn. (32)).

### 2.5.1 镜头的失真是从哪里来的?

镜头失真会影响归一化投影坐标 $x$，即在图像到传感器的转换之前（定义通过内参）。在我们研究实际的失真模型之前，我们定义了一个一般的失真函数扭曲：${\mathbb{R}^2} \to {\mathbb{R}^2}$，它将映射非扭曲的2D坐标：$\bf{x}$ 到扭曲的坐标：$\tilde x$（同样是在标准化的投影平面上）：

> 畸变发生在相机坐标系（3D，光心为原点）到图像坐标系（2D，光心为原点）的转换过程中，即影响到 $x$. 

$$
\tilde x = warp\left( {x,k} \right)
$$

其中：$k$ 是畸变的参数。有了这个定义，我们就可以在方程式（23）中重新表述投影过程以包括镜头失真为：
$$
{\bf{u}} = A'\hom \left[ {\tilde x} \right] = A'\hom \left[ {warp\left( {x,k} \right)} \right]
$$

$$
= A'\hom \left[ {warp\left( {{{\hom }^{ - 1}}\left[ {W \cdot \hom \left( X \right)} \right],k} \right)} \right]
$$

$$
 = P\left( {A,k,W,X} \right)
$$

在下面，我们将描述如何扭曲函数在Eqns（26）-（27）)被指定和计算。

### 2.5.2 径向畸变模型

径向模型最常用于校正几何透镜变形。通过径向畸变，我们理解位移仅限于从图像中心发出的径向线；径向位移的量（向内或向外）仅是半径的函数。在归一化投影平面中，光轴与成像面的交点：$x_c=(0,0)$，假设它是镜头畸变中心。投影点 ${\bf{x}}=(x,y)^T$ 径向距离 $r_i$ 可以被简单计算：
$$
{r_i} = \left\| {{x_i} - {x_c}} \right\| = \left\| {{x_i}} \right\| = \sqrt {x_i^2 + y_i^2}
$$
对于某个点$x_i$ 它的径向畸变只依赖它的径向距离：$r_i$，因此，失真模型可以由一个单变量函数来指定：$D(r,k)$，因此径向扭曲：
$$
\tilde r = {f_{rad}}\left( r \right) = r \cdot \left[ {1 + D\left( {r,k} \right)} \right]
$$
因此（公式25），扭曲投影点是${\tilde x_i} = warp\left( {{x_i},k} \right)$：
$$
warp\left( {{x_i},k} \right) = {x_i} \cdot \left[ {1 + D\left( {\left\| {{x_i}} \right\|,k} \right)} \right]
$$
方程 $D(r,k)$ 指定了（正向或者逆向）的 偏离对于一个给定的径向直径：$r$，基于多项式函数的一种简单而有效的径向失真模型：
$$
D\left( {r,k} \right) = {k_0} \cdot {r^2} + {k_1} \cdot {r^4} = k\left( {\begin{array}{*{20}{c}}
  {{r^2}} \\ 
  {{r^4}} 
\end{array}} \right)
$$
未知的系数 $k=(k_0,k_1)$，如图1所示，如果 $k_0=k_1=0$，那么畸变为0。

| ![image-20220506115431599](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220506115431599.png) |
| :----------------------------------------------------------: |
|        图1 径向畸变 $k=(k_0,k_1)=(-0.2286,0.190335)$         |

| <img src="https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220506115836533.png" alt="image-20220506115836533" style="zoom:67%;" /> |
| :----------------------------------------------------------: |
|                          图2 公式30                          |



## 2.6 投影过程总结

综上所述，以下步骤模拟了完整的投影过程（见图3）：

| ![image-20220506120050225](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220506120050225.png) |
| :----------------------------------------------------------: |
| 图3 投影链的摘要（从右到左）。<br />在图（c）中3D点 $\bf{X}$（相机坐标系）被投影到“理想 $f=1$” 平面，得到归一化坐标：$x=(x,y)^T$。<br />径向镜头畸变 （b） 映射到 $\bf{x}$ 到 $\tilde x = {\left( {\tilde x,\tilde y} \right)^T}$. <br />仿射变换通过内参转换（矩阵 $A$），最终产生观察到的传感器图像坐标：${\bf{u}}=(u,v)^T$ 在图（a）. |

1. **世界-相机转换**：给定一个点 $[X,Y,Z]^T$，以3D世界坐标系表示，它的3D相机坐标通过视角转换 $\bf{W}$ ：
   $$
   \left[ {\begin{array}{*{20}{c}}
     {X'} \\ 
     {Y'} \\ 
     {Z'} \\ 
     1 
   \end{array}} \right] = \underbrace {\left[ {\begin{array}{*{20}{c}}
     {{r_{11}}}&{{r_{12}}}&{{r_{13}}}&{{t_x}} \\ 
     {{r_{21}}}&{{r_{22}}}&{{r_{23}}}&{{t_y}} \\ 
     {{r_{31}}}&{{r_{32}}}&{{r_{33}}}&{{t_z}} \\ 
     0&0&0&1 
   \end{array}} \right]}_W\left[ {\begin{array}{*{20}{c}}
     X \\ 
     Y \\ 
     Z \\ 
     1 
   \end{array}} \right]
   $$
   在齐次坐标下，或者简单的：
   $$
   X'=W \cdot X
   $$

2. **投影到“归一化”（理想）图像平面上**：从3D点 $X'=(x',Y',Z')^T$ 透视投影到连续的规范化的2D坐标${\bf{x}}=(x,y)^T$ 在图像平面定义：
   $$
   \left[ {\begin{array}{*{20}{c}}
     x \\ 
     y 
   \end{array}} \right] = \frac{1}{{Z'}}\left[ {\begin{array}{*{20}{c}}
     {X'} \\ 
     {Y'} 
   \end{array}} \right] = {\hom ^{ - 1}}\left( {\left[ {\begin{array}{*{20}{c}}
     1&0&0&0 \\ 
     0&1&0&0 \\ 
     0&0&1&0 
   \end{array}} \right]\left[ {\begin{array}{*{20}{c}}
     {X'} \\ 
     {Y'} \\ 
     {Z'} \\ 
     1 
   \end{array}} \right]} \right) = \hat P\left( {W,X} \right)
   $$
   这相当于一个焦距为f=1的理想针孔投影

3. **镜头径向畸变**：归一化的二维投影坐标 ${\bf{x}}=(x,y)^T$ 服从非线性径向畸变，原始中心为：$x_c=[0,0]^T$，表达式如下：
   $$
   \tilde x = wrap\left( {x,k} \right)or\left[ {\begin{array}{*{20}{c}}
     {\tilde x} \\ 
     {\tilde y} 
   \end{array}} \right] = \left[ {\begin{array}{*{20}{c}}
     x \\ 
     y 
   \end{array}} \right]\left[ {1 + D\left( {r,k} \right)} \right]
   $$
   其中 $r=\sqrt{x^2 + y ^ 2}$，并且 $D(r,k)$ 如公式32定义。$\tilde x =(\tilde x, \tilde y) ^T$ 是镜头畸变后的2D坐标，仍然在归一化图像平面。

4. 仿射二维变换到传感器坐标：归一化的投影点最终被映射到缩放和倾斜的传感器坐标（公式13），通过仿射变换：
   $$
   u = A' \cdot \hom \left[ {\tilde x} \right]or\left[ {\begin{array}{*{20}{c}}
     u \\ 
     v 
   \end{array}} \right] = \left[ {\begin{array}{*{20}{c}}
     \alpha &\gamma &{{u_0}} \\ 
     0&\beta &{{v_0}} 
   \end{array}} \right]\left[ {\begin{array}{*{20}{c}}
     {\tilde x} \\ 
     {\tilde y} \\ 
     1 
   \end{array}} \right]
   $$
   其中：$\alpha, \beta,\gamma,u_0,v_0$ 是相机内参（如公式15和22所示）。

# 03 基于平面的自标定

Zhang[14,15]流行的相机校准方法使用了一些（至少两个）平面校准模式的视图，称为“模型”或“目标”，其布局和度量尺寸是精确已知的。校准程序大致如下：

1. 移动相机、或者标定板（全部）到不同视角，获得图像：$I_0,...,I_{M-1}$
2. 对于每个图像 $I_i(i=0,...,M-1)$，$N$ 个传感器点（特征点）${\dot u_{i,0}},...,{\dot u_{i,j}}$ 提取，假设：$1:1$ 这些点和模型平面
3. 从观测的点，相关的单应性矩阵：$H_0,...,H_{M-1}$（来自模型点和观察到的二维图像点的线性映射），通过每个视角 $i$ 进行估计
4. 从同一个单应性：$H_i$，五个内参参数 $(\alpha,\beta,\gamma,u_0,v_0)$ 估计得到封闭解（线性），忽略任意镜头畸变。$M\ge3$ 视角提供一个唯一解（达到一个不确定的比例因子）。如果假设传感器平面没有倾斜（即γ=0，这是一个合理的假设），那么N张=2张图像就足够了。更多的视角通常会导致更准确的结果(见Sec3.3).
5. 一旦知道了相机的内参，就会计算出每个相机视图i的外部三维参数$R_i，t_i$(见Sec3.4).
6. 径向畸真参数k0，k1通过线性最小二乘最小化估计（Sec3.5）
7. 最后，利用估计的参数值作为初始猜测，对所有M个视图进行非线性优化，对所有参数进行优化（Sec 3.6）

下面将更详细地解释这些步骤（有关表1中的完整描述和符号列表，请参见[14]）。

## 3.1 标定模型和观察视角

标定模型包含 $N$ 个参考点：$X_0,...,X_N-1$，这些参考点的3D位置是已知的。这些点假设在$XY$平面，即它们的$Z$轴为0.

我们假设该平面被 $M$ 个不同视角（即图片）拍摄，使用 $i=0,...,M-1$ 代表第 $i$ 个模型，从每副相机图片，我们可以得到观察的传感器点：
$$
{\dot u_{i,j}} \in {\mathbb{R}^2}
$$
其中：视角编号$i=0,...,M-1$，点编号 $j=0,...,N-1$。对于每个观察到的点：${\dot u_{i,j}}$ 必须对应模型点：$X_j$。于是模型点 $X_j$ 和图像点 ${\dot u_{i,j}}$ 必须以相同的顺序提供。满足这个条件是必要的，因为否则校准将提供无效的结果。

## 3.2 Step1: 对每个视角计算单应性矩阵

使用公式17，观察到的 ${\dot u_{i,j}} = \left( {{{\dot u}_{i,j}},{{\dot v}_{i,j}}} \right)$ 映射（单应性），相应的3D点：$X_J$ 能够表示为：
$$
s\left[ {\begin{array}{*{20}{c}}
  {{{\dot u}_{i,j}}} \\ 
  {{{\dot v}_{i,j}}} \\ 
  1 
\end{array}} \right] = A\left[ {{R_i}|{t_i}} \right]\left[ {\begin{array}{*{20}{c}}
  {{X_j}} \\ 
  {{Y_j}} \\ 
  {{Z_j}} \\ 
  1 
\end{array}} \right]
$$
或者写作：
$$
s \cdot {\dot {\underset{\raise0.3em\hbox{$\smash{\scriptscriptstyle-}$}}{u} }_{i,j}} = A\left[ {{R_i}|{t_i}} \right]{X{\underset{\raise0.3em\hbox{$\smash{\scriptscriptstyle-}$}}{u} }_j}
$$
$\underset{\raise0.3em\hbox{$\smash{\scriptscriptstyle-}$}}{u}, \underset{\raise0.3em\hbox{$\smash{\scriptscriptstyle-}$}}{X} $ 代表齐次坐标，其中：
$$
A = \left[ {\begin{array}{*{20}{c}}
  \alpha &\gamma &{{u_0}} \\ 
  0&\beta &{{v_0}} \\ 
  0&0&1 
\end{array}} \right]
$$
代表相机内参（对每个视角都相同），$s$ 是任意的，非零尺度系数。$R_i,t_i$ 是3D旋转矩阵和平移矩阵不同视角 $i$。

由于模型点：$X_j$ 被假设到 $XY$ 平面上世界坐标系（即 $Z_j=0$ 对所有 $j$），可以写作：
$$
s\left[ {\begin{array}{*{20}{c}}
  {{{\dot u}_{i,j}}} \\ 
  {{{\dot v}_{i,j}}} \\ 
  1 
\end{array}} \right] = A\left[ {{r_{i,0}},{r_{i,1}},{r_{i,2}},{t_i}} \right]\left[ {\begin{array}{*{20}{c}}
  {{X_j}} \\ 
  {{Y_j}} \\ 
  0 \\ 
  1 
\end{array}} \right] = A\left[ {{r_{i,0}},{r_{i,1}},{t_i}} \right]\left[ {\begin{array}{*{20}{c}}
  {{X_j}} \\ 
  {{Y_j}} \\ 
  1 
\end{array}} \right]
$$
其中：$r_{i,0},r_{i,1},r_{i,2}$ 代表 $R_i$ 的三列。注意：$Z_j=0$ 令第三行向量：$(r_{i,2})$ 多余，因此它在公式42中被消除。它等于2D单应性映射：
$$
s\left[ {\begin{array}{*{20}{c}}
  {{u_{i,j}}} \\ 
  {{v_{i,j}}} \\ 
  1 
\end{array}} \right] = {H_i}\left[ {\begin{array}{*{20}{c}}
  {{X_j}} \\ 
  {{Y_j}} \\ 
  1 
\end{array}} \right]
$$
其中：$s$ 是一个非零尺度系数，$H_i=\lambda A (R_i | t_i)$ 是一个 $3 \times 3$ 单应性矩阵（$\lambda$ 是任意系数，可以被忽略，当我们使用齐次坐标系）。矩阵 $H_i$ 由 $h_{i,0},h_{i,1},h_{i,2}$ 是：
$$
{H_i} = \left[ {{h_{i,0}},{h_{i,1}},{h_{i,2}}} \right] = \lambda A\left[ {{r_{i,0}},{r_{i,1}},{t_i}} \right]
$$
于是一个集合 $(u_{i,j},v_{i,j})^T$ 和 $(X_j,Y_j) ^T$ 计算一个单应性矩阵 $H_i$。

### 3.2.1 用直接线性变换的单应性估计(DLT)

在估计单应性映射的几种方法中，DLT是最简单（CH4），它也经常被用于张正友标定法原始部署。我们假设两个对应2D点序列：模型点 $X=(x_0,....,x_{N-1})$，点：$X_j=(X_j,Y_j)$ 和它相应传感器点（观察）传感器点：$\dot U = \left( {{{\dot u}_0},...,{{\dot u}_{N - 1}}} \right)$，每个 $u_j=(u_j,v_j)^T$，相关的单应性转换，如下（使用齐次坐标）：
$$
{u_j} = H \cdot {X_j}
$$
或者：
$$
\left[ {\begin{array}{*{20}{c}}
  {{u_j}} \\ 
  {{v_j}} \\ 
  {{w_j}} 
\end{array}} \right] = \left[ {\begin{array}{*{20}{c}}
  {{H_{0,0}}}&{{H_{0,1}}}&{{H_{0,2}}} \\ 
  {{H_{1,0}}}&{{H_{1,1}}}&{{H_{1,2}}} \\ 
  {{H_{2,0}}}&{{H_{2,1}}}&{{H_{2,2}}} 
\end{array}} \right]\left[ {\begin{array}{*{20}{c}}
  {{X_j}} \\ 
  {{Y_j}} \\ 
  {{Z_j}} 
\end{array}} \right]
$$
对于 $j=0,...,N-1$。不失去一般下，我们设置齐次坐标：$Z_j=1$（于是$X_j=X_j$）:

![image-20220507145953002](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507145953002.png)



在笛卡尔坐标系，一对非线性方程如下：

![image-20220507150202430](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507150202430.png)

可以重写排为：

![image-20220507150333125](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507150333125.png)

最终：

![image-20220507150353732](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507150353732.png)

这是一对齐次方程（因为右手边为零），在未知系数中是线性的：$H_{r,c}$ （虽然，由于混合项，仍然是非线性的w.r.t.坐标的坐标）。通过将未知同源矩阵H的9个元素收集到向量中：

![image-20220507150645010](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507150645010.png)

“公式”（52）和（53）可以用该形式编写：

![image-20220507150720899](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507150720899.png)

对于每个点：$\left( {{{\dot u}_j},{{\dot v}_j}} \right) \leftrightarrow \left( {{X_j},{Y_j}} \right)$ 。因此，N个点对，假设是由相同的单应性矩阵 $H$，产生 2N 个单应性方程组：

![image-20220507151204826](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507151204826.png)

或者用矩阵标记：

![image-20220507151523966](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507151523966.png)

其中：$M$ 是 $2 N \times 9$ 矩阵（所有的元素都是已知的常数），$\bf{h}$ 九个未知参数，$0$  是长度为2N的零向量。

---

**求解齐次线性方程组**：而Eqn（57）看起来很类似于一个普通的线性方程组的形式：$M x = b$，它不能以通常的方式来解决（没有额外的约束），因为它总是有 $h=0$ 作为一个平凡解。但是该方程可以通过奇异值分解获得（章节6.11，4.5.3）矩阵M，通过分解矩阵 $M (尺寸 2N \times 9)$ 到三个矩阵$U,S,V$的乘积：

![image-20220507154809676](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507154809676.png)

其中：$U$ 是一个 $2 N \times 2N$（该案例） 矩阵，$S$ 是 $2N \times 9$ 的对角矩阵，$V$ 是 $9 \times 9$ 的矩阵。M的具体分解是这样的：

> 定理：令 $A \in \mathbb{R} {^{m \times n}}$，则存在正交矩阵 $U \in {\mathbb{R}^{m \times m}}$ 和 $V \in {\mathbb{R}^{n \times n}}$，使得：
> $$
> A=U \sum V ^T
> $$
> 其中：
> $$
> \sum {}  = \left[ {\begin{array}{*{20}{c}}
>   {\sum\nolimits_1 {} }&0 \\ 
>   0&0 
> \end{array}} \right]
> $$
> 且 $\sum _1 =diag(\sigma_1, \sigma_2,...,\sigma_r)$，其对角元素按照顺序：
> $$
> {\sigma _1} \geqslant {\sigma _2} \geqslant  \cdots  \geqslant {\sigma _r} > 0,r = rank\left( A \right)
> $$



![image-20220507155339507](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507155339507.png)

对角矩阵 $S$，$s_0,...,s_8$ 称为矩阵 $M$ 的奇异值。每个奇异值 $s_i$ 是否具有关联的列向量 $u_i$ 在 $U$（称为 $M$ 的左奇异向量），行向量 $v_i^T$ 在 $V^T$（即 $V$ 的行向量 $v_i$），称为右奇异向量。因此Eqn（59）同样可以写成：

![image-20220507160203860](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507160203860.png)

其中 $U$ 是由 $(u_0,...,u_{2N-1})$ 组成的，$V$ 是由 $v_0^T,...,v_8^T$。

---

> 在总体最小二乘法中，矩阵方程为：
> $$
> \left( {A + E} \right)x = b + e
> $$
> 其中：$A_{m \times n}$ 数据矩阵，$E_{m \times n}$ 扰动矩阵，$x$ 是$n \times 1$，$b_{m \times 1}$ 数据向量，$e_{m \times 1}$ 扰动向量。
>
> 显然上式可以写作：
> $$
> \left( {\left[ { - b,A} \right] + \left[ { - e,E} \right]} \right)\left[ {\begin{array}{*{20}{c}}
>   1 \\ 
>   x 
> \end{array}} \right] = 0
> $$
> 记为：
> $$
> (B+D)z=0
> $$
> 其中：
>
> - $B=[-b,A]$：增广矩阵，$m \times (n + 1)$
> - $D=[-e,E]$：扰动矩阵，$m \times (n + 1)$
> - $z = \left[ {\begin{array}{*{20}{c}}
>     1 \\ 
>     x 
>   \end{array}} \right]$：$(n+1) \times 1$​
>
> > 情况1：$\sigma_n$ 明显大于 $\sigma_{n + 1}$，即最小的奇异值只有1个（奇异值是从大到小排列的）。
>
> 令 $m \times (n + 1)$ 的增广矩阵 $B$ 的奇异值为：
> $$
> B = U \sum V^H
> $$
> 并且其奇异值按照顺序：$\sigma_1 \ge \sigma_2 \ge ... \ge \sigma_{n+1}$ 排列，与这些奇异值对应的右奇异向量为：$v_1,v_2,...,v_{n+1}$。根据上面的分析，总体最下二乘解为：$z= v_{n+1}$。于是最小二乘解：
> $$
> {x_{TLS}} = \frac{1}{{v\left( {1,n + 1} \right)}}\left[ {\begin{array}{*{20}{c}}
>   {v\left( {2,n + 1} \right)} \\ 
>    \vdots  \\ 
>   {v\left( {n + 1,n + 1} \right)} 
> \end{array}} \right]
> $$
> 其中：$ v (i,n+1)$ 是 $V$ 的第 $n + 1 $ 列的第 $i$ 个元素。

公式57的解，即未知参数向量：$h$，最终通过右奇异值：$v_k$ 与最小的奇异值相关联的$s_k= min(s_0,...,s_8)$：

![image-20220507163136968](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220507163136968.png)

如果相应的单应性变换对应点集是精确的，则：$s_k = 0$。这通常是由4个对应的点对计算出单应性的情况，这是求解8个自由度所需的最小数。

---

如果涉及超过4个点对，则计算公式（57）是超定的（这是通常的情况）。在这里，$s_k$ 代表单应性的残差。当然，如果所有拟合点都是精确的话，$s_k$ 也应该是零。在超定系统的情况下，所得到的解使方程式（57）最小化，在最小二乘意义上：

![image-20220509213247055](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220509213247055.png)

在许多SVD部署中，$S$ 中的奇异值被排列为：非单调递增（即 $s_i \ge s_{i+1} $），因此，在我们的例子中，$s_8$ 是最小值，并且：$v_8$ （$V$ 的最后一列）是相应的解。例如，在JAVA的ACM部署中：

![image-20220509215028175](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220509215028175.png)

注意 $Eqn（62）$ 中的公式。最小化了一个与几何投影误差不直接相关的代数残差。这通常不会导致问题，因为剩余的错误在最终的整体优化步骤中被消除(见章节3.6)。然而，在这一阶段最小化同型态的投影误差(在Zhang的实现中没有这样做)可能有助于提高最终优化的收敛性。它需要非线性优化，对于上述解可以作为一个很好的初始猜测(见章节3.2.3)

### 3.2.2 归一化输入数据

为了提高数值计算的稳定性，我们最好对2D点：$\chi $ 和 $U$，在执行单应性估计之前。

归一化通过转换每个点集通过乘以一个 $3 \times 3$ 的归一化矩阵 $N$ （在齐次坐标下），因此转化的点集中心在原点，并且标准直径：

![image-20220510194902486](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510194902486.png)

可以看附录 Sec B 如何计算归一化矩阵 $N_X,N_U$。单应性估计（相似于公式45），执行获得归一化点：$X',U'$，通过计算矩阵 $H'$ （最小二乘场景）：

![image-20220510195400681](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510195400681.png)

对于 $j=0,...,N-1$，通过替代 ${u_j}$，可以得到：

![image-20220510195721010](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510195721010.png)

通过消去 $x_j$ 两边，可以得到：

![image-20220510195939201](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510195939201.png)

于是归一化后的单应性矩阵：

![image-20220510201405992](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510201405992.png)



### 3.2.3 对单应性矩阵进行非线性优化

如上所述，用DLT方法获得的同源性估计通常不会最小化传感器图像中的投影误差。在这一步中，通过最小化投影误差，对单个视图的估计单应性矩阵 $H$ 进行了数值细化。最小化投影误差需要非线性优化，这通常是通过迭代方案来实现的，如附录的Sec E中描述的黎文堡-马卡特(LM)方法。







# Section B 归一化2D点

一般来说，规范化一个给定的二维点集 $X=(X_0,...,X_{N-1} )$ 是通过移动和缩放所有点来完成的，这样变换集的质心与原点对齐，其直径具有预定义的大小：

![image-20220510204342053](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510204342053.png)

通过 $x_j'=N_x x_j$，使用单应性坐标，即：

![image-20220510204443237](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510204443237.png)

$N_x$ 是针对点集 $X$ 的一个专门的归一化矩阵，通常结构：

![image-20220510204703369](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510204703369.png)

其中：$\bar x = {\left( {\bar x,\bar y} \right)^T} = \frac{1}{N}\sum\nolimits_{j = 0}^{N - 1} {{x_j}} $ 是原始点集的中心。通过计算尺度系数，在文献中可以找到几种方法，其中两种方法如下所示：

**方法1：**这种方法对点集（沿两个轴均匀地）进行缩放，使点到原点的平均距离等于 $\sqrt{2}$，不应用旋转。在这种情况下：

![image-20220510210054036](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510210054036.png)

**方法2：**这里的缩放在x和y方向上不均匀地应用，这样沿两个轴的方差都被归一化，即：

![image-20220510210125028](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510210125028.png)

![image-20220510210257334](https://flyman-cjb.oss-cn-hangzhou.aliyuncs.com/image-20220510210257334.png)

其中：$\sigma ^ 2$ 代表方差。

请注意，对于不与坐标轴对齐的具有高偏心度的点集，上述方法没有一个是最优的。可能有更好的方法来归一化这些点集，其中还包括旋转，以使主导方向与坐标轴对齐，例如，通过主成分分析(PCA)。

