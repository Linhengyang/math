# 综述: 线性方程组/线性空间/线性映射
贯穿全章节的计算手段主要是:
* 初等行变换(elementary row operation)
* 矩阵向量乘积(matrix-vector product)
* 矩阵乘法(matrix-matrix product)
  
这块内容从**线性方程组(矩阵方程)**讲起，判断是否**有解/解集结构**，到用**线性组合(线性相关/无关)**的方式理解**矩阵方程**，然后两个分支
* 分支1：从**线性映射**重新理解**矩阵**，包括**矩阵乘法(复合映射)**和**可逆矩阵**。这里重新用**初等矩阵左乘**重新理解了初等行变换，通过解答**矩阵方程**解集的结构理解**矩阵**的**列空间**和**零空间**，最后引出**LU分解**
* 分支2：从**线性子空间**理解矩阵的**列空间**/**零空间**并描述，然后引出**基**/**维度**/**坐标**的概念，并推出矩阵的**秩**。之后，引入完整的**线性空间**，并给出**线性空间**的**基**/**维度**/**坐标**，最后说明了同一个向量在不同的**基**下如何变换**坐标**

## 分支1综述
分支1脉络中, 主线是矩阵方程 $A\vec{x}=\vec{b}$ 的求解，其中 $A \in \mathbb{R^{m,n}}$。
记系数矩阵 $A$ 的echelon form(阶梯型)是 $U$。  
求解过程中, 初等行变换（左乘初等矩阵）的操作，将 $A$ 转化 $U$，甚至 $A$ 的reduced echelon form(简约阶梯型)形式。  
初等行变换非常重要，在这个操作变换下，矩阵有两个不变量：
1. 矩阵的秩 $rankA$
2. 矩阵的列向量组的序

用初等矩阵 $E_i$ 左乘来理解初等行变换, 可以得到：
$$E_p\ \cdots\ E_1\ A=U$$
即 $A=({E_p\ \cdots\ E_1})^{-1}U$。记 $P={E_p\ \cdots\ E_1}$，则 $A=P^{-1}\ U$。

对于线性方程组 $$A\vec{x}=\vec{b}$$ 来说, 它等价于 $$P\ A\vec{x}=P\ \vec{b}$$, 即 $$U\vec{x}=P\ \vec{b}$$  

---
第一步：判读是否有解  
对比 $U$ 的非零行数量 $k$ 和 $P\ \vec{b}$ 的非零行数量
* 当 $k$ 小于 $P\ \vec{b}$ 的非零行数量时，无解
* 当 $k$ 大于等于 $P\ \vec{b}$ 的非零行数量时，有解

第二步：解集结构  
对比 $U$ 的非零行数量 $k$ 和 $n$
* 当 $k=n$ 时，说明矩阵 $A$ 列向量组 $\mathcal{A}$ 线性无关，它作为 $\mathbb{R^n}$ 的一个基（basis），对于任何 $\vec{b} \in \mathbb{R^{n}}$，在此 basis 下有唯一表示，即 $A\vec{x}=\vec{b}$ 只有唯一解，即 $\vec{b}$ 在 $A$ 的列向量组 $\mathcal{A}$ 下的坐标 $\displaystyle \left[\vec{b}\right]_\mathcal{A}\in\mathbb{R^n}$。从 $U\ \vec{x}=P\ \vec{b}$ 也可以看出，echelon form $U$ 的每一列都是主元所在列，所以 $\vec{x}$ 只有唯一解。
* 当 $k \lt n$ 时，因为矩阵的行秩等于列秩，所以 $U$ 的主元所在列数量也是 $k$。echelon form $U$ 的主元所在列组成了 $A$ 的列向量组 $\mathcal{A}$ 的一个basis，记为 $\mathcal{B}$。而 $A$ 的列向量组 $\mathcal{A}$ 在自由元所在列可以构成一个basis $\mathcal{B}$，即 $span\displaystyle \left(\mathcal{A}\right) \equiv span \left(\mathcal{B}\right)$，是 $\mathbb{R^m}$ 的一个 $dim = k$ 的子空间。有解的意思即指 $\vec{b} \in span\left(\mathcal{A}\right)$。

### 矩阵方程解集
#### 齐次  
首先考虑齐次线性方程组 $$A\vec{x}=\vec{0}$$
echelon form下是 $$U\vec{x}=\vec{0}$$
想象下 $U\vec{x}=\vec{0}$ 的形式，因为 $U$ 是echelon form, 再考虑它的reduced echelon form，此时每行/每列都最多只有一个主元，所以每一行可以写成：
```math
x_{pivot\_q} + \sum_{j=1}^{n-k} coef_{j}\ x_{free\_j} = 0, q \in (1,\cdots,k)
```
移项，得：
```math
x_{pivot\_q} = \sum_{j=1}^{n-k} coef_{j}\ x_{free\_j}, q \in (1,\cdots,k)
```
即：
```math
\vec{x}_{pivot} = C\vec{x}_{free}, C \in \mathbb{R^{k, n-k}}
```
那么，此时考虑解向量
```math
\vec{x} = 
\begin{bmatrix}
\vec{x}_{pivot}\\\vec{x}_{free}
\end{bmatrix} = 
\begin{bmatrix}
C\ \vec{x}_{pivot}\\I\ \vec{x}_{free}
\end{bmatrix} = 
\begin{bmatrix}
C\\I
\end{bmatrix}\ \vec{x}_{free}
```
记
```math
\begin{bmatrix}
C\\I
\end{bmatrix} = Q, \ \ Q \in \mathbb{R^{n, n-k}}
```
考虑独热部分，所以 $rank(Q) = n-k$，同时 $\vec{x}_{free}$ 取遍 $\mathbb{R^{n-k}}$，所以此时解集 $set$ (又称 $A$ 的零空间 $Nul(A)$ )是由 $Q$ 的列向量组 $\mathcal{Q}$ 张成的线性空间，即 $\mathbb{R}^{n}$ 的一个 $dim = n-k$的subspace。

#### 非齐次
特解+系数矩阵的零空间，即为非齐次线性方程组的解集。

---
### LU分解  
当矩阵 $A$ 转化为 echelon form $U$ 的初等行变换操作中，没有「两行互换」时，初等矩阵 $E_p,\cdots,E_1$ 都是下三角单位矩阵（主对角线元素都是1），从而 $P = E_p\ \cdots\ E_1$ 以及 $P^{-1}$都是下三角单位矩阵，将其记作 $L$，即得到
```math
A = LU
```
其中 $L$ 是下三角单位矩阵，而 $U$ 是 $A$ 的echelon form。

## 分支2综述
分支2脉络中，主线是**子空间**的**基**。**基**的大小就是子空间的**维度**，向量 $\vec{x} \in subspace\ H \subset \mathbb{R^m}$ 被基 $\mathcal{B}$ 唯一表示时的权重，就是向量 $\vec{x}$ 在这个基下的**坐标**$\left[\vec{x}\right]_\mathcal{B} \in \mathbb{R^p}, if\ dimH = p$。  
  
考虑映射 $T:H \rightarrow \mathbb{R^p}, dimH=p, \mathcal{B}\ is\ a\ basis\ of\ H$，那么有
```math
T(\vec{x}) = \left[\vec{x}\right]_{\mathcal{B}}
```
即将一个 $H$ 中的向量映射到 基basis下的坐标向量。易证这是一个保持**加法**和**数乘**运算的**同构映射**。所以 $H$ “像是” $\mathbb{R^p}$，尽管 $H$ 中的元素的分量数目可能是超过 $p$的。  
  
矩阵 $A$ 的**列空间**和**零空间**的**维度**和**基**，都已经知道怎么计算。  
  
一个子空间 $H$ 的任一组**基**都是能互相表示(等价)的，而且能够张成 $H$ 的 $dimH$ 个向量，就是 $H$ 的一组基。
  
### 线性空间
完整的**线性空间** $V$ 的定义：
* 一个非空集合 $V$，两个运算（加法和数乘）
* 十条法则  
  
完整的**子空间** $H$ 的定义:
* 包含 $0$ 向量
* 对 $V$ 的加法和数乘封闭  

由此也引入了向量组张成的子空间 $span(\mathcal{B})$。在这里把矩阵 $A$ 的零空间Null Space $Nul(A)$ 和列空间 Col Space $Col(A)$重新理解一遍（其实没什么变化）。  

但是如果从映射的视角出发：  
考虑映射 $T:\mathbb{R^n} \rightarrow \mathbb{R^m},\ with\ matrix\ A \in \mathbb{R^{m,n}}$，那么
1. 它的零空间 $Nul(A)$是 $domain\ \mathbb{R^n}$ 的一个子空间，代表这个子空间中的向量在映射 $T$ 下都被映射到了 $codomain\ \mathbb{R^m}$ 的 $\vec{0}$。
2. 它的列空间 $Col(A)$是 $codomain\ \mathbb{R^m}$ 的一个子空间，可以记作 $span(\left[\vec{a}_1, \vec{a}_2, \cdots, \vec{a}_n\right])$，代表映射T的值域 $range$ 。  
  
通过探讨映射 $T$ 是否是单射，可以得出：
1. $rank\ A = n$
2. $A$ 的列向量组线性无关
3. $A\vec{x} = \vec{0}$ 只有0解  
  
这三个相互等价的条件时，映射 $T$ 是单射，即 $Nul(A) = set(\vec{0})$ 。
  
通过探讨映射 $T$ 是否是满射，可以得出：
1. $rank\ A = m$
2. $A$ 的列向量组张成 $\mathbb{R^m}$（记住 $dim\ \mathbb{R^m} = m$）
  
这两个相互等价的条件时，映射 $T$ 是满射，即 $Col(A) = \mathbb{R^m}$ 。
  
矩阵 $A \in \mathbb{R^{m,n}}$ 的**行空间Row Space**：行向量组张成的空间，是 $\mathbb{R^n}$ 的一个子空间。矩阵 $A$ 的echelon form $U$ 的非零行是它的一个**基**（注意这里是 $U$ 的非零行，跟列空间不同，列空间的基是从 $A$ 中选 $U$ 的主元所在列）。  

#### 线性映射
如果将线性映射的概念从**数值向量空间**之间的映射，推广到**线性空间**之间的映射 $T$，那么相应地，
* 把 $T$ 的“零空间” $Nul(T)$ 定义为 $T$ 的**kernel**，是**domain**的子空间：映射 $T$ 把kernel中的元素 $\vec{v}$ 映射到codomain的 $\vec{0}$ 元素
* 把 $T$ 的“列空间” $Col(T)$ 定义为 $T$ 的**range**，是**codomain**的子空间：对于range中的元素 $\vec{u}$，能从domain中找到对应元素 $\vec{v}$，使得 $T(\vec{v})=\vec{u}$
  
#### 基/坐标/维度/秩
对于无限维线性空间，线性代数中不作过多研究，仅需知道它的维度定义为无限。 
  
对于**有限维线性空间**，即可以由有限个向量张成的**线性空间**下，**线性组合**、**线性相关/无关**、**张成span**、**极大线性无关组**的定义都如出一辙，从而**基**和**坐标**的定义和性质也如出一辙：  
  
考虑线性空间 $V$ 的一个基
```math
\mathcal{B} = \left[\vec{b}_1,\cdots,\vec{b}_n\right]
```
对任意 $\vec{x} \in V$，相对这个基都有唯一表示
```math
\vec{x} = \sum_{i=1}^{n}c_i\ \vec{b}_i
```
把权重按序表示成列向量 $\left[c_1,\cdots,c_n\right]^T$，即是 $\vec{x}$ 相对基 $base\ \mathcal{B}$ 的坐标，写作
$\left[\vec{x}\right]_\mathcal{B} \in \mathbb{R^n}$。  
  
坐标映射coordinates-mapping
```math
T:V \rightarrow \mathbb{R^n},\ T(\vec{x})=\left[\vec{x}\right]_\mathcal{B}
```
是一个双射线性映射，也是一个关于加法和数乘的**同构**映射，线性空间 $V$ 和 $\mathbb{R^n}$ 同构。同构意味着坐标映射保持加法和数乘运算，即意味着保持线性相关/无关的关系。

---

当 $V = \mathbb{R^n}$ 时，根据「坐标是一个列向量，**基**以**坐标**为权重作**线性组合**得到原向量」的原则，可得此时这是一个matrix-vector product运算：
```math
\vec{x} = P_{\mathcal{B}}\ \left[\vec{x}\right]_{\mathcal{B}},\ P_{\mathcal{B}} = \left[\vec{b}_1,\cdots,\vec{b}_n\right]
```
而向量 $\vec{x}$ 自身可以看作其在**标准基**下的坐标向量，也就是说，基 $\mathcal{B}$ 作为列向量组组成的矩阵 $P_{\mathcal{B}}$ 成为了两个基之间的**坐标转移矩阵(change-of-coordinates matrix)**。在这里，矩阵 $P_{\mathcal{B}}$ 是从 $\mathcal{B}$ 到标准基的坐标转移矩阵，而  $P_{\mathcal{B}}^{-1}$ 是从标准基到 $\mathcal{B}$ 的坐标转移矩阵，因为有：
```math
\left[\vec{x}\right]_{\mathcal{B}} = P_{\mathcal{B}}^{-1}\ \vec{x}
```
   

#### 坐标映射(坐标变换)
前面已经讲述过，当**线性空间** $V, dimV = n$ 就是 $\mathbb{R^n}$ 时，对于它的一个**基** $\mathcal{B}$ ，有
```math
\vec{x} = P_{\mathcal{B}}\ \left[\vec{x}\right]_{\mathcal{B}}
```
这里 $P_{\mathcal{B}}$ 就是 $\mathcal{B}$ 作为列向量组组成的矩阵。然而对于一般的线性空间 $V$，很可能没办法把基 $\mathcal{B}$ 写成矩阵的形式。  
对于一般的线性空间 $V, dimV = n$，考虑它的两个基 $\mathcal{B}$ 和 $\mathcal{C}$，那么对于 $V$ 中的任何一个 $\vec{x} \in V$，它在两个基下有不同的坐标，分别是
```math
\left[\vec{x}\right]_\mathcal{B},\ \left[\vec{x}\right]_\mathcal{C}\ \in\mathbb{R^n}
```
从 $\mathcal{B}$ 下的坐标向量到 $\mathcal{C}$ 下的坐标向量，存在**唯一**的**坐标转移映射** $T,\ matrix\ P_\mathcal{C \leftarrow B} \in \mathbb{R^n}$，使得下式成立：
```math
\left[\vec{x}\right]_\mathcal{C} = P_\mathcal{C \leftarrow B}\ \left[\vec{x}\right]_\mathcal{B}
```
坐标转移矩阵 $P_\mathcal{C \leftarrow B}$ 的列向量，是基 $\mathcal{B}$ 里的向量分别按序在另一个基 $\mathcal{C}$ 下的坐标向量，即
```math
Consider\ \mathcal{B}\ as\ set(\vec{b}_1,\cdots,\vec{b}_n),\ P_\mathcal{C \leftarrow B}=
\begin{bmatrix}
\left[\vec{b}_1\right]_\mathcal{C},\ \cdots\ \left[\vec{b}_n\right]_\mathcal{C}
\end{bmatrix} \in \mathbb{R^n}
```
坐标转移映射都是可逆的，即**坐标转移矩阵**都是**可逆**的，有 $P_\mathcal{B \leftarrow C} = P_\mathcal{C \leftarrow B}^{-1}$ 。
  
---

#### 矩阵(方阵)相似
**坐标**概念的建立，有一个很大的意义在于，原**线性空间** $V,\ dimV = n$ 中的向量可能不方便用数字表示，无法参与进一步的分析（比如研究 $V$ 到自身的线性映射）。但是如果去研究与之**同构**的 $\mathbb{R^n}$，即**坐标向量**和**坐标空间**，就方便了。特别地，如果研究从 $V$ 到 $V$ 的线性映射，就可以转而研究从坐标空间 $\mathbb{R^n}$ 到 $\mathbb{R^n}$ 的线性映射，即一个形状为(n,n)的方阵。  
  
甚至有时候，一个**基** $\mathcal{B}$ 得到的坐标空间 $\mathbb{R^n}\_{\mathcal{B}}$ 可能还「不够好」，我们会换一个**基** $\mathcal{C}$，即换一个坐标空间 $\mathbb{R^n}\_{\mathcal{C}}$。  
  
更具体地说，一个线性映射 $T:\mathbb{R^n}\_{\mathcal{B}}\rightarrow\mathbb{R^n}\_{\mathcal{B}},\ matrix\ T\in \mathbb{R^{n,n}}$，与另一个线性映射 $S:\mathbb{R^n}\_{\mathcal{C}}\rightarrow\mathbb{R^n}\_{\mathcal{C}},\ matrix\ S\in \mathbb{R^{n,n}}$，是 $V$ 的两个不同基下的**坐标空间**里各自的线性映射，但其实它们都是 $V$ 到 $V$ 上的一个线性映射在不同**同构**的坐标空间里的“翻版”，此谓两个映射**相似**，或者说两个**矩阵相似**。  
  
用数学来表示，即：  
```math
V(dimV = n,\ base\ \mathcal{B}\ \&\ \mathcal{C}),\ V\cong\mathbb{R^n}_{\mathcal{B}},\ V\cong\mathbb{R^n}_{\mathcal{C}}
```
现在，有一个 $V$ 到 $V$ 的线性映射
```math
F:V\rightarrow V,\ \ F(\vec{x})=\vec{y},\ \vec{x}\in V,\ \vec{y} \in V
```
考虑 $F$ 在两个坐标空间的“翻版”，即如下两个线性映射：
```math
\begin{cases}
F_\mathcal{B}:\mathbb{R^n}_{\mathcal{B}}\rightarrow\mathbb{R^n}_{\mathcal{B}},\ matrix\ T\in \mathbb{R^{n,n}},\ T\left[\vec{x}\right]_\mathcal{B} = \left[\vec{y}\right]_\mathcal{B} \\
F_\mathcal{C}:\mathbb{R^n}_{\mathcal{C}}\rightarrow\mathbb{R^n}_{\mathcal{C}},\ matrix\ S\in \mathbb{R^{n,n}},\ S\left[\vec{x}\right]_\mathcal{C} = \left[\vec{y}\right]_\mathcal{C} \\
\end{cases}
```
因为有坐标转移映射
```math
\begin{cases}
\left[\vec{x}\right]_\mathcal{C}=P_\mathcal{C \leftarrow B}\left[\vec{x}\right]_\mathcal{B},\ P_\mathcal{C \leftarrow B}\in\mathbb{R^{n,n}} \\
\left[\vec{y}\right]_\mathcal{B}=P_\mathcal{B \leftarrow C}\left[\vec{y}\right]_\mathcal{C},\ P_\mathcal{B \leftarrow C}\in\mathbb{R^{n,n}} \\
\end{cases}
```
代入，得到：
```math
T\left[\vec{x}\right]_\mathcal{B}=\left[\vec{y}\right]_\mathcal{B}=P_\mathcal{B \leftarrow C}\left[\vec{y}\right]_\mathcal{C}=P_\mathcal{B \leftarrow C}S\left[\vec{x}\right]_\mathcal{C}=P_\mathcal{B \leftarrow C}SP_\mathcal{C \leftarrow B}\left[\vec{x}\right]_\mathcal{B}
```
即
```math
T\left[\vec{x}\right]_\mathcal{B} \equiv P_\mathcal{B \leftarrow C}SP_\mathcal{C \leftarrow B}\left[\vec{x}\right]_\mathcal{B}
```
记 $P_\mathcal{B \leftarrow C} = P$，得到：
```math
T = P\ S\ P^{-1}
```
去掉一切推导归纳来看，对于矩阵 $A,B \in \mathbb{R^{n,n}}$，如果存在**可逆矩阵** $P$，使得 $A = P\ B\ P^{-1}$ 成立，则称矩阵 $A$ 和 $B$ **相似**。