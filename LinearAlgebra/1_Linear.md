# 综述: 线性方程组/线性空间/线性映射
贯穿全章节的计算手段主要是:
* 初等行变换(elementary row operation)
* 矩阵向量乘积(matrix-vector product)
* 矩阵乘法(matrix-matrix product)

正文：
* 铺垫：从**线性方程组(矩阵方程)**讲起，重点是**线性组合**，用它理解**生成子空间(span)**、**矩阵方程**、**矩阵乘法**，然后理解**线性无关**/**线性相关**、**极大线性无关组**的概念
* 主干1：从**线性映射**重新理解**矩阵**，包括**矩阵乘法(复合映射)**和**可逆矩阵**。这里用**初等矩阵左乘**重新理解了初等行变换，然后通过解答**矩阵方程是否有解和解集的结构**理解**矩阵**的**列空间**和**零空间**，最后引出**LU分解**
* 主干2：从**线性子空间**理解矩阵的**列空间**/**零空间**并描述，然后引出**基**/**维度**/**坐标**的概念，并推出矩阵的**秩**。之后，引入完整的**线性空间**，并给出**线性空间**的**基**/**维度**/**坐标**，给出同一个向量在不同的**基**下如何**坐标变化/转移**

后述：
* 理解**线性空间**、**坐标**和**坐标转移**、用**线性映射**理解**矩阵**之后，引出**相似**的概念，并解释为什么要研究矩阵的特征值和特征向量

## 铺垫综述
从线性方程组（矩阵方程）开始铺垫。用矩阵乘以向量（matrix-vector product）的方式描述矩阵方程：
```math
\begin{bmatrix}
a_{11}\ \ a_{12}\ \ \cdots\ \ a_{1n}\\
a_{21}\ \ a_{22}\ \ \cdots\ \ a_{2n}\\
\vdots\\
a_{m1}\ \ a_{m2}\ \ \cdots\ \ a_{mn}\\
\end{bmatrix}
\begin{bmatrix}
x_{1}\\
x_{2}\\
\vdots\\
x_{n}\\
\end{bmatrix}=
\begin{bmatrix}
b_{1}\\
b_{2}\\
\vdots\\
b_{n}\\
\end{bmatrix}
```
即
```math
A\vec{x}=\vec{b}
```
### 线性组合
用**线性组合**的观点去看线性方程组：即 $\vec{b}$ 是 $A$ 的列向量组的线性组合。  
有了**线性组合**的观点，就可以理解**生成/张成子空间**的概念，即向量组合的权重遍历数域。用**生成子空间**的观点去看矩阵方程：即 $\vec{b} \in span(\vec{a_1},\ \vec{a_2}\ \cdots \vec{a_n})$。  
甚至可以用**线性组合**的观点重新看**matrix-vector product**：即矩阵列向量组的线性组合，权重是向量的元素。  
  
**线性无关**和**线性相关**，即：向量组在线性组合成零向量时是否只有权重全0这一种方式。如果是，则称该向量组**线性无关**，否则就是**线性相关**。  
这里有一系列的命题和判断，并引出**极大线性无关组**的概念。即：
* 线性相关的向量组，其中至少有一个向量可以由其他向量线性表出
* 考虑到单个非零向量是线性无关的，那么该向量组一定存在子集，使得子集本身线性无关，但只要多放进一个其他向量就线性相关。这类子集被称为该向量组的**极大线性无关组**
* 向量组的所有**极大线性无关组**都是等价的，即极大线性无关组们可以互相线性表出
  
### 矩阵乘法
用线性组合的观点理解矩阵乘法，即对于：
```math
A \in \mathbb{R}^{m,n},\ \ \ 
B \in \mathbb{R}^{n,p},\ B=\left[\vec{b}_{1},\ \vec{b}_{2},\ \cdots, \vec{b}_{n}\right]\\
```
有
```math
AB = \left[A\vec{b}_{1},\ A\vec{b}_{2},\ \cdots, A\vec{b}_{n}\right]
```
即**左矩阵分别与右矩阵的列向量作matrix-vector product**。即：
```math
Col_j(A\ B)=A\ Col_j(B)
```
同样地，也可以有：
```math
Row_i(A\ B)=Row_i(A)\ B
```

### 线性映射
线性映射是一种保持了**加法**和**数乘**运算的映射，即：  
```math
For\ \vec{x} \in \mathbb{R}^{n}, \mathbb{R}^{n}\ is\ domain,\ \mathbb{R}^{m}\ is\ codomain
```
映射
```math
Mapping\ \ \ T:\mathbb{R^n}\rightarrow\mathbb{R^m},\\
\vec{x}\ \rightarrow\ T(\vec{x}),\\
we\ have\ image\ of\ \vec{x}\ as\ T(\vec{x}),\ \ set\{T(\vec{x}) \vert\ \forall \ \vec{x}\in\mathbb{R}^{n}\}\ as\ range\ of\ T
```
当映射 T 保持加法和数乘运算时，即：  
* $\forall\ \vec{u},\ \vec{v}\ \in\ domain,\ T(\vec{u}+\vec{v})=T(\vec{u})+T(\vec{v})$  
* $\forall\ c\ \in\ \mathbb{R},\ \vec{u}\ \in\ domain,\ T(c\cdot\vec{u})=c\cdot T(\vec{u})$  
  
那么此时 T 是一个**线性映射**。  
  
可以证明，任何一个线性映射 $T:\mathbb{R^n} \rightarrow \mathbb{R^m}$ 一定和一个矩阵 $A\in\mathbb{R}^{m,n}$ 一一对应，而且**映射的像就是矩阵乘以向量元素的结果**。也就是说，矩阵可以理解为线性映射。  

## 主干1综述
用**线性映射**重新理解**矩阵**之后，可以根据映射的性质，很快给出系数矩阵相应的结论，具体来说，考虑线性映射
```math
T:\mathbb{R^n}\rightarrow\mathbb{R^m}\ \ \ with\ \ \ matrix\ \ \ A\ \in\ \mathbb{R}^{m,n},\ A = \left[\vec{a}_1,\ \vec{a}_2,\cdots\vec{a}_n\right]
```
（当前暂时用**向量组的极大线性无关组的向量个数**作为**秩**的定义），那么
* T是满射，等价于 $A$ 的列向量组**张成（span）** $\mathbb{R}^{m}$，等价于 $A$ 的列向量组的秩等于 $m$。
* T是单射，等价于 $A$ 的列向量组**线性无关**，等价于 $A$ 的列向量组的秩等于 $n$。
* T是双射（有可逆映射），等价于 $A$ 的列向量组的秩等于 $n$ 等于 $m$，等价于 $A$ 是方阵，且列向量组线性无关。
  
由此：  
我们可以用**复合映射**来重新理解**矩阵乘法**，这为我们求解矩阵方程 $A\vec{x}=\vec{b}$ 提供了一条思路：如果可以找到一个双射 $T:\mathbb{R^n}\rightarrow\mathbb{R^m}\ \ \ with\ \ \ matrix\ \ \ P\ \in\ \mathbb{R}^{m,n}$ ，此时矩阵方程 $A\vec{x}=\vec{b}$ 和 矩阵方程 $P\ A\vec{x}=P\ \vec{b}$ 等价。只要 $P\ A$ 的结构足够好（比如说是阶梯形矩阵），那么求解 $\vec{x}$ 就是很方便的事情。  
综上，开始需要研究**逆映射**，即**可逆矩阵**和**逆矩阵**。具体地说，为了求解矩阵方程 $A\vec{x}=\vec{b}$，需要找到 $P\ A$ 为阶梯型矩阵的可逆矩阵 $P$。从最简单的可逆矩阵--初等矩阵入手。  
  
### 可逆矩阵/逆映射
对于矩阵 $A \in \mathbb{R^n}$，如果存在矩阵 $B \in \mathbb{R^n}$，使得 $B A = I$，那么称矩阵 $B$ 是矩阵 $A$ 的**逆矩阵**。  
如果矩阵 $A$ 的逆矩阵存在，那么它是唯一的，记为 $A^{-1}$，有 $A^{-1} A = A A^{-1} = I$。  
逆矩阵代表逆映射，意味着原矩阵代表的线性映射可逆，意味着它是双射，即：矩阵 $A$ 可逆 等价于 $rank(A) = n$。  
  
### 最简单的可逆矩阵：初等矩阵
主干1脉络中, 主线是矩阵方程 $A\vec{x}=\vec{b}$ 的通用求解，其中 $A \in \mathbb{R^{m,n}}$。
记系数矩阵 $A$ 的echelon form(阶梯型)是 $U$。  
求解过程中, 初等行变换的操作，将 $A$ 转化 $U$，甚至 $A$ 的reduced echelon form(简约阶梯型)形式。**初等行变换等价于对应的初等矩阵左乘**。  
  
初等行变换非常重要，初等行变换中、后的所有矩阵都和原矩阵**行等价**。在这个操作变换中，矩阵有两个不变量：
1. 矩阵的秩 $rankA$
2. 矩阵的列向量组的序
  
初等行变换有三型，都有其对应的**初等矩阵**(即由单位矩阵作一次相应的初等行变换的结果矩阵)
* replacement:即 $row(i)$ 换成 $row(i) + k * row(j),\ \ j \neq i$
* interchange:即 $row(i)$ 和 $row(j),\ \ j \neq i$ 两行互换
* scale by **r**:即 $row(i)$ 换成 $r * row(i),\ \ r \neq 0$
  
用初等矩阵 $E_i$ 左乘来理解初等行变换, 可以得到：
$$E_p\ \cdots\ E_1\ A=U$$
即 $A=({E_p\ \cdots\ E_1})^{-1}U$。记 $P={E_p\ \cdots\ E_1}$，则 $A=P^{-1}\ U$。

对于线性方程组 $$A\vec{x}=\vec{b}$$ 来说, 它等价于 $$P\ A\vec{x}=P\ \vec{b}$$
即 $$U\vec{x}=P\ \vec{b}$$  

---
### 矩阵方程求解并描述解集
第一步：判读是否有解  
对比 $U$ 的非零行数量 $k$ 和 $P\ \vec{b}$ 的非零行数量
* 当 $k$ 小于 $P\ \vec{b}$ 的非零行数量时，无解
* 当 $k$ 大于等于 $P\ \vec{b}$ 的非零行数量时，有解

第二步：解集结构  
对比 $U$ 的非零行数量 $k$ 和 $n$
* 当 $k \lt n$ 时，因为非零行都在零行上面，所以每一行非零行都有主元，同时有主元的行都是非零行，所以主元数量就是 $k$，所以 $U$ 的主元所在列的数量也是 $k$，echelon form $U$ 的主元所在列指明了 $U$ 的极大线性无关组的位置。考虑到初等行变换不改变列与列之间的线形关系，以及 $U$ 和 $A$ 的列向量组之间的顺序相同，于是这些位置在 $A$ 的列也组成了 $A$ 的列向量组 $\mathcal{A}$ 的一个极大线性无关组，定义为 $A$ 的列向量组 $\mathcal{A}$ 的一个**基**，记为 $\mathcal{B}$。这里 $size\ of\ \mathcal{B} = k$。  
---- 这里，定义向量组的**基**并给出了计算方法，即对矩阵作初等行变换至阶梯型，此时主元所在列的位置，就是原矩阵的列向量组的一个基的位置。  
---- 这里，证明了一个不是很重要的结论：矩阵的列秩等于行秩，良定义为**矩阵的秩**。  
---- 这里，引入概念**矩阵A的列空间Col(A)**，即矩阵 $A$ 的列向量组 $\mathcal{A}$ 生成的空间 $span\displaystyle \left(\mathcal{A}\right)$，是 $\mathbb{R}^{m}$ 的子空间。可以看出 $span\displaystyle \left(\mathcal{A}\right) \equiv span \left(\mathcal{B}\right)$。定义向量组的基 $\mathcal{B}$ 为它们所生成的**子空间的基**，并定义基的向量个数为**生成子空间的维度**，记为 $dim\ span\displaystyle \left(\mathcal{A}\right) = size\ of\ \mathcal{B}$。  
  
    有解，即说明 $\vec{b} \in span\left(\mathcal{A}\right) = span\left(\mathcal{B}\right)$。由于 $\mathcal{B}$ 是一个基，所以 $\vec{b}$ 被 $\mathcal{B}$表示方法是唯一的。而由于 $A$ 的列向量组在自由元所在列，也可以被 $\mathcal{B}$ 线性表示，于是向量 $\vec{b}$ 被 $\mathcal{A}$表示就不是唯一的（有点绕，具体分析看下文解集结构分析就行）。
* 当 $k=n$ 时，说明矩阵 $A$ 列向量组 $\mathcal{A}$ 线性无关，它作为 $span\left(\mathcal{A}\right)$ 的一个基（basis），对于任何 $\vec{b} \in span\left(\mathcal{A}\right)$，在此 basis 下有唯一表示，即 $A\vec{x}=\vec{b}$ 只有唯一解，即 $\vec{b}$ 在 $A$ 的列向量组 $\mathcal{A}$ 下的坐标 $\displaystyle \left[\vec{b}\right]_\mathcal{A}\in\mathbb{R^n}$。从 $U\ \vec{x}=P\ \vec{b}$ 也可以看出，echelon form $U$ 的每一列都是主元所在列，所以 $\vec{x}$ 只有唯一解。

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
C\ \vec{x}_{free}\\I\ \vec{x}_{free}
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
考虑独热部分 $I_{n-k,n-k}$，所以 $rank(Q) = n-k$，同时 $\vec{x}_{free}$ 取遍 $\mathbb{R^{n-k}}$，所以此时解集 $set$ 是由 $Q$ 的列向量组 $\mathcal{Q}$ 张成的线性空间，即 $\mathbb{R}^{n}$ 的一个 $dim = n-k$的subspace。  
---- 这里，引入概念**矩阵A的零空间Nul(A)**，即矩阵 $A$ 对应的矩阵方程 $A\vec{x}=\vec{b}$ 的解集。通过上述证明，可以得出矩阵 $A$ 的零空间 $Nul(A)$ 是 $\mathbb{R}^{n}$ 的一个 $dim = n - rank(A)$ 的subspace，它的基是矩阵方程 $A\vec{x}=\vec{b}$ 的 $n-rank(A)$ 个线性无关的解。

#### 非齐次
特解+系数矩阵的零空间，即为非齐次线性方程组的解集。

---
### LU分解  
矩阵 $A$ 转化为 echelon form $U$ 的初等行变换操作中，因为所有的「replacement」都是用主对角线上的非零元素去消除下方元素，所以对应的初等矩阵都是下三角单位矩阵（主对角线元素都是1）。在这个过程中，若矩阵的主对角线元素（不断变化）都是非零的，那么全程不需要「interchange」操作，则化简过程中所有的初等矩阵 $E_p,\cdots,E_1$ 都是下三角单位矩阵，从而 $P = E_p\ \cdots\ E_1$ 以及 $P^{-1}$ 都是下三角单位矩阵，将 $P^{-1}$ 记作 $L$，即得到
```math
E_p\ \cdots\ E_1 A = PA = U
```
即得：
```math
A = P^{-1}U = LU
```
其中 $L$ 是下三角单位矩阵，而 $U$ 是 $A$ 的echelon form 阶梯型。这就是矩阵 $A$ 的LU分解。从前述论述中可以看出，不是所有矩阵都可以作LU分解，可以LU分解的矩阵需要满足在初等行变换至echelon form的过程中，不能出现「interchange」操作（即主对角线元素上不能出现零）。

## 主干2综述
主干2脉络中，主线是**子空间**的**基**，**基**的大小（size，即所含向量个数）就是子空间的**维度**。用**基**唯一表示时的权重即为**坐标**。  
  
### 线性空间
完整的**线性空间**定义：定义在**域** $F$ 上的非空集合 $V$，：
* 定了**向量加法**，满足
    - $V$ 中存在零元素
    - $V$ 可定义关于向量加法的负元素
    - 向量加法具有结合律
    - 向量加法具有交换率
* 定了**数量乘法**，满足
    - $F$ 的乘法单位元，使得数量乘法映射成为恒等映射
    - 对向量加法满足分配律
    - 对域 $F$ 的加法满足分配律
    - 对域 $F$ 的乘法满足结合律
  
完整的 $V$ 的**子空间** $H$ 的定义:
* 包含 $V$ 的 $0$ 向量  
* $H$ 对 $V$ 的加法和数乘封闭  

由此顺利引入**向量组张成的子空间** $span(\mathcal{B})$ 的概念，因为 $V$ 的向量张成的子集正好满足子空间的定义。  
一个子空间 $H$ 的任一组**基**都是能互相表示(等价)的，而且能够张成 $H$ 的 $dimH$ 个向量，就是 $H$ 的一组基。

### 同构/同构映射
**同构映射**是指**保持结构的双射**，即：  
考虑定义了 $+$ 运算的群 $G$ 和定义了 $\times$ 运算的群 $F$，有 $\vec{a}, \vec{b} \in G, T(\vec{a}), T(\vec{b}) \in F$。若 $one-one \  Mapping \  T:G \rightarrow F$，满足 $T(\vec{a} + \vec{b}) = T(\vec{a}) \times T(\vec{b})$，称 $T$ 是一个同构映射，群 $G$ 和 $F$ 同构，记作 $G \cong F$。  
  
同理，线性空间之间的**同构**指：**同一个域上定义的两个线性空间之间，存在保持向量加法和数乘的双射。**  
  
思考：线性映射也能保持运算，那么线性映射是同构映射吗？
答：不是。同构映射要求该映射是「双射」  
  
线性空间同构的充要条件：
* 域F上两个有限维的线性空间 $G$ 和 $F$，则 $G$ 和 $F$ 同构 $\leftrightarrow$ $G$ 和 $F$ 维数相同。   
  
#### 子空间和坐标空间之间的同构关系
考虑 $dim H = p$，向量 $\vec{x} \in subspace\ H \subset \mathbb{R^m}$ 被其基 $\mathcal{B}$ 唯一表示时的权重，就是向量 $\vec{x}$ 在这个基下的**坐标向量**，简称**坐标** $\left[\vec{x}\right]_\mathcal{B} \in \mathbb{R^p}$。  
  
考虑映射 $T:H \rightarrow \mathbb{R^p}, \mathcal{B}\ is\ a\ basis\ of\ H$，
```math
T(\vec{x}) = \left[\vec{x}\right]_{\mathcal{B}}
```
即将一个 $H$ 中的向量映射到 基basis下的坐标向量。易证这是一个保持**加法**和**数乘**运算的**同构映射**，即 $H \cong \mathbb{R^p}$。形象地理解，有 $H$ “像是” $\mathbb{R^p}$，尽管 $H$ 中的元素的分量数目可能是超过 $p$的。  
    
### 矩阵的零空间和列空间/线性映射的kernel和range
在前述中已经介绍了矩阵 $A$ 的**列空间Col(A)**和**零空间Nul(A)**。  
在这里，引入概念矩阵 $A \in \mathbb{R^{m,n}}$ 的**行空间Row(A)**，即 $A$ 的行向量组张成的空间，是 $\mathbb{R^n}$ 的一个子空间。
矩阵 $A$ 的echelon form $U$ 的非零行是它的一个**基**（注意这里是 $U$ 的非零行，跟列空间不同，列空间的基是从 $A$ 中选 $U$ 的主元所在列。原因也很简单：因为 $U$ 是 $A$ 作初等行变换得到的，所以它们的行向量组互相等价，所以 $U$ 的一个基也是 $A$ 的基）。  

矩阵 $A$ 的**零空间**和**列空间**的**维度**和**基**，都已经知道怎么计算：
* 零空间 Nul(A)：求解矩阵方程 $A\vec{x} = \vec{0}$，等价于 $U\vec{x} = \vec{0}$，把主元和自由元都用自由元表示，表示的系数矩阵的列向量组就是 Nul(A) 的一个基，自由元个数 $n - rank(A)$ 就是 Nul(A) 的维度。  
* 列空间 Col(A)：对 $A$ 作初等行变换至阶梯型矩阵 echelon form $U$，主元所在列（从 $A$ 中找）就是 Col(A) 的一个基，非零行数量等于 $rank(A)$ 就是 Col(A) 的维度。  

如果从**线性映射**的视角出发：  
考虑线性映射 $T:\mathbb{R^n} \rightarrow \mathbb{R^m},\ with\ matrix\ A \in \mathbb{R^{m,n}}$，那么
1. 零空间 $Nul(A)$是 $domain\ \mathbb{R^n}$ 的一个子空间，代表这个子空间中的向量在映射 $T$ 下都被映射到了 $codomain\ \mathbb{R^m}$ 的 $\vec{0}$。
2. 列空间 $Col(A)$是 $codomain\ \mathbb{R^m}$ 的一个子空间，可以记作 $span(\left[\vec{a}_1, \vec{a}_2, \cdots, \vec{a}_n\right])$，代表映射T的在 $codomain\ \mathbb{R^m}$ 的最大可达范围（值域）。  
  
* 通过探讨映射 $T$ 是否是单射，得出：当 $T$ 是单射时，有
    - $rank\ A = n$
    - $A$ 的列向量组线性无关
    - $A\vec{x} = \vec{0}$ 只有0解  
    
  这三个相互等价的条件。此时 $Nul(A) = set(\vec{0})$ 。
  
* 通过探讨映射 $T$ 是否是满射，得出：当 $T$ 是满射时，有
    - $rank\ A = m$
    - $A$ 的列向量组张成 $\mathbb{R^m}$（记住 $dim\ \mathbb{R^m} = m$）  
    
  这两个相互等价的条件。此时 $Col(A) = \mathbb{R^m}$ 。
  
#### kernel和range
如果将线性映射的概念从**数值向量空间**之间的映射，推广到广义**线性空间**之间的映射 $T$，那么相应地，
* 把 $T$ 的“零空间” $Nul(T)$ 定义为 $T$ 的**kernel**，是**domain**的子空间：映射 $T$ 把kernel中的元素 $\vec{v}$ 映射到codomain的 $\vec{0}$ 元素
* 把 $T$ 的“列空间” $Col(T)$ 定义为 $T$ 的**range**，是**codomain**的子空间：对于range中的元素 $\vec{u}$，能从domain中找到对应元素 $\vec{v}$，使得 $T(\vec{v})=\vec{u}$
  

### 基/坐标/维度/秩
对于无限维线性空间，线性代数中不作过多研究，仅需知道它的维度定义为无限。 
  
对于**有限维线性空间**，即可以由有限个向量张成的**线性空间**，关于其向量的**线性组合**、**线性相关/无关**、**张成span**、**极大线性无关组**的定义都和上述如出一辙，从而**基**和**坐标**的定义和性质也相同：即通过极大线性无关组定义向量组的基，然后通过向量组的基定义张成子空间的基。  

#### 基定理 Base theorem
考虑 $\mathbb{R^{n}}$ 的线性子空间 $H$，有 $rank(H) = k$，则有：
* 任一个**线性无关**的向量组 $\left\\{\vec{h_1},\vec{h_2},\cdots, \vec{h_k} \right\\}$ 是 $H$ 的一个**基**  
* 任一个向量组 $\left\\{\vec{h_1},\vec{h_2},\cdots, \vec{h_k} \right\\}$ ，若它可以**张成** $H$，则它是 $H$ 的一个**基**

#### 坐标/坐标映射 coordinates-mapping
考虑线性空间 $V$ 的一个**基**
```math
\mathcal{B} = \left[\vec{b}_1,\cdots,\vec{b}_n\right]
```
对任意 $\vec{x} \in V$，相对这个基都有唯一表示
```math
\vec{x} = \sum_{i=1}^{n}c_i\ \vec{b}_i
```
把权重按序表示成列向量 $\left[c_1,\cdots,c_n\right]^T$，即是 $\vec{x}$ 相对基 $base\ \mathcal{B}$ 的**坐标**，写作
$\left[\vec{x}\right]_\mathcal{B} \in \mathbb{R^n}$。  
  
考虑线性空间 $V$ 到 $\mathbb{R^n}$ 的**坐标映射**
```math
T:V \rightarrow \mathbb{R^n},\ T(\vec{x})=\left[\vec{x}\right]_\mathcal{B}
```
是一个关于向量加法和数乘的**同构映射**，线性空间 $V$ 和 $\mathbb{R^n}$ 同构。同构意味着坐标映射保持加法和数乘运算，即意味着**线性相关/线性无关等关系在同构空间中也是保持的**。这里要提一句，只有同构映射才能保持线性相关/线性无关等性质，因为同构映射保证是双射。普通的线性映射下，线性相关/线性无关等性质并不能保持。  
  
**坐标映射**有点像是**线性组合**的逆运算：后者以**基**线性合成出向量，前者线性分解向量找到权重得到坐标。  


#### 回到n维数值向量空间：看坐标转换及其相关的意义

当 $V = \mathbb{R^n}$ 时，根据「坐标是一个列向量，**基**以**坐标**为权重作**线性组合**得到原向量」的原则，可得此时这是一个matrix-vector product运算：
```math
\vec{x} = P_{\mathcal{B}}\ \left[\vec{x}\right]_{\mathcal{B}},\ P_{\mathcal{B}} = \left[\vec{b}_1,\cdots,\vec{b}_n\right]
```
而向量 $\vec{x}$ 自身可以看作其在**标准基**下的坐标向量，也就是说，基 $\mathcal{B}$ 作为列向量组组成的矩阵 $P_{\mathcal{B}}$ 成为了两个基之间的**坐标转移矩阵(change-of-coordinates matrix)**。在这里，矩阵 $P_{\mathcal{B}}$ 是从 $\mathcal{B}$ 到标准基的坐标转移矩阵，而  $P_{\mathcal{B}}^{-1}$ 是从标准基到 $\mathcal{B}$ 的坐标转移矩阵，因为有：
```math
\left[\vec{x}\right]_{\mathcal{B}} = P_{\mathcal{B}}^{-1}\ \vec{x}
```
   

#### 坐标变换/坐标转移
前面已经讲述过，当**线性空间** $V, dimV = n$ 就是 $\mathbb{R^n}$ 时，对于它的一个**基** $\mathcal{B}$ ，有
```math
\vec{x} = P_{\mathcal{B}}\ \left[\vec{x}\right]_{\mathcal{B}}
```
这里 $P_{\mathcal{B}}$ 就是 $\mathcal{B}$ 作为列向量组组成的矩阵。然而对于一般的线性空间 $V$，很可能没办法把基 $\mathcal{B}$ 写成矩阵的形式。  
对于一般的线性空间 $V, dimV = n$，考虑它的两个基 $\mathcal{B}$ 和 $\mathcal{C}$，那么对于 $V$ 中的任何一个 $\vec{x} \in V$，它在两个基下有不同的坐标，分别是
```math
\left[\vec{x}\right]_\mathcal{B},\ \left[\vec{x}\right]_\mathcal{C}\ \in\mathbb{R^n}
```
从 $\mathcal{B}$ 下的坐标向量到 $\mathcal{C}$ 下的坐标向量，存在**唯一**的**线性映射** $T:\mathbb{R^n} \rightarrow \mathbb{R^n},\ \ with\ \ matrix\ \ P_\mathcal{C \leftarrow B} \in \mathbb{R^n}$，使得有：
```math
\left[\vec{x}\right]_\mathcal{C} = P_\mathcal{C \leftarrow B}\ \left[\vec{x}\right]_\mathcal{B}
```
这个线性映射被称为**坐标变换**或**坐标转移**映射，矩阵 $P_\mathcal{C \leftarrow B}$ 被称为**坐标转移矩阵**，其列向量是基 $\mathcal{B}$ 里的各向量分别按序在另一个基 $\mathcal{C}$ 下的坐标向量，即
```math
Consider\ \mathcal{B}\ as\ set(\vec{b}_1,\cdots,\vec{b}_n),\ P_\mathcal{C \leftarrow B}=
\begin{bmatrix}
\left[\vec{b}_1\right]_\mathcal{C},\ \cdots\ \left[\vec{b}_n\right]_\mathcal{C}
\end{bmatrix} \in \mathbb{R^n}
```
坐标转移映射都是可逆的，即**坐标转移矩阵**都是**可逆**的，显然有 $P_\mathcal{B \leftarrow C} = P_\mathcal{C \leftarrow B}^{-1}$ 。
  

## 后述
### 线性映射为什么是矩阵
前面给出过一个结论，即：任何一个线性映射 $T:\mathbb{R^n} \rightarrow \mathbb{R^m}$ 一定和一个矩阵 $A\in\mathbb{R}^{m,n}$ 一一对应，而且**映射的像就是矩阵乘以向量元素的结果**。现在使用**坐标**的角度，来证明这个结论。注意定义坐标不需要当前这个结论，所以并没有循环论证。  
  
考虑一个从 $V$ 到 $W$ 的线性映射 $T: V \rightarrow W$，有 $\vec{x} \rightarrow T(\vec{x})$。这里 $dimV = n$，基 $\mathcal{B} = \\{\vec{b_1},\vec{b_2}, \cdots, \vec{b_n}\\}$ 是 $V$ 的一个基；还有 $dimW = m$，基 $\mathcal{C} = \\{\vec{c_1},\vec{c_2}, \cdots, \vec{c_m}\\}$ 是 $W$ 的一个基。  
考虑 $x = r_1\vec{b_1}+r_2\vec{b_2}+\cdots+r_n\vec{b_n}$, 其坐标向量为：
```math
\left[\vec{x}\right]_\mathcal{B} = \begin{bmatrix}
r_1\\
r_2\\
\vdots\\
r_n
\end{bmatrix}
```
线性映射 $T$ 的像image有 $T(\vec{x}) = T(r_1\vec{b_1}+r_2\vec{b_2}+\cdots+r_n\vec{b_n}) = r_1T(\vec{b_1})+r_2T(\vec{b_2})+\cdots+r_nT(\vec{b_n})$。考虑其在基 $\mathcal{C}$ 下的坐标向量，由于坐标映射保持加法和数乘，可得：
```math
\left[T(\vec{x})\right]_\mathcal{C} = \left[r_1T(\vec{b_1})+r_2T(\vec{b_2})+\cdots+r_nT(\vec{b_n})\right]_\mathcal{C} = 
r1\left[T(\vec{b_1})\right]_\mathcal{C}+r_2\left[T(\vec{b_2})\right]_\mathcal{C}+\cdots+r_n\left[T(\vec{b_n})\right]_\mathcal{C} = 
\left[\left[T(\vec{b_1})\right]_\mathcal{C}, \left[T(\vec{b_2})\right]_\mathcal{C}, \cdots, \left[T(\vec{b_n})\right]_\mathcal{C}\right] \begin{bmatrix}
r_1\\
r_2\\
\vdots\\
r_n
\end{bmatrix}
```
令：
```math
\left[\left[T(\vec{b_1})\right]_\mathcal{C}, \left[T(\vec{b_2})\right]_\mathcal{C}, \cdots, \left[T(\vec{b_n})\right]_\mathcal{C}\right] = M \in \mathcal{R}^{m,n}
```
可得
```math
\left[T(\vec{x})\right]_\mathcal{C} = M \left[\vec{x}\right]_\mathcal{B}
```
综上，线性映射意味着一个 $n$ 维线性空间（定义域空间）到 $m$ 维线性空间（陪域空间）的映射。
同步地作映射，向量 $\vec{x}$ 在定义域空间的基下的坐标
$\left[\vec{x}\right]_\mathcal{B}$ 左乘一个 $m$ 行 $n$ 列的系数矩阵 $M$ 后，得到像image $T(\vec{x})$ 在陪域空间的基下的坐标，即：
```math
M \left[\vec{x}\right]_\mathcal{B} = \left[ T(\vec{x}) \right]_\mathcal{C}
```
这里对应的系数矩阵是**定义域线性空间的基向量映射到陪域空间后，它们的像image在陪域空间的基下的坐标**。
  
### 为什么变换/方阵有相似
**坐标**概念的建立，有一个很大的意义在于，原**线性空间** $V,\ dimV = n$ 中的向量可能不方便用数字表示，无法参与进一步的分析（比如研究 $V$ 到自身的线性变换）。但是如果去研究与之**同构**的 $\mathbb{R^n}$，即**坐标向量**和**坐标空间**，就方便了。特别地，如果研究从 $V$ 到 $V$ 的线性变换，就可以转而研究从坐标空间 $\mathbb{R^n}$ 到 $\mathbb{R^n}$ 的线性变换，即一个形状为(n,n)的方阵。  
  
甚至有时候，一个**基** $\mathcal{B}$ 得到的坐标空间 $\mathbb{R^n}\_{\mathcal{B}}$ 可能还「不够好」，我们会换一个**基** $\mathcal{C}$，即换一个坐标空间 $\mathbb{R^n}\_{\mathcal{C}}$。  
  
更具体地说，一个线性变换 $T:\mathbb{R^n}\_{\mathcal{B}}\rightarrow\mathbb{R^n}\_{\mathcal{B}},\ matrix\ T\in \mathbb{R^{n,n}}$，与另一个线性变换 $S:\mathbb{R^n}\_{\mathcal{C}}\rightarrow\mathbb{R^n}\_{\mathcal{C}},\ matrix\ S\in \mathbb{R^{n,n}}$，是 $V$ 的两个不同基下的**坐标空间**里各自的线性变换，但其实它们都是 $V$ 到 $V$ 上的一个线性变换在不同基下**同构**的坐标空间里的“翻版”，此谓两个变换**相似**，或者说两个**矩（方）阵相似**。  
  
用数学来表示，即：  
```math
V(dimV = n,\ base\ \mathcal{B}\ \&\ \mathcal{C}),\ V\cong\mathbb{R^n}_{\mathcal{B}},\ V\cong\mathbb{R^n}_{\mathcal{C}}
```
现在，有一个 $V$ 到 $V$ 的线性变换
```math
F:V\rightarrow V,\ \ F(\vec{x})=\vec{y},\ \vec{x}\in V,\ \vec{y} \in V
```
考虑 $F$ 在两个坐标空间的“翻版”，即如下两个线性变换：
```math
F_\mathcal{B}:\mathbb{R^n}_{\mathcal{B}}\rightarrow\mathbb{R^n}_{\mathcal{B}},\ matrix\ T\in \mathbb{R^{n,n}},\ T\left[\vec{x}\right]_\mathcal{B} = \left[\vec{y}\right]_\mathcal{B}
```
```math
F_\mathcal{C}:\mathbb{R^n}_{\mathcal{C}}\rightarrow\mathbb{R^n}_{\mathcal{C}},\ matrix\ S\in \mathbb{R^{n,n}},\ S\left[\vec{x}\right]_\mathcal{C} = \left[\vec{y}\right]_\mathcal{C}
```
因为有坐标转移映射
```math
\left[\vec{x}\right]_\mathcal{C}=P_\mathcal{C \leftarrow B}\left[\vec{x}\right]_\mathcal{B},\ P_\mathcal{C \leftarrow B}\in\mathbb{R^{n,n}}
```
```math
\left[\vec{y}\right]_\mathcal{B}=P_\mathcal{B \leftarrow C}\left[\vec{y}\right]_\mathcal{C},\ P_\mathcal{B \leftarrow C}\in\mathbb{R^{n,n}}
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

### 相似的意义
我们稍微联系一下特征向量和特征值的意义，就可以发现**相似**的意义。 

#### 回顾 
首先复习一下**特征向量**和**特征值**的概念。考虑**坐标空间** 
```math
\mathbb{R^n}_\mathcal{B},\ with\ \ basis\ \ \mathcal{B} = set\left(\vec{b_1}, \cdots,\vec{b_n}\right)
```
在这个坐标空间中，有一个线性变换 $A$
```math
A:\mathbb{R^n}_\mathcal{B}\rightarrow\mathbb{R^n}_\mathcal{B},\ \left[\vec{x}\right]_\mathcal{B}\rightarrow A\left[\vec{x}\right]_\mathcal{B}
```
如果对**非零向量坐标** $\left[\vec{e}\right]_{\mathcal{B}}$ 和实数 $\lambda \in \mathbb{R}$ 有下式成立：
```math
A\left[\vec{e}\right]_{\mathcal{B}} = \lambda \left[\vec{e}\right]_{\mathcal{B}}
```
则称向量坐标 $\left[\vec{e}\right]_{\mathcal{B}}$ 和实数 $\lambda \in \mathbb{R}$ 分别是矩阵A的**特征向量**和**特征值**。  

#### 推导
我们来推导下**同构坐标空间**中，映射(或矩阵) $A$ 的**相似**映射(矩阵) $A'$ 的特征向量和特征值。  
  
考虑同构坐标空间 $\mathbb{R^n}_{\mathcal{C}}$ ,
```math
\mathbb{R^n}_\mathcal{C},\ with\ \ basis\ \ \mathcal{C} = set\left(\vec{c_1}, \cdots,\vec{c_n}\right)
```
写出坐标转移矩阵：
```math
P_{\mathcal{B}\leftarrow \mathcal{C}} = \left[\left[\vec{c_1}\right]_\mathcal{B},\cdots,\left[\vec{c_n}\right]_\mathcal{B}\right]
```
```math
P_{\mathcal{C}\leftarrow \mathcal{B}}  = P_{\mathcal{B}\leftarrow \mathcal{C}}^{-1}
```
利用前面的**相似**的概念，矩阵 $A$ 和 $A'$相似，那么可以得出坐标空间 $\mathbb{R^n}_{\mathcal{C}}$ 中的线性映射 $A'$，有：
```math
A\left[\vec{x}\right]_\mathcal{B} = P_{\mathcal{B}\leftarrow \mathcal{C}}A'P_{\mathcal{B}\leftarrow \mathcal{C}}^{-1}\left[\vec{x}\right]_\mathcal{B}
```
恒成立。把它写成另一个形式
```math
A\left[\vec{x}\right]_\mathcal{B} = P_{\mathcal{B}\leftarrow \mathcal{C}}A'\left[\vec{x}\right]_\mathcal{C}\tag1
```
也就是说，在各自坐标空间下，线性映射 $A$ 和 $A'$ 的结果只差一个坐标转移。  
  
假设已知矩阵 $A$ 的一对特征向量 $\vec{e}$ 和特征值 $\lambda$，即
```math
A \left[\vec{e}\right]_\mathcal{B} = \lambda \left[\vec{e}\right]_\mathcal{B}\ \tag{2.1}
```
(2.1)代入(1)式，得到
```math
\lambda \left[\vec{e}\right]_\mathcal{B} = A \left[\vec{e}\right]_\mathcal{B}=P_{\mathcal{B}\leftarrow \mathcal{C}}A'\left[\vec{e}\right]_\mathcal{C}
```
即
```math
P_{\mathcal{B}\leftarrow \mathcal{C}}^{-1}\lambda \left[\vec{e}\right]_\mathcal{B} = A'\left[\vec{e}\right]_\mathcal{C}
```
即
```math
\lambda P_{\mathcal{B}\leftarrow \mathcal{C}}^{-1} \left[\vec{e}\right]_\mathcal{B} = A'\left[\vec{e}\right]_\mathcal{C}
```
即
```math
A'\left[\vec{e}\right]_\mathcal{C}=\lambda \left[\vec{e}\right]_\mathcal{C},\ \tag{2.2}
```
(2.2)和(2.1)式对比，说明了一个很有意思的道理：「转换坐标空间之后，**原特征向量在新基下的坐标向量，是原矩阵的相似矩阵的特征向量，且对应的特征值不变**」。也就是说，特征向量和特征值是矩阵**本征**的性质，跟坐标系（基）的选取无关。    

重点来了，假如说，现在 $A$ 有 **n 个特征向量 $set\left(\vec{e_1},\cdots,\vec{e_n}\right)$ 且能组成一个基**，即考虑同构坐标空间 $\mathbb{R^n}_{\mathcal{E}}\ \ with\ \ basis\ \ set\left(\vec{e_1},\cdots,\vec{e_n}\right)$ ，将基 $\mathcal{E}$ 代入(2.2)式，这时有
```math
A'\left[\vec{e_i}\right]_\mathcal{E}=\lambda_i \left[\vec{e_i}\right]_\mathcal{E},\ i = 1,\cdots,n\tag{3}
```
而基向量在基下的坐标是可求的，即
```math
\begin{cases}
\left[\vec{e_1}\right]_\mathcal{E} =
\left[
    \begin{matrix}1\\
    \vdots\\
    0
    \end{matrix}\right]
\\
\vdots\\
\left[\vec{e_n}\right]_\mathcal{E} =
\left[
    \begin{matrix}0\\
    \vdots\\
    1
    \end{matrix}\right]
\end{cases}
\tag4
```
将(4)代入(3)式，得到了 n 个 n元一次方程组，共 $n^2$ 个方程，足够求出 $A'$ 。事实上，将 $A'$ 写作列向量组的形式 $\left[\vec{a'_1},\cdots,\vec{a'_n}\right]$，即得：
```math
\left[\vec{a'_1},\cdots,\vec{a'_n}\right]
\left[
    \begin{matrix}0\\
    \vdots\\
    1\\
    \vdots\\
    0
    \end{matrix}
\right]
    \begin{matrix}\\
    \ 
    \\
    i_{th}\ row\\
    \\
    \ 
    \end{matrix}
\ \ \ =\ \ \ 
\lambda_i
\left[
    \begin{matrix}0\\
    \vdots\\
    1\\
    \vdots\\
    0
    \end{matrix}
\right]
    \begin{matrix}\\
    \ 
    \\
    i_{th}\ row\\
    \\
    \ 
    \end{matrix}
,\ \ \ i = 1,\cdots,n
```
即得：
```math
\vec{a'_i} = 
\left[
    \begin{matrix}0\\
    \vdots\\
    \lambda_{i}\\
    \vdots\\
    0
    \end{matrix}
\right]
    \begin{matrix}\\
    \ 
    \\
    i_{th}\ row\\
    \\
    \ 
    \end{matrix}
,\ \ \ i = 1,\cdots,n
```
即 $A' = diag\left(\lambda_1,\cdots,\lambda_n\right)$。回到(0)式，考虑基 $\mathcal{B}$ 是标准基，那么从 $\mathbb{R^n}\_{\mathcal{E}}$ 到 $\mathbb{R^n}\_{\mathcal{B}}$ 的坐标转移矩阵 $P_{\mathcal{B}\leftarrow \mathcal{E}}$，
```math
P_{\mathcal{B}\leftarrow \mathcal{E}} = \left[\left[\vec{e_1}\right]_\mathcal{B},\cdots,\left[\vec{e_n}\right]_\mathcal{B}\right]
```
即
```math
\begin{cases}
P = \left[\vec{e_1},\cdots,\vec{e_n}\right],\ \vec{e_i}\ \ are\ \ eigenvectors\ \ of\ \ A,\ i=1,\cdots,n\\
A = PA'P^{-1},\ \ A' = diag\left(\lambda_1,\cdots\lambda_n\right)
\end{cases}
\tag5
```
  
去掉一切推导归纳总结：如果 矩阵 $A$ 有 **n 个特征向量线性无关**（即能组成一个基），那么在这个**新基的坐标空间**中，矩阵 $A$ 的相似矩阵 $A'$ 是一个**对角矩阵**，对角线上是 n 个特征向量对应的特征值。n个特征向量作为列向量组构成的矩阵 $P$ 是新基到标准基的**坐标转移矩阵**，即上面的(5)式。  
  
#### 启发
到这里，可以看出，如果能知道一个矩阵的 **n 个线性无关的特征向量**和对应的**特征值**， 那么此时，就得到了一个式子，将这个矩阵分解成可逆矩阵和对角矩阵的乘积，又叫**对角化**。所以下一章开始，我们将重点求解一个矩阵的特征向量和特征值，并判断这 n 个特征向量是否线性无关。