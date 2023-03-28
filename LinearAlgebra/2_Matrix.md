# 综述: 行列式/对角化(特征向量和特征值)/特征多项式
贯穿全章节的计算手段主要是:
* 行列式(determinant)
* 求解特征值和特征向量(没有多项式观点)(eigenvalues & eigenvector without polynomial)
    1. 已知矩阵和特征向量，求特征值
    2. 已知矩阵特征值，求特征向量和特征空间
    3. 仅知矩阵 $A$，特征值和特征向量都未知，求解其特征值和特征向量(即求解特征方程)
* 矩阵的对角化判断和计算

这块内容首先介绍了**行列式**，然后介绍了**矩阵**的**特征向量**和**特征值**，然后接下来有两个线性代数这门课程的超级重点：
* 重点1：从上一章Linear最后的**线性映射(矩阵)的相似**开始，配合**特征向量**和**特征值**，引出了**对角化**，并给出了**对角化**的第一个充要条件
* 重点2：引入**特征方程**和**特征多项式**之后，首先介绍**相似**矩阵的**特征向量/特征值/特征方程**，然后切换成代数的视角思考**对角化**的必要条件和充分条件：通过探讨不同特征值的特征向量的线性无关性质，以及**特征值**的**几何重数**和**代数重数**，给出**对角化**的第二个充要条件

## 行列式
行列式的定义：行列式是一个从数域 $K$ 上的n阶矩阵到数域 $K$的映射，即 $F: M_n\left(K\right) \rightarrow K$，对 $A \in M_n\left(K\right)$，记它的行列式为 $det\left(A\right)$，或者是 $|A|$。  
  
### 行列式的完全展开  
对于 $A \in M_n\left(K\right)$，它的行列式的完全展开表达，是 $n!$ 个乘积项的和：
```math
A = \left(a_{ij}\right),\ det\left(A\right) = \sum_{j_1\ j_2 \cdots\ j_n}(-1)^{\tau\left(j_1\ j_2\ \cdots\ j_n\right)}a_{1j_1}\ a_{2j_2}\ \cdots\ a_{nj_n}
```
其中 $j_1 j_2 \cdots j_n$ 代表一个 $1,2,\cdots n$ 的n元排列，总共有 $n!$ 项。符号 $\tau\left(j_1 j_2 \cdots j_n\right)$ 代表该排列的逆序数。  
  
从行列式的完全展开表达式可以看出，**行列式映射对于矩阵 $A$ 的某一列/行（固定其他列/行）是线性的**。即：
1. 如果矩阵的某一列/行乘以一个数 $k$，那么行列式也相应乘以 $k$
2. 如果矩阵的某一列/行可以看作两列/行之和，那么行列式等于两个分拆矩阵的行列式之和。
  
从行列式的完全展开表达式可以看出，**行列式映射对于矩阵 $A$ 的行列转换是不变的**。即：
```math
det\left(A\right) = det\left(A^{T}\right)
```
  
### 行列式的代数余子式展开（cofactor expansion）
对于 $A \in M_n\left(K\right)$，它的行列式的代数余子式展开表达，是 $n$ 个项的和。考虑
```math
A = \left(a_{ij}\right),\ \ C_{ij} = (-1)^{i+j}det(A_{ij})
```
这里 $A_{ij}$ 是矩阵 $A$ 去掉 $i$ 行和 $j$ 列之后的子矩阵，项 $C_{ij}$ 被称作矩阵 $A$ 的**代数余子式**。可以看出**代数余子式展开是一种递归展开**。  
其展开式表达为：
1. 行列式按 $i$ 行展开
```math
det\left(A\right) = a_{i1}C_{i1}+a_{i2}C_{i2}+\cdots+a_{in}C_{in}
```
2. 行列式按 $j$ 列展开
```math
det\left(A\right) = a_{1j}C_{1j}+a_{2j}C_{2j}+\cdots+a_{nj}C_{nj}
```
从行列式的代数余子式展开可以得出，**上/下三角矩阵的行列式等于对角线元素的乘积**。

### 行列式的矩阵乘积展开
通过行列式的完全展开定义，和矩阵乘法定义，可以证明（证明略，无非是整理项和系数）**行列式映射保持矩阵的乘法运算**，即：
```math
A, B \in M_n\left(K\right), \ \ \ det\left(AB\right) = det\left(A\right)det\left(B\right)
```
已经知道，任意一个矩阵 $A$ 都和它的简约阶梯型(reduced echelon form)矩阵 $U$ 行等价，即互相可以用初等行变换转化，即可以通过左乘相应的初等矩阵转化，即：
```math
A = E_p E_{p-1} \cdots E_{1} U
```
从而有
```math
det\left(A\right) = det\left(E_p\right) det\left(E_{p-1}\right) \cdots det\left(E_1\right) det\left(U\right)
```
这里 $U$ 和 $E$ 的行列式都很好求：
```math
det\left(E\right) =
\begin{cases}
1 & if\ \ \ E\ \ \ is\ \ \ row\ \ \ replacement\\
-1 & if\ \ \ E\ \ \ is\ \ \ interchange\\
r & if\ \ \ E\ \ \ is\ \ \ scale\ \ \ by\ \ \ r \neq 0
\end{cases}
```
从行列式的矩阵乘积展开可以得出：
1. **当 $A$ 可逆时，行列式 $det\left(A\right) \neq 0$**，因为可逆矩阵（invertible/nonsingular）一定和单位矩阵 $I$ 行等价。
2. **当 $A$ 不可逆时，行列式 $det\left(A\right) = 0$**，因为不可逆矩阵（singular）的阶梯型矩阵 $U$ 一定有零行。
  
### 行列式的几何意义
一个 $R^n$ 中的集合 $S$，在 $R^n \rightarrow R^n$ 的线性映射 $T$ 下，“体积”变换前后的伸缩系数，即：
```math
\left\{
\ Area\ \ of\ \ T(S)
\right.
\left.
\right\}
= det(T) *
\left\{
\ Area\ \ of\ \ S
\right.
\left.
\right\}
```