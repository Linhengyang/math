# 综述: 行列式/特征向量和特征值/特征方程和特征多项式/相似和对角化
贯穿全章节的计算手段主要是:
* 行列式(determinant)
    1. 行列式的完全展开
    2. 行列式的代数余子式展开（递归展开 or 按一行/列展开）
    3. 行列式的矩阵乘积展开
* 求解特征值和特征向量(一点多项式观点)(eigenvalues & eigenvector)
    1. 已知矩阵和特征向量，求特征值
    2. 已知矩阵特征值，求特征向量和特征空间
    3. 仅知矩阵 $A$，特征值和特征向量都未知，求解其特征值和特征向量(即求解特征方程 or 因式分解特征多项式)
* 矩阵的对角化判断和计算
    1. 对角化的3个充要条件
    2. 计算对角化：求出所有不同特征值，然后求出每个特征值的特征空间的基

正文：  
这块内容首先介绍了**行列式**，然后介绍了**矩阵**的**特征向量**和**特征值**，然后接下来有两个线性代数这门课程的超级重点：
* 重点1：从上一章Linear最后的**线性映射(矩阵)的相似**开始，配合**特征向量**和**特征值**，引出了**对角化**，并给出了**对角化**的第一个充要条件。这块在上一章Linear已经谈过了，不再赘述。
* 重点2：引入**特征方程**和**特征多项式**之后，首先介绍**相似**矩阵的**特征向量/特征值/特征方程**，然后切换成代数的视角思考**对角化**的必要条件和充分条件：通过探讨不同特征值的特征向量的线性无关性质，以及**特征值**的**几何重数**和**代数重数**，给出**对角化**的第二个充要条件。

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

## 特征向量/特征值/特征方程/特征多项式/特征子空间
这一块围绕着**特征向量**和**特征值**，给出很多冠以特征的名词概念。注意这里还不涉及相似和对角化。
  
* **特征向量**：对于数域 $K$ 上的方阵 $A \in M_{n}(K)$，若**非零**向量 $\vec{x}$ 满足 $A\vec{x}=\lambda\vec{x},\ \ \lambda\in K$，则 $\vec{x}$ 是**矩阵 $A$ 关于特征值 $\lambda$ 的特征向量**（eigenvector of $A$ corresponding to $\lambda$）。
* **特征值**：对于数域 $K$ 上的方阵 $A \in M_{n}(K)$，若 $A\vec{x}=\lambda\vec{x},\ \ \lambda\in K$ 存在**非零解** $\vec{x}$，则 $\lambda$ 是矩阵 $A$ 的**特征值**，且 $\vec{x}$ 是**矩阵 $A$ 关于特征值 $\lambda$ 的特征向量**（eigenvector of $A$ corresponding to $\lambda$）。
  
从定义可以看出，特征值和特征向量是成对出现的。给定一个矩阵，求解它的**特征值**和**特征向量**成为了这块的主线。写出求解的矩阵方程，
```math
A\vec{x}=\lambda\vec{x}
```
有非零解，即：
```math
(A-\lambda I)\vec{x}=\vec{0}
```
有非零解。这些非零解是**矩阵 $A$ 关于特征值 $\lambda$ 的所有特征向量**，加上零向量，就组成了：
* 矩阵 $(A-\lambda I)$ 的零空间
* 矩阵方程 $(A-\lambda I)\vec{x}=\vec{0}$ 的解空间

现在它们有一个新的定义：矩阵 $A$ 关于特征值 $\lambda$ 的 **特征（子）空间**。

### 求解特征值和特征向量

#### 特殊
上/下三角矩阵的特征值是对角线元素，特征向量是标准基。  
  

#### 通用
矩阵方程 $(A-\lambda I)\vec{x}=\vec{0}$ 有非零解，等价于 $det(A-\lambda I)=0$，由此得到两个概念：
* **特征方程** $det(A-\lambda I)=0$：一元n次方程，未知元是 $\lambda$，它在数域 $K$ 上的解，给出了矩阵 $A$ 的所有特征值。
* **特征多项式** $det(A-\lambda I)$：一元n阶多项式，未知元是 $\lambda$，在数域 $K$ 上对它作因式分解，所有**一次因式**给出了矩阵 $A$ 的所有特征值。对于特征值 $\lambda_{i}$，**一次因子** $(\lambda - \lambda_{i})$ 在**特征多项式**中的**重数**，被称为**特征值 $\lambda_{i}$ 的代数重数**。

求解特征方程（因式分解特征多项式），可以给出矩阵 $A$ 的所有特征值。  
对每一个特征值 $\lambda_{i}$：  
1. **一次因子** $(\lambda - \lambda_{i})$ 在**特征多项式**中的**重数**，被称为**特征值 $\lambda_{i}$ 的代数重数**。
2. 矩阵方程 $(A-\lambda_{i} I)\vec{x}=\vec{0}$ 的解空间，被称作**矩阵 $A$ 关于特征值 $\lambda_i$ 的特征（子）空间**。解空间除去零向量，就得到了矩阵 $A$ 关于特征值 $\lambda_i$ 的所有**特征向量**。解空间的**维度**，被称为**特征值 $\lambda_{i}$ 的几何重数**。  
  
至此，矩阵的特征值和特征向量的通用解法完毕。通过一些一元n次多项式的知识，可以得到如下结论：  
1. 在复数域上，一元n次多项式一定可以分解成一次因式，所以任何矩阵一定有复数域上的特征值。考虑重数的话，任何矩阵在复数域上都有n个特征值。
2. 在数域K上，考虑重数的话，一元n次多项式最多有n个解，所以任何矩阵在数域K上最多有n个特征值（考虑重数），那么不考虑重数的话，不同的特征值数量只会更少。  
  
总结一下特征值和特征向量的通用解法：
1. **特征值**直接求解于：**特征多项式**在数域上的线形因式分解。在这个过程中，有多少个不同特征值、不同特征值各自的**代数重数**都可以得出。
2. **特征向量**的求解在**特征值**之后，求解于：特征值代入后在数域K上的齐次线性方程组解集，即求解特征空间。这样，不同特征值各自的**几何重数**可以得出。




## 相似和对角化
在上一章Linear中已经讨论了矩阵**相似**的定义和意义。在这里重新复习一下：  
### 相似的定义
对于矩阵 $A,B \in \mathbb{R^{n,n}}$，如果存在**可逆矩阵** $P$，使得 $A = P\ B\ P^{-1}$ 成立，则称矩阵 $A$ 和 $B$ **相似**。  

### 相似矩阵的特征多项式/特征值/特征向量
结合相似和特征值/特征向量的定义，容易得到下面这些结论：
1. **相似的矩阵有相同的特征多项式**。
```math
det(A-\lambda I)=det(PBP^{-1}-\lambda I)=det(P(B-\lambda I)P^{-1})=det(P)det(B-\lambda I)det(P^{-1})=det(B-\lambda I)
```

2. **相似的矩阵有相同的特征值（以及对应的代数重数也相同）**。  
        显然，因为它们有相同的特征多项式，所以其因式分解也相同  

3. 相似的矩阵有相同的特征向量（坐标变化意义下的）。
```math
A\vec{x}=\lambda\vec{x} \leftrightarrow PBP^{-1}\vec{x}=\lambda\vec{x} \leftrightarrow BP^{-1}\vec{x}=\lambda P^{-1}\vec{x} \leftrightarrow 
B\left[\vec{x}\right]_\mathcal{B} = \lambda \left[\vec{x}\right]_\mathcal{B}
```
```math
basis\ \mathcal{B}\ is\ column\ vectors\ of\ P
```

### 对角化的定义
正如上一章Linear最后讨论地那样，如果矩阵 $A$ 有 n 个线性无关的特征向量，那么以这 n 个线性无关的特征向量作为 $\mathbb{R}^{n}$ 的**基**，矩阵 $A$ 的相似矩阵 $D$ 是一个**对角矩阵**，对角线上是 n 个特征向量对应的**特征值**，n个特征向量作为列向量组构成了矩阵 $P$，即：
```math
\begin{cases}
P = \left[\vec{e_1},\cdots,\vec{e_n}\right],\ \vec{e_i}\ \ are\ \ eigenvectors\ \ of\ \ A,\ i=1,\cdots,n\\
A = PDP^{-1},\ \ D = diag\left(\lambda_1,\cdots\lambda_n\right)
\end{cases}
```
此即为**对角化**。  
  
我们得到了一个对角化的充分条件，其实它也是对角化的一个必要条件。因为当 $A = PDP^{-1}$时，令 $P = \left[\vec{p_1},\cdots,\vec{p_n}\right]$。则有：
```math
A\vec{p}_i =
PDP^{-1}\vec{p}_i = 
\left[\vec{p_1},\cdots,\vec{p_n}\right]
\begin{bmatrix}
\lambda_1 & 0 & \cdots & 0\\
0 & \lambda_2 & \cdots & 0\\
\vdots\\
\vdots\\
0 & 0 & \cdots & \lambda_n\\
\end{bmatrix}
\ 
\begin{bmatrix}
0 \\
\vdots \\
1 \\
\vdots \\
0 \\
\end{bmatrix}
\begin{matrix}
\ \\
\ \\
i_{th}\\
\ \\
\ \\
\end{matrix}=
\lambda_i \vec{p}_i
```
那么可以看出， $\vec{p}_i$ 都是矩阵 $A$ 的特征向量。这里
```math
P^{-1}\vec{p}_i=\left[
\begin{matrix}
0\\
\vdots\\
1\\
\vdots\\
0
\end{matrix}
\right]
\begin{matrix}
\ \\
\ \\
i_{th}\\
\ \\
\ 
\end{matrix}
```
是因为考虑
```math
P^{-1}P = I \leftrightarrow P^{-1}\left[\vec{p_1},\cdots,\vec{p_n}\right]=I \leftrightarrow \left[P^{-1}\vec{p_1},\ P^{-1}\vec{p_2},\cdots P^{-1}\vec{p_n}\right] = I \leftrightarrow 
P^{-1}\vec{p}_i=\left[
\begin{matrix}
0\\
\vdots\\
1\\
\vdots\\
0
\end{matrix}
\right]
\begin{matrix}
\ \\
\ \\
i_{th}\\
\ \\
\ 
\end{matrix}
```
#### 矩阵 $A$ 可对角化的充要条件1（从对角化的定义而来）：矩阵 $A$ 有 $n$ 个线性无关的特征向量

### 可对角化的其他充要条件
矩阵可对角化的充要条件1缺乏操作性，仅仅是从对角化定义出发得到的。现在来思考下，如何找到矩阵 $A \in M_n(K)$ 的 $n$ 个线性无关的特征向量。
假如说，矩阵 $A$ 有 $m$ 个不同的特征值，分别是 $\lambda_1,\ \lambda_2,\ \cdots,\ \lambda_m$。  
对于每一个特征值 $\lambda_i$，矩阵 $A$ 关于 $\lambda_i$ 的几何重数是 $r_i$，也就是说，特征值 $\lambda_i$ 可以贡献出 $r_i$ 个线性无关的特征向量。那么首先考虑，把这些**属于不同特征值的线性无关的特征向量**集合，取并集后是否依然线性无关呢？答案是是的（证明放在最后补充）。即：  
  
**属于不同特征值的、线性无关的特征向量集合，它们的并集仍然是线性无关的**  
  
所以，从 $m$ 个不同的特征值中，最多可得到 $r_1+r_2+\cdots+r_m$ 个线性无关的特征向量，所以若这个值大于等于 $n$，那么矩阵 $A$ 有 $n$ 个线性无关的特征向量，可对角化；如果这个值小于 $n$，那么矩阵 $A$ 不可对角化。实际上下文可以知道，这个和值不会大于 $n$。  

#### 矩阵 $A$ 可对角化的充要条件2（从特征子空间的几何而来）：矩阵 $A$ 的不同特征值的几何重数之和等于 $n$  

这里可以立即给出一个矩阵 $A$ 可对角化的充分条件：**若 $A$ 有 $n$ 个不同的特征值，那么 $A$ 可对角化**  
原因很简单，每个不同的特征值，至少能给出一个特征向量。

  
矩阵可对角化的充要条件2已经具有很强的操作性了，它需要求出：
1. 矩阵 $A$ 的所有不同的特征值
2. 矩阵 $A$ 的所有不同特征值的几何重数

那么有时候会想要得到一个“对不同的特征值逐一检查”的方法，使得当矩阵 $A$ 不可对角化时，能提前发现。考虑到
1. **特征值 $\lambda_i$ 的几何重数 小于等于 代数重数**（证明放在最后补充）
2. **不同特征值的代数重数之和 小于等于 $n$**（特征值 $\lambda_i$ 的代数重数，是 $n$ 阶的特征多项式里，一次因式 $(\lambda - \lambda_i)$ 的重数，从而不同特征值的代数重数之和 小于等于 最高阶次 $n$）  
  
那么对于矩阵 $A \in M_n(K)$，可以得到如下的等价结论：
1. 矩阵 $A$ 可对角化
2. 矩阵 $A$ 的不同特征值的几何重数之和 = $n$
3. 矩阵 $A$ 的不同特征值的代数重数之和 = $n$，即矩阵 $A$ 的特征多项式可以在数域 $K$ 上完全因式分解为一次因式的乘积
4. 矩阵 $A$ 的不同特征值的几何重数之和 = 矩阵 $A$ 的不同特征值的代数重数之和 = $n$
5. 对于每一个不同的特征值 $\lambda_i$，矩阵 $A$ 关于 $\lambda_i$ 的几何重数 = 矩阵 $A$ 关于 $\lambda_i$ 的代数重数，且矩阵 $A$ 在数域 $K$ 上有 $n$ 个特征值（考虑重数）
6. 对于每一个不同的特征值 $\lambda_i$，矩阵 $A$ 关于 $\lambda_i$ 的几何重数 = 矩阵 $A$ 关于 $\lambda_i$ 的代数重数，且矩阵 $A$ 的特征多项式可以在数域 $K$ 上完全因式分解为一次因式的乘积。
7. 对于每一个不同的特征值 $\lambda_i$，矩阵 $A$ 关于 $\lambda_i$ 的几何重数 = 矩阵 $A$ 关于 $\lambda_i$ 的代数重数，且矩阵 $A$ 的特征多项式在复数域上的 $n$ 个根，都属于数域 $K$ 。  

#### 矩阵 $A$ 可对角化的充要条件3（从特征多项式的代数而来）：矩阵 $A$ 的每个特征值的几何重数 等于 它的代数重数，且矩阵 $A$ 的特征多项式可以在数域 $K$ 上可完全因式分解为一次因式的乘积

---

### 可对角化的充要条件总结  

矩阵 $A \in M_n(K)$ 可对角化的充要条件有如下三个：
1. 矩阵 $A$ 有 $n$ 个线性无关的特征向量  
2. 矩阵 $A$ 的不同特征值的几何重数之和等于 $n$  
3. 矩阵 $A$ 的每个特征值的几何重数 等于 它的代数重数，且矩阵 $A$ 在数域 $K$ 上有 $n$ 个特征值（考虑重数）

---

### 补充两个证明

#### 矩阵 $A$ 的、属于不同特征值的、线性无关的特征向量集合，它们的并集仍然是线性无关的

只需证明 两个不同特征值 的情况。多个不同特征值的情况用数学归纳法可马上得到。  
  
对于矩阵 $A$，它有两个不同的特征值 $\lambda_1$ 和 $\lambda_2$，因为不能同时为 $0$，考虑 $\lambda_1 \neq 0$ 。特征向量集合 $set\left(\vec{a}_1,\ \vec{a}_2,\ \cdots \vec{a}_p\right)$ 和 $set\left(\vec{b}_1,\ \vec{b}_2,\ \cdots \vec{b}_q\right)$ 分别是属于特征值 $\lambda_1$ 和 $\lambda_2$ 的线性无关的特征向量组成的集合。那么考虑等式
```math
k_1\vec{a}_1+k_2\vec{a}_2+\cdots+k_p\vec{a}_p+
s_1\vec{b}_1+s_2\vec{b}_2+\cdots+s_q\vec{b}_q=\vec{0}
\tag0
```
等价于
```math
k_1 \lambda_1 \vec{a}_1+k_2 \lambda_1 \vec{a}_2+\cdots+k_p \lambda_1 \vec{a}_p+
s_1 \lambda_1 \vec{b}_1+s_2 \lambda_1 \vec{b}_2+\cdots+s_q \lambda_1 \vec{b}_q=\vec{0}
\tag1
```
0式 两边左乘 矩阵 $A$，得到
```math
k_1A\vec{a}_1+k_2A\vec{a}_2+\cdots+k_pA\vec{a}_p+
s_1A\vec{b}_1+s_2A\vec{b}_2+\cdots+s_qA\vec{b}_q=\vec{0}
```
即得
```math
k_1 \lambda_1 \vec{a}_1+k_2 \lambda_1 \vec{a}_2+\cdots+k_p \lambda_1 \vec{a}_p+
s_1 \lambda_2 \vec{b}_1+s_2 \lambda_2 \vec{b}_2+\cdots+s_q \lambda_2 \vec{b}_q=\vec{0}
\tag2
```
2式 - 1式，再根据 $\lambda_1 \neq \lambda_2$，可以马上得出 $s_1=s_2=\cdots =s_q=0$，然后可以马上得出 $k_1=k_2=\cdots =k_p=0$，于是 $set\left(\vec{a}_1,\ \vec{a}_2,\ \cdots \vec{a}_p,\ \vec{b}_1,\ \vec{b}_2,\ \cdots \vec{b}_q\right)$ 是线性无关的。QED.  

#### 矩阵 $A$ 的每个特征值 $\lambda_i$ ，它的几何重数 小于等于 代数重数

对于矩阵 $A$ 和它的一个特征值 $\lambda_i$，其几何重数为 $d$。那也就是说，矩阵 $A$ 关于特征值 $\lambda_i$ 的特征子空间的维数是 $d$。从该特征子空间中取 $d$ 个线性无关的向量 $set\left(\vec{e}_1,\vec{e}_2,\cdots \vec{e}_d\right)$，然后将它扩充为 $\mathbb{R}^d$ 的一个基
```math
\mathcal{B} = set\left(\vec{e}_1,\vec{e}_2,\cdots \vec{e}_d, \vec{b}_1,\vec{b}_2,\cdots \vec{b}_{n-d}\right)
```
那么坐标转移矩阵
```math
P_\mathcal{E \leftarrow B} = \left[\vec{e}_1,\vec{e}_2,\cdots \vec{e}_d, \vec{b}_1,\vec{b}_2,\cdots \vec{b}_{n-d}\right]
```
简写为 $P$。  
  
上述工作的目的是求出 矩阵 $A$ 以 $P$ 为坐标转移矩阵的相似矩阵 $X$，即求出矩阵 $X$ 满足 $A=PXP^{-1}$：  
```math
A = PXP^{-1}
```
等价于
```math
AP = PX
```
即
```math
A\left[\vec{e}_1,\vec{e}_2,\cdots \vec{e}_d, \vec{b}_1,\vec{b}_2,\cdots \vec{b}_{n-d}\right] = PX
```
即
```math
\left[A\vec{e}_1,A\vec{e}_2,\cdots A\vec{e}_d, A\vec{b}_1,A\vec{b}_2,\cdots A\vec{b}_{n-d}\right] = PX
```
即
```math
\left[\lambda_i\vec{e}_1,\lambda_i\vec{e}_2,\cdots \lambda_i\vec{e}_d, A\vec{b}_1,A\vec{b}_2,\cdots A\vec{b}_{n-d}\right] = PX
```
把未知矩阵 $X$ 写成列向量组的形式，即：
```math
X=
\left[\vec{x}_1, \vec{x}_2, \cdots \vec{x}_d, \vec{x}_{d+1},\cdots \vec{x}_n
\right]
```
代入得到：
```math
\left[\lambda_i\vec{e}_1,\lambda_i\vec{e}_2,\cdots \lambda_i\vec{e}_d, A\vec{b}_1,A\vec{b}_2,\cdots A\vec{b}_{n-d}\right]
=
\left[
P\vec{x}_1,\ 
P\vec{x}_2,\ 
\cdots
P\vec{x}_d,\ 
\cdots 
P\vec{x}_n
\right]
```
考虑前 $d$ 列相等，得到：
```math
\lambda_i\vec{e}_j
=
P\vec{x}_j\ \ ,
j=1,2,\cdots d
```
即
```math
\lambda_i\vec{e}_j
=
\left[\vec{e}_1,\vec{e}_2,\cdots \vec{e}_d, \vec{b}_1,\vec{b}_2,\cdots \vec{b}_{n-d}\right]\vec{x}_j\ \ ,
j=1,2,\cdots d
\tag{*}
```
方程（*）只有唯一一个解，即:
```math
\vec{x}_j=
\lambda_i
\left[
\begin{matrix}
0 \\
\vdots\\
1 \\
\vdots\\
0\\
\end{matrix}
\right]
\begin{matrix}
\ \\
\ \\
j_{th}\ row \\
\ \\
\ \\
\end{matrix}
\ \ \ ,j = 1,2,\cdots d
```
这是因为基
```math
set\left(\vec{e}_1,\vec{e}_2,\cdots \vec{e}_d, \vec{b}_1,\vec{b}_2,\cdots \vec{b}_{n-d}\right)
```
线性无关，所以只有一种线性表出 $\lambda_i\vec{e}_j$ 的系数权重。  
这样，就得到了矩阵 $X$ 的部分（前 $d$ 列）表示，即：
```math
X = 
\left[
\left[
\begin{matrix}
\lambda_i \\
0 \\
\vec{0}\\
0 \\
0 \\
\vec{0}\\
0
\end{matrix}
\right],
\left[
\begin{matrix}
0 \\
\lambda_i \\
0 \\
\vec{0}\\
0 \\
\vec{0}\\
0
\end{matrix}
\right],
\cdots,
\left[
\begin{matrix}
0 \\
\vec{0}\\
0 \\
\lambda_i \\
0 \\
\vec{0}\\
0
\end{matrix}
\right]
\begin{matrix}
\ \\
\ \\
\ \\
d_{th} \\
\ \\
\ \\
\end{matrix}
\ \ \ ,
\vec{x}_{d+1},
\cdots
\vec{x}_{n}
\right]
```
整理一下，即：
```math
A = P
\left[
\begin{matrix}
\lambda_i I_d & B \\
0 & C
\end{matrix}
\right]
P^{-1}
```
因为**相似的矩阵有相同的特征多项式**，所以得到：
```math
det(A-\lambda I) = 
det(
\left[
\begin{matrix}
\lambda_i I_d & B \\
0 & C
\end{matrix}
\right]-\lambda I
)=
det(
\left[
\begin{matrix}
(\lambda_i-\lambda) I_d & B \\
0 & C-\lambda I_{n-d}
\end{matrix}
\right]
)=
det(
(\lambda_i-\lambda) I_d
)\ 
det(
C-\lambda I_{n-d}
)=
(\lambda_i-\lambda)^{d}\ det(
C-\lambda I_{n-d}
)
```
所以 $(\lambda_i-\lambda)^{d}$ 是 特征多项式 $det(A-\lambda I)$ 的子式，所以 $\lambda$ 的几何重数 $d \le r$，这里 $r$ 是 $\lambda$ 的代数重数，即一次因式 $(\lambda_i-\lambda)$ 在特征多项式 $det(A-\lambda I)$ 里的最高次幂。QED.  
  
从这个证明也可以看出，矩阵 $A$ 知道多少个特征向量，就能确定多少个对角线元素。容易证明：  
对于一个不能彻底对角化的矩阵 $A$ 来说，如果能确定 $A$ 的 $m$ 个不同的特征值 $\lambda_1, \lambda_2,\cdots \lambda_m$，对应的几何重数 $r_1,r_2,\cdots r_m$，有 $r = n - (r_1+r_2+\cdots+r_m) > 0$，那么 $A$ 的**相似约当标准型**是：
```math
\begin{bmatrix}
\lambda_1 I_{r_1} & 0 & \cdots & 0 & B_1\\
0 & \lambda_2 I_{r_2} & \cdots & 0 & B_2\\
\vdots & \ddots\\
0 & 0 & \cdots & \lambda_m I_{r_m} & B_m\\
\vec{0} & \vec{0} & \cdots & \vec{0} & C\\
\end{bmatrix}
```
矩阵 $C$ 是一个 $r$ 行 $r$ 列 的矩阵，它的行列式就是 $A$ 的特征多项式因式中不能化成一次因式的部分。