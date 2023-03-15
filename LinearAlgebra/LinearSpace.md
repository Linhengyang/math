# 综述: 线性方程组/线性空间/线性映射
贯穿全章节的计算手段主要是:
* 初等行变换(elementary row operation)
* 矩阵向量乘积(matrix-vector product)
* 矩阵乘法(matrix-matrix product)
  
这块内容从**线性方程组**讲起, 判断是否**有解/解集结构**, 到**线性组合(相关/无关)**理解**矩阵**, 然后两个分支
* 分支1: 从**线性映射**重新理解**矩阵**, 包括**矩阵乘法(复合映射)**和**可逆矩阵**. 这里重新用**初等矩阵左乘**重新理解了初等行变换, 并引出**LU分解**.
* 分支2: 从**线性子空间**理解矩阵的**列空间**/**零空间**并描述, 然后引出**基**/**坐标**/**维度**的概念, 并推出矩阵的**秩**.

## 分支1综述
分支1脉络中, 主线是矩阵方程 $A\vec{x}=\vec{b}$ 的求解，其中$A\in\mathbb{R^{m,n}}$。
记系数矩阵 $A$ 的echelon form(阶梯型)是 $U$。  
求解过程中, 初等行变换（左乘初等矩阵）的操作，将 $A$ 转化 $U$，甚至 $A$ 的reduced echelon form(简约阶梯型)形式。  
初等行变换非常重要，在这个操作变换下，矩阵有两个不变量：
1. 矩阵的秩 $rankA$
2. 矩阵的列向量组的序

用初等矩阵 $E_i$ 左乘来理解初等行变换, 可以得到：$$E_p\,\cdots\,E_1\,A=U$$
即 $A=({E_p\,\cdots\,E_1})^{-1}U$。记 $P={E_p\,\cdots\,E_1}$，则 $A=P^{-1}\,U$。

对于线性方程组 $$A\vec{x}=\vec{b}$$ 来说, 它等价于 $$P\,A\vec{x}=P\,\vec{b}$$, 即 $$U\vec{x}=P\,\vec{b}$$  

---
第一步：判读是否有解  
对比 $U$ 的非零行数量 $k$ 和 $P\,\vec{b}$ 的非零行数量
* 当 $k$ 小于 $P\,\vec{b}$ 的非零行数量时，无解
* 当 $k$ 大于等于 $P\,\vec{b}$ 的非零行数量时，有解

第二步：解集结构  
对比 $U$ 的非零行数量 $k$ 和 $n$
* 当 $k=n$ 时，说明矩阵 $A$ 列向量组 $\mathcal{A}$ 线性无关，它作为 $\mathbb{R^n}$ 的一个基（basis），对于任何 $\vec{b} \in \mathbb{R^{n}}$，在此 basis 下有唯一表示，即 $A\,\vec{x}=\vec{b}$ 只有唯一解，即 $\vec{b}$ 在 $A$ 的列向量组 $\mathcal{A}$ 下的坐标 $\displaystyle \left[\vec{b}\right]_\mathcal{A}\in\mathbb{R^n}$。从 $U\,\vec{x}=P\,\vec{b}$ 也可以看出，$U$ 的每一列都是主元所在列，所以 $\vec{x}$ 只有唯一解。
* 当 $k \lt n$ 时，因为矩阵的行秩等于列秩，所以 $U$ 的主元所在列数量也是 $k$。$U$ 的主元所在列组成了 $A$ 的列向量组 $\mathcal{A}$ 的一个basis，记为 $\mathcal{B}$。$A$ 的列向量组 $\mathcal{A}$ 在自由元所在列可以构成一个basis $\mathcal{B}$，$span \left\{\mathcal{A}\right\}=span \left\{\mathcal{B}\right\}$，是 $\mathbb{R^m}$ 的一个 $dim = k$ 的子空间