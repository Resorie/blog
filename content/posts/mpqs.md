+++
title = "冒泡潜水"
date = 2025-12-10
description = "二次筛法 Quadratic Sieve"

[taxonomies]
tags = ["数论", "密码学", "算法"]

[extra]
feature_image = "https://images.unsplash.com/photo-1764591696226-ea4e8d655bc7?crop=entropy&cs=tinysrgb&fit=max&fm=jpg&ixid=M3w4NDIyMzZ8MHwxfHJhbmRvbXx8fHx8fHx8fDE3NjU0MzMzMzl8&ixlib=rb-4.1.0&q=80&w=1080"
pinned = true
+++

填一个[坑](https://www.luogu.com.cn/article/vy3tfjw4)。这个鸽了三年多了。

## Preface

众所周知，Pollard-Rho 能以 `$O(n^{1/4}\log n)$` 的时间复杂度（通过一定优化能做到 `$O(n^{1/4})$`）找到合数 `$n$` 的一个非平凡因子。但事实上，我们有更快的方法做到这件事。

本篇文章将介绍一种亚指数级复杂度分解质因数的算法：二次筛法（Quadratic Sieve, QS）。

注：

- 本篇文章中的复杂度分析将不包含对大整数进行运算的开销。
- 本篇文章偏教程向，会混杂很多比较 trivial 的内容。
- 其实我自己写的代码基本是贺的 zzt 的提交并加以修改。

## Smooth Numbers

数论中，一般认为素数的性质是比较好的。比较自然地，研究合数时，我们会希望研究那些“比较素”的合数。

对这样比较好的合数的刻画，容易会想到，使用素因子较少的合数。但这显然是荒谬的，因为质因数分解并不容易：对于大素数 `$p,q$`，合数 `$n=pq$` 很难被逆向表示为 `$p,q$`。

因而，我们转向另一个方向：研究*素因子较小*的合数，这便引出了 **smooth numbers**。

> 称一个数是 `$y\text{-smooth}$` 的，当且仅当其所有质因子不超过 `$n$`。

容易发现，smooth number 的一个很好的性质是，它们容易被分解。我们只需简单筛出不超过阈值 `$y$` 的素数 `$p_1,p_2,\dots,p_k$`，便能通过简单的试除法（因为 `$y$` 相对较小，这些开销是可以承受的）确定 `$n$` 的分解。

于是，我们会萌生一个想法：对于不 smooth 的 `$n$`，我们能否通过某些方法将其转化为比较 smooth 的数，再对其进行操作？想必，转化之后的操作应当是简单的。

事实上，数论中许多算法，如 Index Calculus 和本篇文章所要介绍的 Quadratic Sieve 都基于这样的转化。

### Smooth Number 的密度

Smooth number 的密度由 Dickman 函数给出：

```math-display
$$
\psi(x,y):=\left|\left\{i\le x\mid i \text{ is } y\text{-smooth}\right\}\right|
$$
```

当 `$x\to\infty$` 时，有：

```math-display
$$
\psi(x,y)\sim x \rho(u)\quad\text{where }x=y^u
$$
```

其中 `$\rho(u)$` 是 Dickman-de Brujin 函数。这表明，当 `$x$` 足够大时，smooth number 的密度不依赖于具体的 `$x,y$`，而只和 `$u=\frac{\ln x}{\ln y}$` 有关。

Dickman 证明了如下估计：

```math-display
$$
\rho(u)\approx u^{-u+o(u)}
$$
```

证明涉及微分学，这里不作详细展开。

## Quadratic Sieve

### Dixon's Random Square Method

我们的算法基于以下思想：若能找到 `$\alpha^2\equiv \beta^2\pmod N$` 且 `$\alpha\not\equiv\pm \beta$`，那么 `$\gcd(\alpha+\beta,N)$` 能够给出一个 `$N$` 的非平凡因子。问题在于，我们如何生成这样的 `$\alpha, \beta$`。

一个比较简单的想法是随机。我们随机这样的一个 `$\alpha_i$` 并计算 `$\gamma_i:=\alpha_i^2\bmod N$`，这是一个模 `$N$` 意义下的二次剩余。我们希望找到 `$\gamma_i$` 的一个非平凡平方根 `$\beta_i$`。假定 `$\gamma_i\perp N$`，那么 `$N$` 的每个质因子对应了 `$x^2\equiv \gamma_i$` 的两个根，因而找到的 `$\beta_i\equiv\pm \alpha_i$` 的概率只有 `$\frac{1}{2^{\omega(N)-1}}\le \frac 12$`。

理想情况下，若 `$\gamma_i$` 是比较 smooth 的，我们可以对其分解 `$\gamma_i=\prod p_j^{e_{i,j}}$`。若 `$\forall j$` 有 `$2\mid e_{i,j}$`，那么我们可以轻易地令 `$\beta_i=\prod p_j^{e_{i,j}/2}$`，这便得到一组 `$\alpha,\beta$`。

当然，我们运气应当不致于这么好，不仅 `$\gamma_i$` 是 smooth 的，而且 `$e_{i,j}$` 都是偶数。前者我们只能寄希望于概率，但是后者我们可以想想办法。考虑选取多个 `$\alpha_i$` 并计算 `$\gamma_i$`，再令向量 `$\boldsymbol v_i=(e_{i,1}\bmod 2,e_{i,2}\bmod 2,\dots,e_{i,k}\bmod 2)\in \mathbb Z_2^n$`。如果我们能找到一些线性相关的 `$\boldsymbol v_i$`，即存在指标集 `$I$` 使得 `$\sum_{i\in I}\boldsymbol v_i=\boldsymbol 0$`，那么可以取

```math-display
$$
\begin{aligned}
\alpha:=&\prod_{i\in I}\alpha_i\\\
\gamma:=& \prod_{i\in I}\gamma_i=\prod p_j^{\sum_{i\in I} e_{i,j}}\\
\beta:=&\prod p_j^{\left(\sum_{i\in I}e_{i,j}\right)/2}
\end{aligned}
$$
```


一般的实现中，可以令 `$\beta_i:=\prod p_j^{\left\lfloor e_{i,j} / 2\right\rfloor}$`，那么就有：

```math-display
$$
\beta=\prod_{i\in I} \beta_i\cdot\prod p_j^{\left(\sum_{i\in I}e_{i,j}\bmod 2\right) / 2}
$$
```

假设选取的 smooth 阈值为 `$B$`，那么我们就是有若干个 `$\pi(B)$` 维的向量。只要我们找到 `$\pi(B)+k$` 个向量，我们就能通过高斯消元找出 `$k$` 组 `$\alpha,\beta$`，每组至少有 `$\frac 12$` 概率成功分解，那么总的分解成功概率就是 `$1-2^{-k}$`。

#### 复杂度分析

~~我被一篇傻逼文献骗了好久~~

显然该方法的复杂度依赖于阈值 `$B$` 的选取。为方便计算，令 `$k = \pi(B),u=\frac{\ln N}{\ln B}$`。

`$\gamma_i=\alpha_i^2\bmod N$` 是模 `$N$` 的二次剩余，它大约在 `$1\sim N$` 中均匀分布，因此我们可以用前面的 smooth number 密度公式估算它 smooth 的概率，即 `$u^{-u}$`。我们需要找到 `$O(k)=O\left(\frac{B}{\ln B}\right)$` 个向量。每次尝试需要 `$k$` 次试除。最后，我们需要进行高斯消元。故可得到：

```math-display
$$
T\approx k^2u^u+k^3
$$
```

这种式子谁爱算谁算，我扔给 AI 了。算出来 `$T$` 最小时有：

```math-display
$$
B\approx \exp(\frac 1{\sqrt2}\sqrt{\ln N\cdot\ln\ln N}) \\
T\approx\exp(\frac 3{\sqrt2}\sqrt{\ln N\cdot\ln\ln N})
$$
```


#### The L-notation

细心的读者容易发现上面的 `$B$` 和 `$T$` 有相似的结构。事实上，这种结构的渐进表示在 smooth number 相关算法中相当常见。

一般地，我们定义 **L-notation**：

```math-display
$$
L_n[\alpha, c]:=\exp\left(\left(c+o(1)\right)\left(\ln n\right)^\alpha\left(\ln \ln n\right)^{1-\alpha}\right)\quad \text{where }c>0, 0\le\alpha\le1
$$
```

该函数的增长量级主要由 `$\alpha$` 控制：

- `$\alpha=1$` 时，`$L_n[1,c]=O\left(n^{c+o(1)}\right)$`（指数级复杂度）
- `$\alpha=0$` 时，`$L_n[0,c]=O\left((\ln n)^{c+o(1)}\right)=O(\text{polylog})$`
- `$\alpha\in(0,1)$` 时，`$L_n[\alpha,c]$` 表示亚指数级复杂度。

于是我们可以重新表示上面结果中的 `$B,T$`：

```math-display
$$
B=L_N\left[\frac 12,\frac1{\sqrt 2}\right]\\
T=L_N\left[\frac 12,\frac3{\sqrt2}\right]
$$
```

读者可以自行尝试证明该算法复杂度的另一种方法：设 `$B=L_N[1/2,c]$` 并带入以上复杂度式子，解出 `$c$`。

#### Optimized Dixon

朴素的 Dixon 方法使用了随机选取的 `$\alpha$`，因而 `$\gamma:=\alpha^2\bmod N$` 的 smoothness 似乎没有保证。我们会希望使用更好的 `$\alpha$` 来优化这一点。

一个比较显然的想法是选择 `$\alpha_i=\left\lceil \sqrt N\right\rceil+i$`。当 `$i$` 不太大时，这样得到的 `$\gamma_i\sim2i\sqrt N$` 比较小，应当会比较 smooth。

假设我们需要 `$M$` 个这样的 `$\alpha_i$`，那么应当有 `$\gamma_i \sim 2M\sqrt N$`，以这个数为上限套用之前的 smooth number 密度公式，计算得到，当 `$B=L_N[1/2,1/2]$` 时，有 `$T\approx L_N[1/2,3/2]$`。

[一个未作太多优化的实现](https://loj.ac/s/2472888)

### QS

以上的算法在 `$N\approx 10^{30}$` 时就略显吃力。我们当然不满足于此。

注意到上面算法的复杂度瓶颈在于，对于每个 `$\alpha_i$`，我们都必须通过试除法确定 `$\gamma_i$` 是否是 smooth 的。由于选取合适的 `$\alpha_i$` 需要大量的尝试，这样反复试除的开销会变得不可接受。考虑如何减少试除次数。

我们都知道求质数可以用筛法，其实求 smooth 数也可以。更进一步，我们还可以筛多项式（因为 `$p\mid f(c)\Leftrightarrow p\mid f(kp+c)$`），比如一般的埃氏筛就是在筛 `$f(i)=i$`（如果跳过偶数可看作在筛 `$f(i)=2i+1$`）。

考察这样一个多项式：`$f(x)=x^2-N$`。当 `$i$` 较小时，有 `$f\left(\left\lceil\sqrt N\right\rceil+i\right)=\alpha_i^2\bmod N$`。我们希望筛一下它来优化上面的 Optimized Dixon 中的试除法。

类似埃氏筛地，我们枚举小质数 `$p$`。若 `$p^k\mid f(i)$` 就将 `$p^k$` 从 `$f(i)$` 的值中除去，筛完一遍变成 `$1$` 的就是 smooth 的 `$f(i)$`。若我们有 `$c$` 是 `$N$` 在模 `$p$` 意义下的二次剩余，那么当且仅当 `$x\equiv \pm c\pmod p$` 时会有 `$p\mid f(x)$`，枚举每个 `$x$` 并筛去 `$p$` 即可。

[一个实现](https://loj.ac/s/2473437)

#### 优化

对 QS 的改进方案中，提升较大的是选择更好的多项式去筛，这点留待后文详叙。这里讲几个小优化。

1. `$f(i)$` 和 `$f(-i)$` 可以一起筛。注意消元的时候要把负号消掉。
2. 注意到 `$p\mid f(i)$` 当且仅当 `$N$` 是模 `$p$` 意义下的二次剩余，因此在筛小素数时可以算一下 Legendre 符号 `$\left(\frac  Np \right)$`，如果不为 `$1$` 这个素数就是没用的。
3. 如果 `$f(i)$` 除掉小质数还剩一个不太大的数（一般的阈值是 `$B^2$`），可以把剩下这个值扔到 `std::unordered_map` 里。如果后面找到的 `$f(j)$` 剩下同样的数，那么可以使用 `$f(i)f(j)$`，这个剩下的值就被平方了，乘到 `$\beta$` 里即可。（Lenstra 等提出，数据范围足够大时，可以尝试将除剩下的数分解为两个不大于 `$B^2$` 的数的乘积，并放到图上求解[^1]）
4. 为了节省空间，我们在筛 `$f$` 时通常无法直接算出 `$\beta_i$` 和 `$\boldsymbol v_i$`，而需要对筛出的 smooth 的 `$f$` 再进行一次试除。这就使得筛 `$f$` 时的一堆大整数除法有点意义不明。可以考虑使筛 `$f$` 的条件变松：我们计算其对数，如果 $\log f(i) $ 减去 $ \sum_{p_j\mid f(i)}\log p_j$ 比较小，那么这个 `$f(i)$` 大概会比较 smooth，再对其进行试除即可。
5. 朴素的高斯消元会成为算法的复杂度瓶颈。注意到 `$\boldsymbol v_i$` 中至多有 `$O(\log N)$` 个数非零，这说明我们最终要进行高斯消元的矩阵是一个稀疏矩阵。可以使用 Block Wiedemann 或 Block Lanczos 算法（Montgomery 给出过一种为素因数分解特化的在 `$GF(2)$` 上的 Block Lanczos[^2]）将消元部分复杂度降至 `$O(k^2\log N)$`。实现中可以使用 `std::bitset` 优化的在线高斯消元，性能不会有明显劣势。

#### 复杂度分析

同 Dixon，我们需要筛出 `$O(k)$` 个数，每个数成功概率为 `$u^{-u}$`，筛每个数的成本为 `$\sum_{p\le B}p^{-1}=O(\log\log B)$`。高斯消元复杂度 `$O(k^2\log N)$`。故有：

```math-display
$$
T\approx ku^u\log\log B+k^2\log N
$$
```

当 `$T$` 最小时，有：

```math-display
$$
B\approx \exp(\frac 12\sqrt{\log N\log\log N})=L_N\left[\frac 12,\frac12\right]\\
T\approx\exp(\sqrt{\log N\log \log N})= L_N\left[\frac 12,1\right]
$$
```


### MPQS

Multiple Polynomial Quadratic Sieve，简称 MPQS ~~（标题冒泡潜水的由来）~~。

朴素 QS 的问题是，当 `$N$` 较大时，我们需要筛的 `$f(i)$` 个数（设为 `$M$`）会较大，而 `$f(i)\sim 2M\sqrt N$` 也会较大，这使得 `$f$` 较 smooth 的猜想不再成立。为了解决这个问题，Montgomery 提出：我们可以筛多个多项式。

具体而言，我们希望用一个线性函数 `$ax+b$` 替代之前的 `$x$`，即考虑 `$(ax+b)^2-N$`。如果 `$b^2\equiv N\pmod a$`，那么 `$(ax+b)^2-n$` 可以写成 `$af(x)=a(ax^2+2bx+c)$` 的形式，其中 `$b^2-ac=N$`。为了方便处理，我们设 `$a$` 是一个 smooth 数的平方（比如 `$p^2$`，也有文献建议取一系列小素数乘积的平方），那么我们只需要去筛 smooth 的 `$f(x)$`。

假设我们只去筛 `$|x|\le M$` 的 `$x$`。与之前同样地，我们希望使 `$|f(x)|$` 较小以使得其较 smooth。根据中学数学知识，我们会希望抛物线对称轴 `$x=-\frac ba$` 靠近所筛区间 `$[-M,M]$` 的中点 `$x=0$`，这要求 `$b$` 较小。另外，我们希望区间端点和中点处 `$|f(x)|$` 值相近，即 `$|f(M)|\approx |f(0)|$`，解出 `$a\approx\frac{\sqrt{2N}}{M}$`。

实现中，我们可以从 `$\frac{(2N)^{1/4}}{M^{1/2}}$` 开始枚举小素数 `$p$`，令 `$a=p^2$`，取 `$b$` 为 `$b^2\equiv N\pmod a$` 的较小解，再算出 `$c$`。然后按照 QS 的方法筛 `$(ax+b)^2-N$` 即可。

MPQS 相比于 QS 复杂度不变，但有较大常数优势。zzq 说有分析认为提升倍数是 `$\frac 12\sqrt{\ln N\ln\ln N}$`，我不知道是怎么算出来的。

[zzt 的实现](https://loj.ac/s/316691) 

### HMPQS

The Hypercube variation of Multiple Polynomial Quadratic Sieve，简称 HMPQS。

Peralta 提出，MPQS 中切换多项式还要重算一遍 `$b,c$` 之类的，比较麻烦。考虑生成一些有联系的多项式来消灭掉这个代价。~~其实和 Hypercube 没啥必然联系，Peralta 为了容易理解随便举的例子~~

我们发扬光大一下 MPQS 中的一个细节。令 `$t$` 为 `$n$` 个素数的乘积（`$t:=\prod_{i=1}^nq_i$`），且 `$t\sim N^{1/4}$`。假设 `$a=t^2$`，`$b$` 为方程 `$b^2\equiv N\pmod {t^2}$` 的一个解。显然这里 `$b$` 有 `$2^n$` 个解。我们希望这样所生成出的 `$2^n$` 个 `$(ax+b)^2$` 有比较好的性质。

首先考虑如何表示 `$b$`。由 CRT，只需解出 `$\alpha_i^2\equiv N\pmod {q_i^2}$`，那么就可以令 `$b=\sum \delta_i\alpha_i\beta_i$`，其中 `$\delta_i=\pm1$`，`$\beta_i$` 是 CRT 算出的系数。`$\alpha_i$` 有两个解，此处我们选择使得 `$\gamma_i:=\alpha_i\beta_i\bmod t^2\le\frac{t^2}2$` 的解。

然后考虑如何遍历所有的 `$b$`。令数列 `$s_1,\dots,s_{2^n}$` 遍历所有的 `$b$` 解。我们希望遍历过程中每次只更改一个 `$\delta_i$`，那么可以有递推关系：

```math-display
$$
s_{i+1}=s_i+2\mu_i\gamma_{k_i}-\omega_it^2
$$
```

其中 `$\mu_i=\pm 1$`，代表 `$\delta_{k_i}$` 的改变。Peralta 的叙述中把这个递推关系理解为 `$n$` 维 Hypercube 上的一个 Hamilton 路径，`$\mu_i$` 代表第 `$k_i$` 维坐标的改变；当然其实还有其他理解方法，比如相信大家都知道格雷码，原理是类似的。`$\omega_i$` 可以理解为对 `$t^2$` 取模的一个修正系数，使得 `$s_{i+1}\in(0,t^2)$`。由之前对 `$\alpha_i$` 的选择，有 `$\omega_i\in\{\pm 1,0\}$`。一般的实现中可以直接取 `$\omega_i\equiv 0$`，此时仍有 `$s_i\in(-nt^2,nt^2)$`。

当然我们并不关心具体值。在筛法中，`$s_i$` 作为多项式的系数，我们关心其模 `$p$` 意义下的值。Peralta 提出，为了方便计算 `$s_i\bmod p$`，我们可以打表一下其差分。令：

```math-display
$$
\Delta_i:=s_{i+1}-s_i=2\mu_i\gamma_{k_i}-\omega_it^2
$$
```

注意到 `$\mu_i$` 有两种取值，`$\omega_i$` 有三种取值，`$\gamma_i$` 有 `$n$` 种取值，那么 `$\Delta_i$` 至多有 `$6n$` 种取值（如果固定 `$\omega_i=0$` 则只有 `$2n$` 种）。这个可以打表算出来。那么转移就只需极小的代价。

现在可以构造我们需要的多项式了。对 MPQS 中的多项式略作变换，除掉 `$a$` 这个平方因子，就得到我们构造的多项式：

```math-display
$$
\begin{aligned}
f_i(x):=&\ (s_i/t+x t)^2\bmod N\\
=&\ \frac{s_i^2-N}{t^2}+2s_ix+x^2t^2
\end{aligned}
$$
```

我们仍须验证这样生成的 `$f_i(x)$` 不太大，以保证其 smooth 的概率。假设我们要筛的区间是 `$x\in[-M,M]$`，我们取 `$t\le\frac{N^{1/4}}{\sqrt M}$`。由于 `$s_i^2<t^4<N$`，通过简单放缩可以得到 `$\frac{s_i^2-N}{t^2}<0$` 且 `$\left|\frac{s_i^2-N}{t^2}\right|\le M\sqrt N$`。二次项 `$x^2t^2\le M^2\left(\frac{N^{1/4}}{\sqrt M}\right)^2=M\sqrt N$`。一次项 `$2s_ix\le 2M\sqrt N$`。因此我们大概有 `$f_i(x)=O(M\sqrt N)$`，这说明其值不太大，也就保证其大概比较 smooth。

在筛法中，令 `$p\mid f_i(x)$`，容易得到 `$x\equiv(-s_i\pm\sqrt N)t^{-2}\pmod p$`。这就得到筛 `$f_i(x)$` 的两个起点：

```math-display
$$
D_i=(-s_i+\sqrt N)t^{-2}\bmod p\\
E_i=(-s_i-\sqrt N)t^{-2}\bmod p=(D_i-2\sqrt Nt^{-2})\bmod p
$$
```

`$s_i\bmod p$` 容易由上述递推算出。如果预先算出 `$\sqrt N\bmod p, t^{-2}\bmod p,2\sqrt Nt^{-2}\bmod p$`，就可以很快算出 `$D_i,E_i$`。这就节省了切换多项式的时间。

还剩下最后一个细节：如果 `$p\mid t$`，那么以上同余式的推导不成立。事实上，当 `$p\mid t$` 时，当且仅当 `$x\equiv -\frac{(s_i^2-N)/t^2}{2s_i}\pmod p$` 时有 `$p\mid f_i(x)$`，即某个 `$x$` 使得 `$p\mid f_i(x)$` 的概率减半为 `$1/p$`。这提示我们，`$t$` 的因数尽量不要选取 smooth check 所使用的小素数，以增大 `$f_i(x)$` smooth 的概率。实现中，我们可以选择素数 `$q\in[B,2B]$` 来生成 `$t$`，这样大约有 `$n\approx\log_B N\approx2\sqrt\frac{\log N}{\log\log N}$`。

~~我并没有写这个的代码~~

-----

References:

- [zzq's blog](https://zhuanlan.zhihu.com/p/106650020)

- [Quadratic Sieve - Wikipedia](https://en.wikipedia.org/wiki/Quadratic_sieve)

- Prime Numbers: A Computational Perspective (Crandall, Pomerance)（zzq 那篇的 reference）

- A Quadratic Sieve on the n-Dimensional Cube (Peralta)

[^1]: Factoring with Two Large Primes (Lenstra, Manasse)

[^2]: A Block Lanczos Algorithm for Finding Dependencies over GF(2) (Montgomery)

