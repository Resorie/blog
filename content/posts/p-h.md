+++
title = "地狱人"
date = 2026-04-23
description = "Pohlig-Hellman 算法"

[taxonomies]
tags = ["数论", "离散对数", "算法"]

[extra]
+++

求解 `$g^x\equiv a\pmod P$`，即 `$\operatorname{ind}_ga\equiv x\pmod{(P-1)}$`。要求 `$P-1$` 较 smooth。

考虑对模数动手。首先我们要明确离散对数在其他模数下的定义。假设模数 `$m\mid(P-1)$`，`$\operatorname{ind}_ga\equiv x_0\pmod m$`。带入原始的指数方程：

<div class="math-display">$$
g^{km+x_0}\equiv a\pmod P
$$</div>

为分离出 `$x_0$`，我们要消灭掉指数上的 `$m$`。由 Fermat 小定理，有 `$g^{P-1}\equiv 1\pmod P$`，且我们有 `$m\mid(P-1)$`，故可以左右同取 `$\frac{P-1}{m}$` 次幂：

<div class="math-display">$$
\left(g^{km+x_0}\right)^{\frac{P-1}{m}}\equiv g^{k(P-1)}\cdot g^{\frac{(P-1)x_0}m}\equiv \left(g^{\frac{P-1}{m}}\right)^{x_0}\equiv a^{\frac{P-1}m}\pmod P
$$</div>

如果令 `$G=g^{\frac{P-1}m},A=a^{\frac{P-1}m}$`，那么我们就得到 `$G^{x_0}\equiv A\pmod P$` 的形式。

这个过程在模意义下也很直观，因为我们可以将整个同余式乘上 `$\frac{P-1}{m}$`：

<div class="math-display">$$
\frac{(P-1)x_0}{m}\equiv \frac{P-1}m\operatorname{ind}_ga\equiv \operatorname{ind}_g a^{\frac{P-1}m}\pmod{(P-1)}
$$</div>

Pohlig-Hellman 算法基于以下观察：若已知 `$x$` 满足 `$\operatorname{ind}_g a\equiv x\pmod {p}$`，那么我们容易求解 `$\operatorname{ind}_ga\equiv x'\pmod {pq}$`。设 `$x'=kp+x(0\le k<q)$`，还原成指数形式有：

<div class="math-display">$$
a^{\frac{P-1}{pq}}\equiv \left(g^{kp+x}\right)^{\frac{P-1}{pq}}\equiv\left(g^{\frac{P-1}{q}}\right)^k\cdot g^{\frac{(P-1)x}{pq}}\pmod P
$$</div>

由这个式子，我们可以用 BSGS 解出 `$k$`，进而得出 `$x'$`。单次递推的复杂度为 `$O(\sqrt q)$`。

利用这个结论，比较自然的想法是：选取一个模数路径 `$m_1\mid m_2\mid \dots\mid m_n=P-1$`，求出 `$\operatorname{ind}_g a\equiv x_1\pmod {m_1}$`，然后依次递推出答案。这当然是可行的，不过稍显复杂。

我们知道求解各种模意义下的方程，经常用到 CRT。在这个算法下，CRT 其实也可以看作是模数 `$p\to pq$` 的一次递推，即先求解模 `$p$` 和 `$q$` 的 case，然后合并。在底数固定的情况下，CRT 合并的开销极低，那么我们大可以将所有 `$p,q$` 互质的情况都用 CRT。于是我们只需要考虑用刚才的递推求解 `$\operatorname{ind}_ga\equiv x\pmod {p^\alpha}$`。

为了减小递推的开销，Pohlig-Hellman 算法直接令模数路径中相邻项之间的商为素因子 `$p$`，即令模数按照 `$p\to p^2\to\dots$` 的路径递推。那么以上的 P-H 递推可以形象地解释为 `$p$` 进制数位的递推：设 `$x$` 的 `$p$` 进制表示为 `$x=x_0+x_1p+\dots+x_{\alpha-1}p^{\alpha - 1}$`。这个对数方程左右同乘上一个 `$p^j$`，就抹去了 `$x$` 在 `$p$` 进制下的较高位：

<div class="math-display">$$
p^j\operatorname{ind}_ga\equiv \sum_{i&lt;\alpha}x_ip^{i+j}\pmod {p^\alpha}
$$</div>

于是等式右边 `$i+j\ge\alpha$` 的项就全部被抹去，得到一个关于 `$x_0,x_1,\dots,x_{\alpha-j-1}$` 的等式。为了能够计算，化为指数形式：

<div class="math-display">$$
A^{p^j}\equiv G^{\sum_{i+j&lt;\alpha}x_ip^{i+j}}\equiv \prod_{i+j&lt;\alpha}G^{x_ip^{i+j}}\pmod P
$$</div>

此时我们可以用 BSGS 算出 `$x_{\alpha-j-1}$`。实际实现中，由于 `$p$` 较小，有时也会直接枚举 `$x_{\alpha-j-1}$` 并验证来求解，以避免 BSGS 较大的常数。解决模 `$p^\alpha$` 的单个 case 的复杂度为 `$O(\sqrt p\alpha)$`，总复杂度为 `$O(\sum (\sqrt{p_i}\alpha_i+\log P))\approx O((\sqrt{\max p}+\omega(P-1))\log P)$`。

回顾 BSGS，我们是对指数做了一个根号分块。那么在 P-H 中，我们通过变换模数，将指数分为了 `$\alpha=O(\log P)$` 块求解，并保证这个分块有易于求解的性质。

更进一步地，模数路径中相邻项之商可以灵活选取。对较大的 `$\alpha$`，可以令模数按照 `$p^{\beta}\to p^{2\beta}\to\dots$` 的路径递推，即将指数看作 `$p^\beta$` 进制数。算法过程与之前类似，此时单个 case 的复杂度为 `$O(p^{\beta/2}\cdot\frac\alpha\beta)$`。

接下来考虑多次询问的情况下，如何预处理以优化单次询问。

- 注意到方程中含 `$x_{\alpha-j-1}$` 的一项中，指数内 `$p$` 的幂次必定被提升至了 `$\alpha-1$`，故在预处理时，可以打这样一个表：算出 `$\forall x\in[0,p)$` 的 `$G^{xp^{\alpha-1}}$`，并构建一个逆映射。这样求解时直接查询即可。

- 设在模 `$p_i^{\alpha_i}$` 的 case 下将指数看作 `$p_i^{\beta_i}$` 进制数，那么该 case 的预处理需要 `$O(p_i^{\beta_i})$` 时间/空间，查询复杂度则为 `$O(\alpha_i/\beta_i+\log P)$`。
- P-H 查询时的复杂度瓶颈之一在于，对 CRT 分离出的每个 `$p^\alpha$` case，我们都要算一遍 `$A=a^{\frac{P-1}{p^\alpha}}$` 及所有 `$A^{p^j}$`，这有时过于昂贵。对较小的 `$p^\alpha$`，我们可以设定阈值 `$S$`，贪心地将若干个 `$p^\alpha$` 分成一组（只要它们的积不超过 `$S$`），并暴力 BSGS 或全部预处理。这样我们就减少了 `$A$` 部分的快速幂次数。如对 `$P=998244353=2^{23}\times7\times17+1$`，可以将 `$7\times17$` 打包成一组。

事实上后两条优化都脱胎于模数路径的选取。原则很简单：合并互质模数时，CRT 优于 P-H 递推；视具体情况，平衡预处理和递推开销。