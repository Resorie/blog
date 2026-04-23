+++
title = "地狱人"
date = 2026-04-23
description = "Pohlig-Hellman 算法"

[taxonomies]
tags = ["数论", "离散对数", "算法"]

[extra]
+++

求解 `$g^x\equiv a\pmod P$`，即 `$\operatorname{ind}_ga\equiv x\pmod{(P-1)}$`。要求 `$P-1$` 较 smooth。

考虑分解 `$P-1$` 然后 CRT，那么我们只需要考虑求解 `$\operatorname{ind}_ga\equiv x\pmod {p^\alpha}$` 的情况。

设 `$x$` 的 `$p$` 进制表示为 `$x=x_0+x_1p+\dots+x_{\alpha-1}p^{\alpha - 1}$`。我们希望从低位向高位递推，因而我们希望有一个操作能够抹去 `$x$` 在 `$p$` 进制下的较高位。考虑在这个对数方程左右同乘上一个 `$p_j$`，那么有：

```math-display
$$
p^j\operatorname{ind}_ga\equiv \sum_{i<\alpha}x_ip^{i+j}\pmod {p^\alpha}
$$
```

那么等式右边 `$i+j\ge\alpha$` 的项就全部被抹去，得到一个关于 `$x_0,x_1,\dots,x_{\alpha-j-1}$` 的等式。

当然，用这个对数的式子，我们是没法实际去算出 `$x_i$` 的。我们要把这个变换模数后的对数式还原成幂次的形式。假设指数 `$x=kp^{\alpha}+x'$`，带入原方程：

```math-display
$$
g^{kp^\alpha+x'}\equiv a\pmod P
$$
```

我们希望消灭 `$kp^\alpha$`。由 Fermat 小定理，有 `$g^{P-1}\equiv 1\pmod P$`，且我们有 `$p^\alpha\mid(P-1)$`，故可以左右同取 `$\frac{P-1}{p^\alpha}$` 次幂：

```math-display
$$
(g^{kp^\alpha+x'})^{\frac{P-1}{p^\alpha}}=g^{k(P-1)}\cdot g^{\frac{(P-1)x'}{p^\alpha}}\equiv (g^{\frac{P-1}{p^\alpha}})^{x'}\equiv a^{\frac{P-1}{p^\alpha}}\pmod P
$$
```

以下设 `$G=g^{\frac{P-1}{p^\alpha}},A=a^{\frac{P-1}{p^\alpha}}$`。我们在该指数模 `$p^\alpha$` 的 case 中，需求解的即 `$G^x\equiv A\pmod P$`。现在我们可以将之前抹去高次系数的对数式子，化为可以计算的形式：

```math-display
$$
A^{p^j}\equiv G^{\sum_{i+j<\alpha}x_ip^{i+j}}\equiv \prod_{i+j<\alpha}G^{x_ip^{i+j}}\pmod P
$$
```

在递推过程中，我们已知 `$x_0,\dots,x_{\alpha-j-2}$`，故可以暴力枚举后验证或用 BSGS 算出 `$x_{\alpha-j-1}$`。解决模 `$p^\alpha$` 的单个 case 的复杂度为 `$O(p\alpha)$`（若用 BSGS 则为 `$O(\sqrt p\alpha)$`，实现中一般会对大的 `$p$` 使用）。总复杂度为 `$O(\sum (\sqrt{p_i}\alpha_i+\log P))\approx O((\sqrt{\max p}+\omega(P-1))\log P)$`。

回顾 BSGS，我们是对指数做了一个根号分块。那么在 P-H 中，我们将指数分为了 `$\alpha=O(\log P)$` 块求解。

更进一步，我们可以灵活选择块长。当 `$\alpha$` 较大时，可以考虑 `$x$` 的 `$p^\beta$` 进制表示，算法过程与之前类似，此时单个 case 的复杂度为 `$O(p^{\beta/2}\cdot\frac{\alpha}{\beta})$`。

接下来考虑多次询问的情况下，如何预处理以优化单次询问。

- 注意到方程中含 `$x_{\alpha-j-1}$` 的一项中，指数内 `$p$` 的幂次必定被提升至了 `$\alpha-1$`，故在预处理时，可以打这样一个表：算出 `$\forall x\in[0,p)$` 的 `$G^{xp^{\alpha-1}}$`，并构建一个逆映射。这样求解时直接查询即可。
- 设在模 `$p_i^{\alpha_i}$` 的 case 下将指数看作 `$p_i^{\beta_i}$` 进制数，那么该 case 的预处理需要 `$O(p_i^{\beta_i})$` 时间/空间，查询复杂度则为 `$O(\alpha_i/\beta_i+\log P)$`。
- P-H 查询时的复杂度瓶颈之一在于，对 CRT 分离出的每个 `$p^\alpha$` case，我们都要算一遍 `$A=a^{\frac{P-1}{p^\alpha}}$` 及所有 `$A^{p^j}$`，这有时过于昂贵。对较小的 `$p^\alpha$`，我们可以设定阈值 `$S$`，贪心地将若干个 `$p^\alpha$` 分成一组（只要它们的积不超过 `$S$`），并暴力 BSGS 或全部预处理。这样我们就减少了 `$A$` 部分的快速幂次数。如对 `$P=998244353=2^{23}\times7\times17+1$`，可以将 `$7\times17$` 打包成一组。

