+++
title = "模法"
date = 2025-07-23
description = "Montgomery Multiplications"

[taxonomies]
tags = ["数论", "算法"]

[extra]
+++

因为滚去 whk 鸽了两年的 Montgomery 取模。

众所周知取模很慢，因为一般的除法很慢。但是除一个 `$2$` 的幂或者对它取模很快。考虑将对固定模数 `$p$` 的取模转化成对一个 `$r=2^k$` 的除法和取模操作。要求 `$p$` 与 `$r$` 互质，具体使用中 `$p$` 通常是奇素数，显然满足要求。

由 `$p$` 与 `$r$` 互质，`$\{xr\mid 0\le x<p\}$` 是一个模 `$p$` 的完全剩余系。考虑用 `$r(x)=xr\bmod p$` 进行一些操作。

为了快速计算 `$r(x)$`，考察一个很邪恶的函数。令 `$r^{-1}$` 为模 `$p$` 意义下 `$r$` 的逆元，定义 `$redc(x)=xr^{-1}\bmod p$`。redc 是 reduce 的缩写，~~Montgomery 的论文里这么写的~~。先计算 `$redc(x)$`。把取模用不定方程替代掉：

```math-display
$$
r\cdot redc(x) + kp=x
\\
redc(x)=\frac{x-kp}{r}
$$
```

然后考察 `$k$`。由 `$r$` 的性质，自然想到对 `$r$` 取模：

```math-display
$$
kp\equiv x\pmod r
\\ k=xp^{-1} \bmod r
$$
```

由于模数固定，`$p^{-1}$` 可以预处理。这样算 `$redc(x)$` 就可以用两次乘法和一些移位搞定了，很快。注意稍微处理一下，保证 `$redc(x)\in[0,p)$`。于是 `$r(x)$` 也很好算：`$r(x)=redc(xr^2)$`。可以预处理一下 `$r^2\bmod p$`。

接着观测一下 `$r(x)$` 的性质。显然 `$r(x\pm y)=r(x)\pm r(y)$`。对乘法有 `$r(xy)=\frac{r(x)r(y)}{r}=redc\left(r(x)r(y)\right)$`。所以可以直接用 `$r(x)$` 代替 `$x$` 进行计算。从 `$r(x)$` 还原时有 `$x=redc(r(x))$`。

需要预处理 `$p^{-1}\bmod r$` 和 `$r^2\bmod p$`。前者可以对 `$f(x)=\dfrac{1}x-p$` 牛迭得到，每次令 `$x'=2x-x^2p$` 即可。后者没啥好办法，直接算吧。

复杂度 `$O\left(\dfrac{\log^2p}{w}\right)-O\left(\dfrac{\log^2p}{w^2}\right)$`，神速。

贴一下 [zzt 在 Loj#6466 的实现](https://loj.ac/s/316691)并且稍微解释一下：

```cpp
struct u256 {
	u128 lo, hi;
	u256() {}
	u256(u128 lo, u128 hi) : lo(lo), hi(hi) {}
	static u256 mul128(u128 a, u128 b) {
		u64 a_hi = a >> 64, a_lo = u64(a);
		u64 b_hi = b >> 64, b_lo = u64(b);
		u128 p01 = u128(a_lo) * b_lo;
		u128 p12 = u128(a_hi) * b_lo + u64(p01 >> 64);
		u64 t_hi = p12 >> 64, t_lo = p12;
		p12 = u128(a_lo) * b_hi + t_lo;
		u128 p23 = u128(a_hi) * b_hi + u64(p12 >> 64) + t_hi;
		return u256(u64(p01) | (p12 << 64), p23);
	}
} ;

struct Mont {
	u128 mod, inv, r2;
	Mont(u128 n) : mod(n) {
		assert(n & 1);
		inv = n;
		for (int i = 0; i < 6; ++i) inv *= 2 - n * inv;
		r2 = -n % n;
		for (int i = 0; i < 4; ++i) if ((r2 <<= 1) >= mod) r2 -= mod;
		for (int i = 0; i < 5; ++i) r2 = mul(r2, r2);
	}
	u128 reduce(u256 x) const {
		u128 y = x.hi - u256::mul128(x.lo * inv, mod).hi;
		return i128(y) < 0 ? y + mod : y;
	}
	u128 reduce(u128 x) const { return reduce(u256(x, 0)); }
	u128 init(u128 n) const { return reduce(u256::mul128(n, r2)); }
	u128 mul(u128 a, u128 b) const { return reduce(u256::mul128(a, b)); }
} mont(1);
```

`u256` 是一个压位 `u128`，`$x=hi\times2^{128}+lo$`。

`Mont` 是 Montgomery 乘法封装成一个结构体。`mod` 是模数 `$p$`，`inv` 是 `$p^{-1}\bmod r$`，`r2` 是 `$r^2\bmod p$`。但是这个 `r2` 的初始化我没看懂。`assert` 用来特判偶数，但是这题数据好像没偶数。

使用时先 `mont=Mont(p)`，每个数都先 `x=mont.init(x)` 一下。乘法直接 `prod=mont.mul(a,b)`。还原的时候 `mont.reduce` 一下自己就行。