+++
title = "PE 随机做题寄录"
date = 2025-12-18
description = "Project Euler 随机做题记录。"

[taxonomies]
tags = ["数学", "数论"]

[extra]
pinned = true
+++

不定期持续（不保证）更新。

## PE169 Sums of Powers of Two

date: 12.17 diff: 50

> [!NOTE]
> 题意：定义 $f(n)$ 为 $n$ 表示为若干个 $2$ 的幂之和的形式，且每个幂出现次数最多两次的方案数。
>
> 已知 $f(10)=5$。求 $f(10^{25})$。

AI 用递推秒了。这里叙述一种我自己的神秘做法。

考虑将 $n$ 的二进制表示从高次幂向低次幂“展开”。具体地说，将 $2^n$ 展开至 $k$ 次表示：
$$
2^n=2^{n-1}+\dots+2\cdot2^k
$$
发现这里的 $k$ 除了不能是原表示中已有的幂以外没有什么限制。于是我们可以安全地令这样的展开穿过某些中间的幂次 $2^m(k<m<n)$。注意展开到 $k$ 次之后所有高于 $k$ 次的幂都不能再展开了，于是可以递归求解。

复杂度 $O(\log n)$。

code：

```python
from functools import lru_cache
from math import log

def bits(N):
    l = [i for i in range(int(log(N) / log(2)) + 2) if N // (2**i) % 2 == 1]
    l.reverse()
    return l

def toN(l):
    return sum(2**i for i in l)

@lru_cache(maxsize=None)
def solve2(N):
    bit = bits(N)
    if len(bit) == 1:
        return bit[0] + 1
    res = bit[-1] + solve2(toN(bit[1:]))
    for ix in range(1, len(bit)):
        res += (bit[ix - 1] - bit[ix] - 1) * solve2(toN(bit[ix:]))
    return res

print(solve2(N))
```

## PE216 The Primality of $2n^2-1$

date: 12.18 diff: 45


> [!NOTE]
> 题意：令 $t(n)=2n^2-1$，求 $2\le n\le N$ 时有多少个 $n$ 使得 $t(n)$ 为素数。
>
> 已知 $N=10000$ 时有 $2202$ 个。求 $N=5\times10^7$ 时有多少。

这种傻逼题怎么能有 45 difficulty 的啊？为什么评论区还有一堆 Miller-Rabin 硬草过去的啊？

直接筛 $t(n)$ 这个多项式即可。$p\mid t(n)\iff n\equiv \pm1/2\pmod p$。

复杂度是常数比较小的 $O(N\log\log N)$。分块筛可以做到更小的时空复杂度。

code:（没怎么优化，空间复杂度很大）

```cpp
#include <bitset>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

using i64 = long long;
using pii = std::pair<i64, i64>;

const int N = 5e7;
const int L = 8e7;

std::vector<i64> primes, rt;
std::bitset<(L >> 1) + 1> isp;
std::bitset<N + 1> tn;

i64 fpow(i64 a, i64 b, i64 m) {
    i64 res = 1;
    for (a %= m; b; b >>= 1, a = a * a % m)
        if (b & 1) res = res * a % m;
    return res;
}

i64 sub(i64 a, i64 b, i64 m) { return ((a - b) % m + m) % m; }

pii mkp(i64 first, i64 second) { return std::make_pair(first, second); }

bool chk_rt(i64 a, i64 p) { return fpow(a, p >> 1, p) == 1; }

// Cipolla
i64 find_rt(i64 a, i64 p) {
    if (a == 0) return 0;
    if ((p & 3) == 3) return fpow(a, (p + 1) >> 2, p);

    static std::mt19937 rd(std::random_device{}());
    i64 r, y;
    do
        r = rd() % p;
    while (chk_rt(y = sub(r * r, a, p), p));

    auto mul = [&](const pii &a, const pii &b) {
        return mkp(
            (a.first * b.first + y * a.second % p * b.second) % p,
            (a.first * b.second + a.second * b.first) % p
        );
    };

    auto ppow = [&](const pii &a, i64 b) {
        pii res = mkp(1, 0), bs = a;
        for (; b; b >>= 1, bs = mul(bs, bs))
            if (b & 1) res = mul(res, bs);
        return res;
    };

    return ppow(mkp(r, 1), (p + 1) >> 1).first;
}

i64 t(i64 n) { return 2 * n * n - 1; }

int main() {
    auto append = [&](i64 p) {
        primes.push_back(p);
        rt.push_back(find_rt((p + 1) >> 1, p));
    };

    for (int i = 3; i <= L; i += 2) { // odd
        if (!isp[i >> 1]) {
            // primes i satisfying 1/2 is a quadratic residue
            if (chk_rt((i + 1) >> 1, i)) append(i);

            for (int j = 3 * i; j <= L; j += (i << 1))
                isp[j >> 1] = true;
        }
    }

    int pcnt = primes.size();
    for (int i = 0; i < pcnt; i++) {
        i64 p = primes[i], c = (rt[i] % p + p) % p;

        // p|t(n) iff n^2 equiv 1/2 mod p iff n equiv pm c mod p
        i64 s1 = c, s2 = p - c;
      
        for (i64 j = 0; j <= N; j += p) {
            if (s1 + j <= N && t(s1 + j) != p) tn[s1 + j] = true;
            if (s2 + j <= N && t(s2 + j) != p) tn[s2 + j] = true;
        }
    }

    int res = 0;
    for (int i = 2; i <= N; i++)
        if (!tn[i]) res++;

    std::cout << res << '\n';
    std::cerr << clock() * 1. / CLOCKS_PER_SEC;
    return 0;
}
```

