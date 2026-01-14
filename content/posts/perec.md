+++
title = "PE 随机做题寄录"
date = 2025-12-18
description = "Project Euler 随机做题记录。"

[taxonomies]
tags = ["数学", "数论"]
+++

不定期持续（不保证）更新。

## PE169 Sums of Powers of Two

date: 12.17 diff: 50

> [!NOTE]
> 题意：定义 `$f(n)$` 为 `$n$` 表示为若干个 `$2$` 的幂之和的形式，且每个幂出现次数最多两次的方案数。
>
> 已知 `$f(10)=5$`。求 `$f(10^{25})$`。

AI 用递推秒了。这里叙述一种我自己的神秘做法。

考虑将 `$n$` 的二进制表示从高次幂向低次幂“展开”。具体地说，将 `$2^n$` 展开至 `$k$` 次表示：

```math-display
$$
2^n=2^{n-1}+\dots+2\cdot2^k
$$
```

发现这里的 `$k$` 除了不能是原表示中已有的幂以外没有什么限制。于是我们可以安全地令这样的展开穿过某些中间的幂次 `$2^m(k<m<n)$`。注意展开到 `$k$` 次之后所有高于 `$k$` 次的幂都不能再展开了，于是可以递归求解。

复杂度 `$O(\log n)$`。

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

## PE216 The Primality of `$2n^2-1$`

date: 12.18 diff: 45


> [!NOTE]
> 题意：令 `$t(n)=2n^2-1$`，求 `$2\le n\le N$` 时有多少个 `$n$` 使得 `$t(n)$` 为素数。
>
> 已知 `$N=10000$` 时有 `$2202$` 个。求 `$N=5\times10^7$` 时有多少。

这种傻逼题怎么能有 45 difficulty 的啊？为什么评论区还有一堆 Miller-Rabin 硬草过去的啊？

直接筛 `$t(n)$` 这个多项式即可。`$p\mid t(n)\iff n\equiv \pm1/2\pmod p$`。

复杂度是常数比较小的 `$O(N\log\log N)$`。分块筛可以做到更小的时空复杂度。

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

## PE433 Steps in Euclid's Algorithm

date: 26.1.13 diff: 65

> [!NOTE]
> 题意：考察欧几里得算法求 `$\gcd(x,y)$` 的流程：
> 1. 令 `$(x_0, y_0)=(x, y)$`；
> 2. 递推计算 `$(x_k, y_k)=(y_{k-1}, x_{k-1}\bmod y_{k-1})$`；
> 3. 最终得到 `$(x_n, y_n)=(\gcd(x,y),0)$`。
> 
> 令 `$E(x,y)=n$`，即使用欧几里得算法求 `$\gcd(x,y)$` 的运算步数。令 `$S(N)=\sum_{1\le x,y\le N} E(x,y)$`。
>
> 求 `$S(5\times 10^6)$`。已知 `$S(1)=1, S(10)=221, S(100)=39826$`。

神秘论文题。reference: <https://link.springer.com/content/pdf/10.1007/978-1-4615-4819-5_7.pdf>

首先显然根据对称性：

```math-display
$$
S(N)=\frac{N(N+1)}2+2\sum_{1\le y<x\le N}E(x,y)
$$
```

论文中有结论，对 `$n>2$`：

```math-display
$$
f(n)=\sum_{1\le y<n}[n\perp y]E(n,y)=\frac32\varphi(n)+2r(n)
$$
```


其中：



```math-display
$$
r(n)=\sum_{x<y,x'<y'}[x\perp y][x'\perp y'][xx'+yy'=n]
$$
```

化简一下 `$f$`，去掉互质的要求：

```math-display
$$
g(n)=\sum_{1\le y<n}E(n,y)=\lfloor 3(n-1)/2\rfloor + \tilde 2r(n)
$$
```

其中：

```math-display
$$
\tilde r(n)=\sum_{x<y,x'<y'}[x\perp y][xx'+yy'=n]
$$
```


回到和式，就变成简单莫反了。


```math-display
$$
\begin{aligned}
\sum_{1\le y<x\le N}E(x,y)
=&\ \sum_{x=2}^Ng(x)\\
=&\ \sum_{x=2}^N\lfloor 3(x-1)/2\rfloor + 2\sum_{x<y,x'<y'}[x\perp y][xx'+yy'\le N]\\
=&\ \sum_{x=2}^N\lfloor 3(x-1)/2\rfloor + 2\sum_{d}\mu(d)\sum_{x<y,x'<y'}[xx'+yy'\le N/d]
\end{aligned}
$$
```


还需要计算 `$\sum_{x<y,x'<y'}[xx'+yy'\le n]$`。枚举 `$P=xx'\le(n-1)/2=L$`，将 `$y>x,y'>x'$` 的限制利用容斥拆开：

```math-display
$$
\begin{aligned}
&\ \sum_{x<y,x'<y'}[xx'+yy'\le n]\\
=&\ \sum_{xx'\le L}\left(\sum_{y,y'}[yy'\le n-xx'](1-[y\le x] - [y'\le x']+[y\le x,y'\le x'])\right)\\
=&\ \sum_{xx'\le L}\left(D(n-xx')-\sum_{y\le x}\lfloor\frac{n-xx'}{y}\rfloor-\sum_{y'\le x'}\lfloor\frac{n-xx'}{y'}\rfloor+xx'\right)\\
=&\ \sum_{P=1}^{L}\left(d(P)D(n-P)+d(p)P-2\sum_{x\mid P}\sum_{y\le x} \lfloor\frac{n-P}{y}\rfloor\right)
\end{aligned}
$$
```

其中 `$d(n)=\sum_{d\mid n}1,D(n)=\sum_{i=1}^nd(i)$`。前面两项容易计算。最后一项拎出来计算：

```math-display
$$
\begin{aligned}
\sum_{P=1}^{L}\sum_{x\mid P}\sum_{y\le x}\lfloor\frac{n-P}y\rfloor
=&\ \sum_{x=1}^L\sum_{k=1}^{L/x}\sum_{y=1}^x\lfloor\frac{n-kx}{y}\rfloor\\
=&\ \sum_{y=1}^L\sum_{k=1}^{L/y}\sum_{x=y}^{L/k}\lfloor\frac{n-kx}{y}\rfloor\\
=&\ \sum_{y=1}^L\sum_{k=1}^{L/y}\sum_{x=0}^{L/k-y}\lfloor\frac{n-k(L/k)+kx}{y}\rfloor
\end{aligned}
$$
```

最内层的一个和式可以类欧 `$O(\log n)$` 计算。枚举外面的两层进行求和，复杂度总共为 `$O(n\log^2 n)$`。

复杂度（设 `$f(n)=O(n\log^2 n)$`）：

```math-display
$$
T(N)\approx O\left(\sum_{n\le \sqrt N}f(n)+f(N/n)\right)=O(N\log^3 N)
$$
```


code：

```cpp
#include <iostream>
#include <vector>

using namespace std;

typedef long long ll;

const ll N = 5e6;

std::vector<int> primes, mu;
std::vector<ll> d, D;

void Init(int n) {
    std::vector<bool> isp(n + 5, true);
    std::vector<int> e(n + 5);
    mu.resize(n + 5), d.resize(n + 5), D.resize(n + 5);

    mu[1] = d[1] = e[1] = 1;

    auto append = [&](int p) {
        primes.push_back(p), mu[p] = -1, d[p] = 2, e[p] = 1;
    };

    for (int i = 2; i <= n; i++) {
        if (isp[i]) append(i);
        for (int p : primes) {
            if (i * p > n) break;
            isp[i * p] = false;

            if (!(i % p)) {
                mu[i * p] = 0;
                d[i * p] = d[i] / (e[i] + 1) * (e[i] + 2);
                e[i * p] = e[i] + 1;
                break;
            }

            mu[i * p] = -mu[i];
            d[i * p] = d[i] * 2, e[i * p] = 1;
        }
    }

    for (int i = 1; i <= n; i++)
        D[i] = D[i - 1] + d[i], mu[i] += mu[i - 1];
}

ll f(ll a, ll b, ll c, ll n) {
    if (!a) return b / c * (n + 1);
    if (a >= c || b >= c)
        return n * (n + 1) / 2 * (a / c) + b / c * (n + 1) +
               f(a % c, b % c, c, n);
    ll m = (a * n + b) / c;
    return n * m - f(c, c - b - 1, a, m - 1);
}

ll S(int n) {
    ll res = 0, T = 0, L = (n - 1) / 2;
    for (int p = 1; p <= (n - 1) / 2; p++)
        res += d[p] * (D[n - p] + p);

    for (int y = 1; y <= L; y++)
        for (int t = L / y, k = 1; k <= t; k++)
            T += f(k, n - L / k * k, y, L / k - y);
    return res - 2 * T;
}

int main() {
    Init(N);

    ll ans = N * (N + 1) / 2;

    for (int i = 2; i <= N; i++)
        ans += 3LL * (i - 1) / 2 * 2;

    for (int r, t, l = 1; l <= N; l = r + 1) {
        t = N / l, r = N / t;
        ans += 4LL * S(t) * (mu[r] - mu[l - 1]);
    }

    std::cout << ans << '\n';

    std::cout << clock() * 1. / CLOCKS_PER_SEC << '\n';
    return 0;
}
```

