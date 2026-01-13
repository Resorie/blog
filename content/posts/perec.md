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
> 题意：考察欧几里得算法求 `$\gcd(a,b)$` 的流程：
> 1. 令 `$(x_n, y_n)=(i, j)$`；
> 2. 递推计算 `$(x_{n-1}, y_{n-1})=(y_n, x_n\bmod y_n)$`；
> 3. 最终得到 `$(x_0, y_0)=(\gcd(a,b),0)$`。
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


`$\sum[xx'+yy'\le n]$` 可以用整除分块 `$O(n\sqrt n)$` 算出。

复杂度：

```math-display
$$
T(N)\approx O\left(\sum_{n\le \sqrt N}(N/n)^{3/2}\right)=O(N^{5/4})
$$
```

code（AI 写的）：

```cpp
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

typedef long long ll;

const int MAXN = 5000005;

// Global arrays for precomputed values
int d_cnt[MAXN]; // d(k): number of divisors
ll D_sum[MAXN];  // D(k): prefix sum of d(k)
int mu[MAXN];    // Mobius function
int primes[MAXN];
int p_cnt = 0;
bool is_prime[MAXN];

// Linear Sieve for Mobius and standard precomputation for Divisors
// time: O(N)
void precompute_all(int n) {
    // 1. Linear Sieve for Mu
    fill(is_prime, is_prime + n + 1, true);
    is_prime[0] = is_prime[1] = false;
    mu[1] = 1;

    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) {
            primes[p_cnt++] = i;
            mu[i] = -1;
        }
        for (int j = 0; j < p_cnt && i * primes[j] <= n; ++j) {
            int p = primes[j];
            is_prime[i * p] = false;
            if (i % p == 0) {
                mu[i * p] = 0;
                break;
            } else {
                mu[i * p] = -mu[i];
            }
        }
    }

    // 2. Compute divisor count d[i]
    // Using harmonic loop is simple and sufficiently fast (O(N log N))
    fill(d_cnt, d_cnt + n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        for (int j = i; j <= n; j += i) {
            d_cnt[j]++;
        }
    }

    // 3. Compute prefix sums D[i]
    D_sum[0] = 0;
    for (int i = 1; i <= n; ++i) {
        D_sum[i] = D_sum[i - 1] + d_cnt[i];
    }
}

// Calculates K(n) = sum_{x<y, x'<y'} [xx' + yy' <= n]
// Logic adapted directly from sol.cpp
ll calc_K(int N) {
    if (N < 2) return 0;

    ll ans = 0;
    int limit = (N - 1) / 2;

    // Term 1: sum d(P) * D(N-P) and Term 2: sum P * d(P)
    for (int P = 1; P <= limit; ++P) {
        ans += (ll)d_cnt[P] * D_sum[N - P];
        ans += (ll)P * d_cnt[P];
    }

    // Term 3: -2 * sum_{P} sum_{x|P} H(x, N-P)
    // Implemented using the sqrt decomposition strategy from sol.cpp
    ll term3 = 0;
    int sqrt_limit = sqrt(limit);

    // Part A: x <= sqrt_limit
    for (int x = 1; x <= sqrt_limit; ++x) {
        int max_k = limit / x;
        for (int k = 1; k <= max_k; ++k) {
            int P = x * k;
            int M = N - P;

            // Calculate H(x, M) naively as x is small
            ll h_val = 0;
            for (int i = 1; i <= x; ++i) {
                h_val += M / i;
            }
            term3 += h_val;
        }
    }

    // Part B: x > sqrt_limit (implies k is small)
    for (int k = 1; k <= limit / (sqrt_limit + 1); ++k) {
        int min_x = sqrt_limit + 1;
        int max_x = limit / k;

        for (int x = min_x; x <= max_x; ++x) {
            int P = x * k;
            int M = N - P;

            // Calculate H(x, M) using integer division blocking (整除分块)
            // This is efficient because M is large
            ll h_val = 0;
            for (int l = 1, r; l <= x; l = r + 1) {
                if (M / l == 0) {
                    r = x;
                } else {
                    r = min(x, M / (M / l));
                    h_val += (ll)(r - l + 1) * (M / l);
                }
            }
            term3 += h_val;
        }
    }

    ans -= 2 * term3;
    return ans;
}

int main() {
    // Optimization for faster IO
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N = 5e6;

    // 1. Precompute globals
    precompute_all(N);

    // 2. Calculate Part 1: sum_{x=2}^N floor(3(x-1)/2)
    // Based on formula derived in pe433.md
    ll sum_part1 = 0;
    for (int x = 2; x <= N; ++x) {
        sum_part1 += (3LL * (x - 1)) / 2;
    }

    // 3. Calculate Part 2: 2 * sum_{d} mu(d) * K(floor(N/d))
    // Based on formula derived in pe433.md
    ll sum_part2 = 0;
    for (int d = 1; d <= N; ++d) {
        if (mu[d] == 0) continue;

        int arg = N / d;
        if (arg < 2) continue; // K(0) = K(1) = 0

        ll k_val = calc_K(arg);

        if (mu[d] == 1)
            sum_part2 += k_val;
        else
            sum_part2 -= k_val;
    }

    // 4. Final Assembly
    // S(N) = N(N+1)/2 + 2 * (Part1 + 2 * Part2)
    // Note: The lower triangle sum E(x,y) = Part1 + 2*Part2
    ll lower_triangle_sum = sum_part1 + 2 * sum_part2;
    ll diagonal_sum = (ll)N * (N + 1) / 2;
    ll total_ans = diagonal_sum + 2 * lower_triangle_sum;

    cout << total_ans << endl;

    cout << clock() * 1. / CLOCKS_PER_SEC << '\n';

    return 0;
}
```





