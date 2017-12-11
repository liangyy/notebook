---
title: "Mendelian randomization and instrumental variable regression"
date: 2017-12-11T13:35:41-06:00
draft : false
author: "Yanyu Liang"
tags: ["mendelian randomization", "causality"]
categories: ["method note"]
---

$$
\newcommand\cov{\text{Cov}}
\newcommand\var{\text{Var}}
\newcommand\iv{\text{IV}}
$$

# Meta data of reading

* **Links**:
  1. <https://www.bauer.uh.edu/rsusmel/phd/ec1-8.pdf>
  2. <http://fmwww.bc.edu/EC-C/F2012/228/EC228.f2012.nn15.pdf>
  3. <http://mathworld.wolfram.com/Lindeberg-FellerCentralLimitTheorem.html>
  4. <http://mathworld.wolfram.com/LindebergCondition.html>
* **Year**: NA
* **DOI**: NA

# The problem and the idea of MR

Suppose we have phenotype $Y$, gene expression $X$, and genotype $Z$. The goal is to see if $X$ and $Y$ has some causal relationship. Since there are some unknown confounders, the residual of $Y \sim X$ is correlated with $X$. Therefore, the OLS estimator of effect size $\hat{\beta}\_{xy}$ is biased.

To account for such drawback, $Z$ is introduced as instrumental variable (IV) since genotype is pre-determined so that there should not be confounders that can affect genotype. Therefore, the residual of $Y \sim Z$ should not be correlated with $Z$. So, OLS estimator is unbiased. The estimator constructed by $\hat{\beta}\_{zy}$ and $\hat{\beta}\_{zx}$ is simply:

$$\begin{aligned}
\hat{\beta}\_{xy} &= \frac{\hat{\beta}\_{zy}}{\hat{\beta}\_{zx}}
\end{aligned}$$

# Consistency and the derivation of the variance

The following is derived from <https://www.bauer.uh.edu/rsusmel/phd/ec1-8.pdf>.

Formally, IV estimator is defined as $b_{\iv} = (Z'X)^{-1} Z'y$. Then, we have:

$$\begin{aligned}
  b_{\iv} &= (Z'X)^{-1} Z'y \cr
  &= (Z'X)^{-1} Z'(X \beta\_{xy} + \epsilon) \cr
  &= \beta\_{xy} + (Z'X)^{-1} Z'\epsilon \cr
  &\xrightarrow{p} \beta\_{xy}
\end{aligned}$$

since $Z$ is independent to error.

The [variance](#ivvar) of $b\_{\iv}$ is inspired by Lindeberg-Feller CLT:

$$\begin{aligned}
  \sqrt{N} (b\_{\iv} - \beta\_{xy}) &= \sqrt{N} (Z'X)^{-1} Z'\epsilon \cr
  &= (Z'X / N)^{-1} \sqrt{N} (Z'\epsilon / N) \cr
  (Z'\epsilon / N) \sqrt{N} &\xrightarrow{d} \mathcal{N}(0, \sigma^2 \var(Z)) \quad\text{, by L-F CLT} \cr
  \sqrt{N} (b\_{\iv} - \beta\_{xy}) &\xrightarrow{d} \mathcal{N}(0, \sigma^2 \cov(Z, X)^{-1} \var(Z) \cov(Z, X)^{-1})
\end{aligned}$$

The last line is kind of heuristic to me but I cannot find a justification for it (anyway ...). Note that $\sigma^2 := \var(\epsilon)$ and it turns out that $\hat{\sigma}^2 = \frac{1}{N} \sum\_i (y\_i - b\_{\iv} x\_i)^2$ is an unbiased estimator of this term. So, $\hat{\var}(b_{\iv}) = \hat{\sigma}^2 \hat{\cov}(Z, X)^{-1} \hat{\var}(Z) \hat{\cov}(Z, X)^{-1}$.

# The general procedure of MR

This section describes the MR procedure as stated in this [paper](https://www.nature.com/articles/ng.3538) (see [post](http://localhost:1313/notebook/posts/zhu-2017-ng/)).

$b\_{\iv}$ can be computed by two-step least squares (2SLS), simply $\hat{b}\_{xy} = \hat{b}\_{zy} / \hat{b}\_{zx}$. By the result of the above section, we have:

$$\begin{aligned}
  \var(\hat{b}\_{xy}) &= \frac{\text{unexplained variance by 2SLS}}{N} \times \frac{\var(Z)}{\cov(Z, X)^{2}} \cr
  &= \frac{\text{unexplained variance by 2SLS}}{N\var(X) \rho^2\_{xz}} \text{, Since } \rho^2\_{xz} := \frac{\cov(X, Z)^2}{\var(Z)\var(X)} \cr
  &= \frac{\var(Y) (1 - R\_{xy}^2)}{N \var(X) R\_{zx}^2}
\end{aligned}$$

, which justifies the result in the [paper](https://www.nature.com/articles/ng.3538) (equation 2).
