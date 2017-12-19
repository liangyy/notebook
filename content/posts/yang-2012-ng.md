---
title: "Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits"
date: 2017-12-16T11:06:05-06:00
draft : false
author: "Yanyu Liang"
tags: ["gwas", 'integrative analysis', 'complex trait', 'linkage disequilibrium', 'conditional analysis']
categories: ["research paper"]
---

# Meta data of reading

* **Journal**: Nature Genetics
* **Year**: 2012
* **DOI**: 10.1038/ng.2213

$$
\newcommand\independent{\perp\\!\\!\\!\\!\perp}
\newcommand\E{\text{E}}
\newcommand\nocr{\nonumber\cr}
\newcommand\cov{\text{Cov}}
\newcommand\var{\text{Var}}
\newcommand\cor{\text{Cor}}
\newcommand\numberthis{\addtocounter{equation}{1}\tag{\theequation}}
$$

# Motivation

The GWAS procedure is to use single-SNP model to test association and select the SNP with strongest signal to represent the genomic region (~ 2Mb) and the genetic variation is computed based on this SNP only. The paper pointed out the underlying assumptions of this procedure:

> Implicit assumptions, often untested, are that the detected association at the top SNP captures the maximum amount of variation in the region by its LD with an unknown causal variant and that other SNPs in the vicinity show association because they are correlated with the top SNP.

They pointed out 2 reasons why this assumption may fail:

1. Suppose there is only one causal variant, a single tagging SNP may not capture all of its variation
2. It is possible that there are more than one causal variant

So, one-SNP-per-locus procedure may underestimate the underlying causal genetic variation. Some studies have performed conditional analysis to find the secondary SNP inside the locus. The paper proposed a systematic approach to perform conditional analysis by combining GWAS meta-analysis and LD correlation from the same population.

# Multi-SNP model and joint effect

The multi-SNP model is:

<div>$$\begin{align}
  \vec{y} &= X\vec{b} + \vec{e} \nonumber
\end{align}$$</div>

, where `$X \in \mathbb{R}^{n \times N}, \vec{b} \in \mathbb{R}^N, \vec{e} \in \mathbb{R}^n$`. The least squares solution is:

<div>$$\begin{align}
  \hat{b} &= (X'X)^{-1} X'y \label{eq:sej} \cr
  \var(\hat{b}) &= \sigma^2_{J} (X'X)^{-1} \label{eq:varj}
\end{align}$$</div>

Note that \eqref{eq:varj} is an N-by-N (co)variance matrix. It is derived as follows:

<div>$$\begin{align}
  \hat{b} &= (X'X)^{-1} X'y = (X'X)^{-1} X'(Xb + e) \nocr
  &= b + (X'X)^{-1} X'e \nocr
  \var(\hat{b}_j) &= \var([(X'X)^{-1} X'e]_j) \text{, where $[\cdot]$ takes the $i$th row}\nocr
  &= \var([(X'X)^{-1}]_j X'e) \nocr
  &= [(X'X)^{-1}_j X'] \odot [(X'X)^{-1}_j X'] \var(e) \label{eq:inter1}
\end{align}$$</div>

Note that in \eqref{eq:inter1}, we have `$[(X'X)^{-1}_j X'] \in \mathbb{R}^{1 \times n}$` and `$\var(e) \in \mathbb{R}^{n \times 1}$` which are vectors. But notice that $e_1, ..., e_n$ are iid. The expression can be simplified as:

<div>$$\begin{align}
  \eqref{eq:inter1} &= t_j X' X t_j' \sigma^2_J \nonumber
\end{align}$$</div>

, where `$t_j := [(X'X)^{-1}]_j$`. Similarly,

<div>$$\begin{align}
  \cov(\hat{b}_i, \hat{b}_j) &= t_i X' X t_j' \sigma^2_J \nonumber
\end{align}$$</div>

So,

<div>$$\begin{align}
  \var(\hat{b}) &= \sigma^2_J \begin{bmatrix} t_1 \cr \vdots \cr t_N \end{bmatrix}
  X'X \begin{bmatrix} t_1 & \cdots & t_N \end{bmatrix} \nocr
  &= \sigma^2_J (X'X)^{-1} (X'X) (X'X)^{-T} \nocr
  &= \sigma^2_J (X'X)^{-1} \text{, with the fact that $X'X$ is symmetric} \nonumber
\end{align}$$</div>

, which is the result of \eqref{eq:varj}.

# Single-SNP model and marginal effect

In practice, only single-SNP model results are available. The single-SNP model is as follow:

<div>$$\begin{align}
  y &= x_j \beta_j + e \nonumber
\end{align}$$</div>

So, the LS estimator is:

<div>$$\begin{align}
  \hat\beta &= D^{-1} X'y \label{eq:se} \cr
  \var(\hat\beta) &= \sigma^2_M D^{-1} \label{eq:varm}
\end{align}$$</div>

, where `$D = diag(\vec{d}), d_i = X_i' X_i$` and $X\_i$ is the $i$th column of $X$. The derivation of \eqref{eq:varm} is as follow:

<div>$$\begin{align}
  \var(\hat\beta_i) &= \var((X_i' X_i)^{-1} X_i' (X_i b_i + e)) \nocr
  &= \var(b_i + (X_i' X_i)^{-1} X_i' e) \nocr
  &= \var((X_i' X_i)^{-1} X_i' e) \nocr
  &= \frac{X_i' \odot X_i'}{(X_i' X_i)^2} \var(e) \nocr
  &= \frac{X_i' X_i}{(X_i' X_i)^2} \sigma^2_M \nocr
  &= \sigma^2_M (X_i' X_i)^{-1} \nocr
\end{align}$$</div>

Here, we treat single-SNP model as the truth. Note that we treat SNPs independently with each other (namely `$\cov(\hat\beta_i, \hat\beta_j) = 0$`). Different SNPs do not have to share the same residual variance, so a more precise expression is `$\var(\hat\beta_i) = \sigma^2_{M(i)} (X_i' X_i)^{-1}$`

# Inferring joint effect from single effect

From \eqref{eq:se}, we have `$X'y = D \hat\beta$`. Under multiple SNP model, the proportion of variance explained by all SNPs is:

<div>$$\begin{align}
  R_J^2 &= \frac{\cov(\hat{y}, y)}{\var(y)} \nocr
  &= \frac{\hat{b}' X' y}{y'y} \nocr
  &= \frac{\hat{b}' D \hat\beta}{y'y} \nonumber
\end{align}$$</div>

Then, we can derive:

<div>$$\begin{align}
  \hat\sigma^2_J &= \frac{(1 - R_J^2) y'y}{n - N} \nocr
  &= \frac{y'y - \hat{b}' D \hat\beta}{n - N}
\end{align}$$</div>

Similarly,

<div>$$\begin{align}
  R_{M(j)}^2 &= \frac{\hat{y}_j' y}{y'y} \nocr
  &= \frac{X_j' \hat\beta_j y}{y'y} \nocr
  &= \frac{\hat\beta_j X_j'y}{y'y} \nocr
  &= \frac{\hat\beta_j D_j \hat\beta_j}{y'y} \nocr
  &= \frac{\hat\beta_j^2 D_j}{y'y} \nocr
  \hat\sigma^2_{M(j)} &= \frac{(1 - R_{M(j)}) y'y}{n - 1} \nocr
  &= \frac{y'y - \hat\beta_j^2 D_j}{n - 1} \nonumber
\end{align}$$</div>

From \eqref{eq:varm}, we have `$\var(\hat\beta_j) = \hat\sigma^2_{M(j)} / D_j$`, so we get:

<div>$$\begin{align}
  y'y = D_j \var(\hat\beta_j) (n - 1) + D_j \hat\beta_j^2 \label{eq:yy}
\end{align}$$</div>

This expression provides a way to obtain $y'y$ with $\hat\beta\_j$ and $\hat{\var}(\hat\beta_j)$ (w/o knowing individual level data $y$). In practice, the paper used the median of inferred $y'y$ obtained from $j = 1, ..., N$.

* Side note

> The reason why the above analysis is performed is to obtain joint model statistic from the single-SNP model without querying individual level data $y$. One useful relation is: $\hat\sigma^2 = \frac{(1 - R^2) y'y}{N - n}$. $R^2$ is computable since it is just the observed proportion of covariance between predictor and response in the overall variance of response. For the single-SNP case, $\hat\sigma^2$ is available via \eqref{eq:varm}.

In meta-analysis, $X$ is not available as well. But since $X'X$ is the (co)variance matrix of SNP genotypes, it can be approximated by the LD score from a matched reference population. The paper used $W$ to denote such population, where `$w_{ij} = -2f_j, 1 - 2 f_j, 2 - 2 f_j$` to denote genotypes: two major alleles, heterozygous, two minor alleles respectively. Under this set up, `$\E(w_j) = 0, \var(w_j) = 2f_j(1 - 2 f_j)$`. Therefore, we have:

<div>$$\begin{align}
  \frac{\sum_i x_{ij} x_{ik}}{\sqrt{\sum_i x_{ij}^2 \sum_i x_{ik}^2}}
  &\approx \frac{\sum_i w_{ij} w_{ik}}{\sqrt{\sum_i w_{ij}^2 \sum_i w_{ik}^2}} \nocr
  \Rightarrow (X'X)_{jk} &= \sum_i x_{ij} x_{ik} \nocr
  &\approx \sum_i w_{ij} w_{ik} \sqrt{\frac{D_j D_k}{D_{W(j)}D_{W(k)}}} := B_{jk} \nocr
  \Rightarrow B &:= D^{1/2}D_W^{-1/2}W'W D_W^{-1/2} D^{1/2} \approx X'X \label{eq:app}
\end{align}$$</div>

where $D, D_W$ is the diagonal matrix with diagonal entries of $X'X, W'W$. $D\_{j}$ is not available without $X$, but it can be approximated by `$2p_j(1 - p_j)n$`. From \eqref{eq:sej}, \eqref{eq:se}, \eqref{eq:app}, we have:

<div>$$\begin{align}
  \hat{b} &= (X'X)^{-1} X'y = (X'X)^{-1} D \hat\beta \nocr
  &\approx B^{-1} D \hat\beta := \tilde{b} \nonumber
\end{align}$$</div>

Similar to \eqref{eq:varj},

<div>$$\begin{align}
  \var(\tilde{b}) &= \sigma^2_J B^{-1}
\end{align}$$</div>

In the paper, distant SNP pairs were assigned zero correlation instead the observed one in $W$. The paper pointed out an additional complexity in practice, that is SNPs may have different effective sample sizes due to imputation failures. Therefore, the paper suggested to estimate the effective sample for each SNP and use the adjusted sample size to compute $B_{jk}$. The procedure is:

  1. From \eqref{eq:yy}, we obtain $y'y$ by taking the median
  2. Obtain $\hat{n}_j$ using \eqref{eq:yy}:
      <div>$$\begin{align}
        \hat{n}_j &= y'y / D_j \var(\hat\beta_j) - \hat\beta_j^2 / \var(\hat\beta_j) + 1 \nonumber
      \end{align}$$</div>
  3. Compute $B\_{jk}$ and $D\_j$ using the adjusted sample size:
      <div>$$\begin{align}
        B_{jk} &= 2 \min(\hat{n}_j, \hat{n}_k) \sqrt{\frac{p_j(1 - p_j)p_k(1 - p_k)}{D\_{W(j)}D\_{W(k)}}} (W'W)\_{jk} \nocr
        W_j &= 2p_j(1 - p_j) \hat{n}_j \nonumber
      \end{align}$$</div>

# Obtain p-value of marginal effect under multi-SNP model

In brief, the above derivation provides a way to infer joint effect distribution from single-SNP model summary statistic. Namely, we get:

<div>$$\begin{align}
  (\tilde{b} - b) \sim \mathcal{N}(0, \var(\tilde{b}) \nonumber
\end{align}$$</div>

The marginal distribution for each SNP's effect size, $\tilde{b}\_i$, is:

<div>$$\begin{align}
  (\tilde{b}_i - b_i) \sim \mathcal(N, \var(\tilde{b}_i)) \nonumber
\end{align}$$</div>

Therefore, we can construct a test as follow:

  * $H_0$: $b_i$ is zero
  * $H_1$: $b_i$ is not zero

Under the null, $\tilde{b}\_i \sim \mathcal{N}(0, \var(\tilde{b}\_i))$. Then `$\mathbb{P}_{H_0}(|b_i| > |\tilde{b}_i| ) = 2(1 - \Phi(|\tilde{b}_i|))$`, which is the marginal effect of $i$th SNP under multi-SNP model.

# Conditional analysis

The logic of this part is not very intuitive to me, but after some struggling, I end up with the following things.

First of all, the conditional analysis takes a two step procedure to estimate $\hat{b}\_2 | \hat{b}\_1$. That is:

  1. Do $y \sim X\_1$ and obtain `$\bar{b}_1 = (X_1'X_1)^{-1} X_1'y$`
  2. Compute $\tilde{y} = y - X\_1 \bar{b}\_1$
  3. Do $\tilde{y} \sim X\_2$ and obtain $\hat{b}\_2 | \hat{b}\_1 = (X\_2' X\_2)^{-1} X\_2' \tilde{y}$, which matches the equation 15 in the text

The variance of $\hat{b}\_2 | \hat{b}\_1$ stuck me for a while because it is unclear how to define $\hat\sigma^2\_C$. It is still not so clear to me but what I get is the following:

<div>$$\begin{align}
  y &= X_1 b_1 + X_2 b_2 + e \label{eq:y}\cr
  \hat{b}_2 | \hat{b}_1 &= M_2^{-1} X_2'y - M_2^{-1} M_{21} M_1^{-1} X_1'y \nocr
  &\text{, where $M_{ij} = X_i' X_j$ and $M_{ii} = M_i$} \nocr
  &= M_2^{-1} M_{21} b_1 + M_2^{-1}M_2b_2 + M_2^{-1}X_2'e \nocr
  &- M_2^{-1}M_{21}M_1^{-1}M_1 b_1 - M_2^{-1}M_{21}M_1^{-1}M_{12}b_2 - M_2^{-1}M_{21}M_1^{-1}X_1'e \nocr
  &= (M_2^{-1}M_2 - M_2^{-1}M_{21}M_1^{-1}M_{12})b_2 \nocr
  &+ (M_2^{-1}X_2' - M_2^{-1}M_{21}M_1^{-1}X_1') e \nocr
  \var(\hat{b}_2 | \hat{b}_1) &= \var((M_2^{-1}X_2' - M_2^{-1}M_{21}M_1^{-1}X_1') e) \nocr
  &= (M_2^{-1}X_2' - M_2^{-1}M_{21}M_1^{-1}X_1')(M_2^{-1}X_2' - M_2^{-1}M_{21}M_1^{-1}X_1')' \sigma^2_J \nocr
  &= M_2^{-1}X_2' (I - X_1 M_1^{-1} X_1')(I - X_1 M_1^{-1} X_1')' X_2M_2^{-1} \sigma^2_J\nocr
  &= M_2^{-1}X_2' (I - 2X_1M_1^{-1}X_1' + X_1M_1^{-1}X_1'X_1M_1^{-1}X_1') X_2M_2^{-1} \sigma^2_J\nocr
  &= M_2^{-1}X_2' (I - X_1M_1^{-1}X_1') X_2M_2^{-1}\sigma^2_J \nocr
  &= [M_2^{-1} - M_2^{-1}M_{21}M_1^{-1}M_{12}M_2^{-1}]\sigma^2_J \nocr
\end{align}$$</div>

So, it seems to me that $\sigma^2\_C$ and $\sigma^2\_J$ are interchangeable if `$X = [X_1, X_2]$`. If such condition is not satisfied, \eqref{eq:y} is not the multi-SNP model, so it is better to denote the residual variance as $\sigma^2\_C$. The paper computed $\hat\sigma^2\_C$ using the following equation:

<div>$$\begin{align}
  \hat\sigma^2_C &= \frac{y'y - \hat{b}_1' D_1 \hat{\beta}_1 - (\hat{b}_2|\hat{b}_1)' D_2\hat\beta_2}{n - N_1 - N_2}
\end{align}$$</div>

Similar to previous derivation, the individual level statistics can be replaced by $D, \beta, B$. Note that `$M_{12}, M_{21}, M_1, M_2$` can be approximated by $B$.  

# Results in brief

The paper performed the analysis using GIANT GWAS and two reference genotype data. If two SNPs were in low LD, the result was similar to singe-SNP model. For positively correlated SNPs (modest LD), single-SNP model tended to overestimate the effect size. While multi-SNP model gave smaller effect size, they still reached genome-wide significance. For negatively correlated SNPs, single-SNP model tended to miss one of the signal because the signal was masked by the other one. Multi-SNP model was more powerful in this case. They found 36 loci with multiple signals with 38 leading SNPs and 49 additional SNPs. The result is robust to the choice of reference sample. Conditional analysis was also performed to identify secondary associations in the loci. The analysis was also applied to case-control study where $y$ is OR instead of quantitative trait.
