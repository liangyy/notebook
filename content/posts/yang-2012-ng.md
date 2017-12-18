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
  \hat{b} &= (X'X)^{-1} X'y \nocr
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

From \eqref{eq:varm}, we have `$\var(\hat\beta_j) = \hat\sigma^2_{M(j)} / D_j$`, so we get `$y'y = D_j \var(\hat\beta_j) (n - 1) + D_j \hat\beta_j^2$`. This expression provides a way to obtain $y'y$ with $\hat\beta\_j$ and $\hat{\var}(\hat\beta_j)$ (w/o knowing individual level data $y$). In practice, the paper used the median of inferred $y'y$ obtained from $j = 1, ..., N$.