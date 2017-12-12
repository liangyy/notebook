---
title: "Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets"
date: 2017-12-08T14:23:10-06:00
draft : false
author: "Yanyu Liang"
tags: ["gwas", 'integrative analysis', 'eqtl', 'target gene', 'mendelian randomization', 'causality']
categories: ["research paper"]
---

$$
\newcommand\independent{\perp\\!\\!\\!\\!\perp}
\newcommand\E{\text{E}}
\newcommand\cov{\text{Cov}}
\newcommand\var{\text{Var}}
$$

# Meta data of reading

* **Journal**: Nature Genetics
* **Year**: 2016
* **DOI**: 10.1038/ng.3538

# Instrumental variable and Mendelian randomization analysis

The paper mentioned it in the background section along with Mendelian randomization (MR). I googled this term and found this [site](http://www.statisticshowto.com/instrumental-variable/). Also, [wiki page](https://en.wikipedia.org/wiki/Instrumental_variables_estimation) is informative.


It appears in regression analysis where $y$ is the responsive variable and $x$ is the explanatory variable. The model is simply $Y = X \beta + \epsilon$. Suppose $Z$ correlates with $X$ and $Z$ is independent to $Y$ given $X$, then $Z$ correlates with $Y$ as well despite such conditional independence.


In some analysis, people sort of take the advantage of $Z$ (in this context $Z$ is called instrumental variable). In regression model, the underlying assumption is that $\epsilon$ is independent to $X$. But it is hard to achieve if there is some other unknown variable $U$ that affect both $X$ and $Y$. In this case, error term "absorbs" the dependencies in $U$ which leads to the violence of $X \independent \epsilon$ (the following graph is an example).

{{<mermaid align="center">}}
graph LR;
	Z((Z)) --- X((X))
    X --- Y((Y))
    U((U)) --- Y
    X --- U
    E((error)) --- Y
{{< /mermaid >}}

In this case, $Y \sim X$ cannot distinguish $U$ and $\epsilon$ but $Y \sim Z$ does not have such problem since the effect of $U$ is captured by $X$ and $\epsilon$ is disentangled.


The idea of MR can be illustrated by the following graph:

{{<mermaid align="center">}}
graph LR;
	Z(variant) --- X(gene expression)
    X --- Y(phenotype)
    U(unknown factors) --- Y
    X --- U
    E(error) --- Y
{{< /mermaid >}}

Here, the unknown factors can be environmental and regression model cannot tease apart error and them. But since variants are fixed, therefore we can use genetic variant as the instrumental variable to infer the dependency between gene expression and phenotype. It is the idea of MR.

# Motivation

MR analysis requires matched genotype and expression profile data and its power is limited by the strength of the dependency between genotype and phenotype along with the one between gene expression and phenotype. So, it is not useful in practice.

This paper tends to use summary statistic of GWAS and eQTL instead, which may not as strong as individual level data, but it is easy to achieve a much larger sample size.

# Method

This paper proposed a summary data-based MR method (SMR). In intuitively, it estimates effect of $Z$ (genotype) on $Y$ (phenotype), $b\_{zy}$ and effect of $Z$ on $X$, $b\_{zx}$ (note that the former is GWAS and the latter is eQTL mapping). Then $b\_{xy}$ is simply $\frac{b\_{zy}}{b\_{zx}}$. This approach, same as MR, estimates the effect of gene expression on phenotype which is only mediated by genetic component (in other word, free of non-genetic confounders). But the **caveat** is that the method cannot distinguish causality and pleiotropic effect (mediated by genetic confounders).

It appears to me that the derivation of MR (or instrumental variable regression) is not intuitive so I leave the note of MR part for another post [here](https://liangyy.github.io/notebook/posts/mendelian-randomization/#ivvar). Note that in MR, the same $Z$ is used for the estimation of $\hat{b}\_{zx}, \hat{b}\_{zy}$. SMR, instead, uses different $Z$, namely $Y \sim Z\_1$ and $X \sim Z\_2$ and $\hat{b}\_{xy} = \hat{b}\_{zy} / \hat{\beta}\_{zx}$. Since $Z\_1 \independent Z\_2$, $\cov(\hat{b}\_{zy}, \hat{\beta}\_{zx}) = 0$. By Delta method, we have:

$$\begin{aligned}
	\var(\hat{b}\_{xy}) &\approx \begin{bmatrix} \frac{1}{\beta\_{zx}} & -\frac{b\_{zy}}{\beta\_{zx}^2} \end{bmatrix} \begin{bmatrix} \var(\hat{b}\_{zy}) & \cov(\hat{b}\_{zy}, \hat{\beta}\_{zx}) \cr
	\cov(\hat{b}\_{zy}, \hat{\beta}\_{zx}) & \var(\hat{\beta}\_{zx}) \end{bmatrix}
	\begin{bmatrix} \frac{1}{\beta\_{zx}} \cr -\frac{b\_{zy}}{\beta\_{zx}^2} \end{bmatrix} \cr
	&= \frac{b\_{zy}^2}{\beta\_{zx}^2} \bigg[ \frac{\var(\hat\beta\_{zx})}{\beta\_{zx}^2} + \frac{\var{\hat{b}\_{zy}}}{b\_{zy}^2} - 2\frac{\cov(\hat\beta\_{zx}, \hat{b}\_{zy})}{\beta\_{zx}b\_{zy}} \bigg]
\end{aligned}$$

Then the $\chi^2$ statistic is:

$$\begin{aligned}
	T\_{\text{SMR}} &= \hat{b}\_{xy} / \var(\hat{b}\_{xy}) \cr
	&\approx \bigg(\frac{\hat{b}\_{zy}}{b\_{zy}}\bigg)^2 \bigg(\frac{\beta\_{zx}}{\hat\beta\_{zx}}\bigg)^2 \bigg/ \bigg[\frac{(\frac{\hat{b}\_{zy}}{b\_{zy}})^2}{z\_{zy}^2} + \frac{(\frac{\hat\beta\_{zx}}{\beta\_{zx}})^2}{z\_{zx}^2}\bigg] \cr
	&\approx \frac{1}{1 / z\_{zy}^2 + 1 / z\_{zx}^2} \cr
	&= \frac{z\_{zy}^2 z\_{zx}^2}{z\_{zy}^2 + z\_{zx}^2}
\end{aligned}$$

, which is described in equation 5.

Furthermore, some data sets only report effect size without variance or z-score. In this case, the variance can be recovered if the allele frequency is known (see the result at supplementary notes page 9 bottom equation). This result is surprising to me, so I scratch the derivation below.

Without loss of generality, let's consider $\hat{b}\_{zy}$. $\hat\beta\_{zx}$ follows the same rule.

* First "=":

<div>$$\begin{aligned}
	\hat{b}_{zy} &= (Z'Z)^{-1} (Z'y) \cr
	&= (Z'Z)^{-1} (Z b_{zy} + \epsilon) \cr
	&= b_{zy} + (Z'Z)^{-1} (Z' \epsilon) \cr
	&= b_{zy} + \frac{\sum_j (z_j - \bar{z}) \epsilon_j}{\sum_i (z_i - \bar{z})^2}
\end{aligned}$$</div>

Here $Z$ is the normalized version of the original $Z$ for simplicity. Since the expression evolves "divid", taking "variance" on both sides seems tedious. Alternatively, we can apply Delta method.

<div>$$\begin{aligned}
	\hat{A} &= \frac{\sum_j (z_j - \bar{z}) \epsilon_j}{n} \cr
	\hat{B} &= \frac{\sum_i (z_i - \bar{z})^2}{n} \cr
	\sqrt{n} (\hat{A} - 0) &\xrightarrow{d} \mathcal{N}(0, \var(z)\var(\epsilon)) \text{, since $z \independent \epsilon$} \cr
	\sqrt{n} (\hat{B} - \var(z)) &\xrightarrow{d} \mathcal{N}(0, \tau^2) \cr
  \sqrt{n} (\frac{\hat{A}}{\hat{B}} - 0) &\xrightarrow{d} \mathcal{N}(0, \nabla g(A, B)' \Sigma \nabla g(A, B)') \cr
	\text{, where } \nabla g(A, B) &= \begin{bmatrix} 1/B \cr 0 \end{bmatrix} \cr
	\Sigma &= \begin{bmatrix} \var(z) \var(\epsilon) & \cov(A, B) \cr
														\cov(A, B) & \tau^2 \end{bmatrix} \cr
	\text{therefore, } \nabla g(A, B)' \Sigma \nabla g(A, B)' &= \frac{\var(z) \var(\epsilon)}{B^2} \cr
	&= \frac{\var(\epsilon)}{\var(z)} \cr
	\text{thus, } \var(\hat{b}_{zy}) &= \var(\frac{\hat{A}}{\hat{B}}) \cr
	&= \frac{\var(\epsilon)}{n \var(z)}
\end{aligned}$$</div>

So that we obtain the first "equal sign": `$SE = \sqrt{\frac{\sigma_\epsilon}{n \var(z)}}$`

* Second "=":

<div>$$\begin{aligned}
	\E(z) &= 0 \times (1 - p)^2 + 1 \times 2p(1 - p) + 2 \times p^2 = 2p\cr
	\E(z^2) &= 0 \times (1 - p)^2 + 1 \times 2p(1 - p) + 4 \times p^2 = 2p(1 + p)\cr
	\var(z) &= \E(z^2) - \E(z)^2 = 2p(1 - p) \cr
	\E(\epsilon) &= \E(y - \tilde{z} b_{zy}) = \E(y) - b_{zy} \E(\tilde{z}) = 0 \text{, since $\E(y) = \E(\tilde{z}) = 0$}\cr
	\E(\epsilon^2) &= \E((y - \tilde{z} b_{zy})^2) = \E(y^2) + b_{zy}^2 \E(\tilde{z}^2) - 2 b_{zy} \E(y \tilde{z})  \cr
	&= \E(y^2) + \E(\tilde{z}^2) b_{zy}^2 - 2 b_{zy} [ b_{zy} \E(\tilde{z}^2) + \E(\tilde{z}) \E(\epsilon) ] \cr
	&= \E(y^2) - 2p (1 + p) b_{zy}^2 \cr
	\sigma_{\epsilon}^2 &= \E(\epsilon^2) - \E(\epsilon)^2 \cr
	&= \E(y^2) - 2p(1 - p) b_{zy}^2 \cr
	&= \var(y) - 2p(1 - p) b_{zy}^2 \cr
	&= 1 - 2p(1 - p) b_{zy}^2 \text{, since $\var(y) = 1$}
\end{aligned}$$</div>

, where $\tilde{z} := z - \E(z)$. Therefore,

<div>$$\begin{aligned}
	SE &= \sqrt{\frac{\sigma_\epsilon}{n \var(z)}} \cr
	&= \sqrt{\frac{1 - 2p(1 - p)b}{2p(1 - p)n}}
\end{aligned}$$</div>

, as the equation stated in supplementary notes page 9 (last one).

# Accounting for linkage

In practice, variants are correlated with each other because of LD. Therefore, the loci that is used for computation may show association if it is correlated with two causal variants which affect transcription and phenotype respectively. Figure 1b illustrated the scenarios.  

{{< figure src="/notebook/images/linkage.png" title="Three possible explanations of an association" >}}

It turns out that it is possible to filter out the "linkage" case. The paper suggested that if there is only one causal variant, the nearby region should have similar association signal. So, they proposed a hypothesis testing procedure where the null hypothesis is that the association signal follows uniform distribution (this corresponds to the non-linkage cases). Those rejected associations were removed. Ideally, causality part is the only biologically interesting one, but this approach cannot tease apart pleiotropic effect.
