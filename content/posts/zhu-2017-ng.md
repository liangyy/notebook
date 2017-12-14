---
title: "Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets"
date: 2017-12-08T14:23:10-06:00
draft : false
author: "Yanyu Liang"
tags: ["gwas", 'integrative analysis', 'eqtl', 'target gene', 'mendelian randomization', 'causality']
categories: ["research paper"]
---

$$
\require{AMSmath}
\newcommand\independent{\perp\\!\\!\\!\\!\perp}
\newcommand\E{\text{E}}
\newcommand\nocr{\nonumber\cr}
\newcommand\cov{\text{Cov}}
\newcommand\var{\text{Var}}
\newcommand\numberthis{\addtocounter{equation}{1}\tag{\theequation}}
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
	&= \frac{b\_{zy}^2}{\beta\_{zx}^2} \bigg[ \frac{\var(\hat\beta\_{zx})}{\beta\_{zx}^2} + \frac{\var(\hat{b}\_{zy})}{b\_{zy}^2} - 2\frac{\cov(\hat\beta\_{zx}, \hat{b}\_{zy})}{\beta\_{zx}b\_{zy}} \bigg]
\end{aligned}$$

Then the $\chi^2$ statistic is:

$$\begin{aligned}
	T\_{\text{SMR}} &= \hat{b}\_{xy} / \var(\hat{b}\_{xy}) \cr
	&\approx \bigg(\frac{\hat{b}\_{zy}}{b\_{zy}}\bigg)^2 \bigg(\frac{\beta\_{zx}}{\hat\beta\_{zx}}\bigg)^2 \bigg/ \bigg[\frac{(\frac{\hat{b}\_{zy}}{b\_{zy}})^2}{z\_{zy}^2} + \frac{(\frac{\hat\beta\_{zx}}{\beta\_{zx}})^2}{z\_{zx}^2}\bigg] \cr
	&\approx \frac{1}{1 / z\_{zy}^2 + 1 / z\_{zx}^2} \cr
	&= \frac{z\_{zy}^2 z\_{zx}^2}{z\_{zy}^2 + z\_{zx}^2}
\end{aligned}$$

, which is described in equation 5.

Furthermore, some data sets only report z-score without effect size. In this case, the inference can still be done but effect size is needed to estimate `$b_{xy}$`. In fact, the effect size can be recovered if the allele frequency is known (see the result at supplementary notes page 9 bottom and page 10 top equations). This result is surprising to me, so I scratch the derivation below.

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

So that we obtain the first "equal sign": `$SE = \sqrt{\frac{\sigma_\epsilon^2}{n \var(z)}}$`

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
	SE &= \sqrt{\frac{\sigma_\epsilon^2}{n \var(z)}} \cr
	&= \sqrt{\frac{1 - 2p(1 - p)b}{2p(1 - p)n}}
\end{aligned}$$</div>

, as the equation stated in supplementary notes page 9 (last one).

With z-score known, we have:

<div>$$\begin{aligned}
	z_{zy} &= \frac{\hat{b}_{zy}}{SE_{zy}} \cr
	\hat{b}_{zy} &= z_{zy} SE_{zy} \cr
	&= z_{zy} \sqrt{\frac{1 - 2p(1 - p)b_{zy}}{2p(1 - p)n}} \cr
	\Rightarrow \hat{b}_{zy}^2 &\approx z_{zy}^2 \frac{1 - M \hat{b}_{zy}^2}{Mn} \text{, where $M = 2p(1 - p)$} \cr
	\Rightarrow \hat{b}_{zy} &\approx z_{zy} / \sqrt{2p (1 - p) (n + z_{zy}^2)}
\end{aligned}$$</div>

, which is the result on page 10 top of supplementary notes.

# Accounting for linkage

In practice, variants are correlated with each other because of LD. Therefore, the loci that is used for computation may show association if it is correlated with two causal variants which affect transcription and phenotype respectively. Figure 1b illustrated the scenarios.  

{{< figure src="/notebook/images/linkage.png" title="Three possible explanations of an association" >}}

It turns out that it is possible to filter out the "linkage" case. The paper suggested that if there is only one causal variant, the nearby region should have similar association signal. So, they proposed a hypothesis testing procedure to get rid of signal resulted from linkage. The following shows how they came up with the distribution under null hypothesis.

Let's say there are $k$ loci in LD to the causal loci which is denoted as loci 0.  For the causal loci and the $i$th locus, we have (consider $y \sim z$):

<div>$$\begin{align}
	y &= b_{yz(0)} z_{(0)} + \epsilon_{(0)} \nonumber\cr
	y &= b_{yz(i)} z_{(i)} + \epsilon_{(i)} \nonumber\cr
	\cov(y, z_{(i)}) &= b_{yz(0)} \cov(z_{(0)}, z_{(i)}) \nonumber\cr
	\cov(y, z_{(i)}) &= b_{yz(i)} \var(z_{(i)}) \nonumber\cr
	& \text{, by the fact that $\epsilon_{(j)} \independent z_{(i)}, j = 0, i$} \nonumber\cr
	\Rightarrow b_{yz(i)} &= b_{yz(0)} \frac{\cov(z_{(0)}, z_{(i)})}{\var(z_{(i)})} \nonumber\cr
	&= b_{yz(0)} r_{0i} \sqrt{\frac{\var(z_{(0)})}{\var(z_{(i)})}} \nonumber\cr
	&= b_{yz(0)} r_{0i} \sqrt{h_0 / h_i} \label{eq:1}
\end{align}$$</div>

, where $h = 2p (1 - p)$ which is precisely $\var(z)$. So equation 7 follows. Intuitively, consider the following situation:

{{<mermaid align="center">}}
graph LR;
	z(z) --- x(x)
    x --- y(y)
{{< /mermaid >}}

In this case, `$b_{zy} = b_{xy} r_{zx} SE_x / SE_y$`. Namely, if we know the dependency between $x$ and $z$, so does $x$ and $y$, then we can derive the dependency between $z$ and $y$. `$r_{zx}$` indicates to what extent the effect size of $x$ on $y$ can be transferred to the one of $z$ on $y$. `$SE_x / SE_y$` is a rescaling term. To match the same scale (namely $y$), $b \times SE$ should of the same scale for $z \sim y$ and $x \sim y$. Therefore, under the null hypothesis, `$\hat{b}_{zy(i)}$` should have the same expected value. To make inference, we need some distribution of $\vec{b}_{zy}$, so we should know the variance/covariance as well.

* Side note:

> I just found that the above derivation has a subtle drawback. It is that suppose we do `$\cov(\cdot, z_{(0)})$` instead and do it carelessly, it will fail. The problem comes from the term `$\cov(\epsilon_{(i)}, z_{(0)})$`. From the derivation above, it seems to me that `$z_{(0)}$` and `$z_{(i)}$` are equivalent and interchangeable but it is really not the case. The graph above has pointed out such difference but in an implicit way. A more explicit illustration is the following:

{{<mermaid align="center">}}
graph LR;
	z(z) --- x(x)
    x --- y(y)
	e0(e0) --> y
	e1(e1) --> x
{{< /mermaid >}}

> The corresponding equations are:

<div>$$\begin{aligned}
	y &= b_0 x + \epsilon_0 \cr
	x &= a z + \epsilon' \cr
	\Rightarrow y &= b_0 (a z + \epsilon') + \epsilon_0 \cr
	&= b_0 a z + b_0 \epsilon' + \epsilon_0 \cr
	&:= b_1 z + \epsilon_1 \text{, where $\epsilon_1 = b_0 \epsilon' + \epsilon_0$}
\end{aligned}$$</div>

> , from which it is easy to see that `$\cov(\epsilon_{(i)}, z_{(0)}) \neq 0$`

From the result of $\eqref{eq:1}$, we have:

<div>$$\begin{align}
	b_{xy(i)} &= \frac{b_{zy(i)}}{\beta_{zx(i)}} \nocr
	&= \frac{b_{zy(0)}r_{i0}\sqrt{h_0 / h_i}}{\beta_{zy(0)}\gamma_{i0}\sqrt{\eta_0 / \eta_i}} \nocr
	&= \frac{b_{zy(0)}}{\beta_{zy(0)}} \label{eq:2}\cr
	&= b_{xy(0)} \nocr
\end{align}$$</div>

, where \eqref{eq:2} follows from the fact that `$r_{0i}$` and `$\gamma_{0i}$` are simply the property of the locus 0 and locus i, so that `$r_{0i} = \gamma_{0i}$`. The same logic follows for $h$ and $\eta$.

This result implies that `$\hat{d}_{i} := \hat{b}_{(i)} - \hat{b}_{(0)}$` has mean zero. Then, under the null, `$(\hat{d}_{1}, ..., \hat{d}_k) \approx \mathcal{N}(0, R)$` (as stated in the text page 9 top left). $R$ is the variance/covariance matrix with entry `$\cov(\hat{d}_i, \hat{d}_j)$`.

<div>$$\begin{align}
  \cov(x - y, z - y) &= \E((x - y)(z - y)) - \E(x - y)\E(z - y) \nocr
	&= \E((x - y)(z - y)) - [\E(x) - \E(y)][\E(z) - \E(y)] \nocr
	&= \E(xz) - \E(x)\E(z) - [\E(yz) - \E(y)\E(z)] \nocr
	&- [\E(xy) - \E(x)\E(y)] + \E(y^2) - \E(y)^2 \nocr
	&= \cov(x, z) - \cov(y, z) - \cov(x, y) + \var(y) \nocr
	\cov(\hat{d}_i, \hat{d}_j) &= \cov(\hat{b}_{(i)}, \hat{b}_{(j)}) - \cov(\hat{b}_{(i)}, \hat{b}_{(0)}) \nocr
	&- \cov(\hat{b}_{(j)}, \hat{b}_{(0)}) + \var(\hat{b}_{(0)}) \label{eq:d}
\end{align}$$</div>

From \eqref{eq:d} (as stated in equation 8), we need to compute `$\cov(\hat{b}_{(i)}, \hat{b}_{(j)}), \forall i, j = 0, 1, ..., k$`. The derivation of this term is sketched in the supplement part 3. The missing part seems unclear for me, so I derive it great details as follow.

The unclear thing is how to get $\E(g(\hat\theta))$ using Delta method. For $\var(g(\hat\theta))$, it is straight forward, but not so for the first order term because $\E(g(\hat\theta)) \approx g(\theta)$ is what we commonly use in practice. But if we take a closer look at Delta method, we will be sure that we can do more than this (see [post](add url) about Delta method). That is to use higher order approximation.

<div>$$\begin{align}
 	g(X) &= g(\mu) + g'(\mu)(X - \mu) + g''(\mu) \frac{(X - \mu)^2}{2!} + o((X - \mu)^2) \nocr
	\E(g(X)) &= g(\mu) + \E(o(X - \mu)) \label{eq:3}\cr
	\E(g(X)) &= g(\mu) + \frac{g''(\mu)}{2} \E((X - \mu)^2) + \E(o((X - \mu)^2)) \nocr
	&\approx g(\mu) + \frac{g''(\mu)}{2} \var(X) \label{eq:4}
\end{align}$$</div>

\eqref{eq:3} is commonly used approximation and \eqref{eq:4} is a more "accurate" one. The corresponding multivariate version is simply:

<div>$$\begin{align}
  g(X) &\approx g(\mu) + \frac{1}{2} \langle H(g(\mu)), \Sigma_X \rangle \nocr
\end{align}$$</div>

, where $H(g(X))$ is the Hessian of $g$ evaluated at $\mu$ and `$\Sigma_X$` is the variance/covariance matrix of $X$ and $\langle \cdot, \cdot \rangle$ is the matrix inner product. Now, we can apply this rule for deriving `$\cov(\hat{b}_{(i)}, \hat{b}_{(j)})$`.

<div>$$\begin{align}
  \cov(\hat{b}_{(i)}, \hat{b}_{(j)}) &= \E(\hat{b}_{(i)} \hat{b}_{(j)}) - \E(\hat{b}_{(i)}) \E(\hat{b}_{(j)}) \nocr
	\E(\hat{b}_{(i)}) &= \E(\hat{b}_{xy(i)}) = \E(\frac{\hat{b}_{zy(i)}}{\hat\beta_{(zx(i))}}) \nocr
	&\approx \frac{b_{zy(i)}}{\beta_{zx(i)}} + \frac{1}{2} \langle
	\begin{bmatrix} 0 & -1 / \beta_{zx(i)}^2 \cr
	-1 / \beta_{zx(i)}^2 & 2b_{zy(i)} / \beta_{zx(i)}^2 \end{bmatrix}
	, \nocr
	& \begin{bmatrix} \var(\hat{b}_{zy(i)}) & \cov(\hat{b}_{zy(i)}, \hat{\beta}_{zy(i)}) \cr
	\cov(\hat{b}_{zy(i)}, \hat{\beta}_{zy(i)}) & \var(\hat{\beta}_{zy(i)}) \end{bmatrix}
	\rangle \nocr
	&= \frac{b_{zy(i)}}{\beta_{zx(i)}} - \frac{\cov(\hat{b}_{zy(i)}, \hat{\beta}_{zy(i)})}{\beta_{zx(i)}^2} + \frac{b_{zy(i)}\var(\hat{\beta}_{zy(i)})}{\beta_{zx(i)}^2} \nocr
	&= \frac{b_{zy(i)}}{\beta_{zx(i)}} \bigg( 1 + \frac{\var(\hat{\beta}_{zy(i)})}{\beta_{zx(i)}} - \frac{\cov(\hat{b}_{zy(i)}, \hat{\beta}_{zy(i)})}{b_{zy(i)}\beta_{zx(i)}} \bigg) \label{eq:6}
\end{align}$$</div>

\eqref{eq:6} is stated in supplement page 10. For an association signal, `$\frac{\hat\beta}{\var(\hat\beta)}$` should be big (namely the $\chi^2$ statistic). Also, `$\cov(\hat{b}_{zy(i)}, \hat{\beta}_{zy(i)}) = 0$`. Therefore,

<div>$$\begin{align}
  \E(\hat{b}_{xy(i)}) &\approx \frac{b_{zy(i)}}{\beta_{zx(i)}} \nonumber
\end{align}$$</div>

This result gives nothing more than the first order approximation. But it does matter for `$\E(\hat{b}_{xy(i)}\hat{b}_{xy(j)})$`. The intuition is that as more and more terms get involved, the first order approximation becomes worse and worse. So, in general, when too many terms involved, you need to be careful. If the (co)variance or Hessian is crazy, maybe it is a good idea to go beyond first order approximation.

It turns out that we need to compute `$\cov(\hat{b}_{zy(i)}, \hat{b}_{zy(j)})$` and `$\cov(\hat\beta_{zx(i)}, \hat\beta_{zx(j)})$`. They are:
