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

This paper proposed a summary data-based MR method (SMR). In intuitively, it estimates effect of $Z$ (genotype) on $Y$ (phenotype), $b\_{zy}$ and effect of $Z$ on $X$, $b\_{zx}$ (note that the former is GWAS and the latter is eQTL mapping). Then $b\_{xy}$ is simply $\frac{b\_{zy}}{b\_{zx}}$. Simulation showed that SMR is equivalent to MR if the confounding variables are non-genetic and causality is mediated by gene expression or both trait and gene expression (? not clear ...). This situation is referred as pleiotropy.

It appears to me that the derivation of MR (or instrumental variable regression) is not intuitive so I leave the note of MR part for another post [here](https://liangyy.github.io/notebook/posts/mendelian-randomization/#ivvar). 
