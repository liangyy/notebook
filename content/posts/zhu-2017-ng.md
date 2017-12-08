---
title: "Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets"
date: 2017-12-08T14:23:10-06:00
draft : false
author: "Yanyu Liang"
tags: ["gwas", 'integrative analysis', 'eqtl', 'target gene']
categories: ["research paper"]
---

# Meta data of reading

* **Journal**:

Nature Genetics

* **Year**:

2016

* **DOI**:

10.1038/ng.3538


# The idea of instrumental variable

The paper mentioned it in the background section along with Mendelian randomization (MR). I googled this term and found this [site](http://www.statisticshowto.com/instrumental-variable/). Also, [wiki page](https://en.wikipedia.org/wiki/Instrumental_variables_estimation) is informative.

It appears in regression analysis where $y$ is the responsive variable and $x$ is the explanatory variable. The model is simply $Y = X \beta + \epsilon$. Suppose $Z$ correlates with $X$ and $Z$ is independent to $Y$ given $X$, then $Z$ correlates with $Y$ as well despite such conditional independence.

In some analysis, people sort of take the advantage of $Z$ (in this context $Z$ is called instrumental variable).


# Motivation
