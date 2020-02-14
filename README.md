# Deep Gaussian Mixture Models with Common Factor Loadings

A model-based method enabling dimensionality reduction, unsupervised clustering, and visualisation.

See example.R for more details on implementation.

This is a culmination of my Thesis project for the Master in Science (Statistics), supervised by Professor Geoffrey McLachlan.
Note, the thesis attached is a preliminary copy. Updates will be made as time progresses.

It is based on the following work, so please reference the following:

https://github.com/suren-rathnayake/deepgmm/blob/master/

https://github.com/suren-rathnayake/EMMIXmfa

[Viroli, C. and McLachlan, G.J. (2019). Deep Gaussian mixture models. Statistics and Computing 29, 43-51.](https://link.springer.com/article/10.1007/s11222-017-9793-z)

[Baek J., McLachlan G. J., Flack L. K. (2010). Mixtures of Factor Analyzers with Common Factor Loadings: Applications to the Clustering and Visualization of High-Dimensional Data. IEEE Transactions on Pattern Analysis and Machine Intelligence 32, 7, 1298-1309.](https://ieeexplore.ieee.org/document/5184847)

It compares well on small datasets, but on larger datasets, tSNE and UMAP could outperform it. I have tested it on a few public single-cell datasets, and so, this could still be useful as a complementary tool.
![](goolam_et_al.png)
