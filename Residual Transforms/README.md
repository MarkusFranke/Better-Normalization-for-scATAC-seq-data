Execute in this order:

0. Overdispersion.ipynb (independent, can be run anytime)
1. metrics.py (compile before importing to metrics.ipynb)
2. Quality_Control.ipynb 
3. Normalizations.ipynb 
4. SCTransform.ipynb 
5. SCTransform_Custom.ipynb 
6. binaryResidualTransform.ipynb 
7. metrics.ipynb

If disk+RAM is less than 100GB, the kernel may crash, deleting intermediate h5ad files may necessary to not run out of space (and cause the kernel crash)
