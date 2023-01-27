Execute in this order:

0. Overdispersion.ipynb (independent, can be run anytime)
1. Quality_Control.ipynb 
2. Normalizations.ipynb 
3. SCTransform.ipynb 
4. SCTransform_Custom.ipynb 
5. binaryResidualTransform.ipynb 
6. metrics.ipynb

If disk+RAM is less than 100GB, the kernel may crash, deleting intermediate h5ad files may necessary to not run out of space (and cause the kernel crash)
