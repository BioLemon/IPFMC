# IPFMC

## Usage of IPFMC

In this section, we will introduce how to use our python package published on PyPI to perform clustering and biological interpretation of cancer multi-omics data. We will show the process using the LUAD cancer datasets as an example.

### Prepare datasets

1. **Omics datasets**

   We apologize for the inconvenience caused by the large size of the cancer omics data file, which prevents us from uploading it to the github repository. Therefore, we provide the source of the omics data used in our experiments and show the downloading method using the omics data of LUAD cancer as an example. Please follow our guidance to download the LUAD cancer datasets.

   All our omics data are derived from the data provided by Duan et al. in their review published in PCB in 2021([Evaluation and comparison of multi-omics data integration methods for cancer subtyping | PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009224)). Their datasets can be downloaded here: [GaoLabXDU/MultiOmicsIntegrationStudy: Experimental evaluation and comparison of multi-omics data integration methods for cancer subtyping (github.com)](https://github.com/GaoLabXDU/MultiOmicsIntegrationStudy/)

2. **Pathway datasets**

   IPFMC is a multi-omics integrated clustering method based on iterative fusion of pathways, and we provide information on the correspondence between processed pathways and genes and the correspondence between pathways and mirnas that we used in our experiments. These files can be found in ‘IPFMC/Sample_Files’ folder.

### Import neccessary packages

```python
import pandas as pd
import numpy as np
from snf import snf
import IPFMC
```
