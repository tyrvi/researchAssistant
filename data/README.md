expected_count
--------------
The maximum likelihood gene expression levels estimated using RSEM, i.e. the expected_count in RSEM’s output. There are 52 data files in this subdirectory, each being a sample-gene matrix of a certain tissue type. These files can be provided to programs such as EBSeq, DESeq, or edgeR for identifying differentially expressed genes. The expression of only protein coding genes was provided here.


unnormalized
--------------
The gene expression levels calculated from fpkm of RSEM’s output. The data matrices in this subdirectory were not the direct output of RSEM. They underwent quantile normalization, but were not corrected for batch effects.


normalized
--------------
The normalized gene expression levels (FPKM) calculated using RSEM. This set of data files not only was quantile normalized, but was corrected for batch effects (using tool ComBat).
