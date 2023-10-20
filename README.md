# Anomaly_Gene_Visualizer

**R Dependencies**
"tidyr","tidyverse", "data.table", "reshape2", "colorspace", "readxl"

**Usage**
RScript Boxplot_Automation.R \
-i  <normalized_data>.csv \
-gs <gene_states>.csv \
-gp <gene_panel>.xlsx \
-gpn "gene panel name" \
-n normal_count \
-t "tumor sample name" \
-p "title for your boxplot" \
-st state_type


**where** :
-i   = normalized or scaled data. It is an output from the secodary analysis tool. <project>-scaled.csv. Where first two columns are ENSEMBL ID and GENE IDs
-gs  = gene state data. It is also an output from the secondary analysis tool. <project>-gene_states.csv
-gp  = gene panel. It should be an excel sheet that has the list of genes that need to be plotted.
-gpn = The excel sheet name that should be used. You can have a same gene panel file with multiple gene lists in different sheets and just change this option as needed.
-n   = Number of normal sample present in the data. These columns will be used to make boxplots.
-t   = Name of the Tumor sample that needs to be plotted on the boxplots. i.e. The tumor column name. Check Examples.docs tool for visualization.
-p   = Provide a Title for your boxplot. This can be kept NULL if needed
-st  = state type. Available options ("ALL", "VERY HIGH and HIGH", "VERY HIGH", "HIGH", "LOW", "VERY LOW", "VERY LOW and LOW"). 
The state type (-st) option will only select genes from the gene panel that are of that state type in the gene state data. Check Example.docs for visualization.

**Example Command**
RScript Boxplot_Automation.R \
-i  test1-scaled.csv \
-gs test1-gene_states.csv \
-gp gene_panels.xlsx \
-gpn Apoptosis
-n 78 \
-t X91_031323_B \
-p Apoptosis-Genes-Aberrantly-Expressed-as-VERY-HIGH-and-HIGH-for-Tumor-Sample-X91_031323_B
-st VERY\ HIGH\ and\ HIGH
