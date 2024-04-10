# Cancer Metastasis Data Repositories

## Breast Cancer:
### GSE167036 (8 datasets):
- **Main path:** `/fs/ess/PAS1475/Maoteng/Metastasis/GSE167036/Patient_subsample`
- **InferCNV path:** `/fs/scratch/PAS1475/dmt/Benchmark/GSE167036/gene_locations.csv`
- **Subpaths:** `Patient1-Patient8`
- **barcodes.csv:** Cell names
- **features.csv:** Gene names
- **label.csv:** Cell labels (CA for Cancer - primary tumor site, LN for Lymph Node - metastatic site)
- **RNA_matrix.mtx:** Gene*Cell expression matrix
- **Marker Genes Dictionary:**
  - B Cells: `CD79A`, `CD3D`, `MS4A1`
  - T Cells: `CD3D`, `CD3E`
  - CD4 T Cells: `CD4`, `CD3D`
  - CD8 T Cells: `CD8A`, `CD3D`
  - Treg Cells: `FOXP3`, `CTLA4`, `TNFRSF18`
  - Cytotoxic T Cells: `CD8A`, `CD8B`, `GZMB`, `PRF1`
  - TAM Cells: `CD68`, `CD163`, `MRC1`
  - NK Cells: `GNLY`, `NKG7`, `XCL2`
  - Epithelial Cells: `EPCAM`, `KRT19`, `KRT7`, `KRT17`
  - Cancer-associated Fibroblasts (CAFs): `ACTA2`, `FAP`, `PDGFRA`
  - PVL Cells: `RGS5`, `ANGPT1`
  - Mast Cells: `TPSD1`, `CST3`, `CPA3`, `KIT`
  - Plasma Cells: `IGHG1`, `MZB1`, `SDC1`, `CD79A`
  - Endothelial Cells: `PECAM1`, `VWF`, `CDH5`
  
### GSE180286 (10 datasets):
- **Main path:** `/fs/ess/PAS1475/Maoteng/Metastasis/GSE180286/All_patient`
- **InferCNV path:** `/fs/scratch/PAS1475/dmt/Benchmark/GSE180286/(patientA2019_2-patientE2020_3)/gene_locations.csv`
- **Subpaths:** `patientA2019_2- patientE2020_3`
- **cell_barcodes.csv:** Cell names
- **gene_symbol.csv:** Gene names
- **cell_label.csv:** Cell labels (P for primary tumor site, M for metastatic site)
- **gene_cell.mtx:** Gene*Cell expression matrix
- **Marker Genes Dictionary:**
  - B Cells: `CD79A`, `MS4A1`
  - T Cells: `CD3D`, `CD3E`
  - CD4 T Cells: `CD4`, `CD3D`
  - CD8 T Cells: `CD8A`, `CD3D`
  - Treg Cells: `FOXP3`, `CTLA4`, `TNFRSF18`
  - Cytotoxic T Cells: `CD8A`, `CD8B`, `GZMB`, `PRF1`
  - Tumor-Associated Macrophages (TAM) Cells: `CD68`, `CD163`, `MRC1`
  - NK Cells: `GNLY`, `NKG7`, `XCL2`
  - Epithelial Cells: `EPCAM`, `KRT19`, `KRT7`, `KRT17`
  - Cancer-associated Fibroblasts (CAFs): `ACTA2`, `FAP`
  - PVL Cells: `RGS5`, `ANGPT1`
  - Mast Cells: `CST3`, `CPA3`, `KIT`
  - Plasma Cells: `IGHG1`, `MZB1`, `SDC1`, `CD79A`
  - Endothelial Cells: `PECAM1`, `CDH5`
  
## Head and Neck Cancer:
### GSE188737 (7 datasets):
- **Main path:** `/fs/ess/PAS1475/Maoteng/Metastasis/GSE167036/Patient_subsample`
- **InferCNV path:** `/fs/scratch/PAS1475/dmt/Benchmark/GSE188737/gene_locations.csv`
- **Subpaths:** `HN237- HN279`
- **cell_barcodes.csv:** Cell names
- **gene_symbol.csv:** Gene names
- **cell_label.csv:** Cell labels (P for primary tumor site, M for metastatic site)
- **gene_cell.mtx:** Gene*Cell expression matrix
- **Marker Genes Dictionary:**
  - Epithelial Cells: `KRT7`, `KRT17`, `EPCAM`, `KRT19`
  - Salivary Cells: `STATH`
  - Plasma Cells: `IGHG1`
  - Cancer-associated Fibroblasts (CAFs): `COL1A2`, `MMP2`, `PDGFRA`
  - NK Cells: `NKG7`, `XCL2`
  - B Cells: `CD79A`
  - T Cells: `CD3D`, `CD3E`
  - CD4 T Cells: `CD4`
  - CD8 T Cells: `CD8A`
  - Treg Cells: `FOXP3`, `CTLA4`, `TNFRSF18`
  - Cytotoxic T Cells: `CD8A`, `CD8B`, `GZMB`, `PRF1`
  - TAM Cells: `CD68`, `CD163`, `MRC1`
  
## Pancreatic Cancer:
### GSE197177 (3 datasets):
- **Main path:** `/fs/ess/PAS1475/Maoteng/Metastasis/GSE197177/Patient_subsample`
- **Subpaths:** `Patient1-Patient3`
- **cell_barcodes.csv:** Cell names
- **gene_symbol.csv:** Gene names
- **cell_label.csv:** Cell labels (P for primary tumor site, M for metastatic site)
- **gene_cell.mtx:** Gene*Cell expression matrix
- **Marker Genes Dictionary:**
  - Ductal Cells: `EPCAM`, `KRT19`, `KRT7`, `EPCAM`
  - T Cells: `CD3D`, `CD3E`
  - NK Cells: `GNLY`, `NKG7`, `KLRD1`
  - B Cells: `CD79A`, `MS4A1`
  - Mast Cells: `TPSAB1`
  - Plasma Cells: `MZB1`
  - Endothelial Cells: `PECAM1`
  - Fibroblasts: `DCN`, `VIM`, `FAP`, `PDPN`
  - Myeloid Cells: `CD68`
  - Acinar Cells: `PRSS1`
  - Endocrine Cells: `CHGA`
  - MKI67+ Cycling Ductal Cells: `MKI67`
  - Cancer-Associated Fibroblasts (CAFs): `ACTA2`, `COL3A1`
  - Lipid-Associated Macrophages (LAMs): `APOE`, `APOC1`, `FABP5`
  - Regulatory T Cells (Treg): `FOXP3`
  - T Follicular Helper Cells (TFH): `CXCL13`
  - CD8+ T Cells: `CD8A`, `CD8B`
