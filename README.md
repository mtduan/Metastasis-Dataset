# Cancer Metastasis Data Repositories

## Breast Cancer:
### GSE167036 (8 datasets):
- **Main path:** `/fs/ess/PAS1475/Maoteng/Metastasis/GSE167036/Patient_subsample`
- **Subpaths:** `Patient1-Patient8`
- **barcodes.csv:** Cell names
- **features.csv:** Gene names
- **label.csv:** Cell labels (CA for Cancer - primary tumor site, LN for Lymph Node - metastatic site)
- **RNA_matrix.mtx:** Gene*Cell expression matrix
- **Marker Genes Dictionary:**
    B Cells: CD79A, CD3D, MS4A1
    T Cells: CD3D, CD3E
    CD4 T Cells: CD4, CD3D
    CD8 T Cells: CD8A, CD3D
    Treg Cells: FOXP3, CTLA4, TNFRSF18
    Cytotoxic T Cells: CD8A, CD8B, GZMB, PRF1
    TAM Cells: CD68, CD163, MRC1
    NK Cells: GNLY, NKG7, XCL2
    Epithelial Cells: EPCAM, KRT19, KRT7, KRT17
    Cancer-associated Fibroblasts (CAFs): ACTA2, FAP, PDGFRA
    PVL Cells: RGS5, ANGPT1
    Mast Cells: TPSD1, CST3, CPA3, KIT
    Plasma Cells: IGHG1, MZB1, SDC1, CD79A
    Endothelial Cells: PECAM1, VWF, CDH5
  
### GSE180286 (10 datasets):
- **Main path:** `/fs/ess/PAS1475/Maoteng/Metastasis/GSE180286/All_patient`
- **Subpaths:** `patientA2019_2- patientE2020_3`
- **cell_barcodes.csv:** Cell names
- **gene_symbol.csv:** Gene names
- **cell_label.csv:** Cell labels (P for primary tumor site, M for metastatic site)
- **gene_cell.mtx:** Gene*Cell expression matrix

## Head and Neck Cancer:
### GSE188737 (7 datasets):
- **Main path:** `/fs/ess/PAS1475/Maoteng/Metastasis/GSE167036/Patient_subsample`
- **Subpaths:** `HN237- HN279`
- **cell_barcodes.csv:** Cell names
- **gene_symbol.csv:** Gene names
- **cell_label.csv:** Cell labels (P for primary tumor site, M for metastatic site)
- **gene_cell.mtx:** Gene*Cell expression matrix

## Pancreatic Cancer:
### GSE197177 (3 datasets):
- **Main path:** `/fs/ess/PAS1475/Maoteng/Metastasis/GSE197177/Patient_subsample`
- **Subpaths:** `Patient1-Patient3`
- **cell_barcodes.csv:** Cell names
- **gene_symbol.csv:** Gene names
- **cell_label.csv:** Cell labels (P for primary tumor site, M for metastatic site)
- **gene_cell.mtx:** Gene*Cell expression matrix
