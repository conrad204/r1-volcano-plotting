# Volcano plotting (plumber)

this repo uses R plumber API to generate a volcano plot based off `example.csv`
currently, there is no implementation to upload custom dataset (todo)

**NOTE: while *padj* and *log₂ fold change* values reflect real experimental results, ensemblIDs have been randomly reassigned to unrelated entries to preserve integrity of unpublished data.**

apart from the original R script, initial commit for plumber API was primarily generated using ai

## quick start and docker build

you'll need docker installed

PowerShell (from project root):

```powershell
docker pull cfeng204/volcano-api:latest
docker run -p 8000:8000 cfeng204/volcano-api:latest
```

then open the index.html file

## uploading csv files

if you want to try out your own dataset, upload any .csv with the following columns: 

1. EnsemblID (ENSMUSG only) OR GeneName
2. log2FoldChange
3. padj
NOTE: extra columns are fine and will be ignored

## example data: 

**for ensemblID:**

| EnsemblID              | baseMean     | FoldChange   | log2FoldChange | lfcSE       | stat        | pvalue     | padj        |
|-----------------------|--------------|--------------|----------------|-------------|-------------|------------|-------------|
| ENSMUSG00000102151.1  | 3661.785238  | -1.042306356 | -0.059779379   | 0.086325349 | -0.692489289| 0.48863011 | 0.753275158 |
| ENSMUSG00000031375.17 | 73.83677932  | 1.02347041   | 0.033469393    | 0.236027536 | 0.141802915 | 0.887235681| 0.957484933 |
| ENSMUSG00000102817.2  | 1.451684654  | -1.082852719 | -0.114837032   | 1.517656734 | -0.075667329| 0.939683769| NA          |
| ENSMUSG00000026739.13 | 12.87340559  | 1.136545422  | 0.184655342    | 0.49498248  | 0.373054299 | 0.709108031| 0.880737776 |
| ENSMUSG00000018501.17 | 283.7034353  | -1.21858013  | -0.285201122   | 0.250009298 | -1.140762059| 0.253968953| 0.543772989 |


**for GeneName:**

| EnsemblID           | GeneName | GeneType       | baseMean     | log2FoldChange | lfcSE       | stat        | pvalue      | padj        |
|---------------------|----------|---------------|--------------|----------------|-------------|-------------|-------------|-------------|
| ENSRNOG00000000001  | Arsj     | protein_coding| 20.65466513  | -0.585591653   | 0.353336125 | -1.657321772| 0.097454438 | 0.268000159 |
| ENSRNOG00000000007  | Gad1     | protein_coding| 5.11998274   | 0.006560455    | 0.666815683 | 0.009838484 | 0.992150152 | NA          |
| ENSRNOG00000000012  | Tcf15    | protein_coding| 73.74778968  | -0.050829959   | 0.289290583 | -0.175705542| 0.860525281 | 0.93504675  |
| ENSRNOG00000000017  | Steap1   | protein_coding| 390.6887143  | 0.543779638    | 0.142461801 | 3.817020647 | 0.000135073 | 0.003668259 |
| ENSRNOG00000000024  | Hebp1    | protein_coding| 198.5599409  | 0.38480401     | 0.16437658  | 2.340990479 | 0.019232658 | 0.099247946 |
