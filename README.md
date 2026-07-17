# 🧬 HPA Bioinformatics Dashboard

An interactive Streamlit dashboard for exploring Human Protein Atlas (HPA) gene expression data. The dashboard allows users to compare transcriptomic expression levels (nTPM) between normal human tissues and cancer cell lines, identify potential biomarkers, and visualize gene expression patterns.

## Features

- 📊 Interactive gene expression analysis
- 🧬 Compare normal tissue vs cancer cell lines
- 📈 Histogram of gene expression
- 🏆 Top biomarker identification
- 🔍 Search and filter by tissue and gene
- 📉 Gene comparison bar chart
- 📋 Biomarker summary table
- ⚡ Fast interactive dashboard built with Streamlit

## Dataset

Data obtained from the Human Protein Atlas (HPA).

Files used:
- rna_tissue_consensus.tsv
- rna_celline.tsv
- rna_celline_description.tsv

## Tech Stack

- Python
- Streamlit
- Pandas
- NumPy
- Plotly
- Git & GitHub

## Dashboard Preview

(Add screenshots here)

## Installation

```bash
git clone https://github.com/yourusername/hpa-bioinformatics-dashboard.git
cd hpa-bioinformatics-dashboard
pip install -r requirements.txt
streamlit run app.py
