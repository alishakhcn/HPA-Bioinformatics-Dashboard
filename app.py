import streamlit as st
import pandas as pd
import plotly.express as px
import numpy as np

from utils import *

# ==========================================================
# PAGE CONFIG
# ==========================================================
st.set_page_config(
    page_title="Human Protein Atlas Dashboard",
    page_icon="🧬",
    layout="wide"
)

# ==========================================================
# LOAD DATA
# ==========================================================
tissue, cancer = load_data()

# ==========================================================
# HEADER
# ==========================================================
st.title("🧬 Human Protein Atlas Gene Expression Dashboard")

st.markdown("""
### Comparative Transcriptomic Analysis using Human Protein Atlas (HPA)

This dashboard enables interactive comparison of transcriptomic gene expression (nTPM)
between **normal human tissue** and **kidney cancer cell lines** using Human Protein Atlas data.
""")

st.divider()

# ==========================================================
# SIDEBAR
# ==========================================================


st.sidebar.header("🧬 Dashboard Filters")

# Keep only the tissues required for the project
allowed_tissues = [
    "breast",
    "colon",
    "kidney",
    "liver",
    "lung",
    "pancreas",
    "prostate",
    "skin"
]

# Filter the dataframe to these tissues only
tissue = tissue[tissue["Tissue"].isin(allowed_tissues)].copy()

selected_tissue = st.sidebar.selectbox(
    "Select Normal Tissue",
    sorted(tissue["Tissue"].unique())
)
normal = tissue[tissue["Tissue"] == selected_tissue]

genes = sorted(normal["Gene name"].dropna().unique())

selected_gene = st.sidebar.selectbox(
    "Select Gene",
    genes
)

top_n = st.sidebar.slider(
    "Number of Top Biomarkers",
    min_value=5,
    max_value=50,
    value=20,
    step=5
)

# ==========================================================
# BIOMARKERS
# ==========================================================
bio = biomarker_table(
    normal,
    cancer
)

# ==========================================================
# KPI CARDS
# ==========================================================
c1, c2, c3, c4 = st.columns(4)

c1.metric(
    "🧬 Total Genes",
    normal["Gene name"].nunique()
)

c2.metric(
    "🧫 Kidney Cell Lines",
    cancer["Cell line"].nunique()
)

c3.metric(
    "📍 Selected Tissue",
    selected_tissue
)

c4.metric(
    "🏆 Candidate Biomarkers",
    bio["Gene name"].nunique()
)

st.divider()
# ==========================================================
# GENE EXPRESSION COMPARISON
# ==========================================================

st.subheader("🧬 Gene Expression Comparison")

# Get selected gene data
normal_gene = normal[
    normal["Gene name"] == selected_gene
]

cancer_gene = cancer[
    cancer["Gene name"] == selected_gene
]


if not normal_gene.empty and not cancer_gene.empty:

    # Normal tissue value (single average value)
    normal_value = normal_gene["nTPM"].mean()

    # Cancer cell line values
    cancer_values = cancer_gene["nTPM"]


    # Create dataframe for plotting
    plot_df = pd.DataFrame({

        "Group": 
        ["Normal Kidney"] + 
        ["Kidney Cancer"] * len(cancer_values),

        "Expression": 
        [normal_value] + 
        list(cancer_values)

    })


    # Plot
    fig = px.strip(
        plot_df,
        x="Group",
        y="Expression",
        color="Group",
        stripmode="overlay",
        title=f"{selected_gene} Expression: Normal vs Kidney Cancer"
    )


    # Add normal average marker
    fig.add_scatter(
        x=["Normal Kidney"],
        y=[normal_value],
        mode="markers",
        marker=dict(
            size=18,
            symbol="diamond",
            color="orange"
        ),
        name="Normal Expression"
    )


    # Add cancer average line
    cancer_avg = cancer_values.mean()

    fig.add_scatter(
        x=["Kidney Cancer"],
        y=[cancer_avg],
        mode="markers",
        marker=dict(
            size=18,
            symbol="diamond",
            color="red"
        ),
        name="Cancer Average"
    )


    fig.update_layout(
        xaxis_title="Sample Type",
        yaxis_title="Expression (nTPM)",
        showlegend=True
    )


    st.plotly_chart(
        fig,
        use_container_width=True
    )


    # Simple interpretation
    st.write("### 📌 Interpretation")

    difference = cancer_avg - normal_value

    if difference > 0:
        st.success(
            f"{selected_gene} shows higher expression in kidney cancer "
            f"(+{difference:.2f} nTPM)"
        )

    else:
        st.info(
            f"{selected_gene} shows lower expression in kidney cancer "
            f"({difference:.2f} nTPM)"
        )


else:
    st.warning(
        "Expression data not available for this gene."
    )
# ==========================================================
# TOP BIOMARKERS
# ==========================================================
st.subheader("🔥 Top Differentially Expressed Biomarkers")

top = top_biomarkers(
    normal,
    cancer,
    top_n
)

fig = px.bar(

    top,

    x="Difference",

    y="Gene name",

    orientation="h",

    color="Fold Change",

    title=f"Top {top_n} Biomarkers"

)

fig.update_layout(
    yaxis=dict(
        autorange="reversed"
    )
)

st.plotly_chart(
    fig,
    use_container_width=True
)

st.divider()

# ==========================================================
# HEATMAP
# ==========================================================
st.subheader("🌡 Heatmap of Top Biomarkers")

heat = heatmap_data(
    cancer,
    normal,
    top_n
)

fig2 = px.imshow(

    heat,

    aspect="auto",

    color_continuous_scale="Viridis",

    labels=dict(
        color="nTPM"
    )

)

st.plotly_chart(
    fig2,
    use_container_width=True
)

st.divider()


# ==========================================================
# HISTOGRAM - Gene Expression Distribution
# ==========================================================

st.subheader("📈 Distribution of Kidney Cancer Gene Expression")

# Remove zero expression values for log scale
cancer_positive = cancer[cancer["nTPM"] > 0].copy()

# Log transform expression values
cancer_positive["log_nTPM"] = np.log2(cancer_positive["nTPM"] + 1)

fig3 = px.histogram(
    cancer_positive,
    x="log_nTPM",
    nbins=50,
    title="Log2 Transformed Gene Expression Distribution",
    labels={
        "log_nTPM": "Log2(nTPM + 1)",
        "count": "Number of Genes"
    },
    color_discrete_sequence=["royalblue"]
)

fig3.update_layout(
    xaxis_title="Log2(nTPM + 1)",
    yaxis_title="Number of Genes",
    bargap=0.1
)

st.plotly_chart(
    fig3,
    use_container_width=True
)

st.divider()
# ==========================================================
# BIOMARKER TABLE
# ==========================================================
st.subheader("🧬 Biomarker Ranking")

st.dataframe(

    bio[
        [
            "Gene name",
            "Normal nTPM",
            "Cancer nTPM",
            "Difference",
            "Fold Change"
        ]
    ].style.background_gradient(

        subset=["Difference"],

        cmap="Blues"

    ),

    use_container_width=True,

    height=500

)

st.download_button(

    "📥 Download Biomarker CSV",

    bio.to_csv(index=False),

    "kidney_biomarkers.csv",

    mime="text/csv"

)

st.divider()

# ==========================================================
# GENE INTERPRETATION
# ==========================================================
stats = gene_statistics(
    normal,
    cancer,
    selected_gene
)

if stats:

    st.subheader("🧠 Biological Interpretation")

    st.info(f"""

**Selected Gene : {selected_gene}**

• Normal Expression : **{stats['Normal']} nTPM**

• Cancer Expression : **{stats['Cancer']} nTPM**

• Fold Change : **{stats['Fold Change']}**

• Status : **{stats['Status']}**

This comparison highlights potential transcriptomic differences between normal tissue and kidney cancer cell lines and may assist in identifying candidate biomarkers for further biological investigation.

""")

st.divider()

# ==========================================================
# DATASET PREVIEW
# ==========================================================
with st.expander("📄 Preview Normal Tissue Dataset"):

    st.dataframe(
        normal.head(10)
    )

with st.expander("📄 Preview Kidney Cancer Cell Line Dataset"):

    st.dataframe(
        cancer.head(10)
    )

st.divider()

# ==========================================================
# ABOUT
# ==========================================================
st.subheader("ℹ️ About This Dashboard")

st.write("""
This dashboard was developed using **Python**, **Streamlit**, **Pandas**, and **Plotly**
to visualize Human Protein Atlas (HPA) transcriptomic data.

It enables users to compare normal tissue expression with kidney cancer cell lines,
explore differentially expressed genes, identify candidate biomarkers,
and interpret transcriptomic expression patterns through interactive visualizations.
""")

st.divider()

st.caption(
"Developed by **Alisha Khan (2409036)** | MSc Big Data Analytics II | Bioinformatics Project | Human Protein Atlas Dashboard"
)