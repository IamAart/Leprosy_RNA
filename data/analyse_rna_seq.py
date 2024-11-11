import pandas as pd


def check_count(df):
    nc_df = df[df["biotype"].isin(NON_CODING_BIOTYPES)]
    c_df = df[df["biotype"].isin(CODING_BIOTYPES)]
    print("Total count: ", df.shape[0])
    print("Coding count: ", c_df.shape[0])
    print("Non Coding count: ", nc_df.shape[0])

    return df.shape[0] - nc_df.shape[0] - c_df.shape[0]

NON_CODING_BIOTYPES = ["processed_transcript", "ribozyme", "unitary_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "transcribed_unprocessed_pseudogene", "antisense", "transcribed_unitary_pseudogene", "polymorphic_pseudogene", "lincRNA", "sense_intronic", "transcribed_processed_pseudogene", "sense_overlapping", "IG_V_pseudogene", "pseudogene", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "snRNA", "miRNA", "misc_RNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "TR_V_pseudogene", "TR_J_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "scRNA", "scaRNA", "vaultRNA", "sRNA", "macro_lncRNA", "non_coding", "IG_pseudogene"]
CODING_BIOTYPES = ["IG_D_gene", "protein_coding", "TR_V_gene", "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene", "TR_D_gene", "TEC"]

paths = ["./DESeq2/All_GENES_DESeq2_RNA_SEQ.csv", "./LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.csv", "./EdgeR/All_GENES_EdgeR_RNA_SEQ.csv"]

for path in paths:
    print(path.split("/")[1])
    df = pd.read_csv(path)
    count = check_count(df)
    print("Remaining: ", count)