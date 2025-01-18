import pandas as pd
import os
from dotenv import load_dotenv
import openpyxl

load_dotenv()
NON_CODING_BIOTYPES = os.getenv("NON_CODING_BIOTYPES").split(",")
CODING_BIOTYPES = os.getenv("CODING_BIOTYPES").split(",")

deseq2 = pd.read_csv("./data/DESeq2/All_GENES_DESeq2_RNA_SEQ.csv")
edgeR = pd.read_csv("./data/EdgeR/All_GENES_EdgeR_RNA_SEQ.csv")
limma = pd.read_csv("./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.csv")

def get_division(df, name):
    not_de = df[df["col"] == "not"]
    not_de_c = not_de[not_de["biotype"].isin(CODING_BIOTYPES)].shape[0]
    not_de_nc = not_de[not_de["biotype"].isin(NON_CODING_BIOTYPES)].shape[0]
    up_de = df[df["col"] == "up"]
    up_de_c = up_de[up_de["biotype"].isin(CODING_BIOTYPES)].shape[0]
    up_de_nc = up_de[up_de["biotype"].isin(NON_CODING_BIOTYPES)].shape[0]
    down_de = df[df["col"] == "down"]
    down_de_c = down_de[down_de["biotype"].isin(CODING_BIOTYPES)].shape[0]
    down_de_nc = down_de[down_de["biotype"].isin(NON_CODING_BIOTYPES)].shape[0]
    total_c = up_de_c + down_de_c
    total_nc = up_de_nc + down_de_nc
    total = df.shape[0]

    
    return pd.DataFrame(
        [[name, up_de_c, up_de_nc, down_de_c, down_de_nc, total_c, total_nc, not_de_c, not_de_nc, total]], 
        columns=["DGE", "UP - CODING", "UP - NC", "DOWN - CODING", "DOWN - NC", "TOTAL DE - CODING", "TOTAL DE - NC", "NOT DE - CODING", "NOT DE - NC", "TOTAL"])

dge_method_dict = {"DESeq2": deseq2, "EdgeR": edgeR, "Limma": limma}
final_dataframe = pd.DataFrame(columns=["DGE", "UP - CODING", "UP - NC", "DOWN - CODING", "DOWN - NC", "TOTAL DE - CODING", "TOTAL DE - NC", "NOT DE - CODING", "NOT DE - NC", "TOTAL"])
for name, df in dge_method_dict.items():
    print(get_division(df, name))
    final_dataframe = pd.concat([final_dataframe, get_division(df, name)], axis=0)
final_dataframe.transpose().to_excel("./data/DEG_DIVISION.xlsx", index=False)