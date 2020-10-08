import pandas as pd
def open_exon_annotations(exons_file):
    # exon annotations:
    exons = pd.read_csv(exons_file, sep="\t", names=["Chr","Start","End","Exon"])
    # TODO: define this list in a better place?
    important_genes = ["ERBB2","MET","KRAS","EGFR","CCNE1","KIT","PTEN"]
    # filter important genes:
    exons = exons[exons.Exon.str.replace("_.+", "").isin(important_genes)].copy()
    exons.Chr = exons.apply(lambda r:"chr%d" % r.Chr, axis=1)
    return exons

# Find exons given a range:
def exons_in_range(exons, ch, start, end):
    ex_list = []
    for idx, r in exons[(exons.Chr==ch)].iterrows():
        if start <= r.End and end >= r.Start:
            ex_list.append(r.Exon)
    return ex_list
