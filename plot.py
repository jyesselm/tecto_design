import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict

sns.set(style="white", context="talk", font_scale=1.2)

def times_used(row):
    spl = row["designs_used"].split(";")
    return len(spl)-1


def get_motif_usage(df):
    motifs = defaultdict(str)
    for i, row in df.iterrows():
        motifs_used = row.motifs_uses.split(";")
        for m_name in motifs_used[:-1]:
            if m_name[0] == "H" or m_name[0:1] == "BP":
                continue
            motifs[m_name] += str(row.design_num) + ";"
    motif_usage_df = pd.DataFrame(motifs.items(), columns=("m_name", "designs_used"))
    mtimes_used = motif_usage_df.apply(lambda x: times_used(x), axis=1)
    motif_usage_df['times_used'] = mtimes_used
    return motif_usage_df


wdir = "20170101_11bp_constructs/"

df = pd.read_csv(wdir+"data/bp11_problem.scores")


sns.distplot(df["design_score"])
motif_usage_df = get_motif_usage(df)
plt.savefig(wdir+"/plots/design_scores.png")
plt.clf()

plt.figure()
g = sns.barplot(x="m_name", y="times_used", data=motif_usage_df)
plt.setp(g.get_xticklabels(), rotation=90)
plt.tight_layout()
plt.savefig(wdir+"/plots/motifs_used.png")
plt.clf()



#plt.show()