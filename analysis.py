import argparse
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import math
from collections import defaultdict

from rnamake import motif_graph, vienna

sns.set(style="white", context="talk", font_scale=1.2)

import file_manager

def times_used(row):
    spl = row["designs_used"].split(";")
    return len(spl)-1


def get_motif_usage(df):
    motifs = defaultdict(str)
    for i, row in df.iterrows():
        motifs_used = row['motifs_uses'].split(";")
        for m_name in motifs_used:
            if m_name[0] == "H" or m_name[0:1] == "BP":
                continue
            motifs[m_name] += str(row.design_num) + ";"
    motif_usage_df = pd.DataFrame(motifs.items(), columns=("m_name", "designs_used"))
    mtimes_used = motif_usage_df.apply(lambda x: times_used(x), axis=1)
    motif_usage_df['times_used'] = mtimes_used
    return motif_usage_df


def plot_design_motif_usage(df, fmanager):
    motif_usage_df = get_motif_usage(df)
    plt.figure()
    g = sns.barplot(x="m_name", y="times_used", data=motif_usage_df)
    plt.setp(g.get_xticklabels(), rotation=90)
    plt.tight_layout()
    plt.savefig(fmanager.design_results_path+"/plots/motif_used.png")
    plt.clf()


def plot_design_distros(df, fmanager):
    plt.figure()
    sns.distplot(df["design_score"])
    plt.savefig(fmanager.design_results_path+"/plots/design_score_distro.png")
    plt.clf()

    plt.figure()
    sns.distplot([len(x.split("&")[1]) for x in df["design_sequence"]])
    plt.savefig(fmanager.design_results_path+"/plots/length_distro.png")
    plt.clf()


def plot_seq_opt_distros(df, fmanager):
    plt.figure()
    sns.distplot(df["opt_score"])
    plt.savefig(fmanager.seq_opt_results_path+"best/plots/opt_score_distro.png")
    plt.clf()


def get_summary_df(df, norm_seq=None):
    new_df = df[df.apply(lambda x:  x['avg_hit_count'] != 0, axis=1)].copy()
    avg_count_wildtype = 1196
    dG_predicted = []
    for i, row in new_df.iterrows():
        prediction = 1.9872041e-3*298*math.log(float(avg_count_wildtype)/float(row['avg_hit_count']))
        dG_predicted.append(prediction)

    new_df.loc[:, 'dG_predicted'] = dG_predicted
    return new_df


def plot_simulation_distros(df, fmanager):
    plt.figure()
    sns.distplot(df["dG_predicted"])
    plt.savefig(fmanager.sim_results_path+"/plots/sim_best_dG_distro.png")
    plt.clf()

    plt.figure()
    sns.regplot(x=df["opt_score"], y=df["dG_predicted"], scatter_kws={"s": 100})
    plt.savefig(fmanager.sim_results_path+"/plots/sim_best_opt_score_vs_dG.png")
    plt.clf()


def plot_simulation_distros_alt(df, fmanager):
    plt.figure()
    sns.distplot(df["dG_predicted"],hist=False)
    plt.savefig(fmanager.sim_results_path+"/plots/sim_best_alt_dG_distro.png")
    plt.clf()

    plt.figure()
    sns.regplot(x=df["opt_score"], y=df["dG_predicted"], scatter_kws={"s": 100})
    plt.savefig(fmanager.sim_results_path+"/plots/sim_best_alt_opt_score_vs_dG.png")
    plt.clf()


def get_top_simulation_pdbs(df, fmanager):
    df_new = df.sort(['dG_predicted'], ascending=[1])
    for i, r in df_new.iterrows():
        if r.dG_predicted > 2:
            break
        d_num = int(r.design_num)
        o_num = int(r.opt_num)
        f = open(fmanager.seq_opt_results_path+"/best/"+str(d_num)+".out")
        lines = f.readlines()
        f.close()
        l = lines[o_num-1]
        mg = motif_graph.MotifGraph(l)
        mg.to_pdb(fmanager.seq_opt_results_path+"/best/pdbs/"+"top."+str(d_num)+"."+str(o_num)+".pdb",
                  renumber=1,
                  close_chain=1)


def plot_final_motif_usage(design_df, df, fmanager):
    motifs_used = []
    for i, r in df.iterrows():
        design_row = design_df.iloc[r.design_num]
        motifs_used.append(design_row['motifs_uses'])

    df.loc[:, "motifs_uses"] = motifs_used
    df = df[df.dG_predicted < 2].copy()

    motif_usage_df = get_motif_usage(df)

    plt.figure()
    g = sns.barplot(x="m_name", y="times_used", data=motif_usage_df)
    plt.setp(g.get_xticklabels(), rotation=90)
    plt.tight_layout()
    plt.savefig(fmanager.sim_results_path+"/plots/motif_used.png")
    plt.clf()


def plot_final_motif_usage_alt(design_df, df, fmanager):
    motifs_used = []
    for i, r in df.iterrows():
        design_row = design_df.iloc[r.design_num]
        motifs_used.append(design_row['motifs_uses'])

    df.loc[:, "motifs_uses"] = motifs_used
    df = df[df.dG_predicted < 2].copy()

    motif_usage_df = get_motif_usage(df)
    motif_usage_df.to_csv(fmanager.results_path+"final/motifs_used.csv",index=False)
    df.to_csv(fmanager.results_path+"final/best_alt.csv",index=False)

    plt.figure()
    g = sns.barplot(x="m_name", y="times_used", data=motif_usage_df)
    plt.setp(g.get_xticklabels(), rotation=90)
    plt.tight_layout()
    plt.savefig(fmanager.sim_results_path+"/plots/motif_used_alt.png")
    plt.clf()




def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', help='which problem working on', required=True)
    parser.add_argument('-wdir', help='working directory', required=True)
    parser.add_argument('-design', help='design constructs', action="store_true")
    parser.add_argument('-seq_opt', help='sequence optimization', action="store_true")
    parser.add_argument('-simulation', help='run tecto simulation', action="store_true")

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    fmanager = file_manager.FileManager(args.c, args.wdir)

    if args.design:
        design_df = pd.read_csv(fmanager.design_score_file)

        plot_design_motif_usage(design_df, fmanager)
        plot_design_distros(design_df, fmanager)

    if args.seq_opt:
        best_seq_opt_df = pd.read_csv(fmanager.seq_opt_best_summary_file)
        plot_seq_opt_distros(best_seq_opt_df, fmanager)

    if args.simulation:
        design_df = pd.read_csv(fmanager.design_score_file)
        sim_df = pd.read_csv(fmanager.sim_results_path+"/best.csv", index_col=False)
        sum_df = get_summary_df(sim_df)
        sum_df.to_csv(fmanager.sim_results_path+"/best_sum.csv", index_col=False)

        #plot_simulation_distros(sum_df, fmanager)
        #plot_final_motif_usage(design_df, sum_df, fmanager)
        #get_top_simulation_pdbs(sum_df, fmanager)
        #get_final_sequences(sum_df, fmanager)

    #design_df = pd.read_csv(fmanager.seq_opt_results_path+"/alt_inputs.csv")
    #sim_df = pd.read_csv(fmanager.sim_results_path+"/best_alt.csv", index_col=False)
    #sum_df = get_summary_df(sim_df)
    #sum_df.to_csv(fmanager.sim_results_path+"/best_alt_sum.csv", index_col=False)
    #plot_simulation_distros_alt(sum_df, fmanager)
    #plot_final_motif_usage_alt(design_df, sum_df, fmanager)








