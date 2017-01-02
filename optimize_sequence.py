import pandas as pd
from rnamake.wrappers import sequence_optimizer_wrapper
from rnamake import util
import glob

wdir = "20170101_11bp_constructs/"


def get_best_opt():
    spw = sequence_optimizer_wrapper.SequenceOptimizerWrapper()

    df = pd.read_csv(wdir+"/data/seq_opt.csv")
    for i,r in df.iterrows():
        print r.num
        end_1 = str(int(r.ni1))+","+str(int(r.ei1))
        end_2 = str(int(r.ni2))+","+str(int(r.ei2))
        opts = {
            'mg' : r.fname,
            'end_1' : end_1,
            'end_2' : end_2,
            'optimizer.return_lowest' : True,
            'optimizer.cutoff': 1.0,
            'n' : 3,
            'out_file' : wdir+"data/sequence_opt_out_best/"+str(int(r.num))+".out",
            'score_file' : wdir+"data/sequence_opt_out_best/"+str(int(r.num))+".scores"}
        spw.run(**opts)

    sc_files = glob.glob(wdir+"data/sequence_opt_out_best/*.scores")
    dfs = []
    for sc_file in sc_files:
        fname = util.filename(sc_file)
        dnum = int(fname[:-7])
        df = pd.read_csv(sc_file)
        df['design_num'] = [dnum, dnum, dnum]
        dfs.append(df)

    df_full = pd.concat(dfs)
    df_full.to_csv(wdir+"data/sequence_opt_best.csv", index=False)

get_best_opt()