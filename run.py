import argparse
import pandas as pd
import os
import shutil
import glob

from rnamake.wrappers import design_rna_wrapper, sequence_optimizer_wrapper
from rnamake.wrappers import simulate_tectos_wrapper
from rnamake import util

import design_constructs, file_manager, chip_only

def initial_design_constructs(fmanager):
    drw = design_rna_wrapper.DesignRNAWrapper()
    drw.setup_from_file(fmanager.inputs_path+"design.input")
    opts = {
        'mg'         : fmanager.start_mg_file,
        'score_file' : fmanager.design_score_file,
        'out_file'   : fmanager.design_out_file
    }
    drw.run(**opts)
    output = drw.get_output()
    print output
    f = open(fmanager.outputs_path + "design.log", "w")
    f.write(output)
    f.close()

    old_files = glob.glob( fmanager.design_results_path + "/pdbs/*.pdb")
    for f in old_files:
        os.remove(f)

    files = glob.glob("design*.pdb")
    for f in files:
        shutil.move(f, fmanager.design_results_path + "/pdbs/")


def get_best_opt(fmanager):
    spw = sequence_optimizer_wrapper.SequenceOptimizerWrapper()
    spw.setup_from_file(fmanager.inputs_path+"seq_opt_best.input")

    df = pd.read_csv(fmanager.seq_opt_input_file)
    for i,r in df.iterrows():
        print r.num
        end_1 = str(int(r.ni1))+","+str(int(r.ei1))
        end_2 = str(int(r.ni2))+","+str(int(r.ei2))
        opts = {
            'mg'         : r.fname,
            'end_1'      : end_1,
            'end_2'      : end_2,
            'out_file'   : fmanager.seq_opt_results_path+"best/"+str(int(r.num))+".out",
            'score_file' : fmanager.seq_opt_results_path+"best/"+str(int(r.num))+".scores"
        }
        spw.run(**opts)

    sc_files = glob.glob(fmanager.seq_opt_results_path+"best/*.scores")
    dfs = []
    for sc_file in sc_files:
        fname = util.filename(sc_file)
        dnum = int(fname[:-7])
        df = pd.read_csv(sc_file)
        df['design_num'] = [dnum, dnum, dnum]
        dfs.append(df)

    df_full = pd.concat(dfs)
    df_full.to_csv(fmanager.seq_opt_best_summary_file, index=False)


def run_simulations(fmanager, c):
    df = pd.read_csv(fmanager.seq_opt_best_summary_file)
    seen = []
    df_new = pd.DataFrame(columns=df.columns.values)
    df_new['avg_hit_count'] = []
    pos = 0

    stw = simulate_tectos_wrapper.SimulateTectosWrapper()
    for i, r in df.iterrows():
        if r.opt_sequence in seen:
            continue
        seen.append(r.opt_sequence)
        #if r.design_num != 16:
        #    continue
        # fix tetraloop ss
        css = r.opt_structure
        css = "(((((((.." + css[9:-10] + "...)))))))"
        opts = {
            'cseq' : r.opt_sequence,
            'css'  : css,
            'fseq' : c.fseq,
            'fss'  : c.fss
        }
        stw.run(**opts)
        avg_hit_count = stw.get_output()
        new_row = r.tolist()
        new_row.append(avg_hit_count)
        df_new.loc[pos] = new_row
        pos += 1
        print r.design_num, avg_hit_count

    df_new.to_csv(fmanager.sim_results_path + "/best.csv", index=False)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', help='which problem working on', required=True)
    parser.add_argument('-wdir', help='working directory', required=True)
    parser.add_argument('-design', help='design constructs', action="store_true")
    parser.add_argument('-setup_seq_opt', help='setup sequence optimization', action="store_true")
    parser.add_argument('-seq_opt', help='sequence optimization', action="store_true")
    parser.add_argument('-simulation', help='run tecto simulation', action="store_true")

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    c = design_constructs.Factory.get(args.c)
    fmanager = file_manager.FileManager(args.c, args.wdir)

    c.write_mg_problem(fmanager.start_mg_file)
    c.mg.to_pdb(fmanager.wdir+"/pdbs/start_mg.pdb", renumber=1, close_chain=1)

    if args.design:
        initial_design_constructs(fmanager)

    if args.setup_seq_opt:
        chip_only.generate_chip_only_mgs(fmanager)

    if args.seq_opt:
        get_best_opt(fmanager)

    if args.simulation:
        run_simulations(fmanager, c)
