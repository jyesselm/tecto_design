import argparse
import pandas as pd
import os
import shutil
import glob

from rnamake.wrappers import design_rna_wrapper

import design_constructs, file_manager

fmanager = None

def initial_design_constructs():
    drw = design_rna_wrapper.DesignRNAWrapper()
    drw.setup_from_file(args.wdir+"/inputs/design.input")
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

    files = glob.glob("design*.pdb")
    for f in files:
        shutil.move(f, fmanager.design_results_path + "/pdbs/")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', help='which problem working on', required=True)
    parser.add_argument('-wdir', help='working directory', required=True)
    parser.add_argument('-design', help='design constructs', action="store_true")

    #parser.add_argument('-s', type=int, help='steps in each simulation')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    c = design_constructs.Factory.get(args.c)
    fmanager = file_manager.FileManager(args.c, args.wdir)

    c.write_mg_problem(fmanager.start_mg_file)

    if args.design:
        initial_design_constructs()
