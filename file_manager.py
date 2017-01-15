import os
import glob

class FileManager(object):
    def __init__(self, problem, wdir):
        self.wdir = wdir
        self.problem  = problem

        self.inputs_path = wdir+"/inputs/"
        self.outputs_path = wdir+"/outputs/"
        self.results_path = wdir+"/results/"
        self.runs_path = wdir+"/runs/"

        self.design_results_path = self.results_path+"construct_designs/"

        if not os.path.isdir(wdir):
            self._setup_dir()

        self.start_mg_file = self.inputs_path + self.problem + ".mg"
        self.design_score_file = self.design_results_path + self.problem + ".scores"
        self.design_out_file = self.design_results_path + self.problem + ".out"
        self.design_log_file = self.outputs_path + "design.log"

        self.seq_opt_results_path = self.results_path+"/sequence_opt/"
        self.seq_opt_input_path = self.results_path + "/sequence_opt/inputs/"
        self.seq_opt_input_file = self.results_path + "/sequence_opt/inputs.csv"
        self.seq_opt_best_summary_file = self.results_path + "/sequence_opt/best.csv"

        self.sim_results_path = self.results_path + "/simulation/"

    def _setup_dir(self):
        os.mkdir(self.wdir)
        os.mkdir(self.inputs_path)
        os.mkdir(self.outputs_path)
        os.mkdir(self.results_path)
        os.mkdir(self.wdir+"/pdbs")
        os.mkdir(self.wdir+"/results/construct_designs")
        os.mkdir(self.wdir+"/results/construct_designs/pdbs")
        os.mkdir(self.wdir+"/results/construct_designs/plots")
        os.mkdir(self.wdir+"/results/sequence_opt/")
        os.mkdir(self.wdir+"/results/sequence_opt/inputs")
        os.mkdir(self.wdir+"/results/sequence_opt/best")
        os.mkdir(self.wdir+"/results/sequence_opt/best/pdbs")
        os.mkdir(self.wdir+"/results/sequence_opt/random")
        os.mkdir(self.wdir+"/results/sequence_opt/plots")
        os.mkdir(self.wdir+"/results/simulation")
        os.mkdir(self.wdir+"/results/simulation/plots")
        os.mkdir(self.wdir+"/runs/")

        #setup default input files
        # for design
        f = open(self.wdir+"/inputs/design.input", "w")
        f.write("-verbose\n")
        f.write("-only_ideal\n")
        f.write("-pdbs\n")
        f.write("-designs 9999\n")
        f.write("-search.accept_score 15\n")
        f.write("-search.max_size 70\n")
        f.close()

        # for best sequence opt
        f = open(self.wdir+"/inputs/seq_opt_best.input", "w")
        f.write("-n 3\n")
        f.write("-optimizer.return_lowest \n")
        f.write("-optimizer.cutoff 1.0 \n")
        f.close()