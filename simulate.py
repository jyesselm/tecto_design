from rnamake.unittests.intergration import simulate_tectos
from rnamake import motif_graph
from rnamake.wrappers import simulate_tectos_wrapper
import pandas as pd


def test_with_py():
    st = simulate_tectos.SimulateTectos()
    st.option('fseq', 'CUAGGAAUCUGGAAGUACACGAGGAAACUCGUGUACUUCCUGUGUCCUAG')
    st.option('fss',  '((((((....(((((((((((((....)))))))))))))....))))))')
    st.option('cseq', 'CUAGGAUAUGGUAUACGCAAGUUCGGGAACGAGAGCGUAUACCUAAGUCCUAG')
    st.option('css',  '(((((((..(((((((((....(((....)))..)))))))))...)))))))')

    st.run()
    #st.mset.to_mst().write_pdbs("full")

def run_best():
    df = pd.read_csv("20170101_11bp_constructs/data/sequence_opt_best.csv")
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
            'fseq' : 'CUAGGAAUCUGGAAGUACACGAGGAAACUCGUGUACUUCCUGUGUCCUAG',
            'fss'  : '((((((....(((((((((((((....)))))))))))))....))))))',
        }
        stw.run(**opts)
        avg_hit_count = stw.get_output()
        new_row = r.tolist()
        new_row.append(avg_hit_count)
        df_new.loc[pos] = new_row
        pos += 1
        print r.design_num, avg_hit_count
        #exit()

    df_new.to_csv("20170101_11bp_constructs/data/sequence_opt_best_sim.csv", index=False)


#test_with_py()
run_best()

"""f = open("20170101_11bp_constructs/data/sequence_opt_out_best/16.out")
lines = f.readlines()
f.close()

mg = motif_graph.MotifGraph(lines[0])
for n in mg:
    print n.data.name, n.data.secondary_structure
mg.write_pdbs()
"""
