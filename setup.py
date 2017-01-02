from rnamake import resource_manager as rm
from rnamake import secondary_structure_parser, motif_type, motif_graph



def _get_m_names_from_seq_and_ss(seq, ss):
        parser = secondary_structure_parser.SecondaryStructureParser()
        mg = parser.parse_to_motif_graph(seq, ss)

        start = 0
        names = []
        for n in mg.graph.nodes:
            if n.data.mtype == motif_type.TWOWAY:
                start = 1
                continue

            if n.data.mtype == motif_type.HAIRPIN:
                break

            if not start:
                continue

            name = n.data.end_ids[0]
            new_name = ""
            for e in name:
                if e == "T":
                    new_name += "U"
                else:
                    new_name += e
            names.append(new_name)
        return names


def get_11bp_problem():
    seq = "CUAGGAAUCUGGAAGUACACGAGGAAACUCGUGUACUUCCUGUGUCCUAG"
    ss  = "((((((....(((((((((((((....)))))))))))))....))))))"

    steps = _get_m_names_from_seq_and_ss(seq, ss)
    mg = motif_graph.MotifGraph()

    m = rm.manager.get_bp_step("GG_LL_CC_RR")
    mg.add_motif(m)
    mg.add_motif(m_name="GGAA_tetraloop", m_end_name="A14-A15")
    m = rm.manager.get_bp_step(steps[1])
    mg.add_motif(m, parent_end_name="A7-A22")

    for i in range(2, len(steps)):
        m = rm.manager.get_bp_step(steps[i])
        mg.add_motif(m)

    mg.add_motif(m_name="GAAA_tetraloop", m_end_name="A149-A154")
    mg.add_motif(rm.manager.get_bp_step('CC_LL_GG_RR'), parent_end_name="A229-A245")
    mg.add_motif(rm.manager.get_bp_step('CU_LL_AG_RR'))
    mg.add_motif(rm.manager.get_bp_step('UA_LL_UA_RR'))
    mg.add_motif(rm.manager.get_bp_step('AG_LL_CU_RR'))

    mg.add_motif(rm.manager.get_bp_step('CG_LL_CG_RR'),
                 parent_index=1, parent_end_name="A1-A6")
    mg.add_motif(rm.manager.get_motif(name='HELIX.IDEAL'))

    return mg


def export_11bp_problem():
    mg = get_11bp_problem()
    #for n in mg:
    #    print n.data.name
    mg.write_pdbs()
    mg.get_node(1).data.mtype = motif_type.TWOWAY
    mg.write_pdbs()
    f = open("bp11_problem.mg", "w")
    f.write(mg.to_str() + "\n")
    f.write(str(13) + " A222-A251\n")
    f.write(str(19) + " A4-B8\n")
    f.close()


def get_10bp_problem():
    seq = "CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG"
    ss  = "((((((....((((((((((((....))))))))))))....))))))"

    steps = _get_m_names_from_seq_and_ss(seq, ss)
    mg = motif_graph.MotifGraph()

    mg.add_motif(m_name="GC=GC")
    mg.add_motif(m_name="GGAA_tetraloop", m_end_name="A14-A15")
    mg.add_motif(m_name=steps[1], parent_end_name="A7-A22")

    for i in range(2, len(steps)):
        mg.add_motif(m_name=steps[i])

    mg.add_motif(m_name="GAAA_tetraloop", m_end_name="A149-A154")
    mg.add_motif(m_name="CG=GC", parent_index=1, parent_end_name="A1-A6")
    mg.add_motif(m_name="GC=AU")

    return mg


if __name__ == '__main__':
    export_11bp_problem()

    exit()

    mg = get_10bp_problem()
    mg.get_node(1).data.mtype = motif_type.TWOWAY
    #mg.write_pdbs()
    #exit()
    f = open("bp10_problem.mg", "w")
    f.write(mg.to_str() + "\n")
    f.write(str(12) + " A222-A251\n")
    f.write(str(14) + " R30-R37\n")
    f.close()
