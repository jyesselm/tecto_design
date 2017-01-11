from rnamake import resource_manager as rm
from rnamake import secondary_structure_parser, motif_type, motif_graph


class Factory(object):
    @staticmethod
    def get(name):
        if name == 'Flow11BP':
            return Flow11BP()
        if name == 'Flow10BP':
            return Flow10BP()
        if name == 'Flow16BP':
            return Flow16BP()
        else:
            raise ValueError("no problem with name: " + name)

class TectoDesignConstruct(object):
    def __init__(self):
        self.fseq = ""
        self.fss = ""
        self.mg = None
        self.name = ""

    def _get_m_names_from_seq_and_ss(self, seq, ss):
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


def convert_to_u(seq):
    new_seq = ""
    for e in seq:
        if e == "T":
            new_seq += "U"
        else:
            new_seq += e
    return new_seq


class Flow16BP(TectoDesignConstruct):
    def __init__(self):
        self.name = "Flow16BP"
        self.fseq = "CUAGGAAUCUGGUACCGAGAAGUACGGAAACGUACUUCUCGGUACCUGUGUCCUAG"
        self.fss  = "((((((....((((((((((((((((....))))))))))))))))....))))))"
        self.start_ni = -1
        self.start_ei = -1
        self.end_ni   = -1
        self.end_ei   = -1
        self.mg   = self._setup_mg()


    def _setup_mg(self):

        ref_ss  = '((((((....((((((((((((....))))))))))))....))))))'
        ref_seq = 'CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG'

        seq = ref_seq[11:22]+"+"+ref_seq[26:37]
        db  = ref_ss[11:22]+"+"+ref_ss[26:37]

        new_seq = convert_to_u(seq)

        mg = motif_graph.MotifGraph()

        m = rm.manager.get_bp_step("GG_LL_CC_RR")
        mg.add_motif(m)
        mg.add_motif(m_name="GGAA_tetraloop", m_end_name="A14-A15")

        parser = secondary_structure_parser.SecondaryStructureParser()
        g = parser.parse_to_motif_graph(new_seq, db)
        for i, n in enumerate(g.graph):
            m = rm.manager.get_bp_step(n.data.end_ids[0])
            if i == 0:
                mg.add_motif(m, parent_end_name="A7-A22")
            else:
                mg.add_motif(m)
        for i, n in enumerate(g.graph):
            if i > 5:
                break
            m = rm.manager.get_bp_step(n.data.end_ids[0])
            mg.add_motif(m)

        name1 = mg.last_node().data.ends[1].res1.name
        name2 = mg.last_node().data.ends[1].res2.name
        last = name2 + "G_LL_C" + name1 + "_RR"
        m = rm.manager.get_bp_step(last)
        mg.add_motif(m)

        pos = mg.add_motif(m_name="GAAA_tetraloop", m_end_name="A149-A154")
        self.start_ni = pos
        self.start_ei = mg.get_node(pos).data.get_end_index("A222-A251")
        mg.add_motif(rm.manager.get_bp_step('CC_LL_GG_RR'), parent_end_name="A229-A245")
        mg.add_motif(rm.manager.get_bp_step('CU_LL_AG_RR'))
        mg.add_motif(rm.manager.get_bp_step('UA_LL_UA_RR'))
        mg.add_motif(rm.manager.get_bp_step('AG_LL_CU_RR'))

        #mg.add_motif(rm.manager.get_bp_step('CG_LL_CG_RR'),
        #         parent_index=1, parent_end_name="A1-A6")
        pos = mg.add_motif(rm.manager.get_motif(name='HELIX.IDEAL'),
                            parent_index=1, parent_end_name="A1-A6")
        self.end_ni = pos
        self.end_ei = 1

        seq_spl = mg.sequence().split("&")
        ss_spl  = mg.dot_bracket().split("&")

        self.fseq = "CUAGGAAUCUG" + seq_spl[0][11:-11] + "CUGUGUCCUAG"
        self.fss  = "((((((....(" + ss_spl[0][11:-11]  + ")....))))))"

        return mg

    def write_mg_problem(self, fname="bp11_problem.mg"):
        start_name = "A222-A251"
        end_name = self.mg.get_node(self.end_ni).data.ends[self.end_ei].name()

        f = open(fname, "w")
        f.write(self.mg.to_str() + "\n")
        f.write(str(self.start_ni) + " " + start_name + "\n")
        f.write(str(self.end_ni) + " " + end_name)
        f.close()


class Flow12BP(TectoDesignConstruct):
    def __init__(self):
        self.name = "Flow12BP"
        self.fseq = "CUAGGAAUCUGGAAGUACACGAGGAAACUCGUGUACUUCCUGUGUCCUAG"
        self.fss  = "((((((....(((((((((((((....)))))))))))))....))))))"
        self.start_ni = -1
        self.start_ei = -1
        self.end_ni   = -1
        self.end_ei   = -1
        self.mg   = self._setup_mg()


    def _setup_mg(self):
        steps = self._get_m_names_from_seq_and_ss(self.fseq, self.fss)
        mg = motif_graph.MotifGraph()

        m = rm.manager.get_bp_step("GG_LL_CC_RR")
        mg.add_motif(m)
        mg.add_motif(m_name="GGAA_tetraloop", m_end_name="A14-A15")
        m = rm.manager.get_bp_step(steps[1])
        mg.add_motif(m, parent_end_name="A7-A22")

        for i in range(2, len(steps)):
            m = rm.manager.get_bp_step(steps[i])
            mg.add_motif(m)

        pos = mg.add_motif(m_name="GAAA_tetraloop", m_end_name="A149-A154")
        self.start_ni = pos
        self.start_ei = mg.get_node(pos).data.get_end_index("A222-A251")
        mg.add_motif(rm.manager.get_bp_step('CC_LL_GG_RR'), parent_end_name="A229-A245")
        mg.add_motif(rm.manager.get_bp_step('CU_LL_AG_RR'))
        mg.add_motif(rm.manager.get_bp_step('UA_LL_UA_RR'))
        mg.add_motif(rm.manager.get_bp_step('AG_LL_CU_RR'))

        #mg.add_motif(rm.manager.get_bp_step('CG_LL_CG_RR'),
        #         parent_index=1, parent_end_name="A1-A6")
        pos = mg.add_motif(rm.manager.get_motif(name='HELIX.IDEAL'),
                            parent_index=1, parent_end_name="A1-A6")
        self.end_ni = pos
        self.end_ei = 1

        return mg

    def write_mg_problem(self, fname="bp11_problem.mg"):
        start_name = "A222-A251"
        end_name = self.mg.get_node(self.end_ni).data.ends[self.end_ei].name()

        f = open(fname, "w")
        f.write(self.mg.to_str() + "\n")
        f.write(str(self.start_ni) + " " + start_name + "\n")
        f.write(str(self.end_ni) + " " + end_name)
        f.close()


class Flow11BP(TectoDesignConstruct):
    def __init__(self):
        self.name = "Flow11BP"
        self.fseq = "CUAGGAAUCUGGAAGUACACGAGGAAACUCGUGUACUUCCUGUGUCCUAG"
        self.fss  = "((((((....(((((((((((((....)))))))))))))....))))))"
        self.start_ni = -1
        self.start_ei = -1
        self.end_ni   = -1
        self.end_ei   = -1
        self.mg   = self._setup_mg()


    def _setup_mg(self):
        steps = self._get_m_names_from_seq_and_ss(self.fseq, self.fss)
        mg = motif_graph.MotifGraph()

        m = rm.manager.get_bp_step("GG_LL_CC_RR")
        mg.add_motif(m)
        mg.add_motif(m_name="GGAA_tetraloop", m_end_name="A14-A15")
        m = rm.manager.get_bp_step(steps[1])
        mg.add_motif(m, parent_end_name="A7-A22")

        for i in range(2, len(steps)):
            m = rm.manager.get_bp_step(steps[i])
            mg.add_motif(m)

        pos = mg.add_motif(m_name="GAAA_tetraloop", m_end_name="A149-A154")
        self.start_ni = pos
        self.start_ei = mg.get_node(pos).data.get_end_index("A222-A251")
        mg.add_motif(rm.manager.get_bp_step('CC_LL_GG_RR'), parent_end_name="A229-A245")
        mg.add_motif(rm.manager.get_bp_step('CU_LL_AG_RR'))
        mg.add_motif(rm.manager.get_bp_step('UA_LL_UA_RR'))
        mg.add_motif(rm.manager.get_bp_step('AG_LL_CU_RR'))

        #mg.add_motif(rm.manager.get_bp_step('CG_LL_CG_RR'),
        #         parent_index=1, parent_end_name="A1-A6")
        pos = mg.add_motif(rm.manager.get_motif(name='HELIX.IDEAL'),
                            parent_index=1, parent_end_name="A1-A6")
        self.end_ni = pos
        self.end_ei = 1

        return mg

    def write_mg_problem(self, fname="bp11_problem.mg"):
        start_name = "A222-A251"
        end_name = self.mg.get_node(self.end_ni).data.ends[self.end_ei].name()

        f = open(fname, "w")
        f.write(self.mg.to_str() + "\n")
        f.write(str(self.start_ni) + " " + start_name + "\n")
        f.write(str(self.end_ni) + " " + end_name)
        f.close()


class Flow10BP(TectoDesignConstruct):
    def __init__(self):
        self.name = "Flow10BP"
        self.fseq = "CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG"
        self.fss  = "((((((....((((((((((((....))))))))))))....))))))"
        self.start_ni = -1
        self.start_ei = -1
        self.end_ni   = -1
        self.end_ei   = -1
        self.mg   = self._setup_mg()


    def _setup_mg(self):
        steps = self._get_m_names_from_seq_and_ss(self.fseq, self.fss)
        mg = motif_graph.MotifGraph()

        m = rm.manager.get_bp_step("GG_LL_CC_RR")
        mg.add_motif(m)
        mg.add_motif(m_name="GGAA_tetraloop", m_end_name="A14-A15")
        m = rm.manager.get_bp_step(steps[1])
        mg.add_motif(m, parent_end_name="A7-A22")

        for i in range(2, len(steps)):
            m = rm.manager.get_bp_step(steps[i])
            mg.add_motif(m)

        pos = mg.add_motif(m_name="GAAA_tetraloop", m_end_name="A149-A154")
        self.start_ni = pos
        self.start_ei = mg.get_node(pos).data.get_end_index("A222-A251")
        mg.add_motif(rm.manager.get_bp_step('CC_LL_GG_RR'), parent_end_name="A229-A245")
        mg.add_motif(rm.manager.get_bp_step('CU_LL_AG_RR'))
        mg.add_motif(rm.manager.get_bp_step('UA_LL_UA_RR'))
        mg.add_motif(rm.manager.get_bp_step('AG_LL_CU_RR'))

        #mg.add_motif(rm.manager.get_bp_step('CG_LL_CG_RR'),
        #         parent_index=1, parent_end_name="A1-A6")
        #pos = mg.add_motif(rm.manager.get_motif(name='HELIX.IDEAL'))
        pos = mg.add_motif(rm.manager.get_motif(name='HELIX.IDEAL'),
                            parent_index=1, parent_end_name="A1-A6")
        self.end_ni = pos

        self.end_ni = pos
        self.end_ei = 1

        return mg

    def write_mg_problem(self, fname="bp11_problem.mg"):
        start_name = "A222-A251"
        end_name = self.mg.get_node(self.end_ni).data.ends[self.end_ei].name()

        f = open(fname, "w")
        f.write(self.mg.to_str() + "\n")
        f.write(str(self.start_ni) + " " + start_name + "\n")
        f.write(str(self.end_ni) + " " + end_name)
        f.close()
