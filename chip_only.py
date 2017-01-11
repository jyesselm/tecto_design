from rnamake import motif_graph, motif_factory, motif
from rnamake import resource_manager as rm
import pandas as pd

import file_manager


def get_gaaa_receptor(m):
    chains = m.chains()
    res = []
    bps = []
    for c in chains:
        if c.first().num == 149:
            continue
        res.extend(c.residues)
    for bp in m.basepairs:
        if bp.res1 in res and bp.res2 in res:
            bps.append(bp)

    gaaa_just_receptor = motif_factory.factory.motif_from_res(res, bps)
    gaaa_just_receptor.name = "gaaa_receptor"
    rm.manager.add_motif(motif=gaaa_just_receptor)
    return gaaa_just_receptor


def get_ggaa_tetraloop(m):
    chains = m.chains()
    res = []
    bps = []
    for c in chains:
        if c.first().num != 1:
            continue
        res.extend(c.residues)
    for bp in m.basepairs:
        if bp.res1 in res and bp.res2 in res:
            bps.append(bp)

    ggaa_tetraloop = motif_factory.factory.motif_from_res(res, bps)
    ggaa_tetraloop.name = "ggaa_tetraloop"
    ggaa_tetraloop.block_end_add = -1
    return ggaa_tetraloop


def generate_mg_with_only_chip(mg):
    get_gaaa_receptor(mg.get_node(m_name="GAAA_tetraloop").data)
    ggaa_tetraloop = get_ggaa_tetraloop(mg.get_node(m_name="GGAA_tetraloop").data)

    # get first basepair step at 5' end of chip sequence
    gaaa_org_node = mg.get_node(m_name="GAAA_tetraloop")
    end_i = gaaa_org_node.data.get_end_index(name="A229-A245")
    current = gaaa_org_node.connections[end_i].partner(gaaa_org_node.index)
    while current is not None:
        if current.connections[1] is None:
            break
        current = current.connections[1].partner(current.index)

    new_mg = motif_graph.MotifGraph()
    new_mg.set_options('sterics', 0)

    new_m = motif_graph.flip_alignment(current.data, 1)
    new_mg.add_motif(new_m)

    pos = current.connections[0].partner(current.index).index
    while pos != gaaa_org_node.index:
        new_mg.add_motif(m_name=mg.get_node(pos).data.name,
                         m_end_name=mg.get_node(pos).data.ends[1].name())
        pos -= 1
    new_mg.add_motif(rm.manager.get_motif(name="gaaa_receptor", end_name="A229-A245"))

    c = mg.get_node(m_name="GAAA_tetraloop").connections[1]
    start = c.partner(mg.get_node(m_name="GAAA_tetraloop").index)
    current = start

    while 1:
        c = current.connections[1]
        if c.end_index_1 == c.end_index_2:
            break

        if len(new_mg) == 1:
            new_mg.add_motif(current.data, parent_end_name="A222-A251")
        else:
            new_mg.add_motif(current.data)
        current = current.connections[1].partner(current.index)

    last_pos = new_mg.last_node().index

    new_mg.add_motif(ggaa_tetraloop, orphan=1)
    #new_mg.add_motif(rm.manager.get_bp_step('CG_LL_CG_RR'))
    last_pos_2 = new_mg.add_motif(rm.manager.get_motif(name='HELIX.IDEAL'))

    new_mg.add_connection(last_pos, last_pos_2)
    connect_ni = [last_pos, last_pos_2]

    return new_mg, connect_ni


def generate_chip_only_mgs(fmanager):
    f = open(fmanager.design_out_file)
    lines = f.readlines()
    f.close()

    df = pd.DataFrame(columns="num fname ni1 ei1 ni2 ei2".split())

    pos = 0

    # based on visual inspection
    #exclude = [3, 5]

    for i, l in enumerate(lines):
        print i
        #if i in exclude:
        #    continue
        mg = motif_graph.MotifGraph(mg_str=l)
        try:
            new_mg, connect_ni = generate_mg_with_only_chip(mg)
        except:
            continue
        fname = fmanager.seq_opt_input_path + str(i)+".mg"
        f = open(fname, "w")
        f.write(new_mg.to_str() + "\n")
        f.close()

        df.loc[pos] = [i, fname, connect_ni[0], 1, connect_ni[1], 1]
        pos += 1

    df.to_csv(fmanager.seq_opt_input_file, index=False)





