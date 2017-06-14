import sys
import pandas as pd
import numpy as np
import xmltemplates as xml

def b_from_c(c,omega):
    "c is in uF b in siemens"
    return omega*c*1e-6

def main(fname,nodes_file,edges_file):
    nodes = pd.read_csv(nodes_file,index_col=0)
    edges = pd.read_csv(edges_file,index_col=0)
    #TOP = xml.Topology('feeder_topology.xml')
    #EQ  = xml.Equipment('feeder_equipment.xml')
    TOP = xml.Topology(path=fname)
    EQ  = xml.Equipment(fid=TOP.file)
    f = 50. #[Hz]
    omega = 2*np.pi*f
    
    conductor_lib = {}
    for u in [3,10,20]:
        conductor_lib[u] = pd.read_csv('../cables_%dkV.csv' %(u))
        conductor_lib[u].index += 1
    
    
    ###### Create base voltages ###########
    EQ.base_voltage(unoms=nodes['unom'].unique())
    
    sub_counter = 0
    for i in nodes.index:
        idbase = '_%d_feeder_%d' %(i,nodes.loc[i].fid)
        TOP.top_node(id='_top'+idbase, basekv_id= EQ.kvid[nodes.loc[i].unom])
        if i == 1:
            ###### Create slack bus ###########
            EQ.external_network(id='_external_network' + idbase, sub_id = sub_counter)
            sub_counter += 1
            EQ.terminal(id='_external_network_terminal' + idbase, sequence_num=1, equipment_id='_external_network' + idbase)

            TOP.terminal_connect(terminal_id='_external_network_terminal'+idbase,top_node_id='_top'+idbase,status=True)
        
        if nodes.loc[i].p != 0:
            EQ.load_object(id='_load'+idbase, p=nodes.loc[i].p, q=nodes.loc[i].q)
            EQ.terminal(id='_load_terminal'+idbase, sequence_num=1, equipment_id='_load'+idbase)
            TOP.terminal_connect(terminal_id='_load_terminal'+idbase,top_node_id='_top'+idbase,status=True)
    
    EQ.operational_limit_type(id='_abs_limit_type')
    for i in edges.index:
        f   = int(edges.loc[i].f)
        t   = int(edges.loc[i].t)
        fid = int(edges.loc[i].fid)
        idbase='_%i_feeder_%d' %(i,fid)
        if edges.loc[i].funom == edges.loc[i].tunom:
            u   = edges.loc[i].funom
            cid = edges.loc[i].cable_id
            length = edges.loc[i].length
            EQ.conductor(id='_branch'+idbase,
                    r =length*conductor_lib[u].loc[cid].r,  x =length*conductor_lib[u].loc[cid].x, bch =length*b_from_c(conductor_lib[u].loc[cid].c,omega),
                    r0=length*conductor_lib[u].loc[cid].r0, x0=length*conductor_lib[u].loc[cid].x0, b0ch=length*b_from_c(conductor_lib[u].loc[cid].c0,omega),
                    length=edges.loc[i].length,basekv_id=EQ.kvid[u])

            EQ.terminal(id='_branch_terminal_f'+idbase, sequence_num=1, equipment_id='_branch'+idbase)
            EQ.terminal(id='_branch_terminal_t'+idbase, sequence_num=2, equipment_id='_branch'+idbase)
            
            EQ.operational_limit_set(id='_branch_ratings'+idbase, equipment_id='_branch'+idbase)
            EQ.current_limit(id='_branch_current_limit_'+idbase, val=conductor_lib[u].loc[cid].inom, set_id='_branch_ratings'+idbase, type_id='_abs_limit_type')
        else:
            uf   = int(edges.loc[i].funom)
            ut   = int(edges.loc[i].tunom)
            if uf > ut:
                EQ.power_transformer_end(id='_xfrm_node%d-node%d_wf' %(f,t) + idbase, snom=100, unom=uf, terminal='_branch_terminal_f' + idbase, winding='Yn', clocknum=1)
                EQ.power_transformer_end(id='_xfrm_node%d-node%d_wt' %(f,t) + idbase, snom=100, unom=ut, terminal='_branch_terminal_t' + idbase, winding='D', clocknum=1)
            else:
                EQ.power_transformer_end(id='_xfrm_node%d-node%d_wf' %(f,t) + idbase, snom=100, unom=uf, terminal='_branch_terminal_f' + idbase, winding='D', clocknum=1)
                EQ.power_transformer_end(id='_xfrm_node%d-node%d_wt' %(f,t) + idbase, snom=100, unom=ut, terminal='_branch_terminal_t' + idbase, winding='Yn', clocknum=1)

            EQ.terminal(id='_branch_terminal_f'+idbase, sequence_num=1, equipment_id='_xfrm_node%d-node%d_wf' %(f,t) + idbase)
            EQ.terminal(id='_branch_terminal_t'+idbase, sequence_num=2, equipment_id='_xfrm_node%d-node%d_wt' %(f,t) + idbase)

            EQ.power_transformer(id='_xfrm_node%d-node%d' %(f,t) + idbase, 
                    wf_id='_xfrm_node%d-node%d_wf' %(f,t) + idbase, wt_id='_xfrm_node%d-node%d_wt' %(f,t) + idbase)

        TOP.terminal_connect(terminal_id='_branch_terminal_f' + idbase, top_node_id='_top_%d_feeder_%d' %(f,fid), status=True)
        TOP.terminal_connect(terminal_id='_branch_terminal_t' + idbase, top_node_id='_top_%d_feeder_%d' %(t,fid), status=True)

    EQ.close()
    TOP.file.close()

if __name__ == '__main__':
    main(*sys.argv[1:])
