from mmcif.io.PdbxReader import PdbxReader
import pynmrstar
import sys
import os
import numpy
import plotly.express as px
from operator import itemgetter

class BayesianValidation(object):
    def __init__(self,pdbid):
        cif_data = self.get_coordinates(pdbid)
        self.calculate_distance_matrix(cif_data)
    @staticmethod
    def get_pdb(pdb_id):
        cmd = 'wget https://files.rcsb.org/download/{}.cif -O ../data/cif/{}.cif'.format(pdb_id, pdb_id)
        os.system(cmd)

    @staticmethod
    def get_bmrb(bmrb_id):
        cmd = 'wget http://rest.bmrb.io/bmrb/{}/nmr-star3 -O ./data/star/{}.str'.format(bmrb_id, bmrb_id)
        os.system(cmd)

    def calculate_distance_matrix(self,cif_data):
        seq_id = sorted([i for i in cif_data[1].keys()],key=itemgetter(1, 0))
        mm=[]
        for m in cif_data.keys():
            d=[]
            for i in seq_id:
                d.append([])
                for j in seq_id:
                    d[-1].append(round(self.get_distance(cif_data[m][i],cif_data[m][j]),4))
            mm.append(d)
        meam_mat=[]
        sd_mat=[]
        for i in range(len(mm[0][0])):
            meam_mat.append([])
            sd_mat.append([])
            for j in range(len(mm[0][0])):
                d1=[]
                for k in range(len(mm)):
                    d1.append(mm[k][i][j])
                meam_mat[-1].append(round(numpy.mean(d1),4))
                sd_mat[-1].append(round(numpy.std(d1),4))
        fig1=px.imshow(meam_mat)
        fig2=px.imshow(sd_mat)
        fig1.show()
        fig2.show()


    @staticmethod
    def get_distance(c1, c2):
        """
        Calculates the distance between two coordinate points
        :param c1: array of x,y,z
        :param c2: array of x,y,z
        :return: distance between two ponts
        """
        return numpy.linalg.norm(c1 - c2)

    def get_coordinates(self, pdbid, use_auth_tag=True,atom='CA'):
        """
        Extract coordinate information from cif file as a dictionary
        {model_id : {(seq_id,chain_id,res_id,atom_id) : array[x,y,x],...},...}
        :param cif_file: Input coordinate file
        :return: dictionary
        """
        if not os.path.isdir('../data'):
            os.system('mkdir ../data')
        if not os.path.isdir('../data/cif'):
            os.system('mkdir ../data/cif')
        cif_file = '../data/cif/{}.cif'.format(pdbid)
        if not os.path.isfile(cif_file):
            self.get_pdb(pdbid)
        cif_data = []
        ifh = open(cif_file, 'r')
        pRd = PdbxReader(ifh)
        pRd.read(cif_data)
        ifh.close()
        c0 = cif_data[0]
        atom_site = c0.getObj('atom_site')
        max_models = int(atom_site.getValue('pdbx_PDB_model_num', -1))
        col_names = atom_site.getAttributeList()
        model_id = col_names.index('pdbx_PDB_model_num')
        x_id = col_names.index('Cartn_x')
        y_id = col_names.index('Cartn_y')
        z_id = col_names.index('Cartn_z')
        atom_id = col_names.index('label_atom_id')
        comp_id = col_names.index('label_comp_id')
        asym_id = col_names.index('label_asym_id')
        entity_id = col_names.index('label_entity_id')
        seq_id = col_names.index('label_seq_id')
        icode_id = col_names.index('pdbx_PDB_ins_code')
        alt_id = col_names.index('label_alt_id')
        aut_seq_id = col_names.index('auth_seq_id')
        aut_asym_id = col_names.index('auth_asym_id')
        aut_atom_id = col_names.index('auth_atom_id')
        aut_comp_id = col_names.index('auth_comp_id')
        pdb_models = {}
        atom_ids = {}
        for model in range(1, max_models + 1):
            pdb = {}
            aid = {}
            for dat in atom_site.getRowList():
                if int(dat[model_id]) == model:
                    if use_auth_tag:
                        if dat[aut_atom_id] == atom:
                            aid[(int(dat[aut_seq_id]), dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                                (dat[entity_id], dat[asym_id], dat[comp_id], dat[seq_id], dat[aut_seq_id],
                                 dat[alt_id], dat[icode_id], dat[aut_asym_id])
                            pdb[(int(dat[aut_seq_id]), dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                                numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
                    else:
                        if dat[atom_id] == atom:
                            aid[(int(dat[seq_id]), dat[asym_id], dat[comp_id], dat[atom_id])] = \
                                (dat[entity_id], dat[asym_id], dat[comp_id], dat[seq_id], dat[aut_seq_id],
                                 dat[alt_id], dat[icode_id], dat[aut_asym_id])
                            pdb[int((dat[seq_id]), dat[asym_id], dat[comp_id], dat[atom_id])] = \
                                numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
            pdb_models[model] = pdb
            atom_ids[model] = aid
        return pdb_models


    def read_cif(self,cif_file):
        pass
    def distance_matrix(self,coordinates):
        pass

if __name__ == "__main__":
    p=BayesianValidation('2K2E')