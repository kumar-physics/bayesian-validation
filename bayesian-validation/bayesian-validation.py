from mmcif.io.PdbxReader import PdbxReader
import pynmrstar
import sys
import os
import numpy
import plotly.express as px
from operator import itemgetter

class BayesianValidation(object):
    def __init__(self,xray,nmr,outdir='../data/output'):
        xray_data = self.get_coordinates(xray)
        nmr_data =  self.get_coordinates(nmr)
        x_seq,n_seq=self.find_offset(xray_data,nmr_data)

        if not os.path.isdir('{}/{}-{}'.format(outdir, xray,nmr)):
            os.system('mkdir {}/{}-{}'.format(outdir, xray,nmr))
        odir='{}/{}-{}'.format(outdir,xray,nmr)
        self.write_seq_info(x_seq, n_seq, xray,nmr,odir)
        xray_mean=self.calculate_distance_matrix(xray_data,xray,odir,x_seq)
        nmr_mean=self.calculate_distance_matrix(nmr_data, nmr, odir, n_seq)
        self.cal_difference(xray_mean,nmr_mean,odir)

    def write_seq_info(self,x_seq,n_seq,x,n,outdir):
        out_file='{}/{}-{}_seq_alignment.txt'.format(outdir,x,n)
        fo=open(out_file,'w')
        for i in range(len(x_seq)):
            fo.write('{}-{}-{},{}-{}-{}\n'.format(x_seq[i][0],x_seq[i][1],x_seq[i][2],
                                                n_seq[i][0],n_seq[i][1],n_seq[i][2]))
        fo.close()


    def cal_difference(self,xray_mean,nmr_mean,outdir):
        diff_mat=[]
        out_mat = '{}/diff_mat.txt'.format(outdir)
        fo=open(out_mat,'w')
        for i in range(len(xray_mean)):
            diff_mat.append([])
            for j in range(len(xray_mean[i])):
                diff=abs(xray_mean[i][j]-nmr_mean[i][j])
                diff_mat[-1].append(diff)
                if j < len(xray_mean)-1:
                    fo.write('{},'.format(diff))
                else:
                    fo.write('{}\n'.format(diff))
        out_html='{}/diff_mat.html'.format(outdir)
        fo.close()
        fig=px.imshow(diff_mat)
        fig.write_html(out_html)


    def find_offset(self,xray_data,nmr_data):
        xray_seq = sorted(set([(i[0],i[1],i[2]) for i in xray_data[1].keys() if i[1]=='A']),key=itemgetter(1, 0))
        nmr_seq = sorted(set([(i[0],i[1],i[2]) for i in nmr_data[1].keys()if i[1]=='A']),key=itemgetter(1, 0))
        n=int(len(nmr_seq)/2)
        if len(xray_seq) < len(nmr_seq):
            short_seq=xray_seq
            short_seq2=xray_data[1].keys()
            long_seq = nmr_seq
            long_seq2 = nmr_data[1].keys()
            both_equal=False
            short='xray'
        elif len(xray_seq) > len(nmr_seq):
            short_seq=nmr_seq
            short_seq2 = nmr_data[1].keys()
            long_seq=xray_seq
            long_seq2 = xray_data[1].keys()
            both_equal=False
            short='nmr'
        else:
            short_seq = nmr_seq
            short_seq2 = nmr_data.keys()
            long_seq = xray_seq
            long_seq2 = xray_data.keys()
            both_equal=True
            short=None
        seq_len=9999999
        for i in range(-n,n):
            shifted_short=[(j[0]+i,j[1],j[2]) for j in short_seq]
            union_seq=set(list(shifted_short)+list(long_seq))
            if len(union_seq)<seq_len:
                seq_len=len(union_seq)
                l_k=[]
                s_k=[]
                for k in shifted_short:
                    if (k[0],k[1],k[2],'CA') in long_seq2 \
                            and (k[0],k[1],k[2],'C') in long_seq2 \
                            and (k[0],k[1],k[2],'N') in long_seq2 \
                            and (k[0]-i,k[1],k[2],'CA') in short_seq2 \
                            and (k[0]-i,k[1],k[2],'C') in short_seq2 \
                            and (k[0]-i,k[1],k[2],'N') in short_seq2:
                        l_k.append(k)
                        s_k.append((k[0]-i,k[1],k[2]))
                offset=i
        if both_equal:
            nmr_offset=offset
            xray_offset=0
            n_seq=s_k
            x_seq=l_k
        else:
            if short=='nmr':
                nmr_offset = offset
                xray_offset = 0
                n_seq = s_k
                x_seq = l_k
            else:
                nmr_offset = 0
                xray_offset = offset
                n_seq = l_k
                x_seq = s_k
        return x_seq,n_seq



    @staticmethod
    def get_pdb(pdb_id):
        cmd = 'wget https://files.rcsb.org/download/{}.cif -O ../data/cif/{}.cif'.format(pdb_id, pdb_id)
        os.system(cmd)

    @staticmethod
    def get_bmrb(bmrb_id):
        cmd = 'wget http://rest.bmrb.io/bmrb/{}/nmr-star3 -O ./data/star/{}.str'.format(bmrb_id, bmrb_id)
        os.system(cmd)

    def calculate_distance_matrix(self,cif_data,pdbid,outdir,seq):
        print ('calculating distance matrix for {}'.format(pdbid))
        seq_id = []
        for i in seq:
            for a in ['CA','C','N']:
                seq_id.append((i[0],i[1],i[2],a))
        mm=[]
        for m in cif_data.keys():
            d=[]
            if not os.path.isdir('{}/{}'.format(outdir,pdbid)):
                os.system('mkdir {}/{}'.format(outdir,pdbid))
            fo=open('{}/{}/{}_{}.txt'.format(outdir,pdbid,pdbid,m),'w')
            for i in seq_id:
                d.append([])
                for j in seq_id:
                    dist=round(self.get_distance(cif_data[m][i],cif_data[m][j]),4)
                    d[-1].append(dist)
                    if j!= seq_id[-1]:
                        fo.write('{},'.format(dist))
                    else:
                        fo.write('{}\n'.format(dist))
            fo.close()

            mm.append(d)
        meam_mat=[]
        sd_mat=[]
        fo1 = open('{}/{}/{}_mean.txt'.format(outdir,pdbid,pdbid),'w')
        fo2 = open('{}/{}/{}_std.txt'.format(outdir,pdbid,pdbid),'w')
        htmlout1='{}/{}/{}_mean.html'.format(outdir,pdbid,pdbid)
        htmlout2='{}/{}/{}_std.html'.format(outdir,pdbid,pdbid)
        for i in range(len(mm[0][0])):
            meam_mat.append([])
            sd_mat.append([])
            for j in range(len(mm[0][0])):
                d1=[]
                for k in range(len(mm)):
                    d1.append(mm[k][i][j])
                mn=round(numpy.mean(d1),4)
                sd=round(numpy.std(d1),4)
                meam_mat[-1].append(mn)
                sd_mat[-1].append(sd)
                if j < len(mm[0][0])-1:
                    fo1.write('{},'.format(mn))
                    fo2.write('{},'.format(sd))
                else:
                    fo1.write('{}\n'.format(mn))
                    fo2.write('{}\n'.format(sd))
        fig1=px.imshow(meam_mat)
        fig2=px.imshow(sd_mat)
        fig1.write_html(htmlout1)
        fig2.write_html(htmlout2)
        #fig1.show()
        #fig2.show()
        fo1.close()
        fo2.close()
        return meam_mat

    @staticmethod
    def get_distance(c1, c2):
        """
        Calculates the distance between two coordinate points
        :param c1: array of x,y,z
        :param c2: array of x,y,z
        :return: distance between two ponts
        """
        return numpy.linalg.norm(c1 - c2)

    def get_coordinates(self, pdbid, use_auth_tag=False,atom=['CA','C','N']):
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
        entity_poly = c0.getObj('entity_poly')
        chains = entity_poly.getValue("pdbx_strand_id")
        max_models = int(atom_site.getValue('pdbx_PDB_model_num', -1))
        col_names = atom_site.getAttributeList()
        model_id = col_names.index('pdbx_PDB_model_num')
        group_pdb=col_names.index('group_PDB')
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
                    if dat[group_pdb]=='ATOM':
                        if use_auth_tag:
                            if dat[aut_atom_id] in atom:
                                aid[(int(dat[aut_seq_id]), dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                                    (dat[entity_id], dat[asym_id], dat[comp_id], dat[seq_id], dat[aut_seq_id],
                                     dat[alt_id], dat[icode_id], dat[aut_asym_id])
                                pdb[(int(dat[aut_seq_id]), dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                                    numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
                        else:
                            if dat[atom_id] in atom:
                                #print (atom, dat)
                                aid[(int(dat[seq_id]), dat[asym_id], dat[comp_id], dat[atom_id])] = \
                                    (dat[entity_id], dat[asym_id], dat[comp_id], dat[seq_id], dat[aut_seq_id],
                                     dat[alt_id], dat[icode_id], dat[aut_asym_id])
                                pdb[(int(dat[seq_id]), dat[asym_id], dat[comp_id], dat[atom_id])] = \
                                    numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
            pdb_models[model] = pdb
            atom_ids[model] = aid
        return pdb_models


if __name__ == "__main__":
    #p=BayesianValidation('4FPW','CA')
    # f=open('../data/xray_nmr_pair_single_chain.csv','r').read().split("\n")[:-1]
    # for l in f:
    #     pair=l.split(",")
    #     xray=pair[0]
    #     nmr=pair[1]
    #     BayesianValidation(xray,nmr)
    BayesianValidation('3IDU','2KL6')

