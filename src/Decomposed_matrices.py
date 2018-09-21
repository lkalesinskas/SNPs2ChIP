import numpy as np
import pandas as pd

class Decomposed_matrices:
    '''
    Read three matrices stored in csv.gz files and compute
    square cosine score and contribution score.
    
    Here are some explanation (see more details for the paper).
    
    Let's write our decomposition as X = UDV' 
    where X is input feature matrix, 
    D is diagonal singular value matrix, 
    U is left singular vector matrix (on assay space), 
    V is right singular vector matrix (on genomic bin space), 
    and `'` denotes the transposition of the matrix.
    
    Genomic bin squared cosine score is defined as L2-normalized 
    version of the matrix product (VD) so that any given slice 
    for a given genomic bin has Euclidian norm of 1. 
    
    The interpretation of the score is it represents the relative 
    importance of the component given a genomic bin.    
    '''
    def read_decomposed_matrix_file(self, filename, compression=None):
        '''
        This function reads the csv.gz file
        '''
        if((compression is None) and (len(filename) > 3) and (filename[-3:] == '.gz')):
            compression='gzip'
        df = pd.read_csv(filename, compression=compression)
        mat = df.iloc[:, 1:].as_matrix()
        idx = df.iloc[:, 0].as_matrix()
        return mat, idx    

    def __init__(self, d_file, u_file, v_file):
        '''
        Given the file names of the three csv.gz filess, 
        load them on memory and compute the scores.
        '''
        # read three matrices D, U, V
        d_mat_temp, d_idx = self.read_decomposed_matrix_file(d_file)
        d_vec = d_mat_temp[:, 0]
        u_mat, u_idx = self.read_decomposed_matrix_file(u_file)
        v_mat, v_idx = self.read_decomposed_matrix_file(v_file)
        
        # Compute matrix produces, UD and VD
        u_dot_d = np.dot(u_mat, np.diag(d_vec))
        v_dot_d = np.dot(v_mat, np.diag(d_vec))
        
        # Compute the following four normalized matrices.
        # - v_dot_d_find_pcs: genomic bin --> which PCs? 
        #    genomic bin squared contribution score.
        # - u_dot_d_fine_pcs: assay       --> which PCs? 
        #    assay squared contribution score.
        # - v_dot_d_find_loci: PC --> which genomic bins? 
        #    genomic bin contribution score.
        # - u_dot_d_find_assays: PC --> which assays? 
        #    assay contribution score        
        self.v_dot_d_find_pcs    = (v_dot_d ** 2 ) / (np.sum(v_dot_d ** 2, axis = 1)[:,np.newaxis])
        self.u_dot_d_find_pcs    = (u_dot_d ** 2 ) / (np.sum(u_dot_d ** 2, axis = 1)[:,np.newaxis])
        self.v_dot_d_find_loci   = (v_dot_d ** 2 ) / (np.sum(v_dot_d ** 2, axis = 0)[np.newaxis, :])
        self.u_dot_d_find_assays = (u_dot_d ** 2 ) / (np.sum(u_dot_d ** 2, axis = 0)[np.newaxis, :])
            
    def __find_generic(self, find_matrix, query_idx, topk=5):
        '''
        Generic form of "find_X_given_Y" function.
        '''
        return_idxs = np.argsort(-find_matrix[query_idx, :])[:topk]
        scores = find_matrix[query_idx, return_idxs]
        return return_idxs, scores
    
    def find_pcs_given_locus(self, locus_idx, topk=5):
        '''
        Given a specified locus_idx, return the top k PCs (default 5) and
        the corresponding their relative importance (genomic bin squared contribution score).
        '''
        return self.__find_generic(self.v_dot_d_find_pcs, locus_idx, topk)
    
    def find_pcs_given_assay(self, assay_idx, topk=5):
        '''
        Given a specified assay_idx, return the top k PCs (default 5) and
        the corresponding their relative importance (assay squared contribution score).
        '''
        return self.__find_generic(self.u_dot_d_find_pcs, assay_idx, topk)

    def find_loci_given_pc(self, pc_idx, topk=5):
        '''
        Given a specified PC index, return the top k genomic loci (default 5) and
        the corresponding their relative importance (genomic bin contribution score).
        '''
        return self.__find_generic(np.transpose(self.v_dot_d_find_loci), pc_idx, topk)
        
    def find_assays_given_pc(self, pc_idx, topk=5):
        '''
        Given a specified PC index, return the top k assays (default 5) and
        the corresponding their relative importance (assay contribution score).
        '''
        return self.__find_generic(np.transpose(self.u_dot_d_find_assays), pc_idx, topk)

    # multipe SNPs/assays mode
    def __find_given_list_generic(self, find_matrix, query_idxs, weights = None, topk=5):
        '''
        Generic form of "find_pcs_given_X_list" function.
        '''
        find_matrix_slice = find_matrix[query_idxs, :]
        if(weights is None):
            find_vec = np.sum(find_matrix_slice, axis = 0) / len(query_idxs)
        else:
            weights_vec = np.array(weights)[:, np.newaxis] / np.sum(weights)
            find_vec = np.sum(weights_vec * find_matrix_slice, axis = 0)            
        return_idxs = np.argsort(-find_vec)[:topk]
        scores = find_vec[return_idxs]
        return return_idxs, scores
    
    def find_pcs_given_loci_list(self, loci_idxs, weights = None, topk=5):
        '''
        Given indices of genomic loci as a list, return the top k PCs (default 5) and
        the corresponding their relative importance.
        One can also give weights for the list of loci.        
        '''
        return self.__find_given_list_generic(self.v_dot_d_find_pcs, loci_idxs, weights, topk)
        
    def find_pcs_given_assay_list(self, assay_idxs, weights = None, topk=5):
        '''
        Given indices of assays as a list, return the top k PCs (default 5) and
        the corresponding their relative importance.
        One can also give weights for the list of assays.
        '''
        return self.__find_given_list_generic(self.u_dot_d_find_pcs, assay_idxs, weights, topk)
