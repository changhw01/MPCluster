import numpy as np
from DisjointSets import DisjointSets

def search_first_true(flags):
    ''' Find the first index of True '''
    
    n_flag = len(flags)
    i = 0
    while i < n_flag:
        if flags[i]:
            break
        i+=1
    
    if i == n_flag:
        return None
    else:
        return i

def adjust_rank(rep,r):
    ''' Adjust the ranks of the representation '''
    n_clst, r_pre = rep.shape
    
    rep_new = None
    if r_pre == r:
        rep_new = rep
    elif r_pre > r:
        rep_new = rep[:,0:r]
    else:
        rep_new = np.hstack(rep,np.zeros((n_clst,r-r_pre)))
        
    return rep_new
        
class MPCA:
    ''' Perform multi-persistent clustering analysis from given distance matrix
        and density estimate.

    Parameters
    ----------
    distM: float, numpy 2d-array
        The pairwise distance function.
    density: float, numpy 1d-array
        The density estimate of each data point.


    '''
    def __init__(self, distM, density):
        self.distM = distM
        self.density = density
        self.n_data = len(density)


    def fit(self, thres_dsty, thres_dist):
        ''' Perform multi-persistent analysis using the given thresholds.

        Parameters:
        ----------
        thres_dsty: float, numpy 1d-array
            The density thresholds for multi-persistent analysis.
            It should be sorted from small to large.
        thres_dist: float, numpy 1d-array
            The distance thresholds for multi-persistent analysis.
            It should be sorted from small to large.
        n_h: int
            The number of density thresholds.
        n_r: int
            The number of distance thresholds.
        spM: int, numpy 2d-array
            A n_h-by-n_r array for recording the sub-partitions.
        sM: int, numpy 2d-array
            A n_h-by-n_r array for recording the corresponding cluster size.
        subptts: list of bool array
        size_levelset: int, numpy array            
           
        Note:
        ----------
        h: for density thresholds
        r: for distance thresholds
        ptt: partition
        
        '''

        # Set basic parameters 
        n_h = len(thres_dsty)
        n_r = len(thres_dist)
        spM = np.zeros((n_h, n_r))
        sM = np.zeros((n_h, n_r))
        subptts = []
        size_levelset = np.array([self.n_data] * n_h)
        rank_level = np.zeros(n_h)
        
        # Sort the data by density
        elm2idx = np.argsort(self.density)
        elm2idx = elm2idx[-1::-1]
        
        self.distM = self.distM[elm2idx, :]
        self.distM = self.distM[:, elm2idx]
        self.density = self.density[elm2idx]
        
        # Compute the sizes of all level sets
        for ih in xrange(n_h):
            for i in xrange(self.n_data):
                if self.density[i] < thres_dsty[ih]:
                    size_levelset[ih] = i
                    break
                 
        # Main Loop
        size_ls_pre = 0
        iptt_next = 0
        
        for ih in range(n_h - 1, -1, -1):
            
            #---- if no new data appear ------------------------------------------ #
            if size_ls_pre == size_levelset[ih]:
                spM[ih, :] = spM[ih + 1, :]
                sM[ih, :] = sM[ih + 1, :]
                rank_level[ih] = rank_level[ih + 1]
                continue
            
            #---- Process ------------------------------------------------------ #
            size_ls_cur = size_levelset[ih]
            
            for ir in range(n_r):
                print('ih=%d, ir=%d' % (ih, ir))
                
                if ir == 1:
                    # -- Compute new bases
                    ds = self.nbr_clustering(range(size_ls_cur), thres_dist[ir])
                    sM[ih, ir] = ds.n_subset
                    
                    rank_cur = size_ls_cur
                    rank_level[ih] = rank_cur
                    
                    # -- Compute the representation
                    subptt_cur = np.zeros((ds.n_subset, rank_cur), dtype=np.bool)
                    clsts = ds.get_clsts()
                    initial_elms = clsts.keys()
                    for i in range(ds.n_subset):
                        subptt_cur[i, clsts[initial_elms[i]]] = True
                        
                else:
                    # -- Compute merge graph
                    mG, initial_elms = self.get_merge_graph(ds, thres_dist[ir])
                    
                    # -- Clusters in the merge graph
                    if not mG.any():
                        spM[ih, ir] = spM[ih, ir - 1]
                        sM[ih, ir] = ds.n_subset
                        continue  # continue to next ir
                    else:
                        merge_pairs = np.transpose(np.nonzero(mG))
                        for i, j in merge_pairs:
                            ds.merge_byelm(initial_elms[i], initial_elms[j])
                        sM[ih, ir] = ds.n_subset
                        
                        subptt_cur = np.zeros((ds.n_subset, rank_cur), dtype=np.bool)
                        clsts = ds.get_clsts()
                        initial_elms = clsts.keys()
                        for i in range(ds.n_subset):
                            subptt_cur[i, clsts[initial_elms[i]]] = True
                            
                # -- Check with the subptt below            
                # need to check all subptt of the same size unless an equivlent one been found
                ih_to_test = ih + 1
                while ih_to_test < n_h:
                    # find the subpatt of the same size
                    ir_to_test = search_first_true(sM[ih_to_test, :] == ds.n_subset)                    
                    if ir_to_test is not None:
                        subptt_cur_adj = adjust_rank(subptt_cur,rank_level[ih_to_test])
                        subptt_pre_adj = adjust_rank(subptts[spM[ih_to_test,ir_to_test]],rank_level[ih_to_test])
                        if not np.logical_xor(subptt_cur_adj,subptt_pre_adj).any():
                            spM[ih,ir] = spM[ih_to_test,ir_to_test]
                            # need to update the previous subptt
                            subptts[spM[ih_to_test,ir_to_test]] = subptt_cur
                            break
                    ih_to_test += 1
                    
                # if still not assigned, must be a new one
                if spM[ih,ir]==0:
                    subptts.append(subptt_cur)
                    spM[ih,ir] = iptt_next
                    iptt_next += 1
                    
                #-- Update
                #subptt_pre = subptt_cur
                
            #---- Update levelset size ------------------------------------#
            size_ls_pre = size_ls_cur
        
        
        ## Run through the subptts to get the life time of each basis
        blife = np.zeros(self.n_data)
        for ih in range(n_h-1,-1,-1):
            for ir in range(n_r):
                sptt = subptts[spM[ih,ir]]
                for i in range(sM[ih,ir]):
                    idx = search_first_true(sptt[i,:]>0)
                    blife[idx]+=1
        
        
        ## Data to return
        self.n_h = n_h
        self.n_r = n_r
        self.subptts = subptts
        self.rank_level = rank_level
        self.spM = spM
        self.sM = sM
        self.blife = blife
        self.size_levelset = size_levelset
        
        return self
    
    def nbr_clustering(self, I, r):
        ''' Epsilon ball clustering 
        
        Parameters:
        ----------
        I: int, array
            The indices of elements to be clustered
        r: float, scale
            The distance threshold
        '''

        ds = DisjointSets(I)
        n_data = len(I)
        for i in range(n_data):
            for j in range(i + 1, n_data):
                if self.distM[I[i], I[j]] <= r:
                    ds.merge_byelm(I[i], I[j])
                    
        return ds

    def get_merge_graph(self, ds, r):
        ''' Compute the merging information of the current partition under new threshold r
        
        Parameters:
        ----------
        ds: DisjointSets
        r: float, scale
            The new distance threshold
        '''
        
        mG = np.zeros((ds.n_subset, ds.n_subset))
        clsts = ds.get_clsts()
        initial_elms = np.sort(clsts.keys())
        
        for i in range(ds.n_subset):
            clst_1 = clsts[initial_elms[i]]
            for j in range(i + 1, ds.n_subset):
                clst_2 = clsts[initial_elms[j]]
                if (self.distM[[[k] for k in clst_1], clst_2] <= r).any():
                    mG[i, j] = True
                
        return mG, initial_elms
    
