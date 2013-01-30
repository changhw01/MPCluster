import numpy as np

class DisjointSets:
    ''' The famous disjoint set systems supporting UnionFind operation.
        The internal data structure is disjoint trees.

    Member Data: 
        Internal:
        __elm2idx: Convert the original element to the internal index
        __parent: The __parent of the given node idx in the disjoint trees
        __subsets: The current snapshot of current partition of data
        __size_subsets: The size of this subset
        
        Public:
        n_subset: The current number of __subsets
        n_elm: The size of S
        
    Member Method:
        Internal:
        __init__:
        __find__:
        __union__:
        
        Public:
        add_set:
        add_sets:
        in_same_set:
        merge_byelm:
        update_cidx:
        get_clsts:


    Description:
        elm: can be anything, such as strings or even objects

    '''

    def __init__(self, S=None):

        # initialize by a list of elements
        if isinstance(S,list) or isinstance(S,np.ndarray):
            
            self.n_elm = len(S)
            self.n_subset = self.n_elm

            self.__elm2idx = dict([(S[i],i) for i in xrange(self.n_elm)])
            self.__parent = range(self.n_elm)            
            self.__subsets = dict([(i, [S[i]]) for i in xrange(self.n_elm)])
            self.__size_subsets = dict([(i,1) for i in xrange(self.n_elm)])

        # initialize by the number of elements
        elif isinstance(S,int):
            
            self.n_elm = S
            self.n_subset = self.n_elm

            self.__elm2idx = dict([(i,i) for i in xrange(self.n_elm)])
            self.__parent = range(self.n_elm)
            self.__subsets = dict([(i,[i]) for i in xrange(self.n_elm)])
            self.__size_subsets = dict([(i,1) for i in xrange(self.n_elm)])

        # initialize by another DisjointSets
        elif isinstance(S,DisjointSets):
            self.n_elm = S.n_elm
            self.n_subset = S.n_subset
        
            self.__elm2idx = dict(S.__elm2idx)
            self.__parent = list(S.__parent)
            self.__subsets = dict(S.__subsets)
            self.__size_subsets = dict(S.__size_subsets)

        # initialize as empty 
        else:
            self.n_elm = 0
            self.n_subset = 0

            self.__elm2idx = dict()
            self.__parent = []
            self.__subsets = dict()
            self.__size_subsets = dict()

    def __find__(self, idx):
        ''' Return the root index of the given internal element index

            Will flatten the tree at the same time
        '''
        idx_p = self.__parent[idx]
        idx_gp = self.__parent[idx_p]
        while idx_p!=idx_gp:
            idx_p = idx_gp
            idx_gp = self.__parent[idx_gp]
        
        self.__parent[idx] = idx_p

        return idx_p

    def __union__(self, idx1, idx2, criterion="small_index"):
        ''' Merge the clusters containing the two given elm indices '''
        
        ridx1 = self.__find__(idx1); ridx2 = self.__find__(idx2)
        if ridx1 != ridx2:
            # merge ridx2 to ridx1
            if criterion=="size":
                if self.__size_subsets[ridx1]<self.__size_subsets[ridx1]:
                    temp = ridx1; ridx1=ridx2; ridx2=temp
            else:
                if ridx2<ridx1:
                    temp = ridx1; ridx1=ridx2; ridx2=temp

            self.__parent[ridx2] = ridx1
            self.__subsets[ridx1] = self.__subsets[ridx1]+self.__subsets[ridx2]
            del self.__subsets[ridx2]
            self.__size_subsets[ridx1] = self.__size_subsets[ridx1]+self.__size_subsets[ridx2]
            del self.__size_subsets[ridx2]
            self.n_subset -= 1
        
    def merge_byelm(self, elm1, elm2):
        idx1 = self.__elm2idx[elm1] 
        idx2 = self.__elm2idx[elm2]
        self.__union__(idx1, idx2)
    
    def add_set(self,S):
        ''' Add a new set to the collection. 
            The new elements shouldn't exist in the current collection
            
            Return the index of the new set
        '''
       
        # make it iterable
        if not isinstance(S,list):
            S = list(S)

        # parameters
        # 'self.n_elm' is the idx of the first new element in the collection
        # it will also be the new ridx
        ridx_newset = self.n_elm
        size_newset = len(S)

        # test whether any element of S has already existed
        for elm in S:
            if elm in self.__elm2idx:
                print 'Some element already exists!'
                return None

        # update
        self.__parent.extend([ridx_newset]*size_newset)
        self.__subsets[ridx_newset] = S
        self.__size_subsets[ridx_newset] = size_newset
        self.__elm2idx.update([(S[i],self.n_elm+i) for i in xrange(size_newset)])

        self.n_subset += 1
        self.n_elm += size_newset

        return self.__ridx2cidx[ridx_newset]

    def add_sets(self, S_list):
        ''' Add a list of new sets to the collection.
        '''

        for i in range(len(S_list)):
            if not isinstance(S_list[i],list):
                S_list[i] = list(S_list[i])

        # test whether any element exists in more than one subset
        idx_nextelm = self.n_elm
        temp_dict = dict()
        for S in S_list:
            for elm in S:
                if elm in (self.__elm2idx) or (elm in temp_dict):
                    print 'Some element already exists!'
                    return False
                temp_dict[elm] = idx_nextelm
                idx_nextelm += 1

        # update 
        # the idx of the first element in the collection will be the subset idx in
        # the collection
        self.__parent.extend([temp_dict[S[0]] for S in S_list for s in xrange(len(S))])
        self.__elm2idx.update(temp_dict)
        
        rv = range(self.n_elm,self.n_elm+len(S_list))
        for S in S_list:
            self.__subsets[temp_dict[S[0]]] = S
            self.__size_subsets[temp_dict[S[0]]] = len(S)

        self.n_elm = len(self.__elm2idx)
        self.n_subset = len(self.__subsets)

        return rv
    
    def in_same_set(self, elm1, elm2):
        return self.__find__(self.__elm2idx[elm1])==self.__find__(self.__elm2idx[elm2])
    
    def get_csizes(self):
        return dict([(ridx,s) for ridx,s in self.__size_subsets.iteritems()])
    
    def get_clsts(self):
        ''' Return all the clusters '''
        
        return self.__subsets