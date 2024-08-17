import numpy as np
from .edit_distance import edit_similarity

class directed_potential_graph:
    def __init__(self, m:np.array, max_hop:int):
        # extract edges from adjacency matrix
        self.n:int = m.shape[0]
        self.edges = {}
        for i in range(self.n):
            self.edges[i] = []
            for j in range(self.n):
                if m[i,j] > 0:
                    self.edges[i].append(j)
        
        # extract degrees from edges
        degrees = np.zeros(shape=(self.n,), dtype=np.int32)
        for i in range(self.n):
            if i in self.edges:
                degrees[i] = len(self.edges[i])
        
        # compute hop0 ranks
        sorted_i_len_descendingly = np.argsort(degrees)[::-1]
        hop0_ranks = np.zeros(shape=(self.n,), dtype=np.int32)
        for i in range(self.n):
            if i > 0:
                if degrees[sorted_i_len_descendingly[i]] > degrees[sorted_i_len_descendingly[i - 1]]:
                    hop0_ranks[sorted_i_len_descendingly[i]] = i
                else:
                    hop0_ranks[sorted_i_len_descendingly[i]] = hop0_ranks[sorted_i_len_descendingly[i - 1]]
            else:
                hop0_ranks[sorted_i_len_descendingly[i]] = i
        
        # compute rank-id bidirectional mappings
        self.rank_to_ids = {}
        self.id_to_ranks = hop0_ranks
        for i in range(self.n):
            rank = hop0_ranks[i]
            if rank not in self.rank_to_ids:
                self.rank_to_ids[rank] = []
            self.rank_to_ids[rank].append(i)
        
        # init hop and seq
        self.current_hop = 0
        self.current_seq = np.zeros(shape=(self.n,), dtype=np.int32)
        self.current_subseq_root = np.zeros(shape=(self.n,), dtype=np.int32)
        self.current_subseq_pos = np.zeros(shape=(self.n,), dtype=np.int32)
        self.current_subseq_len = np.zeros(shape=(self.n,), dtype=np.int32)
        for i in range(self.n):
            self.current_seq[i] = degrees[sorted_i_len_descendingly[i]]
            self.current_subseq_root[i] = sorted_i_len_descendingly[i]
            self.current_subseq_pos[sorted_i_len_descendingly[i]] = i
            self.current_subseq_len[i] = 1
        
        # expand seq by increasing hop
        for _ in range(max_hop):
            self.inc_hop()

    def inc_hop(self):
        # compute new subseq length and new seq length
        new_subseq_len = np.zeros(shape=(self.n,), dtype=np.int32)
        new_seq_len = 0
        for i in range(self.n):
            new_subseq_len[i] = 1
            for j in self.edges[i]:
                new_subseq_len[i] += self.current_subseq_len[j]
            new_seq_len += new_subseq_len[i]

        # expand seq
        new_seq = np.zeros(shape=(new_seq_len,),dtype=np.int32)
        new_subseq_pos = np.zeros(shape=(self.n,), dtype=np.int32)
        pos:int = 0
        for idx in range(self.n):
            i = self.current_subseq_root[idx]
            new_subseq_pos[i] = pos

            # sort children before making subseq expansion
            subseq_children_id = np.zeros(shape=(len(self.edges[i]),), dtype=np.int32)
            subseq_children_rank = np.zeros(shape=(len(self.edges[i]),), dtype=np.int32)
            for j in range(len(self.edges[i])):
                neighbor = self.edges[i][j]
                subseq_children_id[j] = neighbor
                subseq_children_rank[j] = self.id_to_ranks[neighbor]
            permutation = np.argsort(subseq_children_rank)[::-1]
            sorted_subseq_children_id = np.zeros(shape=(len(self.edges[i]),), dtype=np.int32)
            for j in range(len(self.edges[i])):
                sorted_subseq_children_id[j] = subseq_children_id[permutation[j]]

            # expand subseq according to sorted_subseq_children_id
            new_seq[pos] = len(self.edges[i])
            pos += 1
            for j in range(len(self.edges[i])):
                neighbor = sorted_subseq_children_id[j]
                neighbor_old_subseq_len = self.current_subseq_len[neighbor]
                neighbor_old_subseq_pos = self.current_subseq_pos[neighbor]
                for k in range(neighbor_old_subseq_len):
                    new_seq[pos] = self.current_seq[neighbor_old_subseq_pos + k]
                    pos += 1
            
        # compute new ranks
        new_rank_to_ids = {}
        new_id_to_ranks = np.zeros(shape=(self.n,), dtype=np.int32)
        for old_rank in self.rank_to_ids:
            ids = self.rank_to_ids[old_rank]
            if len(ids) <= 1:
                for id in ids:
                    new_rank_to_ids[old_rank] = [id]
            else:
                for id1 in ids:
                    drank = 0
                    pos1 = new_subseq_pos[id1]
                    len1 = new_subseq_len[id1]
                    for id2 in ids:
                        if id2 == id1:
                            continue
                        pos2 = new_subseq_pos[id2]
                        len2 = new_subseq_len[id2]
                        min_len = min(len1,len2)
                        
                        cmp = 0
                        for j in range(min_len):
                            if new_seq[pos1 + j] > new_seq[pos2 + j]:
                                cmp = 1
                                break
                            elif new_seq[pos1 + j] < new_seq[pos2 + j]:
                                cmp = -1
                                break
                        if cmp == 0:
                            if len1 > len2:
                                cmp = 1
                            elif len1 < len2:
                                cmp = -1
                        
                        if cmp == -1:
                            drank += 1
                        
                    new_rank = old_rank + drank
                    if new_rank not in new_rank_to_ids:
                        new_rank_to_ids[new_rank] = []
                    new_rank_to_ids[new_rank].append(id1)
                    new_id_to_ranks[id1] = new_rank

        # update rank-id mappings    
        self.rank_to_ids = new_rank_to_ids
        self.id_to_ranks = new_id_to_ranks

        # reorder new seq
        self.current_seq = np.zeros(shape=(new_seq_len,),dtype=np.int32)
        self.current_subseq_len = new_subseq_len
        self.current_subseq_pos = np.zeros(shape=(self.n,), dtype=np.int32)
        self.current_subseq_root = np.zeros(shape=(self.n,), dtype=np.int32)
        pos = 0
        i = 0
        for rank in self.rank_to_ids:
            for root in self.rank_to_ids[rank]:
                self.current_subseq_root[i] = root
                self.current_subseq_pos[root] = pos
                i += 1

                src_pos = new_subseq_pos[root]
                for j in range(self.current_subseq_len[root]):
                    self.current_seq[pos] = new_seq[src_pos]
                    pos += 1
                    src_pos += 1
        self.current_hop += 1
        # print(f'n = {self.n}, hop = {self.current_hop}, seqLen = {self.current_seq.shape[0]}')

    def get_similarity(self, other):
        sim = 0.0
        s1 = self.current_seq
        s2 = other.current_seq
        return edit_similarity(s1, s2)