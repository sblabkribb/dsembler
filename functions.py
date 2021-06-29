from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio import SeqIO
import collections
import re
from Bio import Seq
from Bio import Align
import itertools

aligner = Align.PairwiseAligner()

class TypeError(Exception):
    pass

class DnaAssembly:

    # Reading fasta file to give the DNA sequence 
    # input FASTA file
    def read_fasta(self, file): 

        dna_seq = SeqIO.parse(file, "fasta")
        
        for seq_record in dna_seq:
            seq = str(seq_record.seq).upper() # ensuring all data is capitalized for uniformity
        # returns string
        return seq


    # splicing into oligomers
    # input gene sequence (str), target oligomer size(int), and target overlap size(int)
    def oligomer_splice(self, gene_seq, oligomer_size, overlap_size): 
        global rough_oligo_size

        rough_oligo_size = oligomer_size - overlap_size # calculates the size of oligomers without overlaps
        rough_oligo_list = []

        for start_pos in range(0, len(gene_seq), rough_oligo_size):
            rough_oligo = gene_seq[start_pos : start_pos + rough_oligo_size]

            if rough_oligo:
                rough_oligo_list.append(rough_oligo)
        # returns list of oligomer sequences without overlaps
        return rough_oligo_list


    # generate all possible oligomers within a given overlap range
    # input gene sequence (str), rough oligomers (list), target oligomer size (int), and target overlap size (int)
    def possible_oligomers(self, gene_seq, rough_oligo_list, oligomer_size, overlap_size):

        list_of_oligomers,overlap_length = [], []

        low_overlap = int(overlap_size - (overlap_size/4)) # calculates the shortest possible overlap length
        high_overlap = int(overlap_size * 3 / 2) # calculates the longest possible overlap length
        seq = gene_seq[rough_oligo_size:]

        for i in range(len(rough_oligo_list) - 1):
            possible_overlaps, possible_overlaps_len,= [], []

            for overlap in range(low_overlap, high_overlap):
                overlap_seq = seq[:overlap]
                oligomer = rough_oligo_list[i] + overlap_seq
                possible_overlaps.append(oligomer)
                possible_overlaps_len.append(overlap)

            list_of_oligomers.append(possible_overlaps)
            overlap_length.append(possible_overlaps_len)
            seq = seq[rough_oligo_size:]
        # returns two separate lists of lists that 1) contain oligomers with varying overlap lengths 2) their respective overlap lengths
        return list_of_oligomers, overlap_length

    # function for generating alternative complementary sequences
    # input a list of sequences
    def complement_oli(self, data):

        comp_list = []

        for i in data:

            if data.index(i) % 2 != 0: # allow alternate sequences to be complementary 
                i = str(Seq.Seq(i).complement())
                comp_list.append(i)
            else:
                comp_list.append(i)
        # returns a list containing the complementary sequences in order
        return comp_list


    # returns value in list closest to n
    # input a list and a number you want to check
    def closest(self, data, n):
        # returns int
        return data[min(range(len(data)), key=lambda i: abs(data[i] - n))]


    # selecting oligomers that fall in target requirements (Melting temp, GC content, overlap length, GC clamp)
    # input the list of all possible oligomers, and corresponding overlap lengths, target melting temp (int), temperature range (int), o
    def optimal_oligomers(self, list_of_oligomers, optimal_temp, temp_range, overlap_length):

        low_tm = optimal_temp - temp_range # lower melting temperature
        high_tm = optimal_temp + temp_range # higher melting temperature
        low_gc = 20 # lower gc content limit
        high_gc = 80 # higher gc content limit

        optimal_oligomers, optimal_overlap_len = [], [] 

        for oligo_index in range(len(list_of_oligomers)):
            # removes oligomers if there is T at 3' (alternate sequences will be complementary in the output)
            if oligo_index % 2 == 0: 
                t_3_oligo_index = [list_of_oligomers[oligo_index].index(oligo) for oligo in list_of_oligomers[oligo_index] if oligo[-1] == 'T']
                if len(t_3_oligo_index) != len(list_of_oligomers[oligo_index]):
                    for ele in sorted(t_3_oligo_index, reverse=True):
                        del list_of_oligomers[oligo_index][ele]
                        del overlap_length[oligo_index][ele]
                else:
                    for ele in sorted(t_3_oligo_index[1:], reverse=True):
                        del list_of_oligomers[oligo_index][ele]
                        del overlap_length[oligo_index][ele]
            else:
                while list_of_oligomers[oligo_index][0][0] == 'A':
                    list_of_oligomers[oligo_index][:] = [i[1:] for i in list_of_oligomers[oligo_index]]
                    overlap_length[oligo_index - 1][:] = [number - 1 for number in overlap_length[oligo_index - 1]]
                    
        comp_list = []

        # generating alternate complementary sequences to select based on melting temperature and GC
        for oligomer_pool in list_of_oligomers:
            comp_list.append(self.complement_oli(oligomer_pool))
            
        for oligo_num in range(len(list_of_oligomers)):
            tm, oligo, over_len = [], [], []
            for rough_oligo_index in range(len(list_of_oligomers[oligo_num])):
                overlap = comp_list[oligo_num][rough_oligo_index][-overlap_length[oligo_num][rough_oligo_index]:] # extract overlap sequence
                melting_temp = round(mt.Tm_NN(overlap, nn_table=mt.DNA_NN3), 2) # check melting temperature based on Nearest Neighbour equation
                gc = GC(overlap)
                tm.append(melting_temp) # to keep track of all melting temperatures of overlaps
                if low_gc <= gc <= high_gc and low_tm <= melting_temp <= high_tm: 
                    oligo.append(list_of_oligomers[oligo_num][rough_oligo_index])
                    over_len.append(overlap_length[oligo_num][rough_oligo_index])
            if not oligo: # if the list 'oligo' is empty
                # chooses the overlap with the closest melting temperature to the user's input
                if all(temp >= high_tm for temp in tm) == True: 
                    next_best_oligo = self.closest(tm, high_tm)
                else:
                    next_best_oligo = self.closest(tm, low_tm)
                oligo.append(list_of_oligomers[oligo_num][tm.index(next_best_oligo)])
                over_len.append(overlap_length[oligo_num][tm.index(next_best_oligo)])
            optimal_oligomers.append(oligo)
            optimal_overlap_len.append(over_len)
        # returns 2 list of lists with all the oligomers that fit/ almost fit the criteria specified by user
        return optimal_oligomers, optimal_overlap_len


    # selects the smallest oligomer from the possible optimal oligomers (cost effective)
    # input the list of lists of optimal oligomers and corresponding overlap lengths
    def final_oligomers(self, optimal_oligomers, optimal_overlap_len, rough_oligo_list):

        final_oligomers, final_overlap_len = [], []

        for oligomer_index in range(len(optimal_oligomers)):
            smallest_ovr = min(optimal_overlap_len[oligomer_index])
            smallest_oli = optimal_oligomers[oligomer_index][optimal_overlap_len[oligomer_index].index(smallest_ovr)]
            final_oligomers.append(smallest_oli)
            final_overlap_len.append(smallest_ovr)
        if len(rough_oligo_list[-1]) > optimal_overlap_len[0][-1]:
            final_oligomers.append(rough_oligo_list[-1])
            final_overlap_len.append(0)
        # returns 2 lists that contain the oligomer sequences and overlap lengths
        return final_oligomers, final_overlap_len


    # generating single sequence with all overlap regions
    # input the list of oligomers and their respective overlap length
    def overlap_list(self, oligomer_list, overlap_length):
        
        overlaps = []
        
        for oligomer in oligomer_list:
            overlap_len = overlap_length[oligomer_list.index(oligomer)]
            overlap = oligomer[-overlap_len:]
            overlaps.append(overlap)
        # returns a list of overlap sequences
        return overlaps


    # check overlap alignment of a given list
    # input a list of sequences
    def overlap_alignment(self, sequences):

        alignment_index = []
        
        for a, b in itertools.combinations(sequences, 2):
            align = aligner.align(a, b)
            align_ratio = align.score / max(len(a), len(b))
            if align_ratio > 0.7:
                index = [sequences.index(a), sequences.index(b)]
                alignment_index.append(index)
        # return a list of list containing the indexes of high levels of alignments between two sequences
        return alignment_index


    # make clusters based on alignment score
    def make_clusters(self, data, overlap, cluster_size, cluster_range):
        if not len(data) in range(cluster_size):
            num = []
            l_cluster = cluster_size - cluster_range
            h_cluster = cluster_size + cluster_range + 1
            for i in range(l_cluster, h_cluster, 2):
                num.append(i)
            clusters = []
            h = data
            o = overlap
            while h:
                f = []
                for i in num:
                    j = h[:i]
                    ov = o[:i]
                    y = self.overlap_list(j, ov)
                    d = self.overlap_alignment(y)
                    f.append(d)
                    if f.count(f[0]) == len(f):
                        j = h[:cluster_size]
                        h = h[cluster_size:]
                        o = o[cluster_size:]
                        clusters.append(j)
                    elif any(f.count(element) > 1 for element in f):
                        t = min(f, key=lambda x: len(x))
                        if f.count(t) > 1:
                            c = [num[index] for index, element in enumerate(f) if element == t]
                            p = num[min(range(len(c)), key=lambda i: abs(c[i] - cluster_size))]
                        else:
                            p = num[f.index(t)]
                        j = h[:p]
                        h = h[p:]
                        o = o[p:]
                        clusters.append(j)
        else:
            clusters = [data]
        return clusters


    # orientates sequences in 5 to 3 format
    # input a list of sequences
    def seq_orientation(self, sequences):

        five_to_three = []
        
        for sequence in sequences:
            if sequences.index(sequence) % 2 != 0:
                rc_sequence = str(Seq.Seq(sequence).reverse_complement())
                five_to_three.append(rc_sequence)
            else:
                five_to_three.append(sequence)
        # list of oligomers with the right orientation
        return five_to_three


    # make clusters based on alignment score
    # input list of oligomers and their respective overlap lengths; the (int) user input of cluster size and cluster range
    def complementary_clusters(self, oligomers, overlap_len, cluster_size, cluster_range):
        
        num = []

        low_cluster = cluster_size - cluster_range
        high_cluster = cluster_size + cluster_range + 1
        
        for cluster_len in range(low_cluster, high_cluster, 2):
            num.append(cluster_len)
        
        if not len(oligomers) in range(high_cluster) :
        
            comp_list, five_to_three = self.complement_oli(oligomers), self.seq_orientation(oligomers)
            
            clusters, comp_clusters, cluster_ovr, cluster_five_two_three = [], [], [], []

            oligos = oligomers
            overlap = overlap_len
        
            while oligos:
                alignments = []

                for cluster_length in num:
                    alignment = self.overlap_alignment(self.overlap_list(oligos[:cluster_length], overlap[:cluster_length]))
                    alignments.append(alignment)

                    if alignments.count(alignments[0]) == len(alignments):
                        oligomers = oligos[:cluster_size]
                        overlap_len = overlap[:cluster_size]
                        
                        oligos = oligos[cluster_size:]
                        overlap = overlap[cluster_size:]

                        complementary_cluster = comp_list[:cluster_size]
                        final_cluster = five_to_three[:cluster_size]
                        comp_list = comp_list[cluster_size:]
                        five_to_three = five_to_three[cluster_size:]

                        clusters.append(oligomers)
                        comp_clusters.append(complementary_cluster)
                        cluster_ovr.append(overlap_len)
                        cluster_five_two_three.append(final_cluster)

                    elif any(alignments.count(element) > 1 for element in alignments):
                        min_align = min(alignments, key=lambda x: len(x))
                        if alignments.count(min_align) > 1:
                            c = [num[index] for index, element in enumerate(alignments) if element == min_align]
                            cluster_dim = num[min(range(len(c)), key=lambda cluster_length: abs(c[cluster_length] - cluster_size))]
                        else:
                            cluster_dim = num[alignments.index(min_align)]

                        oligomers = oligos[:cluster_dim]
                        overlap_len = overlap[:cluster_dim]

                        oligos = oligos[cluster_dim:]
                        overlap = overlap[cluster_dim:]

                        complementary_cluster = comp_list[:cluster_dim]
                        final_cluster = five_to_three[:cluster_dim]
                        comp_list = comp_list[cluster_dim:]
                        five_to_three = five_to_three[cluster_dim:]

                        clusters.append(oligomers)
                        cluster_ovr.append(overlap_len)
                        comp_clusters.append(complementary_cluster)
                        cluster_five_two_three.append(final_cluster)
        # if the cluster size inputted is larger than the final cluster
        else:
            clusters = [oligomers]
            comp_clusters = [self.complement_oli(oligomers)]
            cluster_ovr = [overlap_len]
            cluster_five_two_three = [self.seq_orientation(oligomers)]

        # returns three list of lists containing oligomers and their respective overlap lengths in appropriate clusters sizes
        return clusters, comp_clusters, cluster_ovr, cluster_five_two_three

    
    # generates scores and suggests possible areas of error for each oligomer
    def overlap_score(self, data, clusters, cluster_5_3, overlap, optimal_temp, temp_range):
        score, fault, repeats, repeat_seq, repeat_position = [], [], [], [], []
        for i in range(len(data)):
            s, f, r, rs, rp = [], [], [], [], []
            for x in range(len(data[i])):
                s.append(0)
                f.append([])
                rp.append([])
            score.append(s)
            fault.append(f)
            repeats.append(r)
            repeat_seq.append(rs)
            repeat_position.append(rp)
        cluster_length = []
        if len(data) == 1:
            cluster_length.append(0)
        else:
            for i in range(len(data)):
                cluster_length.append(i)
        for cluster in cluster_length:
            g = []
            for x in range(len(data[cluster])):
                overla = overlap[cluster][x]
                data_to_put = data[cluster][x][-overla:]
                g.append(data_to_put)
                d = "".join(g)
                s = self.repeat_seq(d)
            for i in s:
                for a in i:
                    seq_repeat = d[a:a + 11]
                    match_index = [data[cluster].index(h) for h in data[cluster] if seq_repeat in h]
                if len(match_index) > 1:
                    for l in match_index:
                        if overlap[cluster][l] != 0:
                            score[cluster][l] += 10
                            fault[cluster][l] += "R"
                    repeats[cluster].append(match_index)
                    repeat_seq[cluster].append(seq_repeat)
            for x in range(len(repeat_position[cluster])):
                if overlap[cluster][x] != 0:
                    for y in repeats[cluster]:
                        for u in y:
                            if x == u:
                                repeat_position[cluster][x].append(repeat_seq[cluster][repeats[cluster].index(y)])
        low_temp = optimal_temp- temp_range
        high_temp = optimal_temp + temp_range
        for cluster in range(len(clusters)):
            for oligomer in range(len(clusters[cluster])):
                overlap_size = overlap[cluster][oligomer]
                clust = clusters[cluster][oligomer]
                overlap_mt = round(mt.Tm_NN(clust[-overlap_size:], nn_table=mt.DNA_NN3), 2)
                if overlap_size != 0:
                    if low_temp > overlap_mt:
                        score[cluster][oligomer] += round(abs(low_temp - overlap_mt), 2)
                        fault[cluster][oligomer] += "L"
                    elif overlap_mt > high_temp:
                        score[cluster][oligomer] += round(abs(overlap_mt - high_temp), 2)
                        fault[cluster][oligomer] += "H"
                        if overlap_size > int(overlap_size - (overlap_size/4)):
                            if cluster_5_3[cluster][oligomer][-2] == "T":
                                score[cluster][oligomer] += 1
                                fault[cluster][oligomer] += "T"
                    gc = ["CCC", "GGG", "CCG", "CGC", "GCC", "GGC", "GCG", "CGG"]
                    if any((t in cluster_5_3[cluster][oligomer][-5:]) for t in gc):
                        score[cluster][oligomer] += 1
                        fault[cluster][oligomer] += "G"
        return score, fault, repeat_position


 # finding repeat sequences
    # https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
    def repeat_seq(self, data):
        import pdb;pdb.set_trace()
        repeats = 10
        x = []
        seq_no = len(data) - (repeats - 1)
        seq_list = []
        for j in range(seq_no):
            start_pos = j
            end_pos = j + repeats
            b = data[start_pos:end_pos]
            seq_list.append(b)
        s = [item for item, count in collections.Counter(seq_list).items() if count > 1]
        for repeat in s:
            matches = re.finditer(repeat, data)
            matches_position = [match.start() for match in matches]
            x.append(matches_position)
        return x



# def output_webpage(self, data, overlap):
    #     spaces = []
    #     overshot = []
    #     z = 0
    #     for j in range(len(overlap)):
    #         if j == 0:
    #             t = len(data[j]) - overlap[j]
    #         else:
    #             t = len(data[j]) - overlap[(j - 1)] - overlap[j]
    #         if t < 0:
    #             x = 0 - t
    #             t = 0
    #             z += 1
    #         else:
    #             x = 0
    #             z += 0
    #         spaces.append(t)
    #         overshot.append(x)
    #     even, odd_space, odd_overshoot = [], [], []
    #     odd, even_space, even_overshoot = [], [], []
    #     for j in range(len(data)):
    #         if j % 2 != 0:
    #             odd.append(data[j])
    #             even_space.append(spaces[j])
    #             even_overshoot.append(overshot[j])
    #         else:
    #             even.append(data[j])
    #             odd_space.append(spaces[j])
    #             odd_overshoot.append(overshot[j])
    #     forward_seq = even[0]
    #     reverse_seq = ""
    #     for w in range(len(even)):
    #         if w != range(len(even))[-1]:
    #             if even_overshoot[w] == 0:
    #                 forward_seq = forward_seq + (" " * even_space[w]) + even[w + 1]
    #             else:
    #                 b = even[w][-even_overshoot[w]:]
    #                 forward_seq = forward_seq[:-even_overshoot[w]] + b.swapcase() + even[w + 1][even_overshoot[w]:]
    #         else:
    #             forward_seq = forward_seq + even[w]
    #     for v in range(len(odd)):
    #         if odd_overshoot[v] == 0:
    #             reverse_seq = (reverse_seq + " " * odd_space[v] + odd[v])
    #         else:
    #             c = odd[v - 1][-odd_overshoot[v]:]
    #             reverse_seq = (reverse_seq[:-odd_overshoot[v]] + c.swapcase() + odd[v][odd_overshoot[v]:])
    #     return forward_seq, reverse_seq, z
