from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio import SeqIO
import collections
import re
from Bio import Seq
from Bio import Align
import itertools

aligner = Align.PairwiseAligner()

# about
'''
We use five parent classes: Sequence Processing, Oligomers, OligomerSelection, and Cluster. Each class works with a specific type of data and/or give a similar output.
Sequence Processing: Performs specific methods on a DNA sequences present in list or nested list. 
                     Includes methods like generating complementary and 5'->3' oligomers, extracting overlap sequences, identify sequence similarity and repeats between two sequences in a list.
Oligomers: Has two subclasses with the same methods (but different working): LinearAssembly and CircularAssembly. 
          Both have two methods that generate oligomers: oligomer_splice and possible_oligomers.
OligomerSelection: This class is to select for appropriate oligomers based on the user's parameters. It has a child class of Score to characterize oligomer selections.
                   Oligomers can be selected based on the absence of Thymine on the 3' end, melting temperature and GC content of the overlap region, and smallest oligomer in a given list of oligomers.
Cluster: Here a given set of oligomers can be divided into clusters based on the sequence similarities within each cluster.
'''
# file architecture
'''
assembly.py
-| Oligomers
    -| LinearAssembly
        -| oligomer_splice(self)
        -| possible_oligomers(self, rough_oligo_list)
    -| CircularAssembly
        -| @staticmethod rough_oligo_size_cir(gene_seq, rough_oligo_size)
        -| oligomer_splice(self)
        -| possible_oligomers(self, rough_oligo_list)
-| OligomerSelection
    -| @staticmethod closest(list, n)
    -| UserParameterSelection
        -| t_3_free(self, list_of_oligomers, list_of_overlap_lengths)
        -| gc_tm_optimal(self, list_of_oligomers, list_of_overlap_lengths)
        -| smallest_oligomers(self, list_of_oligomers, list_of_overlap_lengths, rough_oligo_list)
        -| smallest_oligomers(self, list_of_oligomers, list_of_overlap_lengths, rough_oligo_list, rough_oligomer_size, high_overlap)
    -| Score
        -| overlap_score(self, cluster, complementary_cluster, 5_to_3_cluster, overlap_length)
-| Clusters
    -| complementary_clusters(self, oligomers, overlap_length) 
-| SequenceProcessing
    -| @staticmethod complement_oli(data)
    -| @staticmethod overlap_list(data)
    -| @staticmethod overlap_alignment(sequences)
    -| @staticmethod seq_orientation(sequences)
    -| @staticmethod repeat_seq(data)
'''

class TypeError(Exception):
    pass

# Reading fasta file to give the DNA sequence
class FastaFile:
    # input FASTA file
    def __init__(self, file):
        self.file = file

    def read_fasta(self):

        seq = None

        dna_seq = SeqIO.parse(self.file, "fasta")

        for seq_record in dna_seq:
            seq = str(seq_record.seq).upper()  # ensuring all data is capitalized for uniformity
        # returns a string
        return seq

# processing DNA sequences within lists
class SequenceProcessing:
    # function for generating alternative complementary sequences
    # input a list of sequences
    @staticmethod
    def complement_oli(data):

        comp_list = []

        for i in data:
            if data.index(i) % 2 != 0:  # allow alternate sequences to be complementary
                i = str(Seq.Seq(i).complement())
                comp_list.append(i)
            else:
                comp_list.append(i)
        # returns a list containing the complementary sequences in order
        return comp_list

    # generating single list with all overlap regions
    # input the list of oligomers and their respective overlap length
    @staticmethod
    def overlap_list(oligomer_list, overlap_length):

        overlaps = []

        for oligomer in oligomer_list:
            overlap_len = overlap_length[oligomer_list.index(oligomer)]
            overlap = oligomer[-overlap_len:]
            overlaps.append(overlap)
        # returns a list of overlap sequences
        return overlaps

    # check overlap alignment of a given list
    # input a list of sequences
    @staticmethod
    def overlap_alignment(sequences):

        alignment_index = []

        for a, b in itertools.combinations(sequences, 2):
            align = aligner.align(a, b)
            align_ratio = align.score / max(len(a), len(b))
            if align_ratio > 0.7:
                index = [sequences.index(a), sequences.index(b)]
                alignment_index.append(index)
        # return a list of list containing the indexes of high levels of alignments between two sequences
        return alignment_index

    # orientates sequences in 5 to 3 format
    # input a list of sequences
    @staticmethod
    def seq_orientation(sequences):

        five_to_three = []

        for sequence in sequences:
            if sequences.index(sequence) % 2 != 0:
                rc_sequence = str(Seq.Seq(sequence).reverse_complement())
                five_to_three.append(rc_sequence)
            else:
                five_to_three.append(sequence)
        # list of oligomers with the right orientation
        return five_to_three

    # finding repeat sequences
    # https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
    @staticmethod
    def repeat_seq(data):
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

# generating oligomer sequences
class Oligomers:
    low_gc = 20  # lower gc content limit
    high_gc = 80  # higher gc content limit

    def __init__(self, gene_seq, oligomer_size, overlap_size):
        self.gene_seq = gene_seq
        self.oligomer_size = oligomer_size
        self.overlap_size = overlap_size
        self.rough_oligo_size =  self.oligomer_size - self.overlap_size

# generating oligomer sequences for linear DNA sequences
class LinearAssembly(Oligomers):
    # splicing into oligomers
    # input gene sequence (str), target oligomer size(int), and target overlap size(int)
    def oligomer_splice(self):
          # calculates the size of oligomers without overlaps
        rough_oligo_list = []
        for start_pos in range(0, len(self.gene_seq), self.rough_oligo_size):
            rough_oligo = self.gene_seq[start_pos: start_pos + self.rough_oligo_size]

            if rough_oligo:
                rough_oligo_list.append(rough_oligo)
        # returns list of oligomer sequences without overlaps
        return rough_oligo_list


    def possible_oligomers(self, rough_oligo_list):

        list_of_oligomers, overlap_length = [], []

        low_overlap = int(self.overlap_size - (self.overlap_size / 4))  # calculates the shortest possible overlap length
        high_overlap = int(self.overlap_size * 3 / 2)  # calculates the longest possible overlap length
        seq = self.gene_seq[self.rough_oligo_size:]

        for i in range(len(rough_oligo_list) - 1):
            possible_overlaps, possible_overlaps_len, = [], []

            for overlap in range(low_overlap, high_overlap):
                overlap_seq = seq[:overlap]
                oligomer = rough_oligo_list[i] + overlap_seq
                possible_overlaps.append(oligomer)
                possible_overlaps_len.append(overlap)
            list_of_oligomers.append(possible_overlaps)
            overlap_length.append(possible_overlaps_len)
            seq = seq[self.rough_oligo_size:]
        # returns two separate lists of lists that 1) contain oligomers with varying overlap lengths 2) their respective overlap lengths
        return list_of_oligomers, overlap_length

# generating oligomer sequences for circular DNA sequences
class CircularAssembly(Oligomers):
    def __init__(self, gene_seq, oligomer_size, overlap_size):
        super().__init__(gene_seq, oligomer_size, overlap_size)
        self.rough_oligomer_size = self.rough_oligo_size_cir(self.gene_seq, self.rough_oligo_size)
        self.low_overlap = int(self.overlap_size - (self.overlap_size / 4))  # calculates the shortest possible overlap length
        self.high_overlap = int(self.overlap_size * 4 / 3)  # calculates the longest possible overlap length

    @staticmethod
    def rough_oligo_size_cir(gene_seq, rough_oligo_size):
        # find the even number when divided with the length of the given gene seq is closest to the user's inputed (rough oligomers)
        number = None
        for i in range(1, 100):
            number = len(gene_seq) / (i * 2)
            if number < rough_oligo_size:
                break
        oligomer_s = round(number)  # estimated approximate overlap size as calculated previously

        oligo_size = oligomer_s - 1  # calculates the size of oligomers without overlaps
        return oligo_size

    # splicing into oligomers
    # input gene sequence (str), target oligomer size(int), and target overlap size(int)
    def oligomer_splice(self):
        # calculates the size of oligomers without overlaps
        rough_oligo_list = []
        # generate rough oligomers as in linear
        for start_pos in range(0, len(self.gene_seq), self.rough_oligomer_size):
            rough_oligo = self.gene_seq[start_pos: start_pos + self.rough_oligomer_size]
            if rough_oligo:
                rough_oligo_list.append(rough_oligo)
        # add an overlap region between the last and first oligomer (to make it circular)
        if len(rough_oligo_list[-1]) < self.rough_oligomer_size:
            rough_oligo_list[-1] = rough_oligo_list[-1] + rough_oligo_list[0]
        else:
            rough_oligo_list.append(rough_oligo_list[0])
        # returns list of oligomer sequences without overlaps
        return rough_oligo_list

    def possible_oligomers(self, rough_oligo_list):

        # as we are using an approximate value of the user's specifications, we use 4/3 as the largest overlap possible
        list_of_oligomers, overlap_length = [], []
        sequence = ''.join(rough_oligo_list)  # uses sequence created by the rough oligomers
        seq = sequence[self.rough_oligomer_size:]

        # generates a list of possible oligomers for the given rough oligomers
        for i in range(len(rough_oligo_list) - 1):
            possible_overlaps, possible_overlaps_len, = [], []

            for overlap in range(self.low_overlap, self.high_overlap):
                overlap_seq = seq[:overlap]
                oligomer = rough_oligo_list[i] + overlap_seq
                possible_overlaps.append(oligomer)
                possible_overlaps_len.append(overlap)
            list_of_oligomers.append(possible_overlaps)
            overlap_length.append(possible_overlaps_len)
            seq = seq[self.rough_oligomer_size:]
        # returns two separate lists of lists that 1) contain oligomers with varying overlap lengths 2) their respective overlap lengths
        return list_of_oligomers, overlap_length, self.rough_oligomer_size, self.high_overlap

    #TODO: define insert and backbone regions

# selects oligomers based on user parameters: Melting temp, GC content, overlap length, GC clamp
class OligomerSelection:
    def __init__(self, optimal_temp, temp_range):
        self.optimal_temp = optimal_temp
        self.temp_range = temp_range
        self.low_tm = optimal_temp - temp_range  # lower melting temperature
        self.high_tm = optimal_temp + temp_range  # higher melting temperature
        self.sp = SequenceProcessing()

    # returns value in list closest to n
    # input a list and a number you want to check
    @staticmethod
    def closest(data, n):
        # returns int
        return data[min(range(len(data)), key=lambda i: abs(data[i] - n))]

# identifies oligomers that fit certain user parameters
class UserParameterSelection(OligomerSelection):

    # selects the oligomers from the possible oligomers without T at 3' end (increased annealing)
    # input the list of all possible oligomers, and corresponding overlap lengths
    def t_3_free(self, list_of_oligomers, list_of_overlap_lengths):

        oligomer_list = list_of_oligomers
        overlap_length = list_of_overlap_lengths

        for oligo_index in range(len(oligomer_list)):
            # removes oligomers if there is T at 3' (alternate sequences will be complementary in the output)
            if oligo_index % 2 == 0:
                t_3_oligo_index = [oligomer_list[oligo_index].index(oligo) for oligo in
                                   oligomer_list[oligo_index] if oligo[-1] == 'T']
                if len(t_3_oligo_index) != len(oligomer_list[oligo_index]):
                    for ele in sorted(t_3_oligo_index, reverse=True):
                        del oligomer_list[oligo_index][ele]
                        del overlap_length[oligo_index][ele]
                else:
                    for ele in sorted(t_3_oligo_index[1:], reverse=True):
                        del oligomer_list[oligo_index][ele]
                        del overlap_length[oligo_index][ele]
            else:
                while oligomer_list[oligo_index][0][0] == 'A':
                    oligomer_list[oligo_index][:] = [i[1:] for i in oligomer_list[oligo_index]]
                    overlap_length[oligo_index - 1][:] = [number - 1 for number in overlap_length[oligo_index - 1]]
        # returns two separate lists of lists that 1) contain oligomers with varying overlap lengths 2) their respective overlap lengths
        return oligomer_list, overlap_length

    # selects the oligomers that fit the required GC range and Tm.
    # input the list of all possible oligomers, and corresponding overlap lengths
    def gc_tm_optimal(self, list_of_oligomers, list_of_overlap_lengths):
        optimal_oligomers, optimal_overlap_len = [], []
        comp_list = []
        # generating alternate complementary sequences to select based on melting temperature and GC
        for oligomer_pool in list_of_oligomers:
            comp_list.append(self.sp.complement_oli(oligomer_pool))

        for oligo_num in range(len(list_of_oligomers)):
            tm, oligo, over_len = [], [], []
            for rough_oligo_index in range(len(list_of_oligomers[oligo_num])):
                overlap = comp_list[oligo_num][rough_oligo_index][
                          -list_of_overlap_lengths[oligo_num][rough_oligo_index]:]  # extract overlap sequence
                melting_temp = round(mt.Tm_NN(overlap, nn_table=mt.DNA_NN3),
                                     2)  # check melting temperature based on Nearest Neighbour equation
                gc = GC(overlap)
                tm.append(melting_temp)  # to keep track of all melting temperatures of overlaps
                if Oligomers.low_gc <= gc <= Oligomers.high_gc and self.low_tm <= melting_temp <= self.high_tm:
                    oligo.append(list_of_oligomers[oligo_num][rough_oligo_index])
                    over_len.append(list_of_overlap_lengths[oligo_num][rough_oligo_index])
            if not oligo:  # if the list 'oligo' is empty
                # chooses the overlap with the closest melting temperature to the user's input
                if all(temp >= self.high_tm for temp in tm) == True:
                    next_best_oligo = self.closest(tm, self.high_tm)
                else:
                    next_best_oligo = self.closest(tm, self.low_tm)
                oligo.append(list_of_oligomers[oligo_num][tm.index(next_best_oligo)])
                over_len.append(list_of_overlap_lengths[oligo_num][tm.index(next_best_oligo)])
            optimal_oligomers.append(oligo)
            optimal_overlap_len.append(over_len)
        # returns 2 list of lists with all the oligomers that fit/ almost fit the criteria specified by user
        return optimal_oligomers, optimal_overlap_len

    # selects the smallest oligomer from the possible optimal oligomers (cost effective)
    # input the list of lists of optimal oligomers and corresponding overlap lengths
    def smallest_oligomers(self, optimal_oligomers, optimal_overlap_len, rough_oligo_list):

        smallest_oligomers, smallest_overlap_len = [], []

        for oligomer_index in range(len(optimal_oligomers)):
            smallest_ovr = min(optimal_overlap_len[oligomer_index])
            smallest_oli = optimal_oligomers[oligomer_index][optimal_overlap_len[oligomer_index].index(smallest_ovr)]
            smallest_oligomers.append(smallest_oli)
            smallest_overlap_len.append(smallest_ovr)
        if len(rough_oligo_list[-1]) > optimal_overlap_len[0][-1]:
            smallest_oligomers.append(rough_oligo_list[-1])
            smallest_overlap_len.append(0)
        # returns 2 lists that contain the oligomer sequences and overlap lengths
        return smallest_oligomers, smallest_overlap_len

    # specific to circular assemblies. Select oligomers based on length (smallest)
    def _smallest_oligomers(self, optimal_oligomers, optimal_overlap_len, rough_oligo_list, rough_oligomer_size, high_overlap):

        smallest_oligomers, smallest_overlap_len = [], []

        for oligomer_index in range(len(optimal_oligomers)):
            smallest_ovr = min(optimal_overlap_len[oligomer_index])
            smallest_oli = optimal_oligomers[oligomer_index][
                optimal_overlap_len[oligomer_index].index(smallest_ovr)]
            smallest_oligomers.append(smallest_oli)
            smallest_overlap_len.append(smallest_ovr)
        if not len(rough_oligo_list) % 2 == 0 & len(
                smallest_oligomers[-1]) < rough_oligomer_size + high_overlap:
            r = (rough_oligomer_size + high_overlap) - len(smallest_oligomers[-1])
            smallest_oligomers[-1] = smallest_oligomers[-1] + rough_oligo_list[-1][
                                                              smallest_overlap_len[-1]:(
                                                                          smallest_overlap_len[-1] + r)]
            return smallest_oligomers, smallest_overlap_len

# generates clusters based on sequence similarities
class Cluster:

    def __init__(self, cluster_size, cluster_range):
        self.cluster_size = cluster_size
        self.cluster_range = cluster_range
        self.low_cluster = cluster_size - cluster_range
        self.high_cluster = cluster_size + cluster_range + 1
        self.sp = SequenceProcessing()

    # make clusters based on alignment score
    # input list of oligomers and their respective overlap lengths; the (int) user input of cluster size and cluster range
    def complementary_clusters(self, oligomers, overlap_len):

        num = []

        for cluster_len in range(self.low_cluster, self.high_cluster, 2):
            num.append(cluster_len)

        if not len(oligomers) in range(self.high_cluster):

            comp_list, five_to_three = self.sp.complement_oli(oligomers), self.sp.seq_orientation(oligomers)

            clusters, comp_clusters, cluster_ovr, cluster_five_two_three = [], [], [], []

            oligos = oligomers
            overlap = overlap_len

            while oligos:
                alignments = []

                for cluster_length in num:
                    alignment = self.sp.overlap_alignment(
                        self.sp.overlap_list(oligos[:cluster_length], overlap[:cluster_length]))
                    alignments.append(alignment)

                    if alignments.count(alignments[0]) == len(alignments):
                        oligomers = oligos[:self.cluster_size]
                        overlap_len = overlap[:self.cluster_size]

                        oligos = oligos[self.cluster_size:]
                        overlap = overlap[self.cluster_size:]

                        complementary_cluster = comp_list[:self.cluster_size]
                        final_cluster = five_to_three[:self.cluster_size]
                        comp_list = comp_list[self.cluster_size:]
                        five_to_three = five_to_three[self.cluster_size:]

                        clusters.append(oligomers)
                        comp_clusters.append(complementary_cluster)
                        cluster_ovr.append(overlap_len)
                        cluster_five_two_three.append(final_cluster)

                    elif any(alignments.count(element) > 1 for element in alignments):
                        min_align = min(alignments, key=lambda x: len(x))
                        if alignments.count(min_align) > 1:
                            c = [num[index] for index, element in enumerate(alignments) if element == min_align]
                            cluster_dim = num[min(range(len(c)),
                                                  key=lambda cluster_length: abs(c[cluster_length] - self.cluster_size))]
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
            comp_clusters = [self.sp.complement_oli(oligomers)]
            cluster_ovr = [overlap_len]
            cluster_five_two_three = [self.sp.seq_orientation(oligomers)]

        # returns three list of lists containing oligomers and their respective overlap lengths in appropriate clusters sizes
        return clusters, comp_clusters, cluster_five_two_three, cluster_ovr

# calculates an arbitrary score for each oligomer
class Score(UserParameterSelection):

    # generates scores and suggests possible areas of error for each oligomer
    def overlap_score(self, clusters, comp_clusters, cluster_5_3, overlap):
        score, fault, repeats, repeat_seq, repeat_position = [], [], [], [], []
        for i in range(len(clusters)):
            s, f, r, rs, rp = [], [], [], [], []
            for x in range(len(clusters[i])):
                s.append(0)
                f.append([])
                rp.append([])
            score.append(s)
            fault.append(f)
            repeats.append(r)
            repeat_seq.append(rs)
            repeat_position.append(rp)
        cluster_length = []
        if len(clusters) == 1:
            cluster_length.append(0)
        else:
            for i in range(len(clusters)):
                cluster_length.append(i)
            print(cluster_length)
        for cluster in cluster_length:
            g = []
            for x in range(len(clusters[cluster])):
                overla = overlap[cluster][x]
                data_to_put = clusters[cluster][x][-overla:]
                g.append(data_to_put)
                d = "".join(g)
                s = self.sp.repeat_seq(d)
            for i in s:
                for a in i:
                    seq_repeat = d[a:a + 11]
                    match_index = [clusters[cluster].index(h) for h in clusters[cluster] if seq_repeat in h]
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
        for cluster in range(len(comp_clusters)):
            for oligomer in range(len(comp_clusters[cluster])):
                overlap_size = overlap[cluster][oligomer]
                clust = comp_clusters[cluster][oligomer]
                overlap_mt = round(mt.Tm_NN(clust[-overlap_size:], nn_table=mt.DNA_NN3), 2)
                if overlap_size != 0:
                    if self.low_tm > overlap_mt:
                        score[cluster][oligomer] += round(abs(self.low_tm - overlap_mt), 2)
                        fault[cluster][oligomer] += "L"
                    elif overlap_mt > self.high_tm:
                        score[cluster][oligomer] += round(abs(overlap_mt - self.high_tm), 2)
                        fault[cluster][oligomer] += "H"
                        if overlap_size > int(overlap_size - (overlap_size / 4)):
                            if cluster_5_3[cluster][oligomer][-2] == "T":
                                score[cluster][oligomer] += 1
                                fault[cluster][oligomer] += "T"
                    gc = ["CCC", "GGG", "CCG", "CGC", "GCC", "GGC", "GCG", "CGG"]
                    if any((t in cluster_5_3[cluster][oligomer][-5:]) for t in gc):
                        score[cluster][oligomer] += 1
                        fault[cluster][oligomer] += "G"
        return score, fault, repeat_position

