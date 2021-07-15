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
We use the following parent classes: Sequence, Oligomers, and Clusters. Each class manipulates their objects.
Sequence: Performs specific methods on any given DNA sequences present as string or in lists or nested list. 
Oligomers: Has one child class: OligomerGroups- Acts on groups of oligomers at a time. 
Clusters: Here a given set of oligomers can be divided into clusters based on the sequence similarities within each cluster.
'''
# file architecture
'''
assembly.py
-| Sequence
    -| @staticmethod read_fasta(file)
    -| @staticmethod complement_oli(data)
    -| @staticmethod overlap_list(oligomer, overlap_len, length)
    -| @staticmethod overlap_alignment(sequences)
    -| @staticmethod seq_orientation(sequences)
    -| @staticmethod repeat_seq(data)
-| Oligomers
    -| @staticmethod gc_clamp(data)
    -| @staticmethod t_3_end(data)
    -| @staticmethod ol_length(oligomers)
    -| rough_oligo_size_cir(self)
    -| overlap_tm(self, data)
    -| add_overlap(high_overlap, seq, oligo)
    -| overlap_score(cluster, complementary_cluster, cluster_5_3, overlap_len)
    -| OligomerGroups
        -| @staticmethod closest(data, n)
        -| oligomer_splice(self)
        -| olgiomer_array(self, rough_oligo_list)
        -| t_3_free(self, oligomer_list, overlap_length)
        -| gc_tm_optimal(self, list_of_oligomers, list_of_overlap_lengths)
        -| smallest_oligomers(self, optimal_oligomers, optimal_overlap_len, rough_oligo_list)
-| Clusters
    -| complementary_clusters(self, oligomers, overlap_length) 
    -| recommended_clusters(self, oligomers, overlap_length)
'''

# processing DNA sequences within lists
class Sequence:
    @staticmethod
    # Reading fasta file to give the DNA sequence
    def read_fasta(file):
        '''
        :param file: file path (str)
        :return seq: extracted gene sequence (str)
        '''
        seq = None

        dna_seq = SeqIO.parse(file, "fasta")
        # ensuring all data is capitalized for uniformity
        for seq_record in dna_seq:
            seq = str(seq_record.seq).upper()

        return seq

    # function for generating alternative complementary sequences
    @staticmethod
    def complement_oli(data):
        '''
        :param data: list of (str) sequences
        :return comp_list: list of (str) alternate sequences being complementary
        '''
        comp_list = []

        for i in data:
            if data.index(i) % 2 != 0:  # allow alternate sequences to be complementary
                i = str(Seq.Seq(i).complement())
                comp_list.append(i)
            else:
                comp_list.append(i)

        return comp_list

    # generating single list with all overlap regions
    @staticmethod
    def overlap_list(oligomer_list, overlap_length, length):
        '''
        :param oligomer_list: list of (str) oligomers
        :param overlap_length: list of (int) overlap lengths
        :param length: size of the oligomer list to be extracted (int)
        :return overlaps: list of extracted (str) overlap sequences
        '''

        overlaps = []

        for oligomer in oligomer_list[:length]:
            overlap_len = overlap_length[oligomer_list.index(oligomer)]
            if overlap_len != 0:
                overlap = oligomer[-overlap_len:]
                overlaps.append(overlap)
            else:
                pass

        return overlaps

    # identifies sequences that have high sequence similarities
    @classmethod
    def overlap_alignment(cls, sequences):
        '''
        :param sequences: input a list of sequences
        :return alignment_index: list of list containing the (int) indexes of high levels of alignments between 2 sequences
        '''

        alignment_index = []
        repeats = 10
        x = []
        for i in sequences:
            seq_no = len(i) - (repeats - 1)
            seq_list = []
            for j in range(seq_no):
                start_pos = j
                end_pos = j + repeats
                b = i[start_pos:end_pos]
                seq_list.append(b)
            x.append(seq_list)
        overlap_segments = cls.flatten(x)
        for a, b in itertools.combinations(overlap_segments, 2):
            align = aligner.align(a, b)
            align_ratio = align.score / max(len(a), len(b))
            if align_ratio > 0.7:
                index = [overlap_segments.index(a), overlap_segments.index(b)]
                alignment_index.append(index)

        return alignment_index

    # orientates sequences in 5 to 3 format
    @staticmethod
    def seq_orientation(sequences):
        '''
        :param sequences: list of (str) sequences
        :return five_to_three: list of (str) alternate sequences being reverse complementary
        '''

        five_to_three = []

        for sequence in sequences:
            if sequences.index(sequence) % 2 != 0:
                rc_sequence = str(Seq.Seq(sequence).reverse_complement())
                five_to_three.append(rc_sequence)
            else:
                five_to_three.append(sequence)
        # list of oligomers with the right orientation
        return five_to_three

    @staticmethod
    def flatten(list_of_list):
        return [item for sublist in list_of_list for item in sublist]

    # finding repeat sequences
    # https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
    @classmethod
    def repeat_seq(cls, data):
        '''
        :param data: list of sequences
        :return y: list of list containing the (int) indexes of the repeats
        :return s: list of list containing the (str) repeat sequences
        '''
        x, y = [], []
        repeats = 10
        for i in data:
            seq_no = len(i) - (repeats - 1)
            seq_list = []
            for j in range(seq_no):
                start_pos = j
                end_pos = j + repeats
                b = i[start_pos:end_pos]
                seq_list.append(b)
            x.append(seq_list)
        overlap_segments = cls.flatten(x)
        s = [item for item, count in collections.Counter(overlap_segments).items() if count > 1]
        for repeat in s:
            matches = re.finditer(repeat, "".join(data))
            matches_position = [match.start() for match in matches]
            y.append(matches_position)
        return y, s

# class that acts on individual oligomer sequences
class Oligomers:
    low_gc = 20  # lower gc content limit
    high_gc = 80  # higher gc content limit

    def __init__(self, gene_seq, assembly_type, oligomer_size, overlap_size, optimal_tm, temp_range):
        self.gene_seq = gene_seq
        self.assembly_type = assembly_type
        self.oligomer_size = oligomer_size
        self.rough_oligo_size = oligomer_size - overlap_size
        # uses user input as shortest possible overlap length
        self.low_overlap = overlap_size
        # calculates the longest possible overlap length
        if self.assembly_type == "l":
            self.high_overlap = int(overlap_size * 3 / 2)
        elif self.assembly_type == "c":
            self.high_overlap = int(overlap_size * 4 / 3)
        self.optimal_temp = optimal_tm
        self.temp_range = temp_range
        self.low_tm = optimal_tm - temp_range  # lower melting temperature
        self.high_tm = optimal_tm + temp_range  # higher melting temperature

    # find the rough oligomer size if the assembly required is Circular
    def rough_oligo_size_cir(self):
        '''
        :return: an adjusted (int) rough oligomer length for circular assemblies
        '''
        # find the even number when divided with the length of the given gene seq is closest to the user's inputted (rough oligomers)
        number = None
        for i in range(1, 100):
            number = len(self.gene_seq) / (i * 2)
            if number < self.rough_oligo_size:
                break
        oligomer_s = round(number)  # estimated approximate overlap size as calculated previously

        oligo_size = oligomer_s - 1  # calculates the size of oligomers without overlaps
        return oligo_size

    # calculate overlap melting temperature of any oligomer (overlap) sequence
    # check if it falls within the target specifications
    def overlap_tm(self, data):
        '''
        :param data: sequence (str)
        :return: fault "H" or "L" or True if there is no fault
        :return: difference in temperature from the required (int)
        '''
        melting_temp = round(mt.Tm_NN(data, nn_table=mt.DNA_NN3), 2)
        if self.low_tm <= melting_temp <= self.high_tm:
            return True, int(0)
        elif melting_temp > self.high_tm:
            difference = melting_temp - self.high_tm
            return "H", difference
        elif melting_temp < self.low_tm:
            difference = self.low_tm - melting_temp
            return "L", difference

    # adds the next few bases to an oligomer in an iterative manner
    def add_overlap(self, seq, oligo):
        '''
        :param seq: whole sequence for which oligomers are designed (str)
        :param oligo: sequence of rough oligomer (str)
        :return possible_oligomers: list of (str) an array of oligomers with varying overlaps
        :return possible_overlaps: list of (int) varying overlaps lengths; corresponds to possible_oligomers
        '''
        possible_overlaps, possible_overlaps_len, = [], []
        for overlap in range(self.low_overlap, self.high_overlap):
            overlap_seq = seq[:overlap]
            oligomer = oligo + overlap_seq
            possible_overlaps.append(oligomer)
            possible_overlaps_len.append(overlap)
        return possible_overlaps, possible_overlaps_len


    # check for the presence of a GC clamp at the 3' end
    @staticmethod
    def gc_clamp(data):
        '''
        :param data: sequence (str)
        :return: boolean for the presence of GC clamp
        '''
        gc = ["CCC", "GGG", "CCG", "CGC", "GCC", "GGC", "GCG", "CGG"]
        if any((t in data[-5:]) for t in gc):
            return True
        else:
            return False

    # checks for the presence of a GC clamp at the 3' end
    @staticmethod
    def t_3_end(data):
        '''
        :param data: sequence (str)
        :return: boolean for the presence of T at the 3' end
        '''
        if data[-1] == "T":
            return True
        else:
            return False

    # calculates the oligomer lengths
    @staticmethod
    def ol_lengths(oligomer):
        oliogmer_length = len(oligomer)
        return oliogmer_length

    # generates scores and suggests possible areas of error for each oligomer
    def overlap_score(self, clusters, comp_clusters, cluster_5_3, overlap):
        '''
        :param clusters: list of lists of (str) oligomers which are directly spliced from the original sequence
        :param comp_clusters: list of lists of (str) with alternate complementary oligomers
        :param cluster_5_3: list of lists of (str) oligomers with reverse complementary oligomers
        :param overlap: list of lists of (int) overlap lengths corresponding to the above lists
        :return score: list of lists of (int) each overlap's score
        :return fault: list of lists of (str) each overlap's faults
        :return oligomer_repeats: list of lists of (str) the 10bp repeats that occur within each individual oligomers
        :return repeats: list of lists of (str) the 10bp repeats that occur between oligomer overlaps within a cluster
        '''
        # generate empty lists with the same shape as the clusters generated
        score, fault, repeats, repeat_seq, oligomer_repeats, repeat_position = [], [], [], [], [], []
        for i in range(len(clusters)):
            s, f, r, rs, rp = [], [], [], [], []
            for x in range(len(clusters[i])):
                # stores overlap score
                s.append(0)
                # stores overlap faults
                f.append([])
                # stores sequence repeats
                rp.append([])
            score.append(s)
            fault.append(f)
            repeats.append(r)
            repeat_seq.append(rs)
            oligomer_repeats.append(rp)
            repeat_position.append(rp)

        cluster_length = []
        if len(clusters) == 1:
            cluster_length.append(0)
        else:
            for i in range(len(clusters)):
                cluster_length.append(i)

        for cluster in cluster_length:
            # recognizes repeats within each oligomer
            for x in range(len(clusters[cluster])):
                seq_repeat_index, seq_repeat = Sequence.repeat_seq([clusters[cluster][x]])
                for count in seq_repeat:
                    score[cluster][x] += 1
                    fault[cluster][x] += "r"
                oligomer_repeats[cluster][x] =  [seq.lower() for seq in seq_repeat]
            # recognizes repeats between overlaps in a cluster
            data = Sequence.overlap_list(clusters[cluster], overlap[cluster], len(clusters[cluster]) -1)
            s, o = Sequence.repeat_seq(data)
            for i in o:
                for a in i:
                    seq_repeat = a
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
                    for y in repeat_seq[cluster]:
                        if x == y.upper():
                            repeat_position[cluster][x].append(repeat_seq[cluster][repeats[cluster].index(y)])

        for cluster in range(len(comp_clusters)):
            for oligomer in range(len(comp_clusters[cluster])):
                overlap_size = overlap[cluster][oligomer]
                clust = comp_clusters[cluster][oligomer]
                # checks if the overlap fits the specified criteria
                if overlap_size != 0:
                    comp_overlap = clust[-overlap_size:]
                    overlap_tm, difference = self.overlap_tm(comp_overlap)
                    # temperature
                    if overlap_tm != True:
                        score[cluster][oligomer] += difference
                        fault[cluster][oligomer] += overlap_tm
                    if overlap_tm == "H":
                        if overlap_size > int(overlap_size - (overlap_size / 4)):
                            # Thymine at 3' end
                            if cluster_5_3[cluster][oligomer][-2] == "T":
                                score[cluster][oligomer] += 1
                                fault[cluster][oligomer] += "T"
                    # GC clamp
                    if self.gc_clamp(cluster_5_3[cluster][oligomer]) == True:
                        score[cluster][oligomer] += 1
                        fault[cluster][oligomer] += "G"

        return score, fault, oligomer_repeats, repeat_position

# class that deals with multiple oligomers at a time
class OligomerGroups(Oligomers):
    def __init__(self, gene_seq, oligomer_size, overlap_size, optimal_tm, temp_range, assembly_type):
        super().__init__(gene_seq, oligomer_size, overlap_size, optimal_tm, temp_range, assembly_type)
        if self.assembly_type == "c":
            self.rough_oligo_size = self.rough_oligo_size_cir()

    # returns value in list closest to n
    # input a list and a number you want to check
    @staticmethod
    def closest(data, n):
        # returns int
        return data[min(range(len(data)), key=lambda i: abs(data[i] - n))]

    # splicing into oligomers for linear assembly
    # input gene sequence (str), target oligomer size(int), and target overlap size(int)
    def oligomer_splice(self):
        '''
        :return: list of (str) oligomer sequences without overlaps
        '''

        # calculates the size of oligomers without overlaps
        rough_oligo_list = []
        for start_pos in range(0, len(self.gene_seq), self.rough_oligo_size):
            rough_oligo = self.gene_seq[start_pos: start_pos + self.rough_oligo_size]
            if rough_oligo:
                rough_oligo_list.append(rough_oligo)

        if self.assembly_type == "c":
            # add an overlap region between the last and first oligomer (to make it circular)
            if len(rough_oligo_list[-1]) < self.rough_oligo_size:
                rough_oligo_list[-1] = rough_oligo_list[-1] + rough_oligo_list[0]
            else:
                rough_oligo_list.append(rough_oligo_list[0])

        return rough_oligo_list

    # generates an array of oligomers with varied lengths of overlaps (linear assembly)
    def oligomer_array(self, rough_oligo_list):
        '''
        :param rough_oligo_list: list of (str) oligomer sequences without overlaps
        :return list_of_oligomers:  list of lists containing (str) oligomers with varying overlaps
        :return overlap_length: list of lists containing (int) the corresponding overlap lengths
        '''
        list_of_oligomers, overlap_length = [], []

        sequence = ''.join(rough_oligo_list)  # uses sequence created by the rough oligomers
        seq = sequence[self.rough_oligo_size:]

        for i in range(len(rough_oligo_list) - 1):
            possible_overlaps, possible_overlaps_len = self.add_overlap(seq, rough_oligo_list[i])
            list_of_oligomers.append(possible_overlaps)
            overlap_length.append(possible_overlaps_len)
            seq = seq[self.rough_oligo_size:]

        return list_of_oligomers, overlap_length

    # selects oligomers if they don't have T at 3' end
    def t_3_free(self, oligomer_list, overlap_length):
        '''
        :param oligomer_list: list of lists containing (str) oligomers with varying overlaps
        :param overlap_length: list of lists containing (int) the corresponding overlap lengths
        :return oligomer_list:  list of lists containing (str) oligomers without T at 3'
        :return overlap_length: list of lists containing (int) the corresponding overlap lengths
        '''
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
    def gc_tm_optimal(self, list_of_oligomers, list_of_overlap_lengths):
        '''
        :param list_of_oligomers: list of lists containing (str) oligomers with varying overlaps
        :param list_of_overlap_lengths: list of lists containing (int) the corresponding overlap lengths
        :return optimal_oligomers:  list of lists containing (str) oligomers within the required GC and Tm
        :return optimal_overlap_len: list of lists containing (int) the corresponding overlap lengths
        '''
        optimal_oligomers, optimal_overlap_len = [], []
        comp_list = []
        # generating alternate complementary sequences to select based on melting temperature and GC
        for oligomer_pool in list_of_oligomers:
            comp_list.append(Sequence.complement_oli(oligomer_pool))

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
    def smallest_oligomers(self, optimal_oligomers, optimal_overlap_len, rough_oligo_list):
        '''
        :param optimal_oligomers: list of lists containing (str) oligomers with varying overlaps
        :param optimal_overlap_len: list of lists containing (int) the corresponding overlap lengths
        :param rough_oligo_list: list of (str) rough oligomers (without overlaps)
        :return smallest_oligomers:  list of (str) oligomers that were smallest from the inputted list
        :return smallest_overlap_len: list (int) corresponding overlap lengths
        '''
        smallest_oligomers, smallest_overlap_len = [], []

        for oligomer_index in range(len(optimal_oligomers)):
            smallest_ovr = min(optimal_overlap_len[oligomer_index])
            smallest_oli = optimal_oligomers[oligomer_index][optimal_overlap_len[oligomer_index].index(smallest_ovr)]
            smallest_oligomers.append(smallest_oli)
            smallest_overlap_len.append(smallest_ovr)
        if self.assembly_type == "l":
            if len(rough_oligo_list[-1]) > optimal_overlap_len[0][-1]:
                smallest_oligomers.append(rough_oligo_list[-1])
                smallest_overlap_len.append(0)
            smallest_overlap_len[-1] = 0 # linear assembly final oligomers don't have overlaps
        else:
            if not len(rough_oligo_list) % 2 == 0 & len(
                    smallest_oligomers[-1]) < self.rough_oligo_size + self.high_overlap:
                r = (self.rough_oligo_size + self.high_overlap) - len(smallest_oligomers[-1])
                smallest_oligomers[-1] = smallest_oligomers[-1] + rough_oligo_list[-1][
                                                                  smallest_overlap_len[-1]:(
                                                                     smallest_overlap_len[-1] + r)]

        return smallest_oligomers, smallest_overlap_len

# class that works with oligomer clusters
class Clusters:
    # make clusters based on alignment score
    # input list of oligomers and their respective overlap lengths; the (int) user input of cluster size and cluster range
    def complementary_clusters(self, oligomers, overlap_len):
        '''
        :param oligomers: list of (str) oligomers
        :param overlap_len: list of (int) corresponding overlap lengths
        :return clusters: list of lists of (str) oligomers which are directly spliced from the original sequence
        :return comp_clusters: list of lists of (str) with alternate complementary oligomers
        :return cluster_5_3: list of lists of (str) oligomers with reverse complementary oligomers
        :return overlap: list of lists of (int) overlap lengths corresponding to the above lists
        '''

        comp_list, five_to_three = Sequence.complement_oli(oligomers), Sequence.seq_orientation(oligomers)

        clusters, comp_clusters, cluster_ovr, cluster_five_two_three = [], [], [], []

        oligos = oligomers
        overlap = overlap_len
        while len(oligos) > 1:
            for cluster_length in range(1, 30):
                alignment, seq_repeat = Sequence.repeat_seq(
                    Sequence.overlap_list(oligos[:cluster_length], overlap[:cluster_length], len(oligos[:cluster_length])))
                if len(alignment) > 0:
                    break

            if cluster_length > 1:
                index = cluster_length -1
            else:
                index = cluster_length

            oligomer = oligos[:index]

            overlap_length = overlap[:index]
            complementary_cluster = comp_list[:index]
            final_cluster = five_to_three[:index]

            oligos = oligos[index:]
            overlap = overlap[index:]
            comp_list = comp_list[index:]
            five_to_three = five_to_three[index:]

            clusters.append(oligomer)
            cluster_ovr.append(overlap_length)
            comp_clusters.append(complementary_cluster)
            cluster_five_two_three.append(final_cluster)

        if len(oligos) == 1:
            clusters.append(oligos)
            cluster_ovr.append(overlap)
            comp_clusters.append(comp_list)
            cluster_five_two_three.append(five_to_three)

        return clusters, comp_clusters, cluster_five_two_three, cluster_ovr
