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

    # Reading fasta file to give the DNA sequence and adding 'N' if the DNA sequence is not exactly divisible by 3
    def read_fasta(self, data):
        dna_seq = SeqIO.parse(data, "fasta")
        for seq_record in dna_seq:
            seq = str(seq_record.seq).upper()
            len_data = len(seq)
            if len_data / 3 != 0:
                seq += 'N'
        return seq

    # splicing into oligomers
    def oligomer_splice(self, data, oligomer_size, overlap_size):
        oligomer_step = oligomer_size - overlap_size
        dna_list = []
        for start_pos in range(0, len(data), oligomer_step):
            data_to_put = data[start_pos:start_pos+oligomer_step]
            if data_to_put:
                dna_list.append(data_to_put)
        return dna_list

    # generate all possible oligomers within a given overlap range
    def possible_oligomers(self, gene_seq, data, oligomer_size, overlap_size):
        list_of_oligomers,overlap_length = [], []
        overlap_region = oligomer_size - overlap_size
        low_overlap = int(overlap_size - (overlap_size/4))
        high_overlap = overlap_region + 1
        seq = gene_seq[overlap_region:]
        for i in range(len(data) - 1):
            possible_overlaps, possible_overlaps_len,= [], []
            for overlap in range(low_overlap, high_overlap):
                overlap_seq = seq[:overlap]
                overlap_len = overlap
                oligomer = data[i] + overlap_seq
                possible_overlaps.append(oligomer)
                possible_overlaps_len.append(overlap_len)
            list_of_oligomers.append(possible_overlaps)
            overlap_length.append(possible_overlaps_len)
            seq = seq[overlap_region:]
        return list_of_oligomers, overlap_length

    # function for generating alternative complementary sequences
    def complement_oli(self, data):
        comp_list = []
        for i in data:
            if data.index(i) % 2 != 0:
                i = str(Seq.Seq(i).complement())
                comp_list.append(i)
            else:
                comp_list.append(i)
        return comp_list

    # returns value in list closest to n
    def closest(self, data, n):
        return data[min(range(len(data)), key=lambda i: abs(data[i] - n))]

    # selecting oligomers that fall in target requirements (Melting temp, GC content, overlap length, GC clamp)
    def optimal_oligomers(self, list_of_oligomers, optimal_temp, temp_range, overlap_length):
        low_tm = optimal_temp - temp_range
        high_tm = optimal_temp + temp_range
        low_gc = 20
        high_gc = 80
        optimal_oligomers, optimal_overlap_len = [], []
        for c in range(len(list_of_oligomers)):
            if c % 2 == 0:
                f = [list_of_oligomers[c].index(i) for i in list_of_oligomers[c] if i[-1] == 'T']
                if len(f) != len(list_of_oligomers[c]):
                    for ele in sorted(f, reverse=True):
                        del list_of_oligomers[c][ele]
                        del overlap_length[c][ele]
                else:
                    for ele in sorted(f[1:], reverse=True):
                        del list_of_oligomers[c][ele]
                        del overlap_length[c][ele]
            else:
                while list_of_oligomers[c][0][0] == 'A':
                    list_of_oligomers[c][:] = [i[1:] for i in list_of_oligomers[c]]
                    overlap_length[c - 1][:] = [number - 1 for number in overlap_length[c - 1]]
        comp_list = []
        for a in list_of_oligomers:
            comp_list.append(self.complement_oli(a))
        for i in range(len(list_of_oligomers)):
            tm, oli, over_len = [], [], []
            for x in range(len(list_of_oligomers[i])):
                overlap = comp_list[i][x][-overlap_length[i][x]:]
                melting_temp = round(mt.Tm_NN(overlap, nn_table=mt.DNA_NN3), 2)
                gc = GC(overlap)
                tm.append(melting_temp)
                if low_gc <= gc <= high_gc and low_tm <= melting_temp <= high_tm:
                    oli.append(list_of_oligomers[i][x])
                    over_len.append(overlap_length[i][x])
            if not oli:
                if all(melt >= high_tm for melt in tm) == True:
                    v = self.closest(tm, high_tm)
                else:
                    v = self.closest(tm, low_tm)
                oli.append(list_of_oligomers[i][tm.index(v)])
                over_len.append(overlap_length[i][tm.index(v)])
            optimal_oligomers.append(oli)
            optimal_overlap_len.append(over_len)
        return optimal_oligomers, optimal_overlap_len

    # selects the smallest oligomer from the possible optimal oligomers (cost effective)
    def final_oligomers(self, optimal_oligomers, optimal_overlap_len):
        final_oligomers = []
        final_overlaps = []
        for i in range(len(optimal_oligomers)):
            smallest_ovr = min(optimal_overlap_len[i])
            smallest_oli = optimal_oligomers[i][optimal_overlap_len[i].index(smallest_ovr)]
            final_oligomers.append(smallest_oli)
            final_overlaps.append(smallest_ovr)
        return final_oligomers, final_overlaps

    # generating single sequence with all overlap regions
    def overlap_list(self, data, overlap_len):
        overlaps = []
        for i in data:
            overlap = overlap_len[data.index(i)]
            data_to_put = i[-overlap:]
            overlaps.append(data_to_put)
        return overlaps

    # check overlap alignment of a given list
    def overlap_alignment(self, data):
        k = []
        for a, b in itertools.combinations(data, 2):
            align = aligner.align(a, b)
            align_ratio = align.score / max(len(a), len(b))
            if align_ratio > 0.7:
                index = [data.index(a), data.index(b)]
                k.append(index)
        return k

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
    def seq_orientation(self, data):
        five_to_three = []
        for i in data:
            if data.index(i) % 2 != 0:
                j = str(Seq.Seq(i).reverse_complement())
                five_to_three.append(j)
            else:
                five_to_three.append(i)
        return five_to_three

    def complementary_clusters(self, data, overlap, cluster_size, cluster_range):
        num = []
        l_cluster = cluster_size - cluster_range
        h_cluster = cluster_size + cluster_range + 1
        for i in range(l_cluster, h_cluster,2):
            num.append(i)
        if not len(data) in range(cluster_size) :
            comp_list, five_to_three = self.complement_oli(data), self.seq_orientation(data)
            clusters = []
            cluster_ovr = []
            cluster_five_two_three = []
            h = data
            o = overlap
            while h:
                f = []
                for i in num:
                    d = self.overlap_alignment(self.overlap_list(h[:i], o[:i]))
                    f.append(d)
                    if f.count(f[0]) == len(f):
                        h = h[cluster_size:]
                        r = o[:cluster_size]
                        o = o[cluster_size:]
                        k = comp_list[:cluster_size]
                        b = five_to_three[:cluster_size]
                        comp_list = comp_list[cluster_size:]
                        five_to_three = five_to_three[cluster_size:]
                        clusters.append(k)
                        cluster_ovr.append(r)
                        cluster_five_two_three.append(b)
                    elif any(f.count(element) > 1 for element in f):
                        t = min(f, key=lambda x: len(x))
                        if f.count(t) > 1:
                            c = [num[index] for index, element in enumerate(f) if element == t]
                            p = num[min(range(len(c)), key=lambda i: abs(c[i] - cluster_size))]
                        else:
                            p = num[f.index(t)]
                        h = h[p:]
                        r = o[:p]
                        o = o[p:]
                        k = comp_list[:p]
                        b = five_to_three[:p]
                        comp_list = comp_list[p:]
                        five_to_three = five_to_three[p:]
                        cluster_ovr.append(r)
                        clusters.append(k)
                        cluster_five_two_three.append(b)
        else:
            clusters = [data]
            cluster_ovr = [overlap]
            cluster_five_two_three = [self.seq_orientation(data)]
        return clusters, cluster_ovr, cluster_five_two_three

    # Visualizing DNA sequences
    def output_webpage(self, data, overlap):
        spaces = []
        overshot = []
        z = 0
        for j in range(len(overlap)):
            if j == 0:
                t = len(data[j]) - overlap[j]
            else:
                t = len(data[j]) - overlap[(j - 1)] - overlap[j]
            if t < 0:
                x = 0 - t
                t = 0
                z += 1
            else:
                x = 0
                z += 0
            spaces.append(t)
            overshot.append(x)
        even, odd_space, odd_overshoot = [], [], []
        odd, even_space, even_overshoot = [], [], []
        for j in range(len(data)):
            if j % 2 != 0:
                odd.append(data[j])
                even_space.append(spaces[j])
                even_overshoot.append(overshot[j])
            else:
                even.append(data[j])
                odd_space.append(spaces[j])
                odd_overshoot.append(overshot[j])
        forward_seq = even[0]
        reverse_seq = ""
        for w in range(len(even)):
            if w != range(len(even))[-1]:
                if even_overshoot[w] == 0:
                    forward_seq = forward_seq + (" " * even_space[w]) + even[w + 1]
                else:
                    b = even[w][-even_overshoot[w]:]
                    forward_seq = forward_seq[:-even_overshoot[w]] + b.swapcase() + even[w + 1][even_overshoot[w]:]
            else:
                forward_seq = forward_seq + even[w]
        for v in range(len(odd)):
            if odd_overshoot[v] == 0:
                reverse_seq = (reverse_seq + " " * odd_space[v] + odd[v])
            else:
                c = odd[v - 1][-odd_overshoot[v]:]
                reverse_seq = (reverse_seq[:-odd_overshoot[v]] + c.swapcase() + odd[v][odd_overshoot[v]:])
        return forward_seq, reverse_seq, z

    # finding repeat sequences
    # https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
    def repeat_seq(self, data):
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
                        score[cluster][l] += 10
                        fault[cluster][l] += "R"
                    repeats[cluster].append(match_index)
                    repeat_seq[cluster].append(seq_repeat)
            for x in range(len(repeat_position[cluster])):
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
