from assembly import OligomerGroups, Clusters
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Seq
import argparse
import time
import xlsxwriter
import re

parser = argparse.ArgumentParser()

# Input parameters and initialize variables
class Input:
    # initialize variables
    def __init__(self, file_name, oligomer_size, overlap_size, optimal_temp, temp_range, user, seq_orientation):
        self.gene_seq = re.sub(r"\s+", "", str(file_name)).upper()
        #self.gene_seq = FastaFile(file_name).read_fasta()
        self.oligomer_size = int(oligomer_size)
        self.overlap_size = int(overlap_size)
        self.optimal_temp = round(float(optimal_temp), 4)
        self.temp_range = float(temp_range)
        self.user = str(user)
        self.seq_orientation = str(seq_orientation)

class Assembly(Input):

    def __init__(self, file_name, oligomer_size, overlap_size, optimal_temp, temp_range, seq_orientation, user):
        super().__init__(file_name, oligomer_size, overlap_size, optimal_temp, temp_range, seq_orientation, user)
        self.og = OligomerGroups(self.gene_seq, self.seq_orientation, self.oligomer_size, self.overlap_size, self.optimal_temp, self.temp_range)
        self.c = Clusters()

    def oligomer_design(self):

        rough_list = self.og.oligomer_splice()
        possible_oligomers, possible_overlap_len = self.og.oligomer_array(rough_list)
        t_oligomers, t_overlaps = self.og.t_3_free(possible_oligomers, possible_overlap_len)
        g_oligomers, g_overlaps = self.og.gc_tm_optimal(t_oligomers, t_overlaps)
        final_oligomers, final_overlaps = self.og.smallest_oligomers(g_oligomers, g_overlaps, rough_list)

        clusters, comp_clusters, cluster_five_two_three, cluster_ovr = self.c.complementary_clusters(final_oligomers, final_overlaps)
        score, fault, oligomer_repeats, repeats = self.og.overlap_score(clusters, comp_clusters, cluster_five_two_three, cluster_ovr)

        return comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, oligomer_repeats

    def output(self, comp_clusters, cluster_five_to_three, overlap_cluster, score, fault, repeats):
        workbook = xlsxwriter.Workbook(
            f'/app/output/oligomers_{self.user}_{self.oligomer_size}_{self.overlap_size}_{int(self.optimal_temp)}_{self.seq_orientation}.xlsx')
        worksheet = workbook.add_worksheet()
        row = 1
        col = 0
        fasta_records = []
        for cluster in range(len(cluster_five_to_three)):
            for oligomer in range(len(cluster_five_to_three[cluster])):
                overlap_len = overlap_cluster[cluster][oligomer]
                oligo_len = len(cluster_five_to_three[cluster][oligomer])
                oligo = cluster_five_to_three[cluster][oligomer]
                clust = comp_clusters[cluster][oligomer]
                if overlap_len != 0:
                    overlap_mt = round(mt.Tm_NN(clust[-overlap_len:], nn_table=mt.DNA_NN3), 2)
                else:
                    overlap_mt = 0
                oligo_len = len(oligo)
                record = SeqRecord(Seq.Seq(oligo), f'Oligomer_{cluster + 1}.{oligomer + 1}',
                                   description=f'Overlap Length: {overlap_len}, Oligomer Length: {oligo_len}, Overlap Melting Temperature: {overlap_mt}')
                fasta_records.append(record)
                worksheet.write(0, 0, 'Cluster Number')
                worksheet.write(0, 1, 'Oligomer Number')
                worksheet.write(0, 2, "Sequence (5' to 3')")
                worksheet.write(0, 3, "Overlap length")
                worksheet.write(0, 4, "Oligomer length")
                worksheet.write(0, 5, 'Melting Temperature of overlap')
                worksheet.write(0, 6, 'Overlap Score')
                worksheet.write(0, 7, 'Cause of errors')
                worksheet.write(0, 8, 'Repeat Sequences within the oligomer')
                worksheet.write(row, col, f'Cluster {cluster + 1}')
                worksheet.write(row, col + 1, f'oligomer {oligomer + 1}')
                worksheet.write(row, col + 2, oligo)
                worksheet.write(row, col + 3, overlap_len)
                worksheet.write(row, col + 4, oligo_len)
                worksheet.write(row, col + 5, overlap_mt)
                worksheet.write(row, col + 6, score[cluster][oligomer])
                worksheet.write(row, col + 7, "".join([i for i in fault[cluster][oligomer]]))
                worksheet.write(row, col + 8, ", ".join([i for i in repeats[cluster][oligomer]]))
                row += 1
        workbook.close()

        SeqIO.write(fasta_records,
                    f'/app/output/oligomers_{self.user}_{self.oligomer_size}_{self.overlap_size}_{int(self.optimal_temp)}_{self.seq_orientation}.fasta',
                    "fasta")
        return cluster_five_to_three



