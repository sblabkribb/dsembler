from functions import DnaAssembly as da
from Bio import Seq
from Bio.SeqUtils import MeltingTemp as mt
import xlsxwriter
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import shutil
import zipfile
import tempfile

da = da()
class DnaAssemblyDesigner:
    def design_oligomers(self, gene_seq, oligomer_size, overlap_size, optimal_temp, temp_range, cluster_size, cluster_range, user):
        
        rough_oligomer = da.oligomer_splice(gene_seq, oligomer_size, overlap_size)
        list_of_oligomers, overlap_length = da.possible_oligomers(gene_seq, rough_oligomer, oligomer_size, overlap_size)
        optimal_oligomers, optimal_overlap_len = da.optimal_oligomers(list_of_oligomers, optimal_temp, temp_range, overlap_length)
        final_oligomers, final_overlaps = da.final_oligomers(optimal_oligomers, optimal_overlap_len, rough_oligomer)
        clust, clusters, overlap_cluster, cluster_five_to_three = da.complementary_clusters(final_oligomers, final_overlaps, cluster_size, cluster_range)

        score, fault, repeats_position =da.overlap_score(clust, clusters, cluster_five_to_three, overlap_cluster, optimal_temp, temp_range)
        
        workbook = xlsxwriter.Workbook(f'/app/output/oligomers_{user}.xlsx')
        worksheet = workbook.add_worksheet()
        row = 1
        col = 0
        workbook = xlsxwriter.Workbook(f'/app/output/oligomers_{user}.xlsx')
        worksheet = workbook.add_worksheet()
        row = 1
        col = 0
        fasta_records = []
        for cluster in range(len(cluster_five_to_three)):
            for oligomer in range(len(cluster_five_to_three[cluster])):
                overlap_len = overlap_cluster[cluster][oligomer]
                oligo_len = len(cluster_five_to_three[cluster][oligomer])
                oligo = cluster_five_to_three[cluster][oligomer]
                clust = clusters[cluster][oligomer]
                if overlap_len != 0:
                    overlap_mt = round(mt.Tm_NN(clust[-overlap_len:], nn_table=mt.DNA_NN3), 2)
                else:
                    overlap_mt = 0
                oligo_len = len(oligo)
                record = SeqRecord(Seq.Seq(oligo), f'Cluster_{cluster + 1}, Oligomer_{oligomer + 1}', description=f'Overlap Length: {overlap_len}, Oligomer Length: {oligo_len}, Overlap Melting Temperature: {overlap_mt}')
                fasta_records.append(record)
                worksheet.write(0, 0, 'Cluster Number')
                worksheet.write(0, 1, 'Oligomer Number')
                worksheet.write(0, 2, "Sequence (5' to 3')")
                worksheet.write(0, 3, "Overlap length")
                worksheet.write(0, 4, "Oligomer length")
                worksheet.write(0, 5, 'Melting Temperature of overlap')
                worksheet.write(0, 6, 'Overlap Score')
                worksheet.write(0, 7, 'Cause of errors')
                worksheet.write(0, 8, 'Repeat Sequences')
                worksheet.write(row, col, f'Cluster {cluster + 1}')
                worksheet.write(row, col + 1, f'oligomer {oligomer + 1}')
                worksheet.write(row, col + 2, oligo)
                worksheet.write(row, col + 3, overlap_len)
                worksheet.write(row, col + 4, oligo_len)
                worksheet.write(row, col + 5, overlap_mt)
                worksheet.write(row, col + 6, score[cluster][oligomer])
                worksheet.write(row, col + 7, "".join([i for i in fault[cluster][oligomer]]))
                worksheet.write(row, col + 8, ", ".join([i for i in repeats_position[cluster][oligomer]]))
                row += 1
        workbook.close()
        
        SeqIO.write(fasta_records, f'/app/output/oligomers_{user}.fasta', "fasta")
        return cluster_five_to_three