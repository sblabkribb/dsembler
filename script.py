from assembly import Sequence, OligomerGroups, Clusters
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Seq
import argparse
import time
import xlsxwriter

parser = argparse.ArgumentParser()

# Input parameters and initialize variables
class Input:
    # initialize variables
    def __init__(self, file_name, oligomer_size, overlap_size, optimal_temp, temp_range, cluster_size, cluster_range, seq_orientation):
        #self.gene_seq = str(file_name).upper()
        self.gene_seq = Sequence.read_fasta(file_name)
        self.oligomer_size = int(oligomer_size)
        self.overlap_size = int(overlap_size)
        self.optimal_temp = round(float(optimal_temp), 4)
        self.cluster_size = int(cluster_size)
        self.temp_range = float(temp_range)
        self.cluster_range = int(cluster_range)
        self.user = round(time.time())
        self.seq_orientation = str(seq_orientation)

    # Obtain values (use python script.py -h to open help)
    @classmethod
    def get_inputs(cls):
        # use the argsparse package to perform command line tasks conveniently
        parser.add_argument('-f', '--file', type=str, required=True, help="Fasta format input file")
        parser.add_argument('-ol', '--oligomer_length', type=int, required=True, help="Maximum length of oligomers")
        parser.add_argument('-ov', '--overlap_length', type=int, required=True, help="Target Overlap length")
        parser.add_argument('-tm', '--optimal_tm', type=float, help="Target Melting Temperature for each overlap", default= 56)
        parser.add_argument('-c', '--cluster_size', type=int, help="Number of Oligomers in one cluster", default=20)
        parser.add_argument('-tmr', '--temp_range', type=float, help="Range in which the overlap Tm is acceptable (+-C)", default=2.5)
        parser.add_argument('-cr', '--cluster_range', type=int, help="Range of Number of Oligomers in one cluster (+- oligomers)", default=4)
        parser.add_argument('-so', '--seq_orientation', type=str, help="Is the sequence linear or circular?", default="l")
        args= parser.parse_args()
        # return inputs to the __init__() function to initialize
        return cls(args.file, args.oligomer_length, args.overlap_length, args.optimal_tm, args.temp_range, args.cluster_size, args.cluster_range, args.seq_orientation)

class Assembly(Input):

    def __init__(self, file_name, oligomer_size, overlap_size, optimal_temp, temp_range, cluster_size, cluster_range, seq_orientation):
        super().__init__(file_name, oligomer_size, overlap_size, optimal_temp, temp_range, cluster_size, cluster_range, seq_orientation)
        self.og = OligomerGroups(self.gene_seq, self.seq_orientation, self.oligomer_size, self.overlap_size, self.optimal_temp, self.temp_range)
        self.c = Clusters(self.cluster_size, self.cluster_range)

    def oligomer_design(self):

        rough_list = self.og.oligomer_splice()
        possible_oligomers, possible_overlap_len = self.og.oligomer_array(rough_list)
        t_oligomers, t_overlaps = self.og.t_3_free(possible_oligomers, possible_overlap_len)
        g_oligomers, g_overlaps = self.og.gc_tm_optimal(t_oligomers, t_overlaps)
        final_oligomers, final_overlaps = self.og.smallest_oligomers(g_oligomers, g_overlaps, rough_list)

        clusters, comp_clusters, cluster_five_two_three, cluster_ovr = self.c.complementary_clusters(final_oligomers, final_overlaps)
        score, fault, repeats = self.og.overlap_score(clusters, comp_clusters, cluster_five_two_three, cluster_ovr)
        try:        
            self.c.recommended_clusters(final_oligomers, final_overlaps)
        except:
            pass

        return comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats

    def output(self, comp_clusters, cluster_five_to_three, overlap_cluster, score, fault, repeats):
        workbook = xlsxwriter.Workbook(
            f'/app/output_script/oligomers_{self.user}_{self.oligomer_size}_{self.overlap_size}_{int(self.optimal_temp)}_{self.cluster_size}_{self.seq_orientation}.xlsx')
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
                worksheet.write(0, 8, 'Repeat Sequences')
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
                    f'/app/output_script/oligomers_{self.user}_{self.oligomer_size}_{self.overlap_size}_{int(self.optimal_temp)}_{self.cluster_size}_{self.seq_orientation}.fasta',
                    "fasta")
        print(
            f'oligomers_{self.user}_{self.oligomer_size}_{self.overlap_size}_{int(self.optimal_temp)}_{self.cluster_size}_{self.seq_orientation} as fasta and xlsx files')
        return cluster_five_to_three

a = Assembly.get_inputs()

comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats = a.oligomer_design()

a.output(comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats)


