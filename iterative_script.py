# import all the required packages and classes
from assembly import Sequence, OligomerGroups, Clusters
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Seq
import argparse
import time
import xlsxwriter
from itertools import chain
import decimal
import os
from progressbar import ProgressBar, Percentage, Bar, ETA

# initiate the parser
# can't import classes from script.py as it would run the script.py file as well
parsers = argparse.ArgumentParser()

# create a range between two floats
def drange(x, y, jump):
  while x < y:
    yield float(x)
    x += decimal.Decimal(jump)

# input file name and sequence orientation; fixed variables
def get_filename():
    parsers.add_argument('-f', '--file', type=str, required=True, help="Fasta format input file")
    parsers.add_argument('-so', '--seq_orientation', type=str, help="Is the sequence linear or circular?",
                         default="l")
    args = parsers.parse_args()
    return args.file, args.seq_orientation

# same Assembly class as in script.py
class Assembly():

    def __init__(self, file_name, oligomer_size, overlap_size, optimal_temp, seq_orientation):
        self.gene_seq = Sequence.read_fasta(file_name)
        self.oligomer_size = oligomer_size
        self.overlap_size = overlap_size
        self.optimal_temp = optimal_temp
        self.temp_range = 2.5 # fixed variable
        self.seq_orientation = seq_orientation
        self.og = OligomerGroups(self.gene_seq, self.seq_orientation, self.oligomer_size, self.overlap_size, self.optimal_temp, self.temp_range)
        self.c = Clusters()

    def oligomer_design(self):

        rough_list = self.og.oligomer_splice()
        possible_oligomers, possible_overlap_len = self.og.oligomer_array(rough_list)
        t_oligomers, t_overlaps = self.og.t_3_free(possible_oligomers, possible_overlap_len)
        g_oligomers, g_overlaps = self.og.gc_tm_optimal(t_oligomers, t_overlaps)
        final_oligomers, final_overlaps = self.og.smallest_oligomers(g_oligomers, g_overlaps, rough_list)
        clusters, comp_clusters, cluster_five_two_three, cluster_ovr = self.c.complementary_clusters(final_oligomers, final_overlaps)
        score, fault, repeats, r = self.og.overlap_score(clusters, comp_clusters, cluster_five_two_three, cluster_ovr)

        return comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats

    def output(self, comp_clusters, cluster_five_to_three, overlap_cluster, score, fault, repeats):
        workbook = xlsxwriter.Workbook(
            f'oligomers_{self.oligomer_size}_{self.overlap_size}_{self.optimal_temp}_{self.seq_orientation}.xlsx')
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
                    f'oligomers_{self.oligomer_size}_{self.overlap_size}_{self.optimal_temp}_{self.seq_orientation}.fasta',
                    "fasta")
        print(
            f'oligomers_{self.oligomer_size}_{self.overlap_size}_{self.optimal_temp}_{self.seq_orientation} as fasta and xlsx files')
        return cluster_five_to_three

# user input of the file name and seq_orientation
file_name, seq_orientation = get_filename()
# set oligomer, overlap, and tm ranges
oligomer_size = range(50, 100) # input final number as + 1 (python counting)
overlap_size = range(20, 30)
melting_temp = list(drange(56, 57, 0.1))

# generate a dictionary
total = {"Oligomer": [],
         "Overlap": [],
         "Melting_Temp": [],
         "Score": []}

# initiate the progress bar
N = len(oligomer_size) * len(overlap_size) * len(melting_temp)
pbar = ProgressBar(widgets=[Bar('>', '[', ']'), ' ', Percentage(), ' ', ETA()],maxval=N)

print("Identifying the best combination of oligomer size, overlap size, and overlap melting temperature")
# identifies the best combination of oligomer size, overlap size, and overlap tm by iterating over by the number of times as specified earlier
for oligo in pbar(oligomer_size):
    for overlap in overlap_size:
        for temp in melting_temp:
            a = Assembly(str(file_name), oligo, overlap, temp, str(seq_orientation))
            comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats = a.oligomer_design()
            # finds the total score of a suggested assembly (higher the score, more possible causes of errors)
            total_score = sum(list(chain(*score)))
            total["Oligomer"].append(oligo)
            total["Overlap"].append(overlap)
            total["Melting_Temp"].append((temp))
            total["Score"].append(total_score)

# identify the instance with the lowest score
final_score = min(total["Score"])
# pick out all instances with the lowest score
indices = [i for i, x in enumerate(total["Score"]) if x == final_score]

user = round(time.time())
os.mkdir(f'/app/output_script/{user}')
os.chdir(f'/app/output_script/{user}')

# generate the output files for each of the combinations identified earlier
for index in indices:
    final_overlap = total["Overlap"][index]
    final_oligomer = total["Oligomer"][index]
    final_melting_temp = total["Melting_Temp"][index]
    a = Assembly(str(file_name), final_oligomer, final_overlap, final_melting_temp, str(seq_orientation))
    comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats = a.oligomer_design()
    a.output(comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats)

os.chdir('/app')
