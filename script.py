from assembly import LinearAssembly, CircularAssembly, UserParameterSelection, Score, Cluster, FastaFile
from Bio import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse
import time
import xlsxwriter

parser = argparse.ArgumentParser()

# about
'''
We use one parent class of Input which has one subclass OligomerGenerator. Depending upon the use case, LinearOligomers or CircularOligomers is selected, and the appropriate Ouptuts are given
'''
# file architecture
'''
script.py
-| Input
    -| OligomerGenerator
        -| LinearOligomers ----------------------|
            -| design_oligomers(self)           -| Output
        -| CircularOligomers---------------------|  -| output_files(self)
            -| design_oligomers(self)
        -| oligomer_design(self)        
'''

# Input parameters and initialize variables
class Input:
    # initialize variables
    def __init__(self, file_name, oligomer_size, overlap_size, optimal_temp, temp_range, cluster_size, cluster_range):
        # self.gene_seq = str(file_name).upper()
        self.gene_seq = FastaFile(file_name).read_fasta()
        self.oligomer_size = int(oligomer_size)
        self.overlap_size = int(overlap_size)
        self.optimal_temp = round(float(optimal_temp), 4)
        self.cluster_size = int(cluster_size)
        self.temp_range = float(temp_range)
        self.cluster_range = int(cluster_range)
        self.user = round(time.time())

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
        args= parser.parse_args()

        # return inputs to the __init__() function to initialize
        return cls(args.file, args.oligomer_length, args.overlap_length, args.optimal_tm, args.temp_range, args.cluster_size, args.cluster_range)

# Run methods from assembly.py to perform a Linear or Circular assembly
class OligomerGenerator(Input):
    def __init__(self, file_name, oligomer_size, overlap_size, optimal_temp, temp_range, cluster_size, cluster_range):
        super().__init__(file_name, oligomer_size, overlap_size, optimal_temp, temp_range, cluster_size, cluster_range)
        self.ups = UserParameterSelection(self.optimal_temp, self.temp_range)

    # common functions for Linear and Circular assembly are executed using this method
    def oligomer_design(self, final_oligomers, final_overlaps):
        # generate clusters
        c = Cluster(self.cluster_size, self.cluster_range)
        self.clust, self.clusters, self.cluster_five_to_three, self.overlap_cluster = c.complementary_clusters(final_oligomers, final_overlaps)

        # calculate score for each overlap
        s = Score(self.optimal_temp, self.temp_range)
        self.score, self.fault, self.repeats_position = s.overlap_score(self.clust, self.clusters, self.cluster_five_to_three, self.overlap_cluster)

        # return variables required by the user
        return self.clusters, self.cluster_five_to_three, self.overlap_cluster, self.score, self.fault, self.repeats_position

# Initiate oligomer splicing with the method developed for LinearAssembly
class LinearOligomers(OligomerGenerator):
    # splice given sequence and generate possible oligomers
    def design_oligomers(self):

        la = LinearAssembly(self.gene_seq, self.oligomer_size, self.overlap_size)

        rough_oligomer = la.oligomer_splice()
        list_of_oligomers, overlap_length = la.possible_oligomers(rough_oligomer)

        # select optimal oligomers
        t_3_free_oligomers, t_3_free_overlap_len = self.ups.t_3_free(list_of_oligomers, overlap_length)
        optimal_oligomers, optimal_overlap_len = self.ups.gc_tm_optimal(t_3_free_oligomers, t_3_free_overlap_len)
        final_oligomers, final_overlaps = self.ups.smallest_oligomers(optimal_oligomers, optimal_overlap_len, rough_oligomer)

        # returns final oligomer list with their corresponding overlap lengths
        return final_oligomers, final_overlaps

# Initiate oligomer splicing with the method developed for CircularAssembly
class CircularOligomers(OligomerGenerator): #TODO
    # splice given sequence and generate possible oligomers
    def _design_oligomers(self):

        ca = CircularAssembly(self.gene_seq, self.oligomer_size, self.overlap_size)

        rough_oligomer = ca.oligomer_splice()
        list_of_oligomers, overlap_length, rough_oligomer_size, high_overlap = ca.possible_oligomers(rough_oligomer)

        # select optimal oligomers
        t_3_free_oligomers, t_3_free_overlap_len = self.ups.t_3_free(list_of_oligomers, overlap_length)
        optimal_oligomers, optimal_overlap_len = self.ups.gc_tm_optimal(t_3_free_oligomers, t_3_free_overlap_len)
        final_oligomers, final_overlaps = self.ups._smallest_oligomers(optimal_oligomers, optimal_overlap_len, rough_oligomer, rough_oligomer_size, high_overlap)

        # returns final oligomer list with their corresponding overlap lengths
        return final_oligomers, final_overlaps

# Generate files containing the output
class Output(LinearOligomers, CircularOligomers):
    # generates excel and fasta file of the oligomers
    def output_files(self, x):
        workbook = xlsxwriter.Workbook(f'/app/output_script/oligomers_{self.user}_{self.oligomer_size}_{self.overlap_size}_{int(self.optimal_temp)}_{self.cluster_size}_{x}.xlsx')
        worksheet = workbook.add_worksheet()
        row = 1
        col = 0
        fasta_records = []
        for cluster in range(len(self.cluster_five_to_three)):
            for oligomer in range(len(self.cluster_five_to_three[cluster])):
                overlap_len = self.overlap_cluster[cluster][oligomer]
                oligo_len = len(self.cluster_five_to_three[cluster][oligomer])
                oligo = self.cluster_five_to_three[cluster][oligomer]
                clust = self.clusters[cluster][oligomer]
                if overlap_len != 0:
                    overlap_mt = round(mt.Tm_NN(clust[-overlap_len:], nn_table=mt.DNA_NN3), 2)
                else:
                    overlap_mt = 0
                oligo_len = len(oligo)
                record = SeqRecord(Seq.Seq(oligo), f'Oligomer_{cluster + 1}.{oligomer +1}', description=f'Overlap Length: {overlap_len}, Oligomer Length: {oligo_len}, Overlap Melting Temperature: {overlap_mt}')
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
                worksheet.write(row, col + 6, self.score[cluster][oligomer])
                worksheet.write(row, col + 7, "".join([i for i in self.fault[cluster][oligomer]]))
                worksheet.write(row, col + 8, ", ".join([i for i in self.repeats_position[cluster][oligomer]]))
                row += 1
        workbook.close()

        SeqIO.write(fasta_records, f'/app/output_script/oligomers_{self.user}_{self.oligomer_size}_{self.overlap_size}_{int(self.optimal_temp)}_{self.cluster_size}_{x}.fasta', "fasta")
        print(
            f'oligomers_{self.user}_{self.oligomer_size}_{self.overlap_size}_{int(self.optimal_temp)}_{self.cluster_size}_{x} as fasta and xlsx files')
        return self.cluster_five_to_three

o = Output.get_inputs()
x = input("Is the sequence linear or circular (l/c)? ")
if x == 'l':
    final_oligomers, final_overlap_len = o.design_oligomers()
else:
    final_oligomers, final_overlap_len = o._design_oligomers()
o.oligomer_design(final_oligomers, final_overlap_len)
v = o.output_files(x)

# parser.add_argument('-d', '--default_tm_c', help="Change defaults of -tm and -c",
#                             default=False)

# parser.add_argument('-dr', '--default_ranges', help="Change defaults of -tmr and -cr",
#                             default=False)
# if args.default_tm_c != False:
#     cls.set_default_tm_c(parser)
#
# if args.default_ranges != False:
#     cls.set_default_range(parser)
#
# args1 = parser.parse_args()

# @staticmethod
    # def set_default_tm_c(parse):
    #     optimal_tm, cluster_size,  = input(
    #         "Enter Parameters <Optimal Melting Temperature> <Oligomers in a Cluster>: ").split()
    #     parse.set_defaults(optimal_tm=int(optimal_tm))
    #     parse.set_defaults(cluster_size=int(cluster_size))
    #
    # @staticmethod
    # def set_default_range(parse):
    #     temp_range, cluster_range,  = input(
    #         "Enter Parameters <Overlap Melting Temperature Range> <Range of Oligomers in a Cluster>: ").split()
    #     parse.set_defaults(temp_range=int(temp_range))
    #     parse.set_defaults(cluster_range=int(cluster_range))
