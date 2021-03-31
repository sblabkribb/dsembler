# DNA Assembly Designer

DNA Assembly Designer is an easy to use web program for assembling short genomes that selects the best possible oligomer sequences based on the users' target parameters. 

The outputs inclue a textfile that shows the appropriate orientation of the oligomers, an excel file that includes all the sequences in 5' to 3' format as well as a scoring system, and a zip file containing fasta files of all the oligomers. The algorithm is outlined below:

1. Parse FASTA file
2. Collect input parameters
3. Generate base oligomers without overlap sequences based on user target preferences
4. Construct an array of oligomers with all the possible overlaps (reasonable range of overlap size)
5. Check for the best possible oligomer(s) based on the target conditions of overlap melting temperature, GC content, GC clamp
> Note: Melting Temperature of overlaps are calculated based on the Nearest Neighbour Equation and the Sugimoto (1996) thermodynamic table.
6. Finalize oligomers based on the smallest oligomer (optimize based on cost as well)
7. Arrange oligmers in appropriate clusters based on overlap alignment ratios
8. Score each oligomer's overlap
    1. **H**: Overlap melting temperature is higher than the specified target temperature _(score = overlap_temp - temp_range_high)_
    2. **L**: Overlap melting temperature is higher than the specified target temperature _(score = temp_range_low - overlap_temp)_
    3. **R**: There are repeats (10bp) within oligomer overlaps within each cluster _(score = 10)_
    4. **T**: The overlap of each oligomer is not the smallest/largest within the possible oligomer range, and there is 'T' at the 3' end _(score = 1)_
    5. **G**: The last 5 bp of 3' end have 3 consecutive 'G', 'C', or a combination of both. _(score = 1)_
9. Write the respective Text, Excel, and Zip files

To build and run the docker image:

```
$ docker build dsembler -t dsembler

$ docker run --publish 5000:5000 dsembler

```

## Python 3.8.5
### Libraries used:
- Biopython:
This extensive library was used to parse the FASTA files, calculate melting temperature, calculate GC, check Alignment scores between two sequences.
- Itertools:
This library was used to create different combinations while checking for alignments between oligomer overlaps
- Flask:
This library/tool was used to make a simple user interface. Related libraries such as Flask_bootstrap, wtforms, flask_pymongo supported the page as well

> Note: M13 bacteriophage genome will be used as an example in this document

==========================================================================
### Usage
The user interface is simple to navigate through as shown below:
![Screenshot](https://github.com/sblabkribb/oligomer_assembler/blob/main/images/input_page.png "Input Page")

**_1. Input_**

![CropScreenshot](https://github.com/sblabkribb/oligomer_assembler/blob/main/images/gene_seq.png "Gene Sequnce Input")

Users can input the gene sequence as plain text or as a FASTA file. The file/text input must only have characters of DNA nucleotides (no spaces).

![CropScreenshot](https://github.com/sblabkribb/oligomer_assembler/blob/main/images/oligomer_overlap.png "Oligomer and Overlap Size Input")

Next input the target oligomer size and overlap size. The oligomers must be at least 50bp, and overlaps at least 20bp and must be at least 30bp smaller than the target oligomer size.

![CropScreenshot](https://github.com/sblabkribb/oligomer_assembler/blob/main/images/melting_temp.png "Melting Temperature Input")

The target melting temperature is inputted in ºC. An acceptable range is required from the user, i.e ± nºC. The default and minimum range is melting temp ± 2.5ºC. 

![CropScreenshot](https://github.com/sblabkribb/oligomer_assembler/blob/main/images/cluster.png "Cluster Size Input")

The target cluster size is obtained as the number of oligomers the user needs in each cluster. A range is also considered to reduce errors that occur due to repeats between oligomer overlaps within clusters.

**_2. Output_**

![Screenshot](https://github.com/sblabkribb/oligomer_assembler/blob/main/images/m13example.png "M13 Bacteriophage Example")

Two buttons appear on the same page once the target parameters are submitted: one to download the excel sheet output, and the other to download the zip file containing oligomers saved as FASTA files. 

**Excel File**

![CropScreenshot](https://github.com/sblabkribb/oligomer_assembler/blob/main/images/m13example_excel.png "M13 Bacteriophage Example- Excel File")

| Cluster Number| Oligomer Number| Oligomer Sequence (5' to 3') | Overlap Length | Overlap Melting Temperature | Overlap Score | Score Faults | Repeats Sequences|
| ------------- |-------------| -------------| -------------| -------------| -------------| -------------| ------------- |

Overlap score is 0 if there are no possible areas of errors.

**Zip File**

![CropScreenshot](https://github.com/sblabkribb/oligomer_assembler/blob/main/images/m13example_zipfile.png "M13 Bacteriophage Example- Zipfile File")

The oligomer FASTA files can be used for easy visualization and amendments on other DNA visualising software, such as SnapGene.

**_3. The refresh button allows you to remove the existing parameters without reloading the page_**

### Support
### Roadmap

The next step is to develop a database so people can access previous genomes or oligomers
### Authors and Acknowledgements
