# DNA Assembly Designer

DNA Assembly Designer is an easy to use web program for assembling short genomes that selects the best possible oligomer sequences based on the users' target parameters. 

## Installation
Dsembler can be installed via GitHub. Users should have docker installed on their computer to run Dsembler locally. Docker can be downloaded from https://docs.docker.com/get-docker/ [Ensure your BIOS settings are compatible with the docker application]. Run the following commands on your terminal

```
$ git clone https://github.com/sblabkribb/dsembler.git

$ ./docker_build.sh

$ ./docker_run.sh
```
==========================================================================
## Usage

The user interface is simple to navigate through as shown below:
![Screenshot](/images/main_form.png "Input Page")
### Form

#### Input

1. **Gene Sequence:** Users can input the gene sequence as plain text. The Application will not run if non-DNA elements are used (spaces and newlines are removed by the program automatically).

2. **Oligomer & Overlap size:** Next, input the target oligomer size and overlap size. The oligomers must be at least 50bp, and overlaps at least 20bp and must be at least 30bp smaller than the target oligomer size.

3. **Melting Temperature & Temperature range:** The target melting temperature is inputted in ºC. An acceptable range is required from the user, i.e ± nºC. The default and minimum range is melting temp ± 2.5ºC. 

4. **Sequence Orientation: ** Check the box if the Sequence you are using is circular. If left unchecked, the sequence is processed as a linear sequence

5. **Save Data: **If you'd like to save your parameters, check this box. Note, you should have signed up to access previous records.

#### Output

Two buttons appear on the same page once the target parameters are submitted: one to download the excel sheet output, and the other to download the FASTA file.

**Excel File**

![CropScreenshot](https://github.com/sblabkribb/oligomer_assembler/blob/main/images/m13example_excel.png "M13 Bacteriophage Example- Excel File")

The data is filled into the excel as the following columns. 

| Cluster Number| Oligomer Number| Oligomer Sequence (5' to 3') | Overlap Length | Overlap Melting Temperature | Overlap Score | Score Faults | Repeats Sequences|
| ------------- |-------------| -------------| -------------| -------------| -------------| -------------| ------------- |

Score of each oligomer's overlap:
- **H**: Overlap melting temperature is higher than the specified target temperature _(score = overlap_temp - temp_range_high)_
- **L**: Overlap melting temperature is higher than the specified target temperature _(score = temp_range_low - overlap_temp)_
- **R**: There are repeats (10bp) within oligomer overlaps within each cluster _(score = 10)_
- **T**: The overlap of each oligomer is not the smallest/largest within the possible oligomer range, and there is 'T' at the 3' end _(score = 1)_
- **G**: The last 5 bp of 3' end have 3 consecutive 'G', 'C', or a combination of both. _(score = 1)_
Overlap score is 0 if there are no possible areas of errors.

**FASTA File**



The oligomer FASTA file can be used for easy visualization of all oligomers. It can allow for amendments on other DNA visualising software, such as SnapGene.

**_3. The refresh button allows you to remove the existing parameters without reloading the page_**

#### Previous Data

Login             | Previous Work
:-------------------------:|:-------------------------:
![Screenshot](/images/login.png "Login") |  ![Screenshot](/images/previous_data.png "Previous Work")

Users can sign up and login to access their previous data. The database was created using the Flask-SQLite.

## Python 3.8.5
### Libraries used:
- Biopython:
This extensive library was used to parse the FASTA files, calculate melting temperature, calculate GC, check Alignment scores between two sequences.
- Flask:
This library/tool was used to make a simple web-based user interface. Related libraries such as Flask_bootstrap, wtforms, flask_pymongo supported the page as well

### Workflow

The outputs inclue a textfile that shows the appropriate orientation of the oligomers, an excel file that includes all the sequences in 5' to 3' format as well as a scoring system, and a zip file containing fasta files of all the oligomers. The algorithm is outlined below figure:

![Schematic](/images/workflow.png "Algorithm Workflow")

> Note: Melting Temperature of overlaps are calculated based on the Nearest Neighbour Equation and the Sugimoto (1996) thermodynamic table.
To build and run the docker image:

### Support
### Roadmap
### Authors and Acknowledgements
