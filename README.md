# Dsembler

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

Detailed User MAnual can be found at https://github.com/sblabkribb/dsembler/blob/main/manual.pdf

The user interface is simple to navigate through as shown below:
![Screenshot](/images/main_form.png "Input Page")
First enter in the required parameters to generate appropriate oligomers for your gene assembly.
- Gene Sequence
- Oligomer Size
- Overlap Size
- Target Melting Temperature of overlaps 
- Acceptable range of target melting temperature
- Sequence Orientation
- User ID (if saving query)

> Note: Melting Temperature of overlaps are calculated based on the Nearest Neighbour Equation and the Sugimoto (1996) thermodynamic table.

Two buttons will appear on the same page once the target parameters are submitted: 
- Download the Excel file
- Download the FASTA file


#### Previous Data

Login             | Previous Work
:-------------------------:|:-------------------------:
![Screenshot](/images/login.png "Login") |  ![Screenshot](/images/previous_data.png "Previous Work")

Users can sign up and login to access their previous data. The database was created using the Flask-SQLite.

==============================================================================

## Python 3.8.5
### Libraries used:
- Biopython:
This extensive library was used to parse the FASTA files, calculate melting temperature, calculate GC, check Alignment scores between two sequences.
- Flask:
This library/tool was used to make a simple web-based user interface. Related libraries such as Flask_bootstrap, wtforms, flask_SQLite supported the page as well

==============================================================================
### Workflow

The outputs inclue a textfile that shows the appropriate orientation of the oligomers, an excel file that includes all the sequences in 5' to 3' format as well as a scoring system, and a FASTA file containing meta data on all the oligomers. The algorithm is outlined below figure:

![Schematic](/images/workflow.png "Algorithm Workflow")
