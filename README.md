# Dsembler (DNA Assembly Designer)

Dsembler is an easy to use web program for assembling short genomes that selects the best possible oligomer sequences based on the users' target parameters.

## 1. Installation

Dsembler can be installed via GitHub. 

```
$ git clone https://github.com/sblabkribb/dsembler.git
```

### 1.1 Docker (Recommended)
Users should have docker installed on their computer to run Dsembler locally. Docker can be downloaded from <https://docs.docker.com/get-docker/> [Ensure your BIOS settings are compatible with the docker application]. Run the following commands on your terminal. Move to the dsembler directory

#### Build docker image build
```
$ docker build -t dsembler:latest .
```

Or pull the docker image from Dockerhub

```
$ docker pull sblabkribb/dsembler:latest
```

####  Docker run for web version (Linux)

```bash
$ docker run -d --rm -v $(pwd):/app --publish 5000:5000 --name dsembler dsembler
```
#### Docker run for web version (Powershell / CMD)
```powershell
> docker run -d --rm -v %cd%:/app --publish 5000:5000 --name dsembler dsembler
```
Replace `$(pwd)` and `%cd%` with the directory on your local machine you want to store the outputs

### 1.2 Python
Users can directly run dsembler via Python. However, many packages will be required to be installed before a successful run  (check requirements.txt).

#### Flask run for web version Bash
```bash Bash
$ export FLASK_APP=app
$ flask run
```
#### Flask run for web version Powershell
```powershell
> $env:FLASK_APP = "app"
> flask run
```

=======================================================================

## 2. Usage

Detailed User Manual can be found at <https://github.com/sblabkribb/dsembler/blob/main/docs/manual.pdf>


=======================================================================

## 3. Python 3.8.5

### Libraries used:

-   Biopython: This extensive library was used to parse the FASTA files, calculate melting temperature, calculate GC, check Alignment scores between two sequences.
-   Flask: This library/tool was used to make a simple web-based user interface. Related libraries such as Flask_bootstrap, wtforms, flask_SQLite supported the page as well

=======================================================================

## 4. Workflow

The algorithm workflow is outlined below figure:

![Schematic](images/workflow.png "Algorithm Workflow")
