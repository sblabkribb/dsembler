FROM python:3.8-slim-buster

WORKDIR /app

RUN apt-get update

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

RUN apt-get install -y procps vim

<<<<<<< HEAD
RUN mkdir -p /app/output

COPY . .

CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0"]

=======
COPY . .

CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0"]
#CMD [ "/bin/bash"]
>>>>>>> 480ff4c5c4a5310ff5414c1474b635e34c6c908b
