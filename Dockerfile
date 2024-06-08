# Dockerfile for EssenCDNA (final project for BENG/CSE/BIMM 182 @ UC San Diego)
# TODO: decide on what base image to use
FROM python:3

LABEL authors="Adrian Layer, Anthony Vasquez, Yasmin A. Jaber, Omar Halawa"

#TODO: decide on packages and their specific versions
# Ensuring up-to-date pip and importing necessary modules 
RUN pip install --upgrade pip && \
    pip install pandas==1.5.3 && \
    pip install scikit-learn==1.2.1

# Build using "docker build -t <TAG> ."
# Run using "docker run -it --rm <IMAGE ID> bash"
