# Finding Connections Supplement

### Table of Contents  
[Introduction](#introduction)  
[Software Install](#software-install)  
[Recreation](#recreation)  
[Reproduction](#reproduction) 
<a name="headers"/>

## Introduction

This repository contains the supplementary files mentioned in the manuscript: 
**Fain and Dobrovolny, (Submitted), "Finding Connections: A survey of within host viral dynamics modeling techniques"**

to recreate and reproduce the figures and analyses in

## Software Install

To recreate the plots in the figures, _Python_ is needed. 
>[!NOTE]
>The provided plotting codes should work with current releases of _Python_ packages. 

To reproduce the data for the plots, _CUDA_, _Python_, and _ollama_ are needed. C code is compiled using "nvcc". 
>[!NOTE]
>After the necessary software is installed, the codes provided will create the data for each plot

## Recreation

_Python_ is used to make the plots used in the manuscript. The data for each plot in a figure is provided and all that needs to be done is to run the provided code. 

#### To create plots for figures 1 and 11

Run the code [ollama-quote-parser.py](https://github.com/BaylorFain/Finding_Connections_Supplement/blob/main/Figure%201%20and%2011/ollama-quote-parser.py)

#### To create plots for figures 4 through 9

Run the code [Plotter.py](https://github.com/BaylorFain/Finding_Connections_Supplement/blob/main/Figure%204%20through%209/Plotter.py)

## Reproduction

#### To create the data for figures 1 and 11

1. Acquire all of the articles in [articles](https://github.com/BaylorFain/Finding_Connections_Supplement/tree/main/articles)
2. Then run [ollama_reader.py](https://github.com/BaylorFain/Finding_Connections_Supplement/blob/main/Ollama_Code/ollama_reader.py)
3. The output is used to create the plots for figure 1 and 11.

#### To create the data for figures 4 through 9

1. Run codes:
    - [AB.cu](https://github.com/BaylorFain/Finding_Connections_Supplement/blob/main/Model_Codes/AB/AB.cu)
    - [ODE.py](https://github.com/BaylorFain/Finding_Connections_Supplement/blob/main/Model_Codes/ODE.py)
    - [AS.py](https://github.com/BaylorFain/Finding_Connections_Supplement/blob/main/Model_Codes/AS.py)
    - [Gillespie_DV.cu](https://github.com/BaylorFain/Finding_Connections_Supplement/blob/main/Model_Codes/Gillespie_DV.cu)
    - [Gillespie_SV.cu](https://github.com/BaylorFain/Finding_Connections_Supplement/blob/main/Model_Codes/Gillespie_SV.cu)
    - [Tau_Leap_DV.py](https://github.com/BaylorFain/Finding_Connections_Supplement/blob/main/Model_Codes/Tau_Leap_DV.py)
    - [Tau_Leap_SV.py](https://github.com/BaylorFain/Finding_Connections_Supplement/blob/main/Model_Codes/Tau_Leap_SV.py)
2. The codes will have put the data into the directory [Data](https://github.com/BaylorFain/Finding_Connections_Supplement/tree/main/Model_Codes/Data). This directory is used to create the plots for figures 4 through 9.
    









