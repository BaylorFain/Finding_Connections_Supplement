#nohup bash -c 'python ollama_reader.py &'
import numpy as np
import pymupdf
import os
import gc
import ollama
import re
import textwrap
import time
import subprocess

def AllText(doc, ToSplit = False):
    try:
        toc = doc.get_toc()

        all_text = ""
        for page in doc.pages(0, toc[-1][-1], 1):
            text = page.get_text("text")
            all_text = all_text + " " + text

    except IndexError:
        num_pages = doc.page_count

        all_text = ""
        for page in doc.pages(0, num_pages, 1):
            text = page.get_text("text")
            all_text = all_text + " " + text

    all_text = all_text.replace("\xa0"," ")
    all_text = all_text.replace(r"","")
    all_text = all_text.replace("\00","")

    all_text = all_text.replace("\n"," ")
    all_text = all_text.replace("  "," ")
    if ToSplit == True:
        all_text = all_text.replace(". ",".. ")
        all_text = all_text.split(". ")

    return(all_text)
    
def Ollama_Reader(context_text, question_text):
    messages = [
        {
            'role': 'system',
            'content': 'You are a helpful assistant that answers questions based on the provided context. If the answer is not in the context, simply state that you cannot find the answer.'
        },
        {
            'role': 'user',
            'content': f"Context: {context_text}\n\nQuestion: {question_text}"
        }
    ]

    response = ollama.chat(
       model='gemma3:4b-it-qat',
        messages=messages
    )

    return(response['message']['content'])

path = r"articles/"
csvname = "quotes"
with open(csvname, "w+") as outfile:
    "Nothing"

dir_list = sorted(os.listdir(path))
file_num = len(dir_list)

in_words = ["model"]
out_words = ["doi ", "doi.", "doi:", "DOI ", "DOI.", "DOI:", "http", "arXiv", "BioRxiv"]

count_file = 1
for pdfname in dir_list:
    try:
        subprocess.run(['sudo', 'systemctl', 'stop', 'ollama.service'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error stopping Ollama service: {e}")

    try:
        subprocess.run(['sudo', 'systemctl', 'start', 'ollama.service'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error starting Ollama service: {e}")
    time.sleep(0.5)

    os.system('clear')
    print(count_file,"/",file_num)
    print(pdfname)
    year = re.sub(r'\D','',pdfname)[:4]

    doc = pymupdf.open(path+pdfname)

    all_text = AllText(doc)

    filename = r"results/" + pdfname.replace(".pdf", " (summary)")
    with open(filename, "w+") as outfile:
        "Nothing"

    ref_break = 0
    big_response = ""
    chunks = textwrap.wrap(all_text, width = 3000, break_long_words = False, replace_whitespace = False)
    for chunk in chunks:
        if ref_break != 0:
            break
        if any(break_word in [word for word in chunk.split() if word.istitle()] for break_word in ["References", "REFERENCES"]):
            ref_break = 1
        if all(in_word in chunk for in_word in in_words):
            if any(2 < chunk.count(out_word) for out_word in out_words):
                "Nothing"
            else: 
                question = '''What models and associated viruses are found in this text? For each model create a section with the model name, then list a one sentence description of the specific math used to build the model and all associated viruses. Do not include any other text, explanations, or conversational filler.

Example Format:
# [insert model name here]
- [insert used math here]
- [insert virus here]
- [insert virus here]
- [insert virus here]

'''

                material = chunk
                llama_response = Ollama_Reader(chunk, question)
                with open(csvname, "a+") as outfile:
                    lines = llama_response.split("\n")
                    for line in lines:
                        string = '|'+ pdfname +'|'+','+str(year)+','+'|'+ line +'|'
                        print(string, file=outfile)
                    print("\n", file=outfile)
                big_response = big_response + llama_response

    with open(filename, "a+") as outfile:
        print(big_response, file=outfile)
        print("\n", file=outfile)

    with open(csvname, "a") as outfile:
        print('',file=outfile)                
    with open(filename, "a") as outfile:
        print('',file=outfile)

    count_file += 1
    gc.collect()
