# python ollama-quote-parser.py.py

import gc
import re
import textwrap
import time
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def addToDict(arg1, arg2, arg3):
    if arg1 not in arg2:
        arg2[arg1] = arg3
    else:
        arg2[arg1] += arg3


models = [
    "ordinary differential equations?",
    "partial differential equations?|PDEs?",
    "agent-?based|individual-?based|automatas?",
    "stochastic differential equations?|stochastic process|SDEs?",
    "statistical|Bayesian",
]

labels = ["ode", "pde", "ab", "sto", "stat"]

viruses = {

    # Human viruses ####################

    "Adenovirus": [
        "Adenovirus",
        "Adenoviruses",
        "AdV",
        "ADV",
    ],

    "Adenovirus type 5": [
        "Adenovirus type 5",
        "Ad5",
    ],

    "Adeno-associated virus": [
        "Adeno-associated virus",
        "AAV",
    ],

    "Chikungunya virus": [
        "Chikungunya virus",
        "CHIKV",
    ],

    "Cytomegalovirus": [
        "Cytomegalovirus",
        "CMV",
        "Human cytomegalovirus",
        "HCMV",
    ],

    "Murine cytomegalovirus": [
        "Murine cytomegalovirus",
        "mCMV",
    ],

    "Dengue virus": [
        "Dengue",
        "Dengue virus",
        "DENV",
        "DENV-1",
        "DENV-2",
        "DENV-3",
        "DENV-4",
        "Dengue virus serotype 1",
        "Dengue virus serotype 2",
        "Dengue virus serotype 3",
        "Dengue virus serotype 4",
        "American Dengue Virus 2",
        "Dengue virus TCI strain",
    ],

    "Epstein-Barr virus": [
        "Epstein-Barr virus",
        "EBV",
    ],

    "Ebola virus": [
        "Ebola virus",
        "EBOV",
    ],

    "Hepatitis A virus": [
        "Hepatitis A virus",
        "HAV",
    ],

    "Hepatitis B virus": [
        "Hepatitis B virus",
        "HBV",
    ],

    "Hepatitis C virus": [
        "Hepatitis C virus",
        "HCV",
    ],

    "Hepatitis D virus": [
        "Hepatitis D",
        "Hepatitis D virus",
        "Hepatitis Delta virus",
        "HDV",
    ],

    "Herpes simplex virus": [
        "Herpes Simplex Virus",
        "Herpes simplex virus",
        "HSV",
        "HSV-1",
        "HSV-2",
        "Herpes simplex virus 1",
        "Herpes simplex virus 2",
    ],

    "Human bocavirus": [
        "Human bocavirus",
        "hBoV",
    ],

    "Human coronavirus": [
        "Human coronavirus",
    ],

    "Human enterovirus": [
        "Human enterovirus",
        "Enterovirus",
    ],

    "Human herpesvirus 6": [
        "Human herpesvirus 6",
        "HHV-6",
        "HH6",
        "Human herpesvirus 6B",
    ],

    "Human herpesvirus 7": [
        "Human herpesvirus 7",
        "HHV-7",
    ],

    "Human immunodeficiency virus": [
        "Human immunodeficiency virus",
        "HIV",
        "HIV-1",
        "Human immunodeficiency virus type 1",
    ],

    "Human metapneumovirus": [
        "Human metapneumovirus",
        "hMPV",
        "MPV",
    ],

    "Human rhinovirus": [
        "Human rhinovirus",
    ],

    "Human T-cell lymphotropic virus": [
        "Human T-cell lymphotropic virus",
        "HTLV",
        "HTLV-1",
        "HTLV-I",
        "HTLV-2",
        "Human T-cell lymphotropic virus type 1",
        "Human T-cell lymphotropic virus type 2",
    ],

    "Lymphocytic choriomeningitis virus": [
        "Lymphocytic choriomeningitis virus",
        "LCMV",
        "LCMV ARM",
    ],

    "Measles virus": [
        "Measles virus",
        "MeV",
    ],

    "MERS-CoV": [
        "MERS-CoV",
        "Middle East respiratory syndrome coronavirus",
    ],

    "Nipah virus": [
        "Nipah virus",
    ],

    "Norovirus": [
        "Norovirus",
        "Norovirus GII.4",
    ],

    "Parainfluenza virus": [
        "Parainfluenza virus",
        "PIV",
    ],

    "Poliovirus": [
        "Poliovirus",
        "Poliovirus type 1",
        "Mahoney strain",
    ],

    "Rabies virus": [
        "Rabies virus",
    ],

    "Respiratory syncytial virus": [
        "Respiratory syncytial virus",
        "Human respiratory syncytial virus",
        "RSV",
        "HRSV",
    ],

    "Rubella virus": [
        "Rubella virus",
    ],

    "SARS coronavirus": [
        "SARS coronavirus",
        "SARS-CoV",
    ],

    "SARS-CoV-2": [
        "SARS-CoV-2",
        "COVID-19",
        "Alpha SARS-CoV-2",
        "Beta SARS-CoV-2",
        "Delta SARS-CoV-2",
        "Omicron SARS-CoV-2",
        "BA.2.12.1 SARS-CoV-2",
        "BA.4 SARS-CoV-2",
        "BA.5 SARS-CoV-2",
        "F13-E SARS-CoV-2",
    ],

    "Varicella-zoster virus": [
        "Varicella-zoster virus",
        "VZV",
    ],

    "Variola virus": [
        "Variola virus",
    ],

    "Vaccinia virus": [
        "Vaccinia virus",
    ],

    "West Nile virus": [
        "West Nile virus",
        "WNV",
    ],

    "Yellow fever virus": [
        "Yellow fever virus",
        "YFV",
    ],

    "Zika virus": [
        "Zika virus",
        "ZIKV",
    ],


    # Influenza  ####################

    "Influenza A virus": [
        "Influenza A virus",
        "IAV",
        "H1N1 influenza virus",
        "H3N2 influenza virus",
        "H5N1 influenza virus",
        "H7N9 influenza virus",
        "H9N2 avian influenza virus",
        "H3N8 virus",
        "A(H1N1)pdm09",
        "pH1N1 2009 virus",
        "1918 influenza virus",
        "H3N2v virus",
        "Influenza A/Puerto Rico/8/34",
        "PR8",
        "Influenza A/WSN/1933",
        "A/Aichi/2/68",
        "A/Texas/91",
        "A/Hong Kong/123/77",
        "IAV 1–126 TX/98",
    ],

    "Influenza B virus": [
        "Influenza B virus",
        "B/Victoria lineage influenza B virus",
    ],

    "Avian influenza virus": [
        "Avian influenza virus",
    ],

    "Junin arenavirus": [
        "Junin arenavirus",
    ],

    "Pichinde virus": [
        "Pichinde virus",
    ],

    "Morogoro virus": [
        "Morogoro virus",
        "Morogoro arenavirus",
        "MORV",
    ],

    "Simian immunodeficiency virus": [
        "Simian immunodeficiency virus",
        "SIV",
        "SIVmac251",
    ],

    "Simian/human immunodeficiency virus": [
        "Simian human immunodeficiency virus",
        "Simian/human immunodeficiency virus",
        "SHIV",
        "SHIV-P3",
        "SHIV-KS661",
    ],


    # Animal viruses ####################

    "African horsesickness virus": [
        "African horsesickness virus",
        "AHSV",
    ],

    "African swine fever virus": [
        "African swine fever virus",
        "ASFV",
        "ASF",
    ],

    "Bovine respiratory syncytial virus": [
        "Bovine respiratory syncytial virus",
    ],

    "Bovine viral diarrhea virus": [
        "Bovine viral diarrhea virus",
        "BVDV",
    ],

    "Foot-and-mouth disease virus": [
        "Bovine FMD virus",
        "Foot-and-mouth disease virus",
        "FMDV",
    ],

    "Canine distemper virus": [
        "Canine distemper virus",
        "CDV",
    ],

    "Equine encephalitis virus": [
        "Equine encephalitis virus",
    ],

    "Equine infectious anemia virus": [
        "Equine infectious anemia virus",
        "EIAV",
    ],

    "Equine influenza virus": [
        "Equine influenza virus",
    ],

    "Eastern equine encephalitis virus": [
        "Eastern equine encephalitis virus",
        "EEEV",
    ],

    "Hantavirus": [
        "Hantavirus",
    ],

    "Hantaan virus": [
        "Hantaan virus",
    ],

    "Hendra virus": [
        "Hendra virus",
    ],

    "Japanese encephalitis virus": [
        "Japanese encephalitis virus",
        "JEV",
    ],

    "Kyasanur Forest disease virus": [
        "Kyasanur Forest disease virus",
        "KFDV",
    ],

    "Lassa virus": [
        "Lassa virus",
        "LASV",
    ],

    "Marburg virus": [
        "Marburg virus",
        "Marburgvirus",
    ],

    "Marek's disease virus": [
        "Marek's disease virus",
        "MDV",
    ],

    "Mouse hepatitis virus": [
        "Mouse hepatitis virus",
        "MHV-1",
    ],

    "Murine leukemia virus": [
        "Murine leukemia virus",
    ],

    "Murray Valley encephalitis virus": [
        "Murray Valley encephalitis virus",
    ],

    "Porcine circovirus type 2": [
        "Porcine circovirus type 2",
        "PCV2",
    ],

    "Porcine respiratory coronavirus": [
        "Porcine respiratory coronavirus",
        "PRCV",
    ],

    "Porcine reproductive and respiratory syndrome virus": [
        "Porcine reproductive and respiratory syndrome virus",
        "PRRSV",
    ],

    "Porcine rotavirus": [
        "Porcine rotavirus",
    ],

    "Pseudorabies virus": [
        "Pseudorabies virus",
        "PRV",
    ],

    "Rabbit hemorrhagic disease virus": [
        "Rabbit hemorrhagic disease virus",
        "RHDV",
    ],

    "Rift Valley fever virus": [
        "Rift Valley fever virus",
        "RVFV",
    ],

    "Ross River virus": [
        "Ross River virus",
        "RRV",
    ],

    "Rous sarcoma virus": [
        "Rous sarcoma virus",
    ],

    "Seoul virus": [
        "Seoul virus",
    ],

    "Semliki Forest virus": [
        "Semliki Forest virus",
    ],

    "Sindbis virus": [
        "Sindbis virus",
    ],

    "Sin Nombre virus": [
        "Sin Nombre virus",
        "SNV",
    ],

    "St. Louis encephalitis virus": [
        "St. Louis encephalitis virus",
    ],

    "Theiler murine encephalomyelitis virus": [
        "Theiler murine encephalomyelitis virus",
    ],

    "Tick-borne encephalitis virus": [
        "Tick-borne encephalitis virus",
        "TBE",
        "TBEV",
    ],

    "Venezuelan equine encephalitis virus": [
        "Venezuelan equine encephalitis virus",
        "VEEV",
    ],

    "Vesicular stomatitis virus": [
        "Vesicular stomatitis virus",
        "VSV",
    ],

    "Western equine encephalitis virus": [
        "Western equine encephalitis virus",
    ],

    "Bluetongue virus": [
        "Bluetongue virus",
        "Bluetongue virus 8",
        "BTV",
        "BTV8",
    ],

    "Ostreid herpesvirus 1": [
        "Ostreid herpesvirus 1",
        "OsHV-1",
        "Oyster herpesvirus",
    ],

    "Ostreococcus tauri virus 1": [
        "Ostreococcus tauri virus 1",
        "OTV1_139",
    ],

    "Tilapia lake virus": [
        "Tilapia lake virus",
        "TiLV",
    ],

    "Tulane virus": [
        "Tulane virus",
    ],

    "Palyam virus": [
        "Palyam virus",
    ],

    "Dermo virus": [
        "Dermo virus",
    ],

    "MSX virus": [
        "MSX virus",
    ],

    "Brown ring disease virus": [
        "Brown ring disease virus",
    ],


    # Plant viruses ####################

    "Barley yellow dwarf virus": [
        "Barley yellow dwarf virus",
        "Barley and cereal yellow dwarf viruses",
    ],

    "Bahia bark scaling of citrus virus": [
        "Bahia bark scaling of citrus virus",
    ],

    "Broad bean mottle virus": [
        "Broad bean mottle virus",
    ],

    "Cauliflower mosaic virus": [
        "Cauliflower mosaic virus",
        "CaMV",
    ],

    "Capsicum mild mottle virus": [
        "Capsicum mild mottle virus",
    ],

    "Clover yellow vein virus": [
        "Clover yellow vein virus",
    ],

    "Cucumber mosaic virus": [
        "Cucumber mosaic virus",
    ],

    "Obuda mosaic virus": [
        "Obuda mosaic virus",
    ],

    "Pepper mild mottle virus": [
        "Pepper mild mottle virus",
        "PMMoV",
    ],

    "Pepper mild green mosaic virus": [
        "Pepper mild green mosaic virus",
        "Tobacco mild green mosaic virus",
    ],

    "Plum pox virus": [
        "Plum pox virus",
        "PPV",
    ],

    "Potato virus X": [
        "Potato virus X",
        "PVX",
    ],

    "Potato virus Y": [
        "Potato virus Y",
        "PVY",
    ],

    "Potato leafroll virus": [
        "Potato leafroll virus",
    ],

    "Rice tungro virus": [
        "Rice tungro virus",
    ],

    "Satellite tobacco mosaic virus": [
        "Satellite tobacco mosaic virus",
        "STMV",
    ],

    "Satellite tobacco necrosis virus": [
        "Satellite tobacco necrosis virus",
        "STNV",
    ],

    "Southern bean mosaic virus": [
        "Southern bean mosaic virus",
        "SBMV",
    ],

    "Tomato mosaic virus": [
        "Tomato mosaic virus",
    ],

    "Tomato spotted wilt virus": [
        "Tomato spotted wilt virus",
    ],

    "Tomato yellow leaf curl virus": [
        "Tomato yellow leaf curl virus",
        "TYLCV",
    ],

    "Tobacco mosaic virus": [
        "Tobacco mosaic virus",
        "TMV",
    ],

    "Tobacco etch virus": [
        "Tobacco etch virus",
        "TEV",
    ],

    "Turnip crinkle virus": [
        "Turnip crinkle virus",
        "TCV",
    ],

    "Turnip mosaic virus": [
        "Turnip mosaic virus",
        "TuMV",
    ],

    "Tobamovirus": [
        "Tobamovirus",
    ],


    # Insect / fungal viruses ####################

    "Baculovirus": [
        "Baculovirus",
        "Baculoviruses",
    ],

    "Autographa californica nucleopolyhedrovirus": [
        "Autographa californica nucleopolyhedrovirus",
        "AcNPV",
    ],

    "Lymantria dispar multiple nucleopolyhedrovirus": [
        "Gypsy moth baculovirus",
        "Lymantria dispar multiple nucleopolyhedrovirus",
        "LdMNPV",
    ],

    "Spodoptera frugiperda multicapsid nucleopolyhedrovirus": [
        "Spodoptera frugiperda multicapsid nucleopolyhedrovirus",
        "SfMNPV",
    ],

    "Helicoverpa zea S nucleopolyhedrovirus": [
        "Helicoverpa zea S nucleopolyhedrovirus",
    ],

    "Nuclear polyhedrosis virus": [
        "Nuclear polyhedrosis virus",
        "Nucleopolyhedroviruses",
    ],

    "Ascovirus": [
        "Ascovirus",
        "Ascoviruses",
    ],

    "Iridovirus": [
        "Iridovirus",
        "Iridoviruses",
    ],

    "Hypoviruses": [
        "Hypoviruses",
        "Chestnut blight hypoviruses",
    ],


    # Bacteriophages ####################

    "Bacteriophages": [
        "Bacteriophages",
        "Bacteriophage",
        "Phages",
    ],

    "T7 phage": [
        "T7 phage",
        "T7 virus",
    ],

    "P0 phage": [
        "P0 phage",
    ],

    "P1 phage": [
        "P1 phage",
    ],

    "Phage PR": [
        "Phage PR",
    ],

    "Phage PS": [
        "Phage PS",
    ],

    "Myoviruses": [
        "Myoviruses",
    ],

    "Podoviruses": [
        "Podoviruses",
    ],

    "CTXφ": [
        "CTXφ",
    ],

    "Lytic phages": [
        "Lytic phages",
    ],

    "Prophages": [
        "Prophages",
    ],
}


def quote_quanta(file_name, age):
    data = defaultdict(lambda: defaultdict(int))
    extra_data = []
    model_data = []
    virus_data = defaultdict(lambda: defaultdict(int))

    used_labels = []

    with open(file_name) as infile:
        model_flag = 0

        for line in infile:
            try:
                s = line.split("|")[-2]
                paper = line.split("|")[1].split(" - ")[0]

                extra_flag = 0

                if (model_flag == 1) and (line.find("#") == -1):
                    used_labels = []

                    for i in range(len(models)):
                        if re.search(models[i], s, flags=re.IGNORECASE):

                            if (
                                models[i]
                                == "partial differential equations?|PDEs?"
                            ):
                                if re.search("age", s, flags=re.IGNORECASE):
                                    age += 1

                            addToDict(
                                labels[i],
                                data[paper],
                                1,
                            )

                            extra_flag = 1
                            used_labels.append(labels[i])

                    if extra_flag == 0:

                        if re.search("ODEs?", s):
                            label = "ode"
                            addToDict(label, data[paper], 1)
                            used_labels.append(label)

                        elif re.search(
                            "differential equations?",
                            s,
                            flags=re.IGNORECASE,
                        ):
                            label = "de"
                            addToDict(label, data[paper], 1)
                            used_labels.append(label)

                        elif re.search(
                            "model",
                            s,
                            flags=re.IGNORECASE,
                        ):
                            label = "mod"
                            addToDict(label, data[paper], 1)
                            used_labels.append(label)
                            model_data.append(s)

                        else:
                            label = "extra"
                            addToDict(label, data[paper], 1)
                            used_labels.append(label)
                            extra_data.append(s)

                    model_flag = 0

                elif line.find("#") == -1:

                    for virus_label, aliases in viruses.items():

                        for alias in aliases:

                            pattern = rf"(?<!\w){re.escape(alias)}(?!\w)"

                            if re.search(
                                pattern,
                                s,
                                flags=re.IGNORECASE,
                            ):
                                for label in used_labels:
                                    addToDict(
                                        virus_label,
                                        virus_data[label],
                                        1,
                                    )
                                break

            except IndexError:
                pass

            if line.find("#") != -1:
                model_flag = 1

    paper_sum = defaultdict(int)

    for paper in data:
        for key, value in data[paper].items():
            addToDict(key, paper_sum, value)

    paper_sum["ode"] += paper_sum.get("de", 0)

    for k in ["extra", "mod", "de"]:
        paper_sum.pop(k, None)

    categories = list(paper_sum.keys())
    values = list(paper_sum.values())

    data_pairs = sorted(
        zip(values, categories),
        reverse=True,
    )

    sorted_sizes = [pair[0] for pair in data_pairs]
    sorted_labels = [pair[1] for pair in data_pairs]

    df_virus = pd.DataFrame(virus_data).fillna(0).astype(int)

    if "de" in df_virus.columns:
        if "ode" not in df_virus.columns:
            df_virus["ode"] = 0

        df_virus["ode"] += df_virus["de"]

    df_virus = df_virus.reindex(
        columns=labels,
        fill_value=0,
    )


    return (
        np.array(sorted_sizes),
        np.array(sorted_labels),
        age,
        df_virus,
    )


def get_top_viruses(df_virus, n=20):
    """Return the top n viruses ranked by total paper count."""
    
    totals = df_virus.sum(axis=1)

    return totals.sort_values(
        ascending=False
    ).head(n).index

def plot_virus_heatmap(
    df_virus,
    filename="figures/virus_heatmap.png",
):

    if df_virus.empty:
        print("No virus data to plot.")
        return

    top_20 = get_top_viruses(df_virus, n=20)

    df_sorted = df_virus.loc[top_20]

    plt.figure()

    sns.heatmap(
        df_sorted,
        annot=True,
        fmt="d",
        cmap="YlGnBu",
        cbar_kws={"label": "Paper Count"},
        linewidths=0.5,
        annot_kws={"size": size * 0.6}
    )

    plt.xlabel("Model Type")
    plt.ylabel("Virus")

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


def plot_virus_log_heatmap(
    df_virus,
    filename="figures/virus_heatmap.png",
):

    if df_virus.empty:
        print("No virus data to plot.")
        return

    top_20 = get_top_viruses(df_virus, n=20)

    df_sorted = df_virus.loc[top_20]

    plt.figure()

    sns.heatmap(
        np.log10(df_sorted.replace(0, np.nan)),
        cmap="YlGnBu",
        cbar_kws={"label": "Paper Count"},
        linewidths=0.5,
    )

    plt.xlabel("Model Type")
    plt.ylabel("Virus")

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


def plot_virus_stacked_bar(
    df_virus,
    filename="figures/virus_stacked.png",
):
    """Horizontal Stacked Bar Chart"""

    if df_virus.empty:
        print("No virus data to plot.")
        return

    top_20 = get_top_viruses(df_virus, n=20)

    df_sorted = df_virus.loc[top_20].copy()

    df_sorted = df_sorted.loc[
        df_sorted.sum(axis=1)
        .sort_values(ascending=True)
        .index
    ]

    ax = df_sorted.plot(
        kind="barh",
        stacked=True,
        colormap="viridis",
    )

    plt.xlabel("Number of Mentions")
    plt.ylabel("Virus")

    plt.legend(
        title="Model Type",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


def plot_virus_log_stacked_bar(
    df_virus,
    filename="figures/virus_stacked.png",
):
    """Horizontal Stacked Bar Chart"""

    if df_virus.empty:
        print("No virus data to plot.")
        return

    top_20 = get_top_viruses(df_virus, n=20)

    df_sorted = df_virus.loc[top_20].copy()

    df_sorted = df_sorted.loc[
        df_sorted.sum(axis=1)
        .sort_values(ascending=True)
        .index
    ]

    ax = df_sorted.plot(
        kind="barh",
        colormap="viridis",
    )

    ax.set_xscale("log")

    plt.xlabel("Number of Mentions")
    plt.ylabel("Virus")

    plt.legend(
        title="Model Type",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()



# Analysis ollama output ####################

sizes_article, labels_article, age, df_virus_art = quote_quanta(
    "quotes_gemma_article",
    age=0,
)

sizes_review, labels_review, age, df_virus_rev = quote_quanta(
    "quotes_gemma",
    age=0,
)

# Plot model comparison ####################
from brokenaxes import brokenaxes


bax = brokenaxes(
    ylims=((0, 400), (5900, 6000)),
    hspace=0.125,
)


ind_art = np.arange(len(labels_article))
ind_rev = np.zeros(len(labels_article))


for i, item in enumerate(labels_article):
    for j in range(len(labels_review)):
        if item == labels_review[j]:
            ind_rev[j] = i


width = 0.45


bax.bar(
    ind_art,
    sizes_article,
    width,
    label="Non-Review",
)

bax.bar(
    ind_rev + width,
    sizes_review,
    width,
    label="Review",
)

bax.axs[1].set_xticks(
    ind_art + width / 2,
    labels_article,
)

bax.legend()

plt.savefig("figures/combine.png")
plt.close()


# Combine article and review virus datasets ####################

#df_virus_total = (
#    df_virus_art
#    .add(df_virus_rev, fill_value=0)
#    .astype(int)
#)

#df_virus_total.to_csv("df_virus_total.csv")

df_virus_total = pd.read_csv(
    "df_virus_total.csv",
    index_col=0,
)

# Plot Virus comparison ####################
size = 20
axessize = 2

plt.rcParams['xtick.labelsize'] = size*0.9
plt.rcParams['ytick.labelsize'] = size
plt.rcParams['axes.labelsize'] = size
plt.rc('legend', fontsize=size)
plt.rc('legend', title_fontsize=size)
plt.rcParams['figure.figsize'] = [12, 12]
plt.rcParams['xtick.major.size'] = axessize*4
plt.rcParams['xtick.major.width'] = axessize
plt.rcParams['xtick.minor.width'] = axessize
plt.rcParams['ytick.major.width'] = axessize
plt.rcParams['ytick.minor.width'] = axessize

plot_virus_heatmap(df_virus_total,"figures/virus_heatmap.png")
plot_virus_stacked_bar(df_virus_total, "figures/virus_stacked_bar.png")

plot_virus_log_heatmap(df_virus_total,"figures/log_virus_heatmap.png")
plot_virus_log_stacked_bar(df_virus_total, "figures/log_virus_stacked_bar.png")


gc.collect()


