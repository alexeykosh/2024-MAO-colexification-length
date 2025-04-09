import pandas as pd
import os
import numpy as np

lg_codes = {'he': 'Hebrew', 
            'mr': 'Marathi',
            'is': 'Icelandic',
            'tl': 'Tagalog',
            'yo': 'Yoruba',
            'pt': 'Portuguese',
            'ml': 'Malayalam',
            'kn': 'Kannada',
            'ar': 'Arabic',
            'bn': 'Bengali',
            'af': 'Afrikaans',
            'fi': 'Finnish',
            'id': 'Indonesian',
            'et': 'Estonian',
            'en': 'English',
            'fa': 'Farsi',
            'tt': 'Tatar',
            'tr': 'Turkish'}

# Read the tsv files in every subfolder of data/results and combine them into one
tsv_files = []
freq_files = []
for root, dirs, files in os.walk("data/results"):
    for file in files:
        if file.endswith(".tsv"):
            tsv_files.append(os.path.join(root, file))
        if file.endswith(".txt"):
            freq_files.append(os.path.join(root, file))

# Read all tsv files into a single dataframe and add language codes (folder names) as a column
df = pd.DataFrame()
for file in tsv_files:
    lang = file.split("/")[2]
    df_lang = pd.read_csv(file, sep="\t")
    df_lang["lang"] = lang
    df = pd.concat([df, df_lang])

# List of all unique words from Pimentel et al. dataset
words = df["word"].unique()

# Read the frequency files and combine them into one single dataframe
df_f = pd.DataFrame()
for file in freq_files:
    lang = file.split("/")[2]
    df_freq = pd.read_csv(file, sep="\t", header=None, 
                          quoting=3, on_bad_lines='skip')
    if len(df_freq.columns) == 3:
        df_freq.columns = ["index", "word", "count"]
    elif len(df_freq.columns) == 4:
        df_freq.columns = ["index", "word", "lemma", "count"]
    df_freq["lang"] = lang
    df_freq = df_freq.dropna()
    df_freq = df_freq[~df_freq["count"].isin([float("inf"), float("-inf")])]
    df_freq["count"] = df_freq["count"].astype(int)
    df_freq['freq'] = df_freq['count'] / df_freq['count'].sum()
    df_freq = df_freq[["lang", "word", "freq"]]
    df_freq = df_freq[df_freq["word"].isin(words)]
    df_freq["corpus"] = file.split("/")[3].replace(".txt", "")
    df_f = pd.concat([df_f, df_freq])

# Combine the frequency and the ambiguity dataframes
df = pd.merge(df, df_f, on=["lang", "word"], how="left")

# Add word length column in characters
df["word_length"] = df["word"].str.len()

# Drop rows with missing or infinite values
df = df.dropna(subset=['word_length', 'bert_polysemy', 'freq'])
df = df[~df['bert_polysemy'].isin([float('inf'), float('-inf')])]
df = df[~df['freq'].isin([float('inf'), float('-inf')])]

# Expand lang to language names
df["language_name"] = df["lang"].replace(lg_codes)

# Drop the first column
df = df.drop(columns=["Unnamed: 0"])

# Write the final dataframe to a CSV file
df.to_csv("data/combined_ambiguity_freq.csv", index=False)