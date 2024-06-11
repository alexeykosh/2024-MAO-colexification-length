import pandas as pd
import zipfile
import os

# Unzip the lexibank-analysed-v1.0.zip file in the data folder
with zipfile.ZipFile('data/lexibank-analysed-v1.0.zip', 'r') as zip_ref:
    zip_ref.extractall('data/')

# Unzip the forms.csv file in the cldf folder
with zipfile.ZipFile('data/lexibank-lexibank-analysed-a4c0952/cldf/forms.csv.zip', 'r') as zip_ref:
    zip_ref.extractall('data/lexibank-lexibank-analysed-a4c0952/cldf/')

# Load the forms, languages, and concepts csv files
forms = pd.read_csv('data/lexibank-lexibank-analysed-a4c0952/cldf/forms.csv')
languages = pd.read_csv('data/lexibank-lexibank-analysed-a4c0952/cldf/languages.csv')
concepts = pd.read_csv('data/lexibank-lexibank-analysed-a4c0952/cldf/concepts.csv')

# match forms with languages by Language_ID and ID
forms = forms.merge(languages, left_on='Language_ID', right_on='ID')
# match concepts by Parameter_ID and ID
forms = forms.merge(concepts, left_on='Parameter_ID', right_on='ID')
forms['Length'] = forms['Segments']\
    .str.strip()\
    .str.count(' ') + 1
# convert Glottocode column to string
forms['Glottocode'] = forms['Glottocode']\
    .astype(str)
# count the number of meanings per form
forms['Number_of_meanings'] = forms\
    .groupby(['Form', 'Language_ID'])['Concepticon_ID']\
    .transform('nunique')

print('Number of forms:', forms.shape[0])
forms.to_csv('data/forms_total.csv', index=False)
