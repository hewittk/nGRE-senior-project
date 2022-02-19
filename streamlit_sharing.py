import streamlit as st
import pandas as pd

st.write("My First Streamlit Web App")

df = pd.DataFrame({"one": [1, 2, 3], "two": [4, 5, 6], "three": [7, 8, 9]})
st.write(df)

maximum_mutations = st.selectbox(
     'What is the maximum number of mutations in the nGRE consensus sequence that you want to tolerate?',
     ('0', '1', '2', '3'))
print("Maximum mutations in potential nGREs: ", maximum_mutations)


user_input = st.text_area("Insert DNA sequence to parse for nGREs")
user_input
