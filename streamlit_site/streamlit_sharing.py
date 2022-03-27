import nGRE_parse
import custom_element_parse
import streamlit as st

PAGES = {
    "nGRE parsing": nGRE_parse,
    "Custom element parsing": custom_element_parse
}
st.sidebar.title('Navigation')
selection = st.sidebar.radio("Go to", list(PAGES.keys()))
page = PAGES[selection]
page.main()
