import nGRE_parse
import custom_element_parse
import streamlit as st
import about_tool
from PIL import Image

PAGES = {
    "About tool": about_tool,
    "nGRE parsing": nGRE_parse,
    "Custom element parsing": custom_element_parse
}
st.sidebar.title('Navigation')
selection = st.sidebar.radio("Go to", list(PAGES.keys()))
page = PAGES[selection]
page.main()
