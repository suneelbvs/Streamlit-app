import streamlit as st
import leafmap.foliumap as leafmap

st.set_page_config(layout="wide")

# Customize the sidebar
markdown = """
This web app is maintained by [Suneel Kumar BVS](https://www.linkedin.com/in/suneelbommisetty). You can follow me on social media: [GitHub](https://github.com/suneelbvs) | [Twitter](https://twitter.com/suneel_bvs) | [WebSite](sbvs.me).
Streamlit template design inspired from [Qiusheng Wu](https://wetlands.io)
"""

st.sidebar.title("About")
st.sidebar.info(markdown)

# Customize page title
#st.title("Streamlit for Geospatial Applications")

st.markdown(
    """
    This my first multipage  web app template created using [streamlit](https://streamlit.io) and [leafmap](https://leafmap.org).
    """
)

st.header("About this streamlit")

markdown = """
1. rdkit cheatsheet: help you in your day to day cheminformatics tasks. I'm still working on the improvement of the current version. Pls DM your suggestions to [My Twitter handle](https://twitter.com/suneel_bvs)

2. rdkit cheatsheet: I'm working on the jupyter notebook to demonstrate this cheatcodes; will keep you update in couple of weeks. This notebooks will help you to understand the cheatsheet functions better

3. MolGenerator: Its SMILES TO 2D & 3D converter and its basic molecular properties.

4. More Pages - work in progress!
"""

st.markdown(markdown)