'''
pytsfit.ui: a Streamlit front-end for PyTsfit.

Thin shell over the existing pytsfit API. No fitting/plotting logic is
re-implemented here; every function delegates to the core modules
(data / models / tsfitting / output / PyTsfit).

    streamlit run src/pytsfit/ui/app.py
'''
