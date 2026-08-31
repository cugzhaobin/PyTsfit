#!/usr/bin/env python
# ----------------------------------------------------------
# Console entry point for the Streamlit UI.
#
#     pytsfit-web [streamlit options]
#
# Resolves the bundled app.py (wherever the package is installed,
# site-packages or editable) and hands it to Streamlit's own CLI,
# so users never have to type the path.
# ----------------------------------------------------------
def main():
    import sys
    from pathlib import Path

    from streamlit.web import cli as stcli

    app = Path(__file__).with_name('app.py')
    sys.argv = ['streamlit', 'run', str(app), *sys.argv[1:]]
    sys.exit(stcli.main())


if __name__ == '__main__':
    main()
