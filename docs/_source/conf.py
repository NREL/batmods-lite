# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'BATMODS-lite'
copyright = '2023, Corey R. Randall'
author = 'Corey R. Randall'
release = '0.0.1'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.todo',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon', 
              'sphinx.ext.autodoc',
              'autoapi.extension', 
              'myst_nb', 
              'sphinx_design',
              'sphinx_favicon',
              'sphinx_copybutton'
              ]

templates_path = ['_templates']

exclude_patterns = ['build', 
                    'Thumbs.db', 
                    '.DS_Store', 
                    '*.ipynb_checkpoints',
                    '__pycache__'
                    ]

source_suffix = {'.rst': 'restructuredtext',
                 '.ipynb': 'myst-nb',
                 '.myst': 'myst-nb'
                 }


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
# https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/layout.html

html_theme = 'pydata_sphinx_theme'

html_favicon = '_static/favicon.ico'

html_static_path = ['_static']
html_js_files = ['custom.js']
html_css_files = ['custom.css']

html_context = {'default_mode': 'light'}

html_sidebars = {'**': ['search-field.html', 'sidebar-nav-bs.html']}

html_theme_options = {'logo': {'image_light': '_static/light_notag.svg',
                               'image_dark': '_static/dark_notag.svg'},
                      'header_links_before_dropdown': 6, 
                      'navbar_align': 'content',
                      'footer_start': ['copyright.html'],
                      'footer_end': ['sphinx-version.html'],
                      'navbar_persistent': [],
                      'secondary_sidebar_items': ['page-toc.html'],
                      'search_bar_text': 'Search...',
                      'show_prev_next': False,
                      'collapse_navigation': True,
                      'show_toc_level': 0,
                      'pygment_light_style': 'tango'
                      }


# -- Options for napoleon ----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html

napoleon_use_rtype = False
napoleon_custom_sections = ['BATMODS-lite', 
                            'How to use the documentation', 
                            'Viewing documentation using IPython',
                            'Material Properties Package',
                            'Math Module',
                            'Mesh Module',
                            'Plotting Utilities',
                            'Single Particle Model Package',
                            'Pseudo-2D Model Package',
                            'Domains Module',
                            'DAE Module',
                            'Root Functions',
                            'Post-processing Utilities',
                            'Solution Classes'
                            ]


# -- Options for autoapi -----------------------------------------------------
# https://sphinx-autoapi.readthedocs.io/en/latest/reference/config.html

autoapi_type = 'python'
autoapi_ignore = ['*/__pycache__/*']
autoapi_dirs = ['../../bmlite']
autoapi_keep_files = True 
autoapi_root = 'api'
autoapi_member_order = 'groupwise'
autoapi_python_class_content = 'both'
autodoc_typehints = 'none'
autoapi_options = ['members', 'inerited-members', 'undoc-members',
                   'show-module-summary', 'imported-members']


# -- Options for myst --------------------------------------------------------
# https://myst-nb.readthedocs.io/en/latest/configuration.html

nb_execution_timeout = 300
nb_number_source_lines = True
myst_enable_extensions = ['amsmath', 'dollarmath']