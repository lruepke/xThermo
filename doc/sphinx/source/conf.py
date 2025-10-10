# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import sys, os
sys.path.append(os.path.abspath('_extensions'))
path_xThermal_PythonAPI='../../../install/API/python'
if(os.path.exists(path_xThermal_PythonAPI)):
    APIs = os.listdir(path_xThermal_PythonAPI)
    try:
        sys.path.append(os.path.abspath(os.path.abspath(path_xThermal_PythonAPI)+'/'+APIs[0]))
    except:
        pass
# path of the python package xThermal: used to generate python API doc
# sys.path.append(os.path.abspath('../../../Library/install/API/python/cp3.9'))
# import xThermal
# import sphinx_gallery
# -- Project information -----------------------------------------------------

project = 'xThermal'
latex_name='xThermal'
copyright = "Zhikui Guo, Lars Ruepke"
author = 'Zhikui Guo'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = ''


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.mathjax', 
            'jinja',
            'sphinx.ext.ifconfig',
            'sphinx_inline_tabs',
            "sphinx_copybutton",
            # "breathe",
            # "exhale",
            # 'plot_directive',
            'sphinxcontrib.bibtex',
            'sphinx.ext.autodoc', 
            'sphinx.ext.viewcode',
            'sphinx_gallery.gen_gallery',
            'sphinx_toolbox.collapse'
            ]

# config gallery
from sphinx_gallery.scrapers import matplotlib_scraper
from sphinx_gallery.sorting import ExplicitOrder
from sphinx_gallery.sorting import ExampleTitleSortKey
from sphinx_gallery.sorting import FileNameSortKey
class matplotlib_svg_scraper(object):
    def __repr__(self):
        return self.__class__.__name__
    def __call__(self, *args, **kwargs):
        return matplotlib_scraper(*args, format='svg', bbox_inches='tight', **kwargs)
class matplotlib_jpg_scraper(object):
    def __repr__(self):
        return self.__class__.__name__
    def __call__(self, *args, **kwargs):
        return matplotlib_scraper(*args, format='jpg',dpi=300, **kwargs)
examples_dirs = ['gallery/H2O','gallery/H2ONaCl', 'gallery/NaCl']
gallery_dirs = ['gallery_H2O','gallery_H2ONaCl', 'gallery_NaCl']
image_scrapers = ('matplotlib',)
min_reported_time = 0
sphinx_gallery_conf = {
    'backreferences_dir': 'gen_modules/backreferences',
    'doc_module': ('sphinx_gallery', 'numpy'),
    'reference_url': {
        'sphinx_gallery': None,
    },
    'examples_dirs': examples_dirs,
    'gallery_dirs': gallery_dirs,
    'remove_config_comments':False,
    # 'ignore_pattern': r'*help*\.py',
    # 'subsection_order': ExplicitOrder(examples_dirs),
    'image_scrapers': matplotlib_svg_scraper(), #matplotlib_jpg_scraper(),
    'compress_images': ('images', 'thumbnails'),
    # specify the order of examples to be according to filename
    'within_subsection_order': ExampleTitleSortKey,
    # 'expected_failing_examples': ['../examples/no_output/plot_raise.py',
    #                               '../examples/no_output/plot_syntaxerror.py'],
    'min_reported_time': min_reported_time,
    # 'binder': {'org': 'sphinx-gallery',
    #            'repo': 'sphinx-gallery.github.io',
    #            'branch': 'master',
    #            'binderhub_url': 'https://mybinder.org',
    #            'dependencies': './binder/requirements.txt',
    #            'notebooks_dir': 'notebooks',
    #            'use_jupyter_lab': True,
    #            },
    'show_memory': False,
    'junit': os.path.join('sphinx-gallery', 'junit-results.xml'),
    # capture raw HTML or, if not present, __repr__ of last expression in
    # each code block
    'capture_repr': ('_repr_html_', '__repr__'),
    'matplotlib_animations': True,
    'image_srcset': ["2x"],
    'first_notebook_cell': ("# This cell is added by sphinx-gallery\n"
                            "# It can be customized to whatever you like\n"
                            "%matplotlib inline"),
    'last_notebook_cell': "# This is the last cell",
    
}
# configure exhale and breathe for doxygen

# =======================================
source_encoding = 'utf-8-sig'
source_suffix = '.rst'
master_doc = 'index'
templates_path = ['_templates']
plot_basedir='Cookbooks/python'
# 配置字体
# import matplotlib.font_manager as font_manager
# path = '_static/fonts/Arial.ttf'
# prop = font_manager.FontProperties(fname=path)
# plot_rcparams={'font.family':prop.get_name(),'mathtext.fontset':'cm'}
# plot_html_show_source_link=True

# internationalization
language = 'en'
locale_dirs = ['locale/']
gettext_compact = True
gettext_auto_build=True
# Set smartquotes_action to 'qe' to disable Smart Quotes transform of -- and ---
smartquotes_action = 'qe'
# customize OpenFOAM syntax highlight
from sphinx.highlighting import lexers
from pygments_OpenFOAM.foam import OpenFOAMLexer
lexers['foam'] = OpenFOAMLexer(startinline=True)
# default language to highlight source code
# highlight_language = 'foam'
# default language to highlight source code
highlight_language = 'cpp'
pygments_style = 'xcode' #xcode
bibtex_bibfiles = ['manual.bib'] 
# -- Project configuration ------------------------------------------------

# The version shown at the top of the sidebar
version = '1.0'
# The full version shown in the page title
release = '1.0'
# Make the "Edit on GitHub" button link to the correct branch
# Default to master branch if Azure Pipelines environmental variable BUILD_SOURCEBRANCHNAME is not defined
# github_version = os.getenv("BUILD_SOURCEBRANCHNAME", 'master')


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_baseurl='www.scibyte.cn'
html_theme = 'gdal_rtd'
html_theme_path = ["_themes"]
html_theme_options = {
    'sticky_navigation': False,
    'includehidden': False,
    'logo_only' : False,
    'sticky_navigation': True,
    'titles_only': True,
    'display_version': False,
    'prev_next_buttons_location': 'both',
    'style_nav_header_background': '#02007e',
    # 'gitlab_url': 'https://gitlab.com/gmdpapers/hydrothermalfoam'
}

html_context = {
    "menu_links": [
        (
            '<i class="fa fa-stamp"></i> Imprint',
            "https://www.sweos.info/imprint/",
        ),
        (
            '<i class="fa fa-home fa-fw"></i> Homepage',
            "https://hydrothermal-openfoam.gitlab.io/saltwatereos/",
        ),
        (
            '<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4603878.svg">',
            'https://zenodo.org/record/4603878/export/hx#.YRPazi8RrSw'
        ),
        (
            '<i class="fa fa-book fa-fw"></i> License:GPL',
            "https://github.com/zguoch/saltwatereos/blob/master/LICENSE",
        ),
        (
            '<i class="fa fa-envelope fa-fw"></i> Contact',
            "mailto:zguo@geomar.de",
        ),
        (
            '<i class="fa fa-gitlab fa-fw"></i> Bug report',
            "https://gitlab.com/hydrothermal-openfoam/saltwatereos/-/issues",
        ),
        (
            '<i class="fa fa-github fa-fw"></i> Source Code',
            "https://github.com/zguoch/saltwatereos",
        ),
        (
            '<i class="fa fa-github fa-fw"></i> Code documentation',
            "https://www.sweos.info/doxygen/",
        ),
    ],
    'project':project,
    'downloads_url':'hydrothermal-openfoam.gitlab.io/saltwatereos/manual/downloads',
    'latex_main':  latex_name, 
    'pdf_versions': [
        (
            'latest',
            'https://hydrothermal-openfoam.gitlab.io/saltwatereos/downloads'
        ),
        (
            '2.0',
            '#'
        ),
        ]
}

# favicon of the docs
html_favicon = "_static/logo.png"
html_static_path = ['_static']
exclude_patterns = ['doxygen']
html_last_updated_fmt = '%b %d, %Y'
# If true, links to the reST sources are added to the pages.
html_logo = "_static/logo_beta.png" 
html_show_sourcelink = False 
# List of custom CSS files (needs sphinx>=1.8)
html_css_files = ["style.css"]

# Redefine supported_image_types for the HTML builder
from sphinx.builders.html import StandaloneHTMLBuilder
StandaloneHTMLBuilder.supported_image_types = [
  'image/gif', 'image/jpeg', 'image/png', 'image/svg+xml'
]


# -- Options for LaTeX output ---------------------------------------------
latex_engine = 'xelatex'
# latex_logo='_static/latex_logo.pdf'
# 设置公式和图片编号依赖于章节
numfig = True
math_numfig = True
math_eqref_format = '({number})'
# 只对make latex有效
# numfig_format = 'Figure. %s'
numfig_secnum_depth = 1
imgmath_latex = 'dvilualatex'
imgmath_image_format = 'svg'
imgmath_dvipng_args = ['-gamma', '1.5', '-bg', 'Transparent']
# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_name='xThermal'
latex_documents = [
    (master_doc, latex_name+'.tex', '\\textcolor{x_color}{x}Thermo: Thermodynamic properties and EOS of \\textcolor{x_color}{water} and \\textcolor{x_color}{saltwater}',
     '\href{mailto:zguo@geomar.de}{Zhikui Guo}, \href{mailto:lruepke@geomar.de}{Lars Rüpke}\\\\ \href{mailto:joerg.hasenclever@uni-hamburg.de}{Jörg Hasenclever}, \href{mailto:falko.vehling@ifg.uni-kiel.de}{Vehling Falko}', 'manual'),
]
latex_toplevel_sectioning="part" 
latex_elements = {
    'papersize': 'a4paper',
    'utf8extra': '',
    'inputenc': '',
    'babel': r'''\usepackage[english]{babel}''',
    'preamble': r'''\usepackage{ctex} 
\definecolor{EOS_green}{rgb}{0.27,0.49,0.36}
\definecolor{x_color}{rgb}{1, 0, 0}
\newcommand{\xThermal}{\textsc{\textcolor{x_color}{x}}Thermo}
\makeatletter
\let\newauthor\@author
\let\newrelease\py@release
\let\newdate\@date
\makeatother

\renewcommand{\sphinxmaketitle}{
\begin{titlepage}
    \begin{center}
    \resizebox{\textwidth}{!}{
        \colorbox{EOS_green}{
        \fontfamily{\rmdefault}\selectfont \textcolor{white}{\hspace{0.1in}\xThermal{}\hspace{0.1in}}
        }
    }
    \\[12pt]
    {\Large Thermodynamic properties and EOS of \textcolor{x_color}{water} and \textcolor{x_color}{saltwater}}
    \end{center}
    \vfill\vfill
    \begin{center}%
    \Large 
    \href{mailto:zguo@geomar.de}{Zhikui Guo}, \href{mailto:lruepke@geomar.de}{Lars Rüpke}\\ \href{mailto:joerg.hasenclever@uni-hamburg.de}{Jörg Hasenclever}, \href{mailto:falko.vehling@ifg.uni-kiel.de}{Vehling Falko} \\
    \newrelease{} \\
    \newdate{}
    \end{center}

\end{titlepage}
\setcounter{footnote}{0}
\clearpage
}
    ''',
}

# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'xThermal', 'Technical note and manual',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'xThermal', 'Thermodynamic properties and EOS of fluids and solids involved in geo-modeling',
     author, 'Zhikui Guo', 'Thermodynamic properties and EOS of fluids and solids involved in geo-modeling.',
     'Miscellaneous'),
]

def setup(app):
    app.add_css_file("style.css")
    # app.add_js_file("js/custom.js")
    app.add_js_file("https://cdn.jsdelivr.net/npm/clipboard@1/dist/clipboard.min.js")
    app.add_js_file("js/echarts/echarts.js")
    
# new defined cite style
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.plugin import register_plugin
from collections import Counter
import re
import unicodedata

from pybtex.style.labels import BaseLabelStyle

_nonalnum_pattern = re.compile('[^A-Za-z0-9 \-]+', re.UNICODE)

def _strip_accents(s):
    return "".join(
        (c for c in unicodedata.normalize('NFD', s)
            if not unicodedata.combining(c)))

def _strip_nonalnum(parts):
    """Strip all non-alphanumerical characters from a list of strings.

    >>> print(_strip_nonalnum([u"ÅA. B. Testing 12+}[.@~_", u" 3%"]))
    AABTesting123
    """
    s = "".join(parts)
    return _nonalnum_pattern.sub("", _strip_accents(s))

class APALabelStyle(BaseLabelStyle):
    def format_labels(self, sorted_entries):
        labels = [self.format_label(entry) for entry in sorted_entries]
        count = Counter(labels)
        counted = Counter()
        for label in labels:
            if count[label] == 1:
                yield label
            else:
                yield label + chr(ord('a') + counted[label])
                counted.update([label])

    def format_label(self, entry):
        label = "Anonymous"
        if 'author' in entry.persons:
            label = self.format_author_or_editor_names(entry.persons['author'])
        elif 'editor' in entry.persons:
            label = self.format_author_or_editor_names(entry.persons['editor'])
        elif 'organization' in entry.fields:
            label = entry.fields['organization']
            if label.startswith("The "):
                label = label[4:]

        if 'year' in entry.fields:
            return "{}, {}".format(label, entry.fields['year'])
        else:
            return "{}, n.d.".format(label)

    def format_author_or_editor_names(self, persons):
        if (len(persons) == 1):
            return _strip_nonalnum(persons[0].last_names)
        elif (len(persons) == 2):
            return "{} & {}".format(
                _strip_nonalnum(persons[0].last_names),
                _strip_nonalnum(persons[1].last_names))
        else:
            return "{} et al.".format(
                _strip_nonalnum(persons[0].last_names))

class APAStyle(UnsrtStyle):

    default_label_style = APALabelStyle

register_plugin('pybtex.style.formatting', 'apa', APAStyle)