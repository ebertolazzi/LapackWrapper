# -*- coding: utf-8 -*-
import os
from pathlib import Path

# pip install breathe
# pip install exhale
# pip install sphinxcontrib-email
# pip install cloud_sptheme

# -- Project information -----------------------------------------------------
exec(open("../project_common.py").read())

extensions.append('breathe');
extensions.append('exhale');

breathe_projects = {
  "doc_cpp"    : "_doxygen/doc_cpp/xml-cpp",
  "doc_sparse" : "_doxygen/doc_sparse/xml-sparse",
}

breathe_default_project = "doc_sparse"

base_dir = os.path.dirname(os.path.realpath(__file__));

dir_path_cpp = base_dir+"../../../src"
dir_path_cpp = Path(dir_path_cpp).resolve()

dir_path_sparse = base_dir+"../../../src/sparse_tool"
dir_path_sparse = Path(dir_path_sparse).resolve()

dir_path_matlab = base_dir+"../../../toolbox/lib"
dir_path_matlab = Path(dir_path_matlab).resolve()

doxy_filter = "        INPUT_FILTER        = "+str(Path(base_dir+"/doxygen_filter.rb").resolve());

doxygen_common_stdin = """
        EXTRACT_ALL         = YES
        SOURCE_BROWSER      = YES
        EXTRACT_STATIC      = YES
        HIDE_SCOPE_NAMES    = NO
        CALLER_GRAPH        = YES
        GRAPHICAL_HIERARCHY = YES
        HAVE_DOT            = YES
        QUIET               = NO
        GENERATE_TREEVIEW   = YES
        SHORT_NAMES         = YES
        IMAGE_PATH          = ../images
        LATEX_HEADER        = ../templates/latex_headers.tex

        XML_PROGRAMLISTING    = YES
        RECURSIVE             = YES
        FULL_PATH_NAMES       = YES
        ENABLE_PREPROCESSING  = YES
        MACRO_EXPANSION       = YES
        SKIP_FUNCTION_MACROS  = NO
        EXPAND_ONLY_PREDEF    = NO
        INHERIT_DOCS          = YES
        INLINE_INHERITED_MEMB = YES
        EXTRACT_PRIVATE       = YES
        PREDEFINED           += protected=private
        GENERATE_HTML         = NO
"""

doc_cpp = {
    'verboseBuild':          True,
    "rootFileName":          "root.rst",
    "createTreeView":        True,
    "exhaleExecutesDoxygen": True,
    "doxygenStripFromPath":  str(dir_path_cpp),
    "exhaleDoxygenStdin":    '''
        INPUT               = ../../src/lapack_wrapper
        PREDEFINED         += protected=private
        XML_OUTPUT          = xml-cpp
'''+doxy_filter+doxygen_common_stdin,
    "containmentFolder":    os.path.realpath('./api-cpp'),
    "rootFileTitle":        "C++ API",
}

doc_sparse = {
    'verboseBuild':          True,
    "rootFileName":          "root.rst",
    "createTreeView":        True,
    "exhaleExecutesDoxygen": True,
    "doxygenStripFromPath":  str(dir_path_sparse),
    "exhaleDoxygenStdin":    '''
        INPUT               = ../../src/sparse_tool
        PREDEFINED         += protected=private
        XML_OUTPUT          = xml-sparse
'''+doxygen_common_stdin,
    "containmentFolder":    os.path.realpath('./api-sparse'),
    "rootFileTitle":        "SparseTool API",
}

exhale_projects_args = {
  "doc_cpp": doc_cpp,
  "doc_sparse" : doc_sparse
}

cpp_index_common_prefix = [
  'lapack_wrapper::',
  'SparseTool::'
]

# If false, no module index is generated.
html_domain_indices = True

# If false, no index is generated.
html_use_index = True
