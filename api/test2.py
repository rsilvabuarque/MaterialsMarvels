# import papermill as pm

# pm.execute_notebook(
#    'testoutput.ipynb',
#    'temp.ipynb',
# #    parameters=dict(alpha=0.6, ratio=0.1)
# )

import nbformat
from nbconvert import HTMLExporter
from nbconvert.preprocessors import ExecutePreprocessor

notebook = nbformat.read("testoutput.ipynb", as_version=4)

ep = ExecutePreprocessor(timeout=600, kernel_name='python3')

ep.preprocess(notebook, {'metadata': {'path': '.'}})

with open('executed_notebook.ipynb', 'w', encoding='utf-8') as f:
   nbformat.write(notebook, f)

# html_exporter = HTMLExporter(template_name="classic")

# (body, resources) = html_exporter.from_notebook_node(notebook)

# with open("index3.html", "w") as f:
#     f.write(body)