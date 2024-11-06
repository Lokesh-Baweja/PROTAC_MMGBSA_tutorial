import nbformat
from nbconvert import MarkdownExporter

# Load your notebook
with open("PROTAC_MMGBSA_TUTORIAL.ipynb", "r") as f:
    notebook_content = nbformat.read(f, as_version=4)

# Convert to markdown
md_exporter = MarkdownExporter()
markdown, resources = md_exporter.from_notebook_node(notebook_content)

# Save as README.md
with open("README.md", "w") as f:
    f.write(markdown)
