#!/usr/bin/env python3

# field_tagger.py | This script adds a YAML front matter with the tag "field:"
# to the markdown files in the specified directory and subdirectories. The 
# value of the tag is the name of the directory where the file is located.

corresponding_tags = {
    "Bioinformática": "bioinformatics",
    "Biología": "biology",
    "Cáncer": "biology",
    "Economía": "economics",
    "Enfermedades": "biology",
    "Estadística": "statistics"}

import os

def parse_frontmatter