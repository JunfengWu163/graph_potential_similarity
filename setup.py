from setuptools import setup, find_packages

setup(
    name='graph_potential_similarity',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        # List your package dependencies here
        'numpy',
        'scipy',
        'rdkit'
    ],
    author='Junfeng Wu',
    author_email='wujunfeng@vip.163.com',
    description='Graph Potential Similarity is a Python package designed to compute the similarity between two graphs based on their potential functions. This tool is particularly useful for a wide range of applications, including network analysis, pattern recognition, genome sequence analysis, chemical compound structure analysis, and more.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/JunfengWu163/graph_potential_similarity',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
