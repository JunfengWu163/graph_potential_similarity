# graph_potential_similarity
 Graph Potential Similarity is a Python package designed to compute the similarity between two graphs based on their potential functions. This tool is particularly useful for a wide range of applications, including network analysis, pattern recognition, genome sequence analysis, chemical compound structure analysis, and more.

## Installation

You can install the package using pip:

```bash
pip install graph-potential-similarity
```

## Usage
Here's a simple example of how to use the Graph Potential Similarity package:

```python
import graph_potential_similarity as gps

# an example of computing distance from SMILES strings
sim1 = gps.smiles_distance('FC(Cl)(Cl)C(F)(F)OC(F)(F)Cl', 'FC(F)(F)C(Cl)(Cl)OC(F)(F)Cl')
print(f"The similarity between the two chemical compounds is: {sim1}")

# another example, compute distance from DNA sequences
seq1 = "GGACCGACAGGAATTCGCTCCTTAGGACGTTATAGTTACGGCCGCCGTTTACTGGGGCTTCAATTCGCAGCTTCGC"
seq2 = "AGCGAACGCTGGCGGCATGCTTAACACATGCAAGTCGCACGAAGGCTTCGGCCTTAGTGGCGGACGGGTGAGTAAC"
sim2 = gps.double_helix_distance(seq1,seq2)
print(f"The similarity between the two DNA sequences is: {sim2}")
```

## Contributing
We welcome contributions! Please see our CONTRIBUTING.md for more details.

## License
This project is licensed under the GNU General Public License v3 (GPLv3) License - see the LICENSE file for details.

## Acknowledgements
This project was made possible by the following resources and contributions:
- **He, J., Chen, J., Huang, G., Cao, J., Zhang, Z., Zheng, H., ... & Van Zundert, A. (2021).** "A polynomial‐time algorithm for simple undirected graph isomorphism." *Concurrency and Computation: Practice and Experience*, 33(7), 1-1. This paper provided foundational insights and algorithms that were instrumental in the development of this package.
- **RDKit**: An open-source cheminformatics software that was used for chemical compound structure analysis.
- **NumPy**: A fundamental package for scientific computing with Python, which was used extensively for numerical operations.

We are grateful for the contributions of these resources and the support of the open-source community.

## Contact
If you have any questions or need support, please contact Junfeng Wu (Email: wujunfeng@vip.163.com).