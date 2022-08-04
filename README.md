## Contrast Subgraphs Allow Comparing Homogeneous and Heterogeneous Networks Derived from Omics Data

[Tommaso Lanciano](https://phd.uniroma1.it/web/LANCIANO-TOMMASO_nP1661409_EN.aspx) (Sapienza University, Rome), [Aurora Savino](https://humantechnopole.it/it/people/aurora-savino/)  (University of Turin, Turin), [Francesca Porcu](https://www.linkedin.com/in/francesca-porcu-189a13168/?originalSubdomain=it)  (Sapienza University, Rome), [Davide Cittaro](https://research.hsr.it/en/centers/omics-sciences/innovation-lab/davide-cittaro.html)  (San Raffaele Scientific Institute IRCSS, Milan), [Francesco Bonchi](http://www.francescobonchi.com/) (Centai, Turin and Eurecat, Barcelona) and [Paolo Provero](https://biotec.campusnet.unito.it/do/docenti.pl/Alias?paolo.provero#tab-profilo) (University of Turin, Turin and San Raffaele Scientific Institute IRCSS, Milan).

_Biological networks are often used to describe the relationships between relevant entities, in particular genes and proteins, and are a powerful tool for functional genomics. Many important biological problems can be investigated by comparing biological networks between different conditions, or networks obtained with different techniques. We show that contrast subgraphs, a recently introduced technique to identify the most important structural differences between two networks, provide a versatile tool for comparing gene and protein networks of diverse origin. We show in three concrete examples how contrast subgraphs can provide new insight in functional genomics by extracting the gene/protein modules whose connectivity is most altered between two conditions or experimental techniques._


<p align="center">
  <img width="700" height="450" src="https://github.com/tlancian/bio_cs/blob/main/fig.png">
</p>


---

This repository contains the code necessary to implement algorithms described in the article "Contrast Subgraphs Allow Comparing Homogeneous and Heterogeneous Networks Derived from Omics Data".

This package is free for research, academic and non-profit making purposes only. If you use this piece of software for your work and got something published please include the citation reported below. The software may not be sold or redistributed without prior approval. One may make copies of the software for their use provided that the copies, are not sold or distributed, are used under the same terms and conditions. As unestablished research software, this code is provided on an "as is" basis without warranty of any kind, either expressed or implied. The downloading, or executing any part of this software constitutes an implicit agreement to these terms. These terms and conditions are subject to change at any time without prior notice.



### Requirements

* The code has been tested with Python 3.8
* Input graph must be a weighted tab-separated edge-list, following the format: "node1 \t node2 \t weight"


### Execution

Run the following command: 'python main.py [-h] [-C] d'

#### Positional arguments:
  * d           &nbsp;&nbsp;&nbsp;&nbsp;Path to input graph.

#### Optional arguments:
  * -h, --help  
    show the help message and exit	
  * -C       &nbsp;&nbsp;&nbsp;&nbsp;C value for execution of Algorithm 1
  	
#### Examples:
'python main.py toy_graph.tsv 1.4'  


---
  
### Contacts
Mail to [paolo.provero@unito.it](mailto:paolo.provero@unito.it) and [bonchi@centai.eu](mailto:bonchi@centai.eu) for any question.
