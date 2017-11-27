---
layout: post
title:  "DNA Visualization"
date:   2017-09-27 12:00:00 -0500
categories: jekyll update
---

## Outline
### 1. Introduction
### 2. The Codeword Design Problem
    A. The Gibbs energy
    B. Codeword Design in DNA Space

### 3. DNA Visualization
    A. Case of Study
    B. Methodology

### 4. Experimentation
    A. Dataset
    B. Change DNA strands to images
    C. Comparison Between species

### 5. Conclusions and Future work

### References

___

<h3><center>1. Introduction</center></h3>

DNA-based applications such as self-assembly ([Garzon & Phan, 2009](#Garzon & Phan, 2009); [Qian & Winfree, 2011](#Qian & Winfree, 2011); [Seeman, 2006](#Seeman, 2006)), Natural Language Processing ([Neel & Garzon, 2006](#Neel & Garzon, 2006); [Bobba et al., 2006](#Bobba et al., 2006)) and DNA-based memories ([Neel & Garzon, 2003](#Neel & Garzon, 2003)); has stimulated the development and analysis of various techniques aimed at finding large sets of short oligonucleotides with noncrosshybridizing (nxh) properties, namely the so-called **Codeword Design problem**. The problem calls for finding large sets of single DNA strands that do noncrosshybridize to themselves and/or to their complements. Self-assembled nanostructures is usually directly related to the largest ensemble of DNA molecules that satisfy a given set of crosshybridization and noncrosshybrization constraints.

([Garzon, Bobba, & Hyde, 2004](#Garzon, Bobba, & Hyde, 2004)), proposed an approximation to encode Abiotic information in DNA Spaces take into account the crosshybridization and noncrosshybrization properties of the DNA spaces, where the signature(Numerical representation of a sequence) of a string $x$ is represented as a vector $V$ based on how $x$ hybridize in some reaction stringency threshold $\tau$ with the noncrosshybridizing set. Then this vector could be visualized as an 1D or 2D image. This approximation has shown good results for species identification ([Garzon, & Wong, 2011](#Garzon, & Wong, 2011)) and microarray design ([Garzon, & Mainali, 2017](#Garzon, & Wong, 2011)).

<h3><center>2. The Codeword Design Problem</center></h3>

### A. The Gibbs energy
DNA duplex formation of two single strands is determined by the familiar Gibbs energy. which generally depends on physical parameters such as the internal energy ($U$), pressure ($p$), volume ($V$), temperature ($T$), and entropy ($S$) of the environment in which the duplex is formed.

The threshold for duplex formation of single strands can be considered in $−6Kcal/mole$; the more negative the energy, the more stable the duplex formed. Strands above this threshold are called noncrosshybridizing (nxh). Knowledge of the corresponding Gibbs energy landscapes would appear critical to a solution of the
Codeword Design problem.

A full use of even Gibbs energy model forces the examination of exponentially many configurations in the minimization of the energy and thus it is computationally intractable, especially when it comes to solving the Codeword Design problem, where in principle is necessary to deal with $2^{4^n}$ possible subsets of the entire DNE space to search for the code set of the largest size.

For this reason, is necessary to use some approximation to the Gibbs Energy. A measure of hybridization affinity introduced by ([Garzon et al., 1997](#Garzon et al., 1997)), the _**h-distance**_, which provides a computationally efficient approximation of the Gibbs energy based solely on composition and sequence. The hybridization is modeled as a binary operation between two single strands $x$, $y$ that forms a double-stranded molecule in which pairs of nucleotides are joined in a duplex molecule by two or three hydrogen bonds, for an $a$−$t$ or a $c$−$g$ pair, respectively. The _**h-distance**_ between two single strands x, y is defined as follows:

<h4><center>
$$
h(x,y) = \displaystyle {\min_{-n < k < n} \{ |k| + H(x,  \sigma^k(y') ) \} }
$$
</center></h4>

where $\sigma^k(y')$ is the shift of $y'$ by $k$ positions from a perfect alignment with $x$ (right-shift if $k > 0$; left-shift if $k < 0$), $y′$ is the Watson-Crick complement of $y$, and the Hamming distance $H$ measures the number of mismatched base pairs in the overlap of $x$ and $y′$ in the specified frame shift $\sigma^k(y')$.

or example, if:

<h4><center>
	$ x = agc, y = tgg $ (and so $y' = cca$)
</center></h4>

- at shift $ k = -2, \begin{smallmatrix} \hphantom{--} agc \\  cca \end{smallmatrix}$  the distance is: $ 2 + H(a,a) = 2$
- at shift $ k = -1, \begin{smallmatrix} \hphantom{-} agc \\   cca \end{smallmatrix}$  &nbsp;&nbsp; the distance is:  $ 1 + H(ag,ca) = 3$
- at shift $ k = \hphantom{-} 0, \begin{smallmatrix} agc \\   cca \end{smallmatrix}$  &nbsp;&nbsp;&nbsp;&nbsp; the distance is:  $ 0 + H(agc,cca) = 3$
- at shift $ k = \hphantom{-} 1, \begin{smallmatrix} agc  \\  cca \end{smallmatrix}$  &nbsp;&nbsp;&nbsp;&nbsp; the distance is:  $ 1 + H(gc,cc) = 2$
- at shift $ k = \hphantom{-} 2, \begin{smallmatrix} agc  \\  cca \end{smallmatrix}$  &nbsp;&nbsp;&nbsp;&nbsp; the distance is: $ 2 + H(c,c) = 2$

Thus:
<center>
    	$h(agc, tgg) = 2$
</center>

Therefore, the _**h-distance**_ quantifies hybridization affinity closely while preserving the advantages of a metric space structure in DNA space. Then, the work to visualize DNA sequences turn to the analysis under this approximation of the Gibbs energy.


```python
# Calculate the h-distance of two strands
# @param: x DNA strand
# @param: y DNA strand
# @return: The h-distance between x and y
def calculateH_distance(x, y):
    # Calculate all right and left shifts of y
    leftY = ["-" for _ in range(len(x))] + list(y)
    rightY = list(y) + ["-" for _ in range(len(x))]

    # Fix x to Calculate the h-distance based on the shifts
    x_fixed = list(x) + ["-" for _ in range(len(y))]
    x_fixed_1 = ["-" for _ in range(len(y))]+ list(x)

    # Calculate the complementary of a sequence
    def check_complementary(x_i, y_i):
        if x_i=='a' and y_i=='t':
            return True
        elif x_i=='t' and y_i=='a':
            return True
        elif x_i=='g' and y_i=='c':
            return True
        elif x_i=='c' and y_i=='g':
            return True
        else:
            return False

    # all possible shifts
    setPossibleShifts = []

    # Left shift
    setPossibleShifts.append((''.join(leftY), len(x) - sum(check_complementary(ch1, ch2) for ch1, ch2 in zip(x_fixed, leftY))))
    while leftY[0] == '-':
        tmp = [leftY.pop(0)]
        leftY = leftY + tmp
        setPossibleShifts.append((''.join(leftY), len(x) - sum(check_complementary(ch1, ch2) for ch1, ch2 in zip(x_fixed, leftY))))


    # Right shift
    setPossibleShifts.append((''.join(rightY), len(x) - sum(check_complementary(ch1, ch2) for ch1, ch2 in zip(x_fixed_1, rightY))))
    while rightY[len(rightY)-1] == '-':
        tmp = [rightY.pop()]
        rightY =  tmp + rightY
        setPossibleShifts.append((''.join(rightY), len(x) - sum(check_complementary(ch1, ch2) for ch1, ch2 in zip(x_fixed_1, rightY))))

    # Return the min h-distance of x and y
    return min(setPossibleShifts, key = lambda t: t[1])[1]

# Calculate the h-distance of x and y' | x' and y_complement
# @param: x DNA strand
# @param: y DNA strand
# @return: The h-distance between the Poligo X and Y
def realHdistance(x, y, n):
    DNA_prime = complementaryReverse(y)
    h_distance = min(calculateH_distance(x, DNA_prime), calculateH_distance(x, y[::-1]))

    # Normalize between 0 and 1
    return float(h_distance - 0) / float(n-0)

# Calculate the complementary of a sequence
# @param: x DNA strand
# @return: The watson crick complementary reverse of x
def complementaryReverse(x):
    def check_complementary(x_i):
        if x_i=='a':
            return 't'
        elif x_i=='t':
            return 'a'
        elif x_i=='g':
            return 'c'
        else:
            return 'g'

    reverse = []
    for xi in list(x):
        reverse.append(check_complementary(xi))

    return ''.join(reverse[::-1])
```



### B. Codeword Design in DNA Space

Solution to the codeword design problem would be enormously facilitated by adequate knowledge of the geometry in the metric spaces $D_n$ . The metric space $D_n$ has the following properties for all $n ≥ 1$ ([Phan & Garzon, 2009](#Phan & Garzon, 2009)):

1. There are $|P| = 4^{n/2}$ $n$-mers consisting of a single palindromic DNA strand that are their own reverse complements (i.e., $|X| = 1$) for $n$ even, and $0$ for $n$ odd.
2. There are $ |D_n | = \frac{4^n - |P|}{2} $, nonpalindromic $n$-mers.
3. There are $ |D_n | = \frac{4^n + |P|}{2} $, $n$-mers in total.

The precise formulation of Codeword Design can now be given as follows in terms of the _**h-distance**_. Given a set of $n$-mers $S$, i.e., $D_n$, possible with a lot of crosshybridization, and a reaction stringency threshold $ \tau $, find a subset of $S$ with no crosshybridization, i.e., find the largest ($n,\tau$)-code in $S$.

___
#### Codeword Design in DNA Spaces
_**Input**_: A set $S$ of $n$-mers, a threshold $\tau$ and an integer $K$;

_**Output**_: Is there an ($n, \tau $)-code subset of $S$ of cardinality at least $K$, i.e., where every two distinct words are at a distance at least $\tau$ from each other?
___

The search for noncrosshybridizing sets(The Codeword Design Problem) was done using SaEA(**Self-adaptive Evolutionary Algorithm**) ([Prieto, Leon, &; Garzon](#Prieto, Leon, &; Garzon)). SaEA uses an Evolutionary Algorithm which allows to use a special encoding for DNA spaces and DNA genetic operators that take advantage of the DNA space structure. The Hybrid Adaptive Evolutionary Algorithm (HAEA) proposed by ([Gomez, 2004](#Gomez, 2004)), is used by SaEA as a parameter adaptation technique that automatically learns the rates of its genetic operators while the individuals are evolved in an Evolutionary Algorithm.

The noncrosshybridizing sets of $n$-mers (with $n = 4, 6, 8$ and $\tau=n/2$) found with SaEA are quantified by two standard measures in information theory. One is the expected value of hybridizations of a random pmer; another is the Shannon Entropy of the corresponding distribution (i.e., the degree of uncertainty with which a pmer is attached to a noncrosshybridizing pmer in the basis). Under ideal Codeword Design problem, the Shannon Entropy is 0 since each pmer in the DNA-space attaches to exactly one pmer with certainty (probability 1) in the basis. The more uncertainty (noise), the larger the entropy and the lesser the reliability.

| **Length** |  $\tau = 50\%$   |
| :------------ | -------------------: |
| 4-mers |  0.97 / 0.88 |
| 6-mers |  0.92 / 0.56 |
| 8-mers |  0.89 / 0.61|

<center>**Table 1.** Noise quality of various nxh bases founded using SaEA, quantified using the expected number of hybridizations of a random pmer and Shannon Entropy of the corresponding distribution (Expected value/Shannon entropy)</center>

With this degree of soundness of the _**h-distance**_ approximation for the Gibbs energy in place and this formulation of the Codeword Design problem and the solution of it using a Self-adaptive Evolutionary
Algorithm, the discussion of how can use this noncrosshybridizing sets in the DNA Space to visualize a DNA sequence can proceed as follow.

<h3><center>3. DNA Visualization</center></h3>
### A. Case of Study
The approach presented by ([Garzon, Bobba, & Hyde, 2004](#Garzon, Bobba, & Hyde, 2004)) is a clear example of how a DNA sequence can be represented through DNA-based methods, in such a way that more information of the sequence may be captured and represented as a signature in 2 Dimensions.

In this approach, a given a string $x$ (ordinarily much larger than $n$) and $B$ the noncrosshybridizin set of $n$-mers; $x$ is said to be h-dependent of $B$ is there some concatenation $c$ of elements of $B$ such that $x$ will hybridize to $c$ under stringency $\tau$, i.e., $|c,x|≤ \tau$. Shredding $x$ to the corresponding fragments according to the components of $c$ in $B$ leads to the following slightly weaker but more manageable definition.

<center>
	$ Signature(x) = V $
</center>

where $V_i$ is the number of fragments of $x$ that hybridize within threshold $\tau$ of strand $i$ in $B$.

The vector V can be visualized as a 1D signature, or as 2D matrix.

<img src="Images/encodeBiotic.png" width="400">
<center>**Figure 1.** Signatures of two plasmid genome in a 2D matrix representation, on a set on noncrosshybridizing probes.</center>

### B. Methodology

To take a more robust representation of the DNA sequence that just a DNA signature, it will be use a RGB color model where each light are represented as

- **R**ed:   &nbsp;&nbsp;&nbsp;&nbsp; nxh set with $n=8$ and $\tau = n/2$

- **G**reen: &nbsp; nxh set with $n=6$ and $\tau = n/2$

- **B**lue: &nbsp;&nbsp;&nbsp;&nbsp; Leave the light at 0 (The nxh set with $n=4$ and $\tau = n/2$ has to noise)

The representation takes into account not only one noncrosshybridizing set but three noncrosshybridizing sets of different sizes as show in Figure 2.

<img src="Images/rgbRepresentation.png" width="400">
<center>**Figure 2.** RGB color model for a DNA sequence.</center>

where the position $(i,j)$ of a color matrix are given by

<center>$Matrix_{i,j}$ = *h-distance*$(B_k, x_k)$</center>

and $B$ and $x$ are the noncrosshybridizing set and the sequence respectively.

<h3><center>4. Experimentation</center></h3>


```python
#######################
# Import libraries
#######################

# Mathematical function: It provides access to the mathematical functions defined by the C standard.
# https://docs.python.org/2/library/math.html
import math

# The glob module#: Finds all the pathnames matching a specified pattern according to the rules used
# by the Unix shell, although results are returned in arbitrary order.
# https://docs.python.org/2/library/glob.html
import glob

# Random: This module implements pseudo-random number generators for various distributions.
# https://docs.python.org/2/library/random.html
import random

# Matplotlib: is a Python 2D plotting library which produces publication quality figures in a variety
# of hardcopy formats and interactive environments across platforms. Matplotlib can be used in Python scripts,
# the Python and IPython shell, the jupyter notebook, web application servers, and four graphical user interface
# toolkits.
# https://matplotlib.org/
import matplotlib.pyplot as plt

# Numpy: Is the fundamental package for scientific computing with Python. It contains among other things:
#    - a powerful N-dimensional array object
#    - sophisticated (broadcasting) functions
#    - tools for integrating C/C++ and Fortran code
#    - useful linear algebra, Fourier transform, and random number capabilities
# http://www.numpy.org/
import numpy as np

# Biopython: Biopython is a set of freely available tools for biological computation
# written in Python by an international team of developers.
# http://biopython.org/
from Bio import SeqIO
```

### A. Dataset

The sequences are from different species and it was provided by the [GenBank](https://www.ncbi.nlm.nih.gov/genbank/), which is the NIH genetic sequence database, an annotated collection of all publicly available DNA sequences.

The species are:
- Plant
    - [Helianthus annuus mitochondrion](https://en.wikipedia.org/wiki/Helianthus_annuus) (Sunflower)
    - [Hordeum vulgare mitochondrion](https://en.wikipedia.org/wiki/Barley) (Barley)
    - [Triticum aestivum mitochondrion](https://en.wikipedia.org/wiki/Common_wheat) (Wheat)


- Virus
    + [Camelpox virus](https://en.wikipedia.org/wiki/Camelpox) (Disease of camels)
    + [Canarypox virus ](https://en.wikipedia.org/wiki/Canarypox) (Disease of wild and captive birds)
    + [Variola major virus ](https://en.wikipedia.org/wiki/Smallpox) (Smallpox)


- Fungi
    + [Ganoderma lucidum mitochondrion](https://en.wikipedia.org/wiki/Lingzhi_mushroom) (Lingzhi mushroom)
    + [Lentinula edodes mitochondrion](https://en.wikipedia.org/wiki/Shiitake) (Shiitake)
    + [Pleurotus ostreatus mitochondrion](https://en.wikipedia.org/wiki/Pleurotus_ostreatus) (Oyster mushroom)


- Bacterium
    + [Anaplasma phagocytophilum](https://en.wikipedia.org/wiki/Anaplasma_phagocytophilum) (Tick-borne fever and pasture fever)
    + [Neisseria gonorrhoeae](https://en.wikipedia.org/wiki/Neisseria_gonorrhoeae) (Gonorrhea)
    + [Streptococcus pyogenes](https://en.wikipedia.org/wiki/Streptococcus_pyogenes) (Otitis media and Mastitis)


Then, load the DNA sequences of each specie.

#### Plant


```python
# All the plant sequences
Plant_sequences = {}

# Fasta files of Plant
fastaFiles = glob.glob("Dataset/Plant/*.fasta")

# Save in the dictionary all the sequences of fasta (key = name, value = sequence)
for file_i in fastaFiles:
    fasta_sequences = SeqIO.parse(open(file_i),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        Plant_sequences[name] = {"strand" : sequence.lower(), "description" : fasta.description}
```

#### Virus


```python
# All the Virus sequences
Virus_sequences = {}

# Fasta files of Virus
fastaFiles = glob.glob("Dataset/Virus/*.fasta")

# Save in the dictionary all the sequences of fasta (key = name, value = sequence)
for file_i in fastaFiles:
    fasta_sequences = SeqIO.parse(open(file_i),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        Virus_sequences[name] = {"strand" : sequence.lower(), "description" : fasta.description}
```

#### Fungi


```python
# All the Fungi sequences
Fungi_sequences ={}

# Fasta files of Fungi
fastaFiles = glob.glob("Dataset/Fungi/*.fasta")

# Save in the dictionary all the sequences of fasta (key = name, value = sequence)
for file_i in fastaFiles:
    fasta_sequences = SeqIO.parse(open(file_i),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        Fungi_sequences[name] = {"strand" : sequence.lower(), "description" : fasta.description}
```

#### Bacterium


```python
# All the Bacterium sequences
Bacterium_sequences = {}

# Fasta files of Bacterium
fastaFiles = glob.glob("Dataset/Bacterium/*.fasta")

# Save in the dictionary all the sequences of fasta (key = name, value = sequence)
for file_i in fastaFiles:
    fasta_sequences = SeqIO.parse(open(file_i),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        Bacterium_sequences[name] = {"strand" : sequence.lower(), "description" : fasta.description}
```

| **Specie** |  **Name**  |  **id **| **DNA lenght** |
| :------------ | :------------------- | :------------------- | :-------------------: |
    | _**Plant**_ |  *Helianthus annuus mitochondrion <br/> Hordeum vulgare mitochondrion <br/>  Triticum aestivum mitochondrion *| CM007908.1 <br/>KU865690.1 <br/>NC_036024.1 |  301004 <br/> 416675 <br/>452526 |
| _**Virus**_ |  *Camelpox virus <br/>Canarypox virus  <br/>Variola major virus* | NC_003391.1 <br/>NC_005309.1 <br/> L22579.1 |  205719 <br/> 359853 <br/> 186103 |
| _**Fungi**_ |  *Ganoderma lucidum mitochondrion <br/> Lentinula edodes mitochondrion <br/> Pleurotus ostreatus mitochondrion* | NC_021750.1 <br/> NC_018365.1 <br/> EF204913.1 |  60635 <br/> 121394 <br/> 73242 |
| _**Bacterium**_ | *Anaplasma phagocytophilum <br/> Neisseria gonorrhoeae <br/> Streptococcus pyogenes* | CP000235.1 <br/> CP016017.1 <br/> CP003116.1 |  1471282 <br/> 942943 <br/> 1750832 |


<center>**Table 2.** Species description </center>

### B. Change DNA strands to images

#### Shreeding the DNA sequence

First, the DNA sequence $x$ will be shredding into $N$ equal pieces for every image


```python
### Chunck the sequences into n parts
# @param: sequence DNA strand
# @param: n the number of partitions
# @return: The sequence splited in n partitions
def chunkIt(sequence, num):
    avg = len(sequence) / float(num)
    out = []
    last = 0.0

    while last < len(sequence):
        out.append(sequence[int(last):int(last + avg)])
        last += avg

    return out

# Save all the dataset shreeding
CompleteDatasetShreeding = {}

print "x shreeding lenght | Description of x"

# Shreding for Plants
for keyPlant in Plant_sequences:
    ShredingKeyPlant = chunkIt(Plant_sequences[keyPlant]["strand"], 1000)
    print len(ShredingKeyPlant), " \t \t   | ", Plant_sequences[keyPlant]["description"]
    CompleteDatasetShreeding[keyPlant] = {"strandShreeding" : ShredingKeyPlant, \
                                          "description" : Plant_sequences[keyPlant]["description"]}

# Shreding for Virus
for keyVirus in Virus_sequences:
    ShredingKeyVirus = chunkIt(Virus_sequences[keyVirus]["strand"], 1000)
    print len(ShredingKeyVirus), " \t \t   | ", Virus_sequences[keyVirus]["description"]
    CompleteDatasetShreeding[keyVirus] = {"strandShreeding" : ShredingKeyVirus, \
                                          "description" : Virus_sequences[keyVirus]["description"]}

# Shreding for Fungi
for keyFungi in Fungi_sequences:
    ShredingKeyFungi = chunkIt(Fungi_sequences[keyFungi]["strand"], 1000)
    print len(ShredingKeyFungi), " \t \t   | ",  Fungi_sequences[keyFungi]["description"]
    CompleteDatasetShreeding[keyFungi] = {"strandShreeding" : ShredingKeyFungi, \
                                          "description" : Fungi_sequences[keyFungi]["description"]}

# Shreding for Bacterium
for keyBacterium in Bacterium_sequences:
    ShredingKeyBacterium = chunkIt(Bacterium_sequences[keyBacterium]["strand"], 1000)
    print len(ShredingKeyBacterium), " \t \t   | ", Bacterium_sequences[keyBacterium]["description"]
    CompleteDatasetShreeding[keyBacterium] = {"strandShreeding" : ShredingKeyBacterium, \
                                          "description" : Bacterium_sequences[keyBacterium]["description"]}
```

    x shreeding lenght | Description of x
    1000  	 	   |  CM007908.1 Helianthus annuus mitochondrion, complete sequence, whole genome shotgun sequence
    1000  	 	   |  NC_036024.1 Triticum aestivum cultivar Chinese Yumai mitochondrion, complete genome
    1000  	 	   |  KU865690.1 UNVERIFIED: Hordeum vulgare subsp. vulgare mitochondrion sequence
    1000  	 	   |  NC_003391.1 Camelpox virus, complete genome
    1000  	 	   |  NC_005309.1 Canarypox virus, complete genome
    1000  	 	   |  L22579.1 Variola major virus (strain Bangladesh-1975) complete genome
    1000  	 	   |  NC_021750.1 Ganoderma lucidum mitochondrion, complete genome
    1000  	 	   |  EF204913.1 Pleurotus ostreatus mitochondrion, complete genome
    1000  	 	   |  NC_018365.1 Lentinula edodes mitochondrion, complete genome
    1000  	 	   |  CP003116.1 Streptococcus pyogenes MGAS15252, complete genome
    1000  	 	   |  CP016017.1 Neisseria gonorrhoeae strain 34769, complete genome
    1000  	 	   |  CP000235.1 Anaplasma phagocytophilum HZ, complete genome


Now, The Shreeding vector $X$ has the same size for all the sequences. It is necesary to load the noncrosshybridizing sets with $n=6,8$ and $\tau = n/2$ found with SaEA(Self-adaptive Evolutionary Algorithm).


```python
### Load the noncrosshybridizin set
# @param: n DNA space of n-mers
# @return: nxhList noncrosshybridizing set of n-mers
def loadNoncrosshybridizinSet(n):
    file_nxh = open("Noncrosshybridizing/" + str(n) + "pmers" + str(n/2) + "tComplete").readlines()
    nxhList = []
    for nxh in file_nxh:
        nxhList.append(nxh.split())
    return nxhList

# Load the nxh sets withn=6, n=8
# Choose a random noncrosshybridizing set of the list
nxh6_mers = random.choice(loadNoncrosshybridizinSet(6))
nxh8_mers = random.choice(loadNoncrosshybridizinSet(8))

print "The noncrosshybridizing set with n=6"
print nxh6_mers
print
print "The noncrosshybridizing set with n=8"
print nxh8_mers
```

    The noncrosshybridizing set with n=6
    ['gaaagc', 'ctaaaa', 'gtagca', 'atattc', 'cgatga', 'acccca', 'caacca', 'acggaa', 'aacttg', 'agtgcg', 'cgtctc', 'caggcc']

    The noncrosshybridizing set with n=8
    ['aagaggtt', 'aggacttc', 'ctacagca', 'gtccattg', 'caggcgca', 'gttcggcg', 'ccactagg', 'ggtaagca', 'aacgctcc', 'atcacgac', 'aaaccggg', 'acggagac', 'ccagataa', 'gtgggtgg', 'tgattcaa', 'actaatat']


Ensure that the noncrosshybridizing set of $n=6,8$ have the same size


```python
lenOfnxh = min([len(nxh6_mers), len(nxh8_mers)])

# Resize the sets to be all the same size
# nxh4_mers = np.random.choice(nxh4_mers, min([len(nxh4_mers), len(nxh6_mers), len(nxh8_mers)]))
nxh6_mers = np.random.choice(nxh6_mers, lenOfnxh)
nxh8_mers = np.random.choice(nxh8_mers, lenOfnxh)

print "The noncrosshybridizing set with n=6 fixed"
print nxh6_mers, ", with lenght: ", len(nxh6_mers)
print
print "The noncrosshybridizing set with n=8 fixed"
print nxh8_mers, ", with lenght: ", len(nxh8_mers)
```

    The noncrosshybridizing set with n=6 fixed
    ['cgtctc' 'caacca' 'cgatga' 'acggaa' 'caggcc' 'ctaaaa' 'acggaa' 'gaaagc'
     'caacca' 'atattc' 'aacttg' 'aacttg'] , with lenght:  12

    The noncrosshybridizing set with n=8 fixed
    ['ggtaagca' 'ccagataa' 'aagaggtt' 'gttcggcg' 'ccagataa' 'ctacagca'
     'aggacttc' 'atcacgac' 'ctacagca' 'aggacttc' 'gtgggtgg' 'aaaccggg'] , with lenght:  12


#### Parce to RGB

The matrices have to be built to each light in the RGB model.


```python
### Load the noncrosshybridizin set
# @param: Shredding sequence of each DNA sequence
# @return: Each matrix of the lights RGB model
def parseToRGB(ShreddingSequence):

    ############################################
    # Red : nxh set with n=8 and tau = n/2
    ############################################

    # Create a empty matrix of size |nxh set| x |Shredding x|
    RedLight = np.zeros((lenOfnxh, len(ShreddingSequence)))

    # Fill the Red matrix
    for nxh_iR in range(lenOfnxh):
        for x_iR in range(len(ShreddingSequence)):
            RedLight[nxh_iR][x_iR] = realHdistance(nxh8_mers[nxh_iR], ShreddingSequence[x_iR], 8)

    ############################################
    # Green : nxh set with n=6 and tau = n/2
    ############################################

    # Create a empty matrix of size |nxh set| x |Shredding x|
    GreenLight = np.zeros((lenOfnxh, len(ShreddingSequence)))

    # Fill the green matrix
    for nxh_iG in range(lenOfnxh):
        for x_iG in range(len(ShreddingSequence)):
            GreenLight[nxh_iG][x_iG] = realHdistance(nxh6_mers[nxh_iG], ShreddingSequence[x_iG], 6)

    ############################################
    # Blue : Fill of 0s
    ############################################

    # Create a empty matrix of size |nxh set| x |Shredding x|
    BlueLight = np.zeros((lenOfnxh, len(ShreddingSequence)))

    return RedLight, GreenLight, BlueLight

# The dictionary with the lights matrix for each DNA sequence
CompleteDatasetRGBMatrix = {}

for key in CompleteDatasetShreeding:
    RedLight, BlueLight, GreenLight = parseToRGB(CompleteDatasetShreeding[key]["strandShreeding"])

    # print CompleteDatasetShreeding[key]["description"], " | ", RedLight.shape, " | ", BlueLight.shape,  " | ", GreenLight.shape
    CompleteDatasetRGBMatrix[key] = {"description" : CompleteDatasetShreeding[key]["description"], \
                                    "RedLight" : RedLight, "BlueLight": BlueLight, "GreenLight": GreenLight}

```

In this moments the image have the size $|nxh\_set|$ x $|Shredding \_ x|$, To be more friendly to the eye, the image is resized to size ($5*|nxh\_set|$ x $|Shredding\_x|/5$). And then join all the matrix lights into a 3D RGB matrix and save the images.


```python
for key in CompleteDatasetRGBMatrix:
    RedLightReshaped = CompleteDatasetRGBMatrix[key]["RedLight"].reshape((lenOfnxh*5, 1000/5))
    GreenLightReshaped = CompleteDatasetRGBMatrix[key]["BlueLight"].reshape((lenOfnxh*5, 1000/5))
    BlueLightReshaped = CompleteDatasetRGBMatrix[key]["GreenLight"].reshape((lenOfnxh*5, 1000/5))

    rgbArray = np.zeros((lenOfnxh*5, 1000/5, 3), 'uint8')
    rgbArray[..., 0] = RedLightReshaped*256
    rgbArray[..., 1] = GreenLightReshaped*256
    rgbArray[..., 2] = BlueLightReshaped*256

    key = str(key)
    CompleteDatasetRGBMatrix[key]["RGB_matrix"] = rgbArray
    scipy.misc.imsave("Images/Results/" + key + '.jpg', rgbArray)
```

### C. Comparison Between species
| |
|:-|:----------------------------------:|:-----------------------------:|:----------------------------------------:|
|**Plant**| <img src="Images/Results/CM007908.1.jpg" width="400">  _Helianthus annuus mitochondrion_ | <img src="Images/Results/KU865690.1.jpg" width="400"> _Hordeum vulgare mitochondrion_ | <img src="Images/Results/NC_036024.1.jpg" width="400"> _Triticum aestivum mitochondrion_ |
|**Virus**| <img src="Images/Results/NC_003391.1.jpg" width="400">  _Camelpox virus_ | <img src="Images/Results/NC_005309.1.jpg" width="400"> _Canarypox virus_ | <img src="Images/Results/L22579.1.jpg" width="400"> _Variola major virus_ |
|**Fungi**| <img src="Images/Results/NC_021750.1.jpg" width="400">  _Ganoderma lucidum mitochondrion_ | <img src="Images/Results/NC_018365.1.jpg" width="400"> _Lentinula edodes mitochondrion_ | <img src="Images/Results/EF204913.1.jpg" width="400"> _Pleurotus ostreatus mitochondrion_ |
|**Bacterium**| <img src="Images/Results/CP000235.1.jpg" width="400">  _Anaplasma phagocytophilum_ | <img src="Images/Results/NC_005309.1.jpg" width="400"> _Neisseria gonorrhoeae_ | <img src="Images/Results/CP003116.1.jpg" width="400"> _Streptococcus pyogenes_ |
<center>**Table 2.** The images that represent each DNA sequence </center>

<h3><center>5. Conclusions and Future work</center></h3>

<h3><center>References</center></h3>

<a id='Neel & Garzon, 2003'></a>
Neel, A.; Garzon, M. H. (2003). Efficiency and Reliability of Genomic Information Storage and Retrieval in DNA-based Memories with Compaction. IEEE Congress on Evolutionary Computation; 2733-2739.

<a id='Bobba et al., 2006'></a>
Bobba, K. C., Neel, A., Phan, V., & Garzon, M. (2006). Reasoning and talking DNA: Can DNA understand English? In C. Mao, & S. Yokomori (Eds.), Proc.12th International Meeting on DNA Computing (Lecture Notes in Computer Science 4287, pp. 337-349). Springer-Verlag.

<a id='Garzon, Bobba, & Hyde, 2004'></a>
Garzon, M. H., Bobba K., & Hyde B. (2004). Digital Information Encoding on DNA. Aspects of Molecular Computing 2004, 152-166.

<a id='Garzon, & Mainali, 2017'></a>
Garzon, M. H., & Mainali, S. (2017). Towards Reliable Microarray Analysis and Design. In: The 9th International Conference on Bioinformatics and Computational Biology.

<a id='Garzon et al., 1997'></a>
Garzon, M. H., Neathery, P., Deaton, R., Murphy, R. C., Franceschetti, D. R., & Stevens Jr, S. E. (1997). A new metric for DNA computing. Proceedings of the 2nd Genetic Programming Conference, 278–472.

<a id='Garzon & Phan, 2009'></a>
Garzon, M. H., & Phan, V. (2009). Optimal codes for computing and self-assembly. International Journal of Nanotechnology and Molecular Computation, 1(1), 1–17.

<a id='Garzon, & Wong, 2011'></a>
Garzon, M. H., & Wong, T. Y. (2011). DNA chips for species identification and biological phylogenies. Natural Computing, 10(1), 375–389.

<a id='Gomez, 2004'></a>
Gomez, J. (2004). Self-adaptation of operator rates in evolutionary algorithms. In Proceedings of
the Genetic and Evolutionary Computation Conference (GECCO 2004), pages 1162–1173.

<a id='Neel & Garzon, 2006'></a>
Neel, A., & Garzon, M. H. (2006). Semantic Retrieval in DNA-Based Memories with Gibbs Energy Models. Biotechnology Progress, 22(1), 86–90.

<a id='Phan & Garzon, 2009'></a>
Phan, V., & Garzon, M. H. (2009). On codeword design in metric DNA spaces. J. of Natural Computing, 8(3), 571–588. Roman, J. (1995). The theory of error-correcting codes. Springer-Verlag.

<a id='Phan & Garzon, 2009'></a>
Prieto, J., Leon, E., &amp; Garzon, M. H. (2017). Nearly optimal DNA Codeword Design using a
Self-adaptive Evolutionary Algorithm. Unpublished manuscript. Universidad Nacional de
Colombia, Colombia.

<a id='Qian & Winfree, 2011'></a>
Qian, L., & Winfree, E. (2011). Scaling up digital circuit computation with DNA strand displacement cascades. Science, 332, 1196–1201.

<a id='Seeman, 2006'></a>
Seeman, N. C. (2003). DNA in a material world. Nature, 421, 427-431.

