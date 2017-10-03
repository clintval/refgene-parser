<h1 align="center">RefGene-Parser</h1>

<h3 align="center">Installation</h3>

Clone the respository and then install locally:

```bash
$ git clone https://github.com/clintval/refgene-parser.git
$ pip install refgene-parser
```

<h3 align="center">Dependencies</h3>

- None

<h3 align="center">Tutorial</h3>

Iterate over the records in a RefGene file:

```python

from refgene import RefGene

refgene = RefGene('mm10.refGene.txt.gz')

for i, gene in enumerate(refgene):
    if i >= 10: break
    print(gene)
```

```python
Gene("chr2", 11705292, 11733985, "+", name="Il15ra", id="NM_001271498")
Gene("chr7", 142434976, 142440396, "+", name="Syt8", id="NM_001285857")
Gene("chr1", 78424744, 78488897, "-", name="Farsb", id="NM_011811")
Gene("chr11", 62574485, 62600305, "+", name="Trpv2", id="NM_011706")
Gene("chr12", 100199434, 100209824, "+", name="Calm1", id="NM_009790")
Gene("chr5", 30933142, 30945480, "-", name="Cgref1", id="NM_026770")
Gene("chr4", 142084297, 142088101, "+", name="Tmem51os1", id="NR_027137")
Gene("chr10", 77257772, 77259223, "-", name="Gm10941", id="NR_026944")
Gene("chr10", 77706586, 77706986, "+", name="Gm10272", id="NR_026831")
Gene("chr7", 100549116, 100607996, "-", name="Mrpl48", id="NR_003559")
Gene("chr7", 7212995, 7278289, "-", name="Vmn2r29", id="NR_003555")
```


<h5 align="center">Exact match for a gene symbol name</h5>

> Will return the first record matching with `name == Kras`

```python
Kras = refgene.gene_by_name('Kras')
print(Kras)
print(Kras.sam_interval)
print(Kras.num_exons)
```

```python
Gene("chr6", 145216698, 145250231, "-", name="Kras", id="NM_021284")
'chr6:145216698-145250231'
5
```


<h5 align="center">Exact match for a gene ID</h5>

> Will return the first record matching `id == NM_009085`

```python
NM_009085 = refgene.gene_by_id('NM_009085')
```

```python
Gene("chr17", 46243919, 46248045, "-", name="Polr1c", id="NM_009085")
```


<h5 align="center">Fuzzy match for a gene symbol name</h5>

```python
list(refgene.genes_by_name_pattern('^.ras$'))
```

```python
[Gene("chr3", 103058284, 103067914, "+", name="Nras", id="NM_010937"),
 Gene("chr7", 45018006, 45021644, "+", name="Rras", id="NM_009101"),
 Gene("chrX", 7924275, 7928607, "-", name="Eras", id="NM_181548"),
 Gene("chr7", 141189933, 141194004, "-", name="Hras", id="NM_008284"),
 Gene("chr6", 145216698, 145250231, "-", name="Kras", id="NM_021284"),
 Gene("chr9", 99385419, 99436712, "-", name="Mras", id="NM_008624"),
 Gene("chr7", 141189933, 141194004, "-", name="Hras", id="NM_001130444"),
 Gene("chr7", 141189933, 141194004, "-", name="Hras", id="NM_001130443")]
```


<h5 align="center">Fuzzy match for a gene ID</h5>

```python
list(refgene.genes_by_id_pattern('NR_00355*'))
```

```python

[Gene("chr7", 100549116, 100607996, "-", name="Mrpl48", id="NR_003559"),
 Gene("chr7", 7212995, 7278289, "-", name="Vmn2r29", id="NR_003555"),
 Gene("chr15", 82636749, 82642045, "-", name="Cyp2d13", id="NR_003552"),
 Gene("chr3", 92373914, 92375229, "+", name="Sprr2g", id="NR_003548"),
 Gene("chr13", 104173720, 104178466, "-", name="Trappc13", id="NR_003546"),
 Gene("chr5", 30950065, 30955095, "+", name="Abhd1", id="NR_003522"),
 Gene("chr17", 3064317, 3084183, "-", name="Pisd-ps2", id="NR_003519"),
 Gene("chrUn_JH584304", 52673, 59689, "-", name="Pisd-ps3", id="NR_003518"),
 Gene("chr11", 3124020, 3131944, "+", name="Pisd-ps1", id="NR_003517"),
 Gene("chr16", 97536080, 97560901, "+", name="Mx2", id="NR_003508"),
 Gene("chr5", 120812634, 120824160, "+", name="Oas1b", id="NR_003507"),
 Gene("chr5", 10865028, 10870808, "+", name="Gm6455", id="NR_003596"),
 Gene("chr19", 5842295, 5845480, "-", name="Neat1", id="NR_003513"),
 Gene("chr15", 62217540, 62219451, "-", name="H2afy3", id="NR_003523"),
 Gene("chrX_GL456233_random", 268798, 270075, "+", name="Zf12", id="NR_003547"),
 Gene("chr16", 97447034, 97462906, "-", name="Mx1", id="NR_003520"),
 Gene("chr11", 88964665, 88966917, "-", name="Gm15698", id="NR_003564"),
 Gene("chr13", 12614064, 12650395, "-", name="Gpr137b-ps", id="NR_003568")]
```
