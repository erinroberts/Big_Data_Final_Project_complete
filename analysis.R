#This script uses biomaRt to get the full ensemble gene id in the two chromosomes found to have GIMAP sequences (chr 4 and 7), 
# and then identifies the flanking regions around 

#download BiomaRt
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# biomaRt help : http://www.ensembl.info/blog/2015/06/01/biomart-or-how-to-access-the-ensembl-data-from-r/


#get ensembl data for all genes on a chromosome. Only need human chromosomes 4 and 7
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
head(listAttributes(useDataset(dataset = "hsapiens_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org"))), 10)

#the getBM function allows you to build a biomaRt query to get information you want
chr4_genes <- getBM(attributes=c('ensembl_gene_id','gene_biotype', 'hgnc_symbol','chromosome_name','start_position','end_position'), filters = 'chromosome_name', values = "4", mart = ensembl)
head(chr4_genes)
chr7_genes <- getBM(attributes=c('ensembl_gene_id','gene_biotype', 'hgnc_symbol','chromosome_name','start_position','end_position'), filters = 'chromosome_name', values = "7", mart = ensembl)
head(chr7_genes)


#access 100kb on left and right of GIMAP genes on these chromosomes with position from tabular output of blast search, get ensemble_gene_ids from these regions. 
#print the output to the command line with the name of each
grab_100kb_NM_001195138.1_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(4, 394, 100394), mart = ensembl)
grab_100kb_NM_001195138.1_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(4, 1, 394), mart = ensembl)
a = 'grab_100kb_NM_001195138.1_right'
print(a)
print(grab_100kb_NM_001195138.1_right)
b = 'grab_100kb_NM_001195138.1_left'
print(b)
print(grab_100kb_NM_001195138.1_left)

grab_100kb_NM_001199577.1_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 794, 102262), mart = ensembl)
grab_100kb_NM_001199577.1_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 794), mart = ensembl)
c = 'grab_100kb_NM_001199577.1_right'
print(c)
print(grab_100kb_NM_001199577.1_right)
d = 'grab_100kb_NM_001199577.1_left'
print(d)
print(grab_100kb_NM_001199577.1_left)

grab_100kb_NM_001244071.1_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 564), mart = ensembl)
grab_100kb_NM_001244071.1_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 564, 103479), mart = ensembl)
e = 'grab_100kb_NM_001244071.1_left'
print(e)
print(grab_100kb_NM_001244071.1_left)
f = 'grab_100kb_NM_001244071.1_right'
print(f)
print(grab_100kb_NM_001244071.1_right)

grab_100kb_NM_001244072.1_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 564), mart = ensembl)
grab_100kb_NM_001244072.1_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 564, 103913), mart = ensembl)
g = 'grab_100kb_NM_001244072.1_left'
print(g)
print(grab_100kb_NM_001244072.1_left)
h = 'grab_100kb_NM_001244072.1_right'
print(h)
print(grab_100kb_NM_001244072.1_right)

grab_100kb_NM_001303630.1_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 435), mart = ensembl)
grab_100kb_NM_001303630.1_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 435, 101903), mart = ensembl)
j = 'grab_100kb_NM_001303630.1_left'
print(j)
print(grab_100kb_NM_001303630.1_left)
k = 'grab_100kb_NM_001303630.1_right'
print(k)
print(grab_100kb_NM_001303630.1_right)

grab_100kb_NM_015660.2_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 116), mart = ensembl)
grab_100kb_NM_015660.2_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 116, 101443), mart = ensembl)
l = 'grab_100kb_NM_015660.2_left'
print(l)
print(grab_100kb_NM_015660.2_left)
m = 'grab_100kb_NM_015660.2_right'
print(m)
print(grab_100kb_NM_015660.2_right)

grab_100kb_NM_018326.2_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 140), mart = ensembl)
grab_100kb_NM_018326.2_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 140, 101965), mart = ensembl)
o = 'grab_100kb_NM_018326.2_left'
print(o)
print(grab_100kb_NM_018326.2_left)
p = 'grab_100kb_NM_018326.2_right'
print(p)
print(grab_100kb_NM_018326.2_right)

grab_100kb_NM_018384.4_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 409), mart = ensembl)
grab_100kb_NM_018384.4_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 409, 101877), mart = ensembl)
q = 'grab_100kb_NM_018384.4_left'
print(q)
print(grab_100kb_NM_018384.4_left)
r = 'grab_100kb_NM_018384.4_right'
print(r)
print(grab_100kb_NM_018384.4_right)

grab_100kb_NM_024711.5_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 564), mart = ensembl)
grab_100kb_NM_024711.5_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 564, 103703), mart = ensembl)
s = 'grab_100kb_NM_024711.5_left'
print(s)
print(grab_100kb_NM_024711.5_left)
t = 'grab_100kb_NM_024711.5_right'
print(t)
print(grab_100kb_NM_024711.5_right)

grab_100kb_NM_130759.3_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 183), mart = ensembl)
grab_100kb_NM_130759.3_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 183, 104420), mart = ensembl)
u = 'grab_100kb_NM_130759.3_left'
print(u)
print(grab_100kb_NM_130759.3_left)
v = 'grab_100kb_NM_130759.3_right'
print(v)
print(grab_100kb_NM_130759.3_right)

grab_100kb_NM_153236.3_left= getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 87), mart = ensembl)
grab_100kb_NM_153236.3_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 87, 101229), mart = ensembl)
w = 'grab_100kb_NM_153236.3_left'
print(w)
print(grab_100kb_NM_153236.3_left)
x = 'grab_100kb_NM_153236.3_right'
print(x)
print(grab_100kb_NM_153236.3_right)

grab_100kb_NM_175571.3_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 1882), mart = ensembl)
grab_100kb_NM_175571.3_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1882, 104187), mart = ensembl)
y = 'grab_100kb_NM_175571.3_left'
print(y)
print(grab_100kb_NM_175571.3_left)
z = 'grab_100kb_NM_175571.3_right'
print(z)
print(grab_100kb_NM_175571.3_right)

grab_100kb_XM_005249950.4_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 1562), mart = ensembl)
grab_100kb_XM_005249950.4_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1562, 103867), mart = ensembl)
a1 = 'grab_100kb_XM_005249950.4_left'
print(a1)
print(grab_100kb_XM_005249950.4_left)
b1 = 'grab_100kb_XM_005249950.4_right'
print(b1)
print(grab_100kb_XM_005249950.4_right)

grab_100kb_XM_005249951.3_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 1252), mart = ensembl)
grab_100kb_XM_005249951.3_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1252, 103560), mart = ensembl)
c1 = 'grab_100kb_XM_005249951.3_left'
print(c1)
print(grab_100kb_XM_005249951.3_left)
d1 = 'grab_100kb_XM_005249951.3_right'
print(d1)
print(grab_100kb_XM_005249951.3_right)

grab_100kb_XM_005250017.1_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 1, 155), mart = ensembl)
grab_100kb_XM_005250017.1_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(7, 155, 101983), mart = ensembl)
e1 = 'grab_100kb_XM_005250017.1_left'
print(e1)
print(grab_100kb_XM_005250017.1_left)
f1= 'grab_100kb_XM_005250017.1_right'
print(f1)
print(grab_100kb_XM_005250017.1_right)

grab_100kb_XM_011531515.1_left = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(4, 1, 472), mart = ensembl)
grab_100kb_XM_011531515.1_right = getBM(c("ensembl_gene_id"), filters = c("chromosome_name", "start", "end"), values = list(4, 472, 101583), mart = ensembl)
g1 = 'grab_100kb_XM_011531515.1_left'
print(g1)
print(grab_100kb_XM_011531515.1_left)
j1 = 'grab_100kb_XM_011531515.1_right'
print(j1)
print(grab_100kb_XM_011531515.1_right)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]





