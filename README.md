#FASTAcrobat

FASTAcrobat is a web application tool to breakdown sequence data and display information about the sequence
and extract DNA, RNA and protein data.

#### Authors
- Olu Coker
- Alberto Scicali
- K. Jeselle Clark
- Chris Snyder

### How to run locally:
- Download the `app.R` file
- In your favorite terminal:
 - `cd` ~/location/of/app.R
 - `Rscript app.r`
 - A local host ip will be displayed such as `http://127.0.0.1:XXXX`
 - Use your favorite web browser to navigate to that URL
 - Enjoy the app

### How to use:

- Upload a `fasta` formatted file
 - The file must have sequence data in `fasta` format to function, ex.
 ```
     >gi|407317075|gb|JQ064943.1| Homo sapiens haplotype 1-46-13-4 D(4) dopamine receptor (DRD4) gene, DRD4-4R allele, partial cds
    CCTTCCCCCACGCCACCCGCGCCCCGCCTCCCCCAGGACCCCTGCGGCCCCGACTGTGCGCCCCCCGCGC
    CGGGCCTTCCCCGGGGTCCCTGCGGCCCCGACTGTGCGCCCGCCGCGCCCGGCCTCCCCCAGGACCCCTG
    CGGCCCCGACTGTGCGCCCCCCGCGCCCGGCCTCCCCCCGGACCCCTGCGGCTCCAACTGTGCTCC
 ```
- Use the dropdown to select which type of sequence it is. Only `DNA`, `RNA` or `Protein` sequences are accepted
- Once the `fasta` file is uploaded, you're free to explore the data and graphs.
- If the uploaded `fasta` has multiple sequences, a drop down selector will be available to select which sequence to be analyzed.


