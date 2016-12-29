<h1 align="center">
    juliet - Minor Variant Caller
</h1>

<p align="center">
  <img src="img/juliet.png" alt="Logo of Juliet" width="400px"/>
</p>

## Install
Install the unanimity suite and one of the binaries is called `juliet`.

## Input data
*Juliet* operates on aligned CCS records in the BAM format.
Reads should be created with [CCS2](../PBCCS.md).
BAM files have to PacBio-compliant, meaning, cigar `M` is forbidden.
*Juliet* currently does not demultiplex barcoded data,
provide one BAM per barcode.

## Scope
Current scope of *Juliet* is identification of codon-wise variants in coding
regions.

## Output
*Juliet* provides a JSON and HTML file. The JSON file contains for each gene
the variant positions. Each variant position consists of the reference codon,
reference aminoacid, relative aminoacid position in the gene, the mutated codon
and aminoacid, the coverage, possible annotated drug resistance mutations, and
counts of the multiple-sequence alignment of the -3 to +3 context positions.

The HTML page is a 1:1 conversion of the JSON file and contains the identical
information, only human-readable.

## Target configuration
*Juliet* is a multi purpose minor variant caller with preinstalled configurations.
A target configuration may contain multiple coding regions and optional drug
resistance mutation positions.

### Predefined target config
Running on predefined organism:
```
$ juliet -c "<HIV>" data.align.bam -o patientZero
```

<img src="img/juliet_hiv-hiv.png" width="500px">

Currently available configs are: `<HIV>`

### Own target config
Create a JSON file. The root child genes contains a list of coding regions, with
begin and end, the name of the gene, and a list of drug resistent mutations
drms. Each DRM consists of its name and the positions it targets. The drms
field is optional. If provided, the referenceSequence is being used to call
mutations, otherwise it will be tested against the major codon. All indices are
with respect to the provided alignment space, 1-based.
Save following as hiv.json:
```
{
    "genes": [
        {
            "begin": 2253,
            "drms": [
                {
                    "name": "PI",
                    "positions": [ 24,32,46,47,50,54,76,82,84,88,90 ]
                },
                {
                    "name": "PI S",
                    "positions": [ 23,24,30,32,46,47,48,50,53,54,73,76,82,83,84,85,88,90 ]
                }
            ],
            "end": 2550,
            "name": "Protease"
        },
        {
            "begin": 2550,
            "drms": [
                {
                    "name": "NNRTI",
                    "positions": [ 100,101,103,106,138,179,181,190,190,227,230 ]
                },
                {
                    "name": "NNRTI S",
                    "positions": [ 65,67,69,70,74,75,77,115,116,141,151,184,210,215,219 ]
                }
            ],
            "end": 3870,
            "name": "Reverse Transcriptase"
        },
        {
            "begin": 3870,
            "drms": [],
            "end": 4230,
            "name": "RNase"
        }
    ],
    "referenceSequence": "TGGAAGGGCT..."
}
```

Run with own target config:
```
$ juliet -c hiv.json data.align.bam -o patientZero
```

<img src="img/juliet_hiv-own.png" width="500px">


### No target config
If no target config has been specific, it is recommended to at least specify the
region of interest to mark the correct reading frame. The output will be labeled
with unknown as gene name:
```
$ juliet data.align.bam -o patientZero
```

<img src="img/juliet_hiv-unknown.png" width="500px">

# FAQ

### Can I use overlapping regions?
Yes! Each gene is treated separately. Overlapping region, even with different
reading frames are possible.

### Can I call a smaller window from a target config?
Use `--region` to specify the begin-end window to subset the target config.