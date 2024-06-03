# serotyper
Pipeline for Viral Serotyping

```
bowtie2-build --threads 20 Bos_taurus.ARS-UCD1.3.dna.toplevel.fa  Bos_taurus.ARS-UCD1.3.dna.toplevel.fa
nextflow run main.nf --input "data/illuminape/SRR27845*_{1,2}.fastq.gz" --skipTrim true --hostile_ext="--index /data/serot
yper/resources/index/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa" -resume
```

```mermaid
graph TB
A[Fastq Files]  --> B((HiCFlow))
C[MSA]  --> B((serotyper))
D[MSA Metadata]  --> B((serotyper))
E[Summary File] -- optional  --> B((serotyper))
B --> F[SE]


F -- optional --> G[pycoQC]
F -- dehosting --> H[HOSTILE] -- Chimera Read Removal--> I[YACRD] -- Serotyping/Coverage/Bacterial read removal --> J[VirStrain + Alignment + Coverage] -- assembly --> K[Canu/Flye/Both + Polishing + Scaffolding]

K --> BLASTN --> Serotyping

B --> L[PE]

L --> FASTQC
L --> Deshosting --> N[Adapter Trimming] --> J
J -- assembly --> P[Unicycler/SPADES/Both] --> BLASTN 
```

