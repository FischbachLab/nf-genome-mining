# NF-BLAST

A Nextflow script that will run BLAST and using fasta files on S3 and databases available on a shared filesystem (EFS).

```{bash}
aws batch submit-job \
    --profile maf \
    --job-name nf-gm-0707-1 \
    --job-queue default-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-genome-mining,\
"--seedfile","s3://genomics-workflow-core/Results/GenomeMining/IMG/20220707/test_seedfile.csv",\
"--project","IMG",\
"--prefix","20220707"
```
