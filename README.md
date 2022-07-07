# NF-BLAST

A Nextflow pipeline to convert annotated genome data from IMG to AWS Neptune property graph format.

```{bash}
aws batch submit-job \
    --profile maf \
    --job-name nf-gm-0707-2 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-genome-mining,\
"--seedfile","s3://genomics-workflow-core/Results/GenomeMining/IMG/20220707/test_seedfile.csv",\
"--project","IMG",\
"--prefix","20220707"
```
