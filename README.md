# NF-BLAST

A Nextflow pipeline to convert annotated genome data from IMG to AWS Neptune property graph format.

```{bash}
aws batch submit-job \
    --profile maf \
    --job-name nf-gm-test-## \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-genome-mining,\
"--seedfile","s3://genomics-workflow-core/Results/GenomeMining/00_TEST/20220707/test.20220707_seedfile.csv",\
"--project","00_TEST",\
"--prefix","20220707"
```
