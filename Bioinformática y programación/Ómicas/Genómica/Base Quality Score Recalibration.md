Base Quality Score Recalibration (BQSR) is a **two-step process** used to correct systematic bias that affects the assignment of base quality scores by the sequencer. Generally, it is performed using the `BaseRecalibrator` and `ApplyBQSR` tools in [[GATK]].

## 1. Building the error model

Creating the error model is done using `BaseRecalibrator`. The `BaseRecalibrator` requires known sites of variation to distinguish between sequencing errors and true variants. These known sites are provided to the tool using the `--known-sites` option, which typically includes databases like dbSNP or the Mills and 1000G gold standard indels.

```
gatk BaseRecalibrator\
-I NA12878_markdup.bam \
-R ref.fa \ 
-O NA12878_markdup_bqsr.table \ 
--known-sites /fdb/GATK_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz 
```

Known sites can be downloaded at the [GATK Resource Bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0).

For legacy bundles, known sites can be downloaded via ftp:

```
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
```

## 2. Applying the model