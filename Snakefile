"""
Snakemake pipeline for analyzing 16S RNA data using QIIME2 (2024.10)

Usage
-----

An example command is as following:

snakemake --cores 10 --use-conda results/AB/AB+fp-f17-r21+bb-t18+rrf-d10000+cls-gg.zip
"""

import glob
from platform import system

_os = system()

if _os == "Linux":
     qiime_env = "envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
elif _os == "Darwin":
     qiime_env = "envs/qiime2-amplicon-2024.10-py310-osx-conda"

rule export_artifact_2:
     """Export Artifact content to a folder"""
     message:
          "Extracting artifact"
     input:
          "results/{cohort}/{cohort}.qza"
     output:
          directory("temp/{cohort}/{cohort}/")
     conda:
          qiime_env
     shell:
          "qiime tools export "
          "--input-path {input} "
          "--output-path {output}"

rule export_artifact:
     """Export Artifact content to a folder"""
     message:
          "Extracting artifact"
     input:
          "results/{cohort}/{cohort}+{etc}.qza"
     output:
          directory("temp/{cohort}/{cohort}+{etc}")
     conda:
          qiime_env
     shell:
          "qiime tools export "
          "--input-path {input} "
          "--output-path {output}"

def get_allfile_names(wildcards):
     """Get all files from a folder"""
     input_folder = os.path.join("temp", wildcards.cohort, wildcards.id)
     return " ".join([os.path.join(input_folder, x) for x in os.listdir(input_folder) if "fastq" in x])

rule help:
     """Print out the information of rules on screen (docstrings)."""
     run:
          from termcolor import cprint
          for rule in workflow.rules:
               try:
                    if not rule.name.startswith('_'):
                         cprint(rule.name, "red", attrs=["bold"])
                         cprint(rule.docstring, "green")
                         cprint("")
               except:
                    pass

rule qza_fastqc:
     """Run fastqc quality control analysis"""
     message:
          "Applying FASTQC"
     input:
          "temp/{cohort}/{id}"
     output:
          directory("results/{cohort}/quality/{id}/fastqc/")
     threads:
          20
     conda:
          "envs/quality.yml"
     params:
          get_allfile_names
     shell:
          "mkdir {output} && "
          "fastqc -o {output} -f fastq -t {threads} {params}"

rule qza_multiqc:
     """Combines multiple Fastqc reports using MultiQC. This rule works for one folder"""
     message:
          "MultiQC"
     input:
          "results/{cohort}/quality/{id}/fastqc/"
     output:
          directory("results/{cohort}/quality/{id}/multiqc/")
     conda:
          "envs/quality.yml"
     shell:
          "multiqc -o {output} {input}"

rule download_taxonomy:
     """Download taxonomy information from ncbi"""
     output:
          "db/nodes.dmp",
          "db/names.dmp"
     shell:
          """
          cd db
          wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
          tar -xvf taxdump.tar.gz names.dmp nodes.dmp
          """

rule create_taxonpath:
   """Create taxonpath.json file"""
   input:
         "db/nodes.dmp"
   output:
         "db/taxonpath.json",
   shell:
      "python scripts/create_taxonpath.py "
      "-i {input} "
      "-o {output} "

rule create_taxon_names:
   """Create names.json file"""
   input:
         "db/names.dmp"
   output:
         "db/names.json",
   shell:
      "python scripts/names2json.py "
      "-i {input} "
      "-o {output} "

rule dataset_multiqc:
     """Combines multiple Fastqc reports using Multiqc. This rule combines all the FastqC reports of one cohort"""
     message:
          "MultiQC"
     input:
          "results/{cohort}/quality/"
     output:
          directory("results/{cohort}/quality_summary/")
     conda:
          "envs/quality.yml"
     shell:
          "multiqc -dd 2 -d -o {output} {input}"

rule manifest:
     """Create manifest file: utilizing scripts/create_manifest_file.py script"""
     message:
          "Creating manifest file"
     input:
          "data/{cohort}/"
     output:
          "results/{cohort}/{cohort}_manifest.tsv"
     conda:
          "envs/other.yml"
     shell:
          "python scripts/create_manifest_file.py -i {input} -o {output}"

rule import_data:
     """Import data: Import the raw fastq files to Qiime2 artifact with qza extension"""
     message:
          "Import data using manifest file"
     input:
          "results/{cohort}/{cohort}_manifest.tsv"
     output:
          "results/{cohort}/{cohort}.qza"
     message:
          "Import data"
     conda:
          qiime_env
     shell:
          "qiime tools import "
          "--type 'SampleData[PairedEndSequencesWithQuality]' "
          "--input-path {input} "
          "--input-format PairedEndFastqManifestPhred33V2 "
          "--output-path {output} "

rule trim_fastp:
     """Crop primers: Crop primers from both end reads."""
     input:
          qza="results/{cohort}/{id}.qza"
     output:
          "results/{cohort}/{id}+fp-f{len1, \d+}-r{len2, \d+}.qza"
     message:
          "Trimming using fastp"
     conda:
          qiime_env
     shell:
          "python scripts/fastp.py --inputf {input} "
          "--len1 {wildcards.len1} --len2 {wildcards.len2} "
          "--outputf {output}"

rule trim_bbduk:
     """Quality trimming using bbduk"""
     input:
          qza="results/{cohort}/{id}.qza"
     output:
          "results/{cohort}/{id}+bb-t{threshold, \d+}.qza"
     message:
          "Trimming using bbduk"
     conda:
          qiime_env
     shell:
          "python scripts/bbduk.py -i {input} "
          "-q {wildcards.threshold} -o {output}"

rule dada2:
     """Dada2 algorithm"""
     input:
          "results/{cohort}/{id}.qza"
     output:
          table="results/{cohort}/{id}+dd_table.qza",
          stats="results/{cohort}/{id}+dd_stats.qza",
          repseq="results/{cohort}/{id}+dd_seq.qza"
     message:
          "Dada2 analysis"
     threads: 30
     conda:
          qiime_env
     shell:
          "qiime dada2 denoise-paired "
          "--p-trunc-len-f 0 --p-trunc-len-r 0 "
          "--i-demultiplexed-seqs {input} "
          "--o-table {output.table} "
          "--o-representative-sequences {output.repseq} "
          "--o-denoising-stats {output.stats} "
          "--verbose --p-n-threads {threads}"

rule rarefy:
     """rarefy feature table."""
     input:
          "results/{cohort}/{id}_table.qza"
     output:
          "results/{cohort}/{id}+rrf-d{r, \d+}_table.qza"
     conda:
          qiime_env
     shell:
          "qiime feature-table rarefy "
          "--i-table {input} "
          "--p-sampling-depth {wildcards.r} "
          "--o-rarefied-table {output}"

rule dump_seq_file:
     """dump seq.qza for rarefied table"""
     input:
          "results/{cohort}/{id}_seq.qza"
     output:
          "results/{cohort}/{id}+rrf-d{r}_seq.qza"
     shell:
          "cp {input} {output}"

rule plot_dada_stats:
     """Plot dada2 stats results in form of distribution plot of the percentage of remained reads from the original data"""
     input:
          "results/{cohort}/{id}+dd_stats.qza"
     output:
          "results/{cohort}/{id}+dd_stats.jpg"
     conda:
          qiime_env
     shell:
          "python scripts/plot_dada.py --inp {input} --plot {output}"

def get_dada_jpgs(wildcards):
     """Get dada2 plot file names from a folder"""
     input_folder = os.path.join("results", wildcards.cohort)
     ret = [os.path.join(input_folder, x) for x in os.listdir(input_folder) if "+dd_stats.qza" in x]
     ret = [x.replace(".qza", ".jpg") for x in ret]
     return ret

def get_dada_jpgs_comma_separated(wildcards):
     """Get dada2 plot file names from a folder separated by comma"""
     input_folder = os.path.join("results", wildcards.cohort)
     ret = [os.path.join(input_folder, x) for x in os.listdir(input_folder) if "+dd_stats.qza" in x]
     ret = [x.replace(".qza", ".jpg") for x in ret]
     return ",".join(ret)

rule dada_stats_report:
     """Combines multiple dada2 stats plots"""
     input:
          get_dada_jpgs
     output:
          "results/{cohort}/{cohort}_dada_stats.pdf"
     conda:
          "envs/other.yml"
     params:
          get_dada_jpgs_comma_separated
     shell:
          "python scripts/report_stats.py --inp {params} --outp {output}"

rule download_silva_classifier:
     """Download pretrained SILVA taxonomy classifier"""
     output:
          "classifiers/silva-classifier.qza"
     conda:
          "envs/other.yml"
     shell:
          "wget https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza -O classifiers/silva-classifier.qza"

rule download_gg_classifier:
     """Download pretrained Greengenes2 taxonomy classifier"""
     output:
          "classifiers/gg-classifier.qza"
     conda:
          "envs/other.yml"
     shell:
          "wget https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes2/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza -O classifiers/gg-classifier.qza"

rule taxonomy:
     """Assign taxonomy to ASVs"""
     input:
          seq = "results/{cohort}/{id}_seq.qza",
          classifier = "classifiers/{cls}-classifier.qza"
     output:
          taxonomy= "results/{cohort}/{id}+cls-{cls}_taxonomy.qza",
     conda:
          qiime_env
     threads: 30
     message:
          "Assign taxonomy using {wildcards.cls} database"
     shell:
          "qiime feature-classifier classify-sklearn "
          "--i-classifier {input.classifier} "
          "--i-reads {input.seq} --p-n-jobs {threads} "
          "--o-classification {output.taxonomy}"

rule mafft:
     """Creating phylogenetic tree step 1 """
     input:
          "results/{cohort}/{id}_seq.qza"
     output:
          "results/{cohort}/tree/{id}_seqaligned.qza"
     conda:
          qiime_env
     message:
          "Mafft alignment"
     shell:
          "qiime alignment mafft "
          "--i-sequences {input} "
          "--o-alignment {output}"

rule alignment_mask:
     """Creating phylogenic tree step 2
     """
     input:
          "results/{cohort}/tree/{id}_seqaligned.qza"
     output:
          "results/{cohort}/tree/{id}_seqaligned_masked.qza"
     message:
          "Alignment mask"
     conda:
          qiime_env
     shell:
          "qiime alignment mask "
          "--i-alignment {input} "
          "--o-masked-alignment {output}"

rule fasttree:
     """Creating phylogenetic tree step 3
     """
     input:
          "results/{cohort}/tree/{id}_seqaligned_masked.qza"
     output:
          "results/{cohort}/tree/{id}+fasttree.qza"
     conda:
          qiime_env
     message:
          "creating tree"
     shell:
          "qiime phylogeny fasttree "
          "--i-alignment {input} "
          "--o-tree {output}"

rule midpoint_root:
     """Creating phylogenetic tree step 4
     """
     input:
          "results/{cohort}/tree/{id}+fasttree.qza"
     output:
          "results/{cohort}/{id}+fasttree_rooted.qza"
     message:
          "Midpoint rooting the tree"
     conda:
          qiime_env
     shell:
          "qiime phylogeny midpoint-root "
          "--i-tree {input} "
          "--o-rooted-tree {output}"

rule export_tree:
     """Creating phylogenitic tree step 5: Exporting tree from Artifact file
     """
     input:
          "results/{cohort}/{id}+fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}+fasttree.nwk"
     conda:
          qiime_env
     message:
          "Save the tree.nwk file"
     params:
          "results/{cohort}/{id}+fasttree_rooted"
     shell:
          "qiime tools export "
          "--input-path {input} "
          "--output-path {params} && "
          "cp {params}/tree.nwk {output} && "
          "rm -r {params}"

rule make_biom:
     """Create Biom table"""
     input:
          table="results/{cohort}/{id}_table.qza",
          taxonomy="results/{cohort}/{id}+cls-{cls}_taxonomy.qza"
     output:
          "results/{cohort}/{id}+cls-{cls}_asv.biom"
     message:
          "Making biom table {output}"
     conda:
          qiime_env
     shell:
          "python scripts/make_biom.py --tablef {input.table} "
          "--taxonomy {input.taxonomy} "
          "--output {output}"

rule extract_biom:
     """Auxillary rule to extract biom from Artifact (QZA)"""
     input:
          "results/{cohort}/{cohort}+{id}.qza"
     output:
          "results/{cohort}/{cohort}+{id}.biom"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype biom"

rule extract_sequence:
     """Auxillary rule to exract dada2 ASV sequences from (qza) Artifact to tsv file"""
     input:
          "results/{cohort}/{cohort}+{id}+dd_seq.qza"
     output:
          "results/{cohort}/{cohort}+{id}+dd_seq.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_stats:
     """Auxillary rule to exract dada2 stats from (qza) Artifact to tsv file"""
     input:
          "results/{cohort}/{cohort}+{id}+dd_stats.qza"
     output:
          "results/{cohort}/{cohort}+{id}+dd_stats.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule export_phyloseq:
     """Create Phyloseq object RDS file"""
     input:
          biom="results/{cohort}/{id}+cls-{cls}_asv.biom",
          tree="results/{cohort}/{id}+fasttree.nwk"
     output:
          "results/{cohort}/{id}+cls-{cls}+phyloseq.RDS"
     conda:
          "envs/phyloseq.yml"
     shell:
          "Rscript scripts/export_phyloseq.R "
          "--biom {input.biom} --tree {input.tree} "
          "--outp {output}"

rule weighted_unifrac:
     """Computes beta diversity weighted unifrac"""
     input:
          table="results/{cohort}/{id}_table.qza",
          tree="results/{cohort}/{id}+fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}+beta_weightedunifrac.qza"
     conda:
          qiime_env
     shell:
          "qiime diversity-lib weighted-unifrac "
          "--i-table {input.table} --i-phylogeny {input.tree} "
          "--o-distance-matrix {output}"

rule unweighted_unifrac:
     """Computes beta diversity weighted unifrac"""
     input:
          table="results/{cohort}/{id}_table.qza",
          tree="results/{cohort}/{id}+fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}+beta_unweightedunifrac.qza"
     conda:
          qiime_env
     shell:
          "qiime diversity-lib unweighted-unifrac "
          "--i-table {input.table} --i-phylogeny {input.tree} "
          "--o-distance-matrix {output}"

rule extract_taxonomy_tsv:
     """converts taxonomy.qza to tsv file"""
     input:
          "results/{cohort}/{id}_taxonomy.qza"
     output:
          "results/{cohort}/{id}_taxonomy.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_dadatable_tsv:
     """converts dada table.qza to tsv"""
     input:
          "results/{cohort}/{id}_table.qza"
     output:
          "results/{cohort}/{id}_table.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_unifrac_tsv:
     """converts unifrac qza artifact to tsv"""
     input:
          "results/{cohort}/{id}unifrac.qza"
     output:
          "results/{cohort}/{id}unifrac.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype distance"

rule merge_dadatable:
     """Merge dada table from 2 datasets"""
     input:
          f1="results/{cohort1}/{cohort1}+{id}_table.qza",
          f2="results/{cohort2}/{cohort2}+{id}_table.qza"
     output:
          "results/{cohort1}-{cohort2}/{cohort1}-{cohort2}+{id}_table.qza"
     conda:
          qiime_env
     shell:
          "qiime feature-table merge "
          "--i-tables {input.f1} "
          "--i-tables {input.f2} "
          "--o-merged-table {output}"

rule merge_dadaseq:
     """Merge dada seq file from 2 datasets"""
     input:
          f1="results/{cohort1}/{cohort1}+{id}_seq.qza",
          f2="results/{cohort2}/{cohort2}+{id}_seq.qza"
     output:
          "results/{cohort1}-{cohort2}/{cohort1}-{cohort2}+{id}_seq.qza"
     conda:
          qiime_env
     shell:
          "qiime feature-table merge-seqs "
          "--i-data {input.f1} "
          "--i-data {input.f2} "
          "--o-merged-data {output}"

rule merge_taxonomy:
     priority:
          1
     input:
          f1="results/{cohort1}/{cohort1}+{id}_taxonomy.qza",
          f2="results/{cohort2}/{cohort2}+{id}_taxonomy.qza"
     output:
          "results/{cohort1}-{cohort2}/{cohort1}-{cohort2}+{id}_taxonomy.qza"
     conda:
          qiime_env
     shell:
          "qiime feature-table merge-taxa "
          "--i-data {input.f1} "
          "--i-data {input.f2} "
          "--o-merged-data {output}"

rule collapse_tax:
     """collapse taxonomy table to genus level"""
     input:
          table="results/{cohort}/{cohort}+{id}_table.qza",
          tax="results/{cohort}/{cohort}+{id}+{etc}_taxonomy.qza"
     output:
          "results/{cohort}/{cohort}+{id}+{etc}+otu_tax.qza"
     conda:
          qiime_env
     shell:
          "qiime taxa collapse --i-table {input.table} "
          "--p-level 6 "
          "--i-taxonomy {input.tax} "
          "--o-collapsed-table {output}"

rule create_metadata_file:
     input:
          "results/{cohort}/{cohort}_manifest.tsv"
     output:
          "results/{cohort}/{cohort}_metadata.tsv"
     conda:
          "envs/other.yml"
     shell:
          "python scripts/create_metadata_file.py -i {input} -o {output}"

rule core_metrics:
     input:
          table="results/{cohort}/{id}_table.qza",
          metadata="results/{cohort}/{cohort}_metadata.tsv"
     output:
          directory("results/{cohort}/{id}+coremetrics/")
     threads:
          30
     conda:
          qiime_env
     shell:
          "qiime diversity core-metrics "
          "--p-sampling-depth {wildcards.r} "
          "--i-table {input.table} "
          "--m-metadata-file {input.metadata} "
          "--output-dir {output} "
          "--p-n-jobs {threads}"

rule alpha_diversity:
     """compoutes alpha diversity"""
     input:
          "results/{cohort}/{id}+rrf-d{r}_table.qza"
     output:
          "results/{cohort}/{id}+rrf-d{r}+alphadiversity.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/alpha_diversity.py --inp {input} "
          "--outp {output}"

rule bridge_alpha_diversity:
     """bridge non-rarefied sample to rrf-d10000"""
     input:
          "results/{cohort}/{id}+rrf-d10000+alphadiversity.tsv"
     output:
          "results/{cohort}/{id}+alphadiversity.tsv"
     shell:
          "cp {input} {output}"

rule beta_diversity:
     """computes non-phylogenetic beta diversity"""
     input:
          "results/{cohort}/{id}+otu_tax.qza"
     output:
          "results/{cohort}/{id}+beta_braycurtis.tsv",
          "results/{cohort}/{id}+beta_jaccard.tsv",
     params:
          "results/{cohort}/{id}+beta.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/beta_diversity.py --inp {input} --outp {params}"

rule biom_to_tsv:
     """converts biom table to tsv"""
     input:
          "results/{cohort}/{id}.biom"
     output:
          "results/{cohort}/{id}_biom.tsv"
     conda:
          qiime_env
     shell:
          "biom convert -i {input} -o {output} --to-tsv"

rule manta:
     """Produces manta output"""
     input:
          tsv="results/{cohort}/{id}+cls-{cls}+otu_tax_biom.tsv",
          taxonpath="db/taxonpath.json",
          names="db/names.json"
     output:
          full="results/{cohort}/{id}+cls-{cls}+manta.csv",
          tax="results/{cohort}/{id}+cls-{cls}+manta_tax.csv",
          abundant="results/{cohort}/{id}+cls-{cls}+manta_abundant_tax.csv",
          sample="results/{cohort}/{id}+cls-{cls}+manta_sample_ids.csv"
     params:
          db=lambda wildcards: "2" if wildcards.cls=="gg" else "1"
     conda:
          "envs/other.yml"
     shell:
          "python scripts/manta.py "
          "-i {input.tsv} "
          "-o {output.full} -x {output.tax} "
          "-s {output.sample} "
          "-a {output.abundant} "
          "-t {input.taxonpath} -n {input.names} "
          "-d {params.db} "

rule manta_alpha_diversity:
     input:
          "results/{cohort}/{id}+alphadiversity.tsv"
     output:
          "results/{cohort}/{id}+manta_alphadiversity.csv"
     conda:
          "envs/other.yml"
     shell:
          "python scripts/prepare_manta_alpha_diversity.py "
          "-i {input} "
          "-o {output}"

rule summary_manta_1:
     """Produces summarized manta results in zipped file"""
     input:
          microbiota="results/{cohort}/{id}+cls-{cls}+manta.csv",
          taxonomy="results/{cohort}/{id}+cls-{cls}+manta_tax.csv",
          dominant_taxon="results/{cohort}/{id}+cls-{cls}+manta_abundant_tax.csv",
          sample_diversity="results/{cohort}/{id}+manta_alphadiversity.csv",
          sample_ids="results/{cohort}/{id}+cls-{cls}+manta_sample_ids.csv"
     output:
          otaxonomy=temp("results/{cohort}/{id}+cls-{cls}/taxonomy.csv"),
          omicrobiota=temp("results/{cohort}/{id}+cls-{cls}/microbiota.csv"),
          osample_diversity=temp("results/{cohort}/{id}+cls-{cls}/sample_diversity.csv"),
          odominant_taxon=temp("results/{cohort}/{id}+cls-{cls}/dominant_taxon.csv"),
          sample_ids=temp("results/{cohort}/{id}+cls-{cls}/sample.csv")
     conda:
          "envs/other.yml"
     shell:
          "cp {input.microbiota} {output.omicrobiota} && "
          "cp {input.taxonomy} {output.otaxonomy} && "
          "cp {input.sample_diversity} {output.osample_diversity} && "
          "cp {input.dominant_taxon} {output.odominant_taxon} &&"
          "cp {input.sample_ids} {output.sample_ids}"

rule summary_manta:
     """Produces summarized manta results in zipped file"""
     input:
          "results/{cohort}/{id}+cls-{cls}/taxonomy.csv",
          "results/{cohort}/{id}+cls-{cls}/microbiota.csv",
          "results/{cohort}/{id}+cls-{cls}/sample_diversity.csv",
          "results/{cohort}/{id}+cls-{cls}/dominant_taxon.csv",
          "results/{cohort}/{id}+cls-{cls}/sample.csv"
     output:
          "results/{cohort}/{id}+cls-{cls}+manta.zip"
     conda:
          "envs/other.yml"
     shell:
          "zip -j {output} {input}"

rule summary:
     """Produces summarized results in zipped file"""
     input:
          "results/{cohort}/{id}+cls-{cls}_taxonomy.tsv",
          "results/{cohort}/{id}_table.tsv",
          "results/{cohort}/{id}+beta_weightedunifrac.tsv",
          "results/{cohort}/{id}+beta_unweightedunifrac.tsv",
          "results/{cohort}/{id}+cls-{cls}+otu_tax.biom",
          "results/{cohort}/{id}+cls-{cls}+otu_tax_biom.tsv",
          "results/{cohort}/{id}+cls-{cls}+phyloseq.RDS",
          "results/{cohort}/{id}+alphadiversity.tsv",
          "results/{cohort}/{id}+cls-{cls}+manta.csv",
          "results/{cohort}/{id}+cls-{cls}+manta_tax.csv",
          "results/{cohort}/{id}_seq.tsv",
          "results/{cohort}/{id}+cls-{cls}+beta_braycurtis.tsv",
          "results/{cohort}/{id}+cls-{cls}+beta_jaccard.tsv",
          "results/{cohort}/{id}+cls-{cls}+manta_abundant_tax.csv"
     output:
          "results/{cohort}/{id}+cls-{cls}.zip"
     conda:
          "envs/other.yml"
     shell:
          "zip -j {output} {input}"

# Start from here, Specifically for JMD...

def get_finished_cohort():
    files = glob.glob("results/*/*+bb-t18+fp-f17-r21+dd+rrf-d10000+cls-silva.zip")
    return [x[8:-4] for x in files]

rule summarize_all_for_jmd:
     """Specifically for JMD"""
     input:
        expand('results/{cohort}+manta.zip', cohort=get_finished_cohort())

rule create_biom_report:
     """generate Excel return report for the collaborators (based on the original biom output)"""
     input:
          adiv_tsv="results/{cohort}/{id}+alphadiversity.tsv",
          biom_tsv="results/{cohort}/{id}+otu_tax_biom.tsv"
     output:
          xlsx="results/{cohort}/report.{id}+biom.xlsx"
     conda:
          "envs/excel.yml"
     shell:
          "python scripts/create_biom_report.py -a {input.adiv_tsv} -d {input.biom_tsv} -o {output.xlsx}"

rule create_report:
     """generate Excel return report for the collaborators (with taxonomy updated)"""
     input:
          manta="results/{cohort}/{id}+cls-{cls}+manta.csv",
          adiv="results/{cohort}/{id}+alphadiversity.tsv",
          taxonpath="db/taxonpath.json",
          names="db/names.json"
     output:
          "results/{cohort}/report.{id}+cls-{cls}.xlsx"
     conda:
          "envs/excel.yml"
     shell:
          "python scripts/create_report.py summarize "
          "-i {input.manta} "
          "-a {input.adiv} "
          "-p {input.taxonpath} -n {input.names} "
          "-o {output} "

# End at here, Specifically for JMD...

ruleorder: merge_taxonomy > taxonomy > manifest
ruleorder: merge_dadatable > rarefy > manifest
ruleorder: export_phyloseq  > extract_biom
ruleorder: extract_biom > make_biom
ruleorder: alpha_diversity > bridge_alpha_diversity
