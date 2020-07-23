version 1.0

workflow consensus_genome {
    input {
        # Required parameters
        File fastqs_0
        File fastqs_1

        String docker_image_id
        File ercc_fasta
        File kraken2_db_tar_gz
        File primer_bed
        File ref_fasta
        File ref_host

        # Sample name : include in tags and files
        String sample

        # Optional prefix to add to the output filenames
        String prefix = ""

        # Step parameters

        # all of the following seem to be parameters in params.* of the nextflow pipeline;
        # could they be extracted to a separate params file? put in .json input to wdl?
        Boolean trim_adapters = true

        Float   intrahost_min_frac = 0.5
        String  intrahost_ploidy   = "2"
        Boolean intrahost_variants = true

        Float ivarFreqThreshold = 0.6
        Int   ivarQualTreshold  = 20
        Int   minDepth          = 10
        Int   mpileupDepth      = 10000

        String no_reads_quast = false

        # assumes about 20 mutations between 2 random samples
        # (this is an overestimate to increase sensitivity)
        Float bcftoolsCallTheta = 0.0006

        # Dummy values - required by SFN interface
        String s3_wd_uri
        String dag_branch
    }

    call removeHost {
        input:
            prefix = prefix,
            fastqs_0 = fastqs_0,
            fastqs_1 = fastqs_1,
            ref_host = ref_host,
            docker_image_id = docker_image_id
    }

    call quantifyERCCs {
        input:
            prefix = prefix,
            fastqs_0 = removeHost.host_removed_fastqs_0,
            fastqs_1 = removeHost.host_removed_fastqs_1,
            ercc_fasta = ercc_fasta,
            docker_image_id = docker_image_id
    }

    call filterReads {
        input:
            prefix = prefix,
            fastqs_0 = removeHost.host_removed_fastqs_0,
            fastqs_1 = removeHost.host_removed_fastqs_1,
            ref_fasta = ref_fasta,
            kraken2_db_tar_gz = kraken2_db_tar_gz,
            docker_image_id = docker_image_id
    }

    # skip remaining steps if no SC2 reads were found (28 = size of an empty bgzip file)
    if (size(filterReads.filtered_fastqs_0) > 28) {
        if (trim_adapters) {
            call trimReads {
                input:
                    fastqs_0 = filterReads.filtered_fastqs_0,
                    fastqs_1 = filterReads.filtered_fastqs_1,
                    docker_image_id = docker_image_id
            }
        }

        call alignReads {
            input:
                prefix = prefix,
                sample = sample,
                # use trimReads output if we ran it; otherwise fall back to filterReads output
                fastqs_0 = select_first([trimReads.trimmed_fastqs_0, filterReads.filtered_fastqs_0]),
                fastqs_1 = select_first([trimReads.trimmed_fastqs_1, filterReads.filtered_fastqs_1]),
                ref_fasta = ref_fasta,
                docker_image_id = docker_image_id
        }

        call trimPrimers {
            input:
                prefix = prefix,
                alignments = alignReads.alignments,
                primer_bed = primer_bed,
                docker_image_id = docker_image_id
        }

        if (intrahost_variants) {
            call intrahostVariants {
                input:
                    bam = trimPrimers.trimmed_bam_ch,
                    ref_fasta = ref_fasta,
                    intrahost_ploidy = intrahost_ploidy,
                    intrahost_min_frac = intrahost_min_frac,
                    docker_image_id = docker_image_id
            }
        }

        call makeConsensus {
            input:
                prefix = prefix,
                sample = sample,
                bam = trimPrimers.trimmed_bam_ch,
                ivarFreqThreshold = ivarFreqThreshold,
                minDepth = minDepth,
                ivarQualThreshold = ivarQualTreshold,
                mpileupDepth = mpileupDepth,
                docker_image_id = docker_image_id
        }

        call quast {
            input:
                prefix = prefix,
                assembly = makeConsensus.consensus_fa,
                bam = trimPrimers.trimmed_bam_ch,
                # use trimReads output if we ran it; otherwise fall back to filterReads output
                fastqs_0 = select_first([trimReads.trimmed_fastqs_0, filterReads.filtered_fastqs_0]),
                fastqs_1 = select_first([trimReads.trimmed_fastqs_1, filterReads.filtered_fastqs_1]),
                ref_fasta = ref_fasta,
                no_reads_quast = no_reads_quast,
                docker_image_id = docker_image_id
        }

        call realignConsensus {
            input:
                prefix = prefix,
                sample = sample,
                realign_fa = makeConsensus.consensus_fa,
                ref_fasta = ref_fasta,
                docker_image_id = docker_image_id
        }
        call callVariants {
            input:
                prefix = prefix,
                call_variants_bam = trimPrimers.trimmed_bam_ch,
                ref_fasta = ref_fasta,
                bcftoolsCallTheta = bcftoolsCallTheta,
                minDepth = minDepth,
                docker_image_id = docker_image_id
        }

        call computeStats {
            input:
                prefix = prefix,
                sample = sample,
                cleaned_bam = trimPrimers.trimmed_bam_ch,
                assembly = makeConsensus.consensus_fa,
                ercc_stats = quantifyERCCs.ercc_out,
                vcf = callVariants.variants_ch,
                fastqs_0 = fastqs_0,
                fastqs_1 = fastqs_1,
                ref_host = ref_host,
                docker_image_id = docker_image_id
        }
    }

    call zipOutputs {
        input:
            prefix = prefix,
            outputFiles = [
                removeHost.host_removed_fastqs_0,
                removeHost.host_removed_fastqs_1
            ]
    }
}

# TODO: task to validate input

task removeHost {
    # TODO: process errors: no reads left
    input {
        String prefix
        File fastqs_0
        File fastqs_1
        File ref_host

        String docker_image_id
    }

    command <<<
        export CPUS=`nproc --all`
        minimap2 -t $CPUS -ax sr ~{ref_host} ~{fastqs_0} ~{fastqs_1} | \
        samtools view -@ $CPUS -b -f 4 | \
        samtools fastq -@ $CPUS -1 "~{prefix}no_host_1.fq.gz" -2 "~{prefix}no_host_2.fq.gz" -0 /dev/null -s /dev/null -n -c 6 -
    >>>

    output {
        File host_removed_fastqs_0 = "~{prefix}no_host_1.fq.gz"
        File host_removed_fastqs_1 = "~{prefix}no_host_2.fq.gz"
    }

    runtime {
        docker: docker_image_id
    }
}

task quantifyERCCs {
    input {
        String prefix
        File fastqs_0
        File fastqs_1
        File ercc_fasta

        String docker_image_id
    }

    command <<<
        minimap2 -ax sr ~{ercc_fasta} ~{fastqs_0} ~{fastqs_1} | samtools view -bo ercc_mapped.bam
        samtools stats ercc_mapped.bam > "~{prefix}ercc_stats.txt"
    >>>

    output {
        File ercc_out = "~{prefix}ercc_stats.txt"
    }

    runtime {
        docker: docker_image_id
    }
}

task filterReads {
    # TODO: process errors: no reads left
    input {
        String prefix
        String taxid = "2697049"
        File fastqs_0
        File fastqs_1

        File ref_fasta
        File kraken2_db_tar_gz

        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        export TMPDIR=${TMPDIR:-/tmp}
        export CPUS=`nproc --all`

        minimap2 -ax sr -t $CPUS "~{ref_fasta}" ~{fastqs_0} ~{fastqs_1} \
            | samtools sort -@ $CPUS -n -O bam -o "${TMPDIR}/mapped.bam"
        samtools fastq -@ $CPUS -G 12 -1 "${TMPDIR}/paired1.fq.gz" -2 "${TMPDIR}/paired2.fq.gz" \
            -0 /dev/null -s /dev/null -n -c 6 "${TMPDIR}/mapped.bam"

        paired1size=$(stat --printf="%s" "${TMPDIR}/paired1.fq.gz")
        if (( paired1size > 28 )); then
            mkdir "${TMPDIR}/kraken_db"
            tar -xzv -C "${TMPDIR}/kraken_db" -f "~{kraken2_db_tar_gz}"
            kraken2 --db "${TMPDIR}"/kraken_db/* \
                --threads $CPUS \
                --report "${TMPDIR}/~{prefix}kraken2_report.txt" \
                --classified-out "${TMPDIR}/~{prefix}classified#.fq" \
                --output - \
                --memory-mapping --gzip-compressed --paired \
                "${TMPDIR}/paired1.fq.gz" "${TMPDIR}/paired2.fq.gz"

            grep --no-group-separator -A3 "kraken:taxid|~{taxid}" \
                "${TMPDIR}/~{prefix}classified_1.fq" \
                > "${TMPDIR}/~{prefix}covid_1.fq" || [[ \$? == 1 ]]
            grep --no-group-separator -A3 "kraken:taxid|~{taxid}" \
                "${TMPDIR}/~{prefix}classified_2.fq" \
                > "${TMPDIR}/~{prefix}covid_2.fq" || [[ \$? == 1 ]]
            bgzip -@ $CPUS -c "${TMPDIR}/~{prefix}covid_1.fq" > "~{prefix}covid_1.fq.gz"
            bgzip -@ $CPUS -c "${TMPDIR}/~{prefix}covid_2.fq" > "~{prefix}covid_2.fq.gz"
        else
            mv "${TMPDIR}/paired1.fq.gz" "~{prefix}covid_1.fq.gz"
            mv "${TMPDIR}/paired2.fq.gz" "~{prefix}covid_2.fq.gz"
        fi
    >>>

    output {
        File filtered_fastqs_0 = "~{prefix}covid_1.fq.gz"
        File filtered_fastqs_1 = "~{prefix}covid_2.fq.gz"
    }

    runtime {
        docker: docker_image_id
    }
}

task trimReads {
    # TODO: process errors: no reads left
    input {
        File fastqs_0
        File fastqs_1

        String docker_image_id
    }

    command {
        trim_galore --fastqc --paired "~{fastqs_0}" "~{fastqs_1}"
    }

    output {
        File trimmed_fastqs_0 = "~{basename(fastqs_0, '.fq.gz')}_val_1.fq.gz"
        File trimmed_fastqs_1 = "~{basename(fastqs_1, '.fq.gz')}_val_2.fq.gz"
    }

    runtime {
        docker: docker_image_id
    }
}

task alignReads {
    # TODO: process errors: no reads left (unlikely)
    input {
        String prefix
        String sample
        File fastqs_0
        File fastqs_1
        File ref_fasta

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        export CORES=`nproc --all`
        minimap2 -ax sr -t $CORES -R '@RG\tID:~{sample}\tSM:~{sample}' "~{ref_fasta}" "~{fastqs_0}" "~{fastqs_1}" \
            | samtools sort -@ $CORES -O bam -o "~{prefix}aligned_reads.bam"
    >>>

    output {
        File alignments = "~{prefix}aligned_reads.bam"
    }

    runtime {
        docker: docker_image_id
    }
}

task trimPrimers {
    input {
        String prefix
        File alignments
        File primer_bed

        String docker_image_id

        Int samQualThreshold = 20
    }

    command <<<
        samtools view -F4 -q "~{samQualThreshold}" -o ivar.bam "~{alignments}"
        samtools index ivar.bam
        ivar trim -e -i ivar.bam -b "~{primer_bed}" -p ivar.out
        samtools sort -O bam -o "~{prefix}primertrimmed.bam" ivar.out.bam
        samtools index "~{prefix}primertrimmed.bam"
    >>>

    output {
        File trimmed_bam_ch = "~{prefix}primertrimmed.bam"
        File trimmed_bam_bai = "~{prefix}primertrimmed.bam.bai"
    }

    runtime {
        docker: docker_image_id
    }

}

task intrahostVariants {
    input {
        File bam
        File ref_fasta

        String intrahost_ploidy
        Float intrahost_min_frac

        String docker_image_id
    }

    command <<<
        export CPUS=`nproc --all`
        ls "~{bam}" | xargs -I % samtools index %
        samtools faidx "~{ref_fasta}"
        freebayes-parallel <(fasta_generate_regions.py "~{ref_fasta}.fai" 1000) $CPUS \
            --ploidy "~{intrahost_ploidy}" \
            --min-alternate-fraction "~{intrahost_min_frac}" \
            -f "~{ref_fasta}" "~{bam}" |
            bcftools view -g het \
            > "intrahost-variants.ploidy~{intrahost_ploidy}-minfrac~{intrahost_min_frac}.vcf"
    >>>

    output {
        File intrahost_vars = "intrahost-variants.ploidy~{intrahost_ploidy}-minfrac~{intrahost_min_frac}.vcf"
    }

    runtime {
        docker: docker_image_id
    }
}

task makeConsensus {
    input {
        String prefix
        String sample
        File bam

        Float ivarFreqThreshold
        Int minDepth
        Int ivarQualThreshold
        Int mpileupDepth

        String docker_image_id
    }

    command <<<
        samtools index ${bam}
        samtools mpileup -A -d "~{mpileupDepth}" -Q0 "~{bam}" | ivar consensus -q "~{ivarQualThreshold}" -t "~{ivarFreqThreshold}" -m "~{minDepth}" -n N -p "~{prefix}primertrimmed.consensus"
        echo ">""~{sample}" > "~{prefix}consensus.fa"
        seqtk seq -l 50 "~{prefix}primertrimmed.consensus.fa" | tail -n +2 >> "~{prefix}consensus.fa"

    >>>

    output {
        File consensus_fa = "~{prefix}consensus.fa"
    }
    runtime {
        docker: docker_image_id
    }
}

task quast {
    input {
        String prefix
        File assembly   # same as consensus_fa
        File bam
        File fastqs_0
        File fastqs_1
        File ref_fasta

        String no_reads_quast

        Int threads = 4

        String docker_image_id

    }

    command <<<
        export CPUS=`nproc --all`
        a=`cat "~{assembly}" | wc -l`

        if [ $a -ne 0 ]; then
            if [ ~{no_reads_quast} = true ]; then
                quast --min-contig 0 -o quast -r "~{ref_fasta}" -t $CPUS --ref-bam "~{bam}" "~{assembly}"
            else
                quast --min-contig 0 -o quast -r "~{ref_fasta}" -t $CPUS --ref-bam "~{bam}" "~{assembly}" -1 "~{fastqs_0}" -2 "~{fastqs_0}"
            fi
        else
            mkdir quast
            echo "quast folder is empty" > "quast/~{prefix}quast.txt"
        fi
    >>>

    output {
        Array[File] quast_dir = glob("quast/*")
    }

    runtime {
        docker: docker_image_id
    }
}

task realignConsensus {
    input {
        String prefix
        String sample
        File realign_fa   # same as consensus_fa
        File ref_fasta

        String docker_image_id
    }

    command <<<
        minimap2 -ax asm5 -R '@RG\tID:~{sample}\tSM:~{sample}' "~{ref_fasta}" "~{realign_fa}" | samtools sort -O bam -o "~{prefix}realigned.bam"
        samtools index "~{prefix}realigned.bam"
    >>>

    output {
        File realigned_bam = "~{prefix}realigned.bam"
        File realigned_bai = "~{prefix}realigned.bam.bai"
    }

    runtime {
        docker: docker_image_id
    }
}

task callVariants {
    input {
        String prefix
        File call_variants_bam  # same as primertrimmed_bam produced by realignConsensus
        File ref_fasta
        Float bcftoolsCallTheta
        Int minDepth

        String docker_image_id
    }

    command <<<
        bcftools mpileup -a FORMAT/AD -f "~{ref_fasta}" "~{call_variants_bam}" | bcftools call --ploidy 1 -m -P "~{bcftoolsCallTheta}" -v - | bcftools view -i 'DP>=~{minDepth}' > "~{prefix}variants.vcf"
        bgzip "~{prefix}variants.vcf"
        tabix "~{prefix}variants.vcf.gz"
        bcftools stats "~{prefix}variants.vcf.gz" > "~{prefix}bcftools_stats.txt"
    >>>

    output {
        File variants_ch = "~{prefix}variants.vcf.gz"
        File bcftools_stats_ch = "~{prefix}bcftools_stats.txt"
    }

    runtime {
        docker: docker_image_id
    }
}

task computeStats {
    input {
        String prefix
        String sample
        File cleaned_bam
        File assembly
        File ercc_stats
        File vcf

        File fastqs_0
        File fastqs_1
        File ref_host

        String docker_image_id
    }

    command <<<
        samtools index "~{cleaned_bam}"
        samtools stats "~{cleaned_bam}" > "~{prefix}samtools_stats.txt"

        python <<CODE

        import argparse
        import collections
        import json
        import re
        import subprocess
        import pysam
        from Bio import SeqIO
        import numpy as np
        from matplotlib import pyplot as plt
        import seaborn as sns

        stats = {"sample_name": "~{sample}"}

        samfile = pysam.AlignmentFile("~{cleaned_bam}", "rb")
        ref_len, = samfile.lengths
        depths = [0] * ref_len
        for column in samfile.pileup():
            depths[column.reference_pos] = column.nsegments
        depths = np.array(depths)

        stats["depth_avg"] = depths.mean()
        stats["depth_q.25"] = np.quantile(depths, .25)
        stats["depth_q.5"] = np.quantile(depths, .5)
        stats["depth_q.75"] = np.quantile(depths, .75)
        stats["depth_frac_above_10x"] = (depths >= 10).mean()
        stats["depth_frac_above_25x"] = (depths >= 30).mean()
        stats["depth_frac_above_50x"] = (depths >= 30).mean()
        stats["depth_frac_above_100x"] = (depths >= 100).mean()

        ax = sns.lineplot(np.arange(1, ref_len+1), depths)
        ax.set_title("~{sample}")
        ax.set(xlabel="position", ylabel="depth")
        plt.yscale("symlog")
        plt.savefig("~{prefix}depths.png")

        seq, = SeqIO.parse("~{assembly}", "fasta")
        stats["allele_counts"] = dict(collections.Counter(str(seq.seq)))

        fq_lines = subprocess.run(" ".join(["zcat"] + ["~{fastqs_0}", "~{fastqs_1}"]) + " | wc -l", shell=True, stdout=subprocess.PIPE).stdout
        stats["total_reads"] = int(int(fq_lines) / 4)

        with open("~{prefix}samtools_stats.txt") as f:
            sam_stats_re = re.compile(r"SN\s+([^\s].*):\s+(\d+)")
            for line in f:
                matched = sam_stats_re.match(line)
                if matched:
                    if matched.group(1) == "reads mapped":
                        stats["mapped_reads"] = int(matched.group(2))
                    elif matched.group(1) == "reads mapped and paired":
                        stats["mapped_paired"] = int(matched.group(2))
                    elif matched.group(1) == "inward oriented pairs":
                        stats["paired_inward"] = int(matched.group(2)) * 2
                    elif matched.group(1) == "outward oriented pairs":
                        stats["paired_outward"] = int(matched.group(2)) * 2
                    elif matched.group(1) == "pairs with other orientation":
                        stats["paired_other_orientation"] = int(matched.group(2)) * 2

        with open("~{ercc_stats}") as f:
            ercc_stats_re = re.compile(r"SN\s+([^\s].*):\s+(\d+)")
            for line in f:
                matched = ercc_stats_re.match(line)
                if matched:
                    if matched.group(1) == "reads mapped":
                        stats["ercc_mapped_reads"] = int(matched.group(2))
                    elif matched.group(1) == "reads mapped and paired":
                        stats["ercc_mapped_paired"] = int(matched.group(2))

        def countVCF(vcf_file, snpcol, mnpcol, indelcol, statsdict):
            vcf = pysam.VariantFile(vcf_file)
            statsdict[snpcol] = 0
            statsdict[mnpcol] = 0
            statsdict[indelcol] = 0
            for rec in vcf.fetch():
                allele_lens = set([len(a) for a in [rec.ref] + list(rec.alts)])
                if len(allele_lens) > 1:
                    statsdict[indelcol] += 1
                else:
                    l, = allele_lens
                    if l == 1:
                        statsdict[snpcol] += 1
                    else:
                        statsdict[mnpcol] += 1
            return statsdict

        stats = {**stats, **countVCF("~{vcf}", 'ref_snps', 'ref_mnps', 'ref_indels', stats)}

        allele_counts = stats["allele_counts"]
        stats["n_actg"] = sum(v for k, v in allele_counts.items() if k in "ACTGU")
        if "N" in allele_counts.keys():
            stats["n_missing"] = allele_counts["N"]
        else:
            stats["n_missing"] = 0
        if "-" in allele_counts.keys():
            stats["n_gap"] = allele_counts["-"]
        else:
            stats["n_gap"] = 0
        stats["n_ambiguous"] = sum(v for k, v in allele_counts.items() if k not in "ACTGUN-")

        with open("~{prefix}stats.json", "w") as f:
            json.dump(stats, f, indent=2)

        CODE
    >>>

    output {
        File sam_stats = "~{prefix}samtools_stats.txt"
        File depths_fig = "~{prefix}depths.png"
        File output_stats = "~{prefix}stats.json"
    }

    runtime {
        docker: docker_image_id
    }
}

task zipOutputs {
    input {
        String prefix
        Array[File] outputFiles
    }

    command <<<
        mkdir "${TMPDIR}/outputs"
        cp ${sep=' ' outputFiles} "${TMPDIR}/outputs"
        tar -czvf "~{prefix}outputs.tar.gz" -C "${TMPDIR}/outputs" .
    >>>

    output {
        File outputs = "~{prefix}outputs.tar.gz"
    }
}