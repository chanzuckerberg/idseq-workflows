version 1.0
import "host_filter.wdl" as stage1
import "non_host_alignment.wdl" as stage2

workflow idseq_short_read_mngs {
    input {
        String docker_image_id
        File non_host_gsnap_index
        File non_host_rapsearch2_index
        String s3_wd_uri = ""
    }
    call stage1.idseq_host_filter as host_filter {
        input:
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri
    }

    call stage2.idseq_non_host_alignment as non_host_alignment {
        input:
        host_filter_out_gsnap_filter_1_fa = host_filter.gsnap_filter_out_gsnap_filter_1_fa,
        host_filter_out_gsnap_filter_2_fa = host_filter.gsnap_filter_out_gsnap_filter_2_fa,
        host_filter_out_gsnap_filter_merged_fa = host_filter.gsnap_filter_out_gsnap_filter_merged_fa,
        cdhitdup_cluster_sizes_cdhitdup_cluster_sizes_tsv = host_filter.cdhitdup_out_cdhitdup_cluster_sizes_tsv,
        cdhitdup_out_dedup1_fa_clstr = host_filter.cdhitdup_out_dedup1_fa_clstr,
        cdhitdup_out_dedup1_fa = host_filter.cdhitdup_out_dedup1_fa,
        local_gsnap_index = non_host_gsnap_index,
        local_rapsearch2_index = non_host_rapsearch2_index,
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri
    }
}
