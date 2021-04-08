from ScanToolsLite import scantools # run_scan.py must be located in ScanToolsLite directory

working_dir = '/gpfs01/home/mbzmch/genome_scans/brassica_fruticulosa_test' # location of associated inputs & outputs
vcf_dir = f'{working_dir}/test.filtered.F4' # location of vcf containgin sample names &/or to be split
table_dir = f'{working_dir}/VCF_test.filtered.F4_DP8.M0.5' # location of split vcf (created by scantools)
scratch_path = f'{working_dir}/scratch' # temporary storage

ref_path = f'{working_dir}/ref.gen/fruticulosa.fasta' # location of reference genome (N.B. GATK requires a ".fasta" or ".fa" suffix)
fai_path = f'{working_dir}/ref.gen/fruticulosa.fasta.fai' # location of reference genome index
gff_path = f'{working_dir}/ref.gen/models_aed0.5.gff' # location of genome annotation
gatk_path = '/gpfs01/home/mbzmch/.conda/envs/gatk3_env/gatk3/GenomeAnalysisTK.jar'

# create scantools object
project = scantools(wrk_dir=working_dir, vcf_dir=vcf_dir)

project.combinePops(to_combine=['ES','RO','BL','TS'], new_pop='coastal')
project.combinePops(to_combine=['TR','Pi','PA'], new_pop='inland') # SC REMOVED

# split vcfs by population
project.splitVCFs(ref_path=ref_path, vcf_dir=vcf_dir, min_dp=8, mffg=0.5, pops=['coastal','inland'], gatk_path=None, conda_env='gatk3_env')

# calculate within population metrics
downsample_size  = min([size for pop,size in project.sizes.items() if pop in ['coastal','inland'] ]) - 1
project.calcwpm(table_dir=table_dir, window_size=10000, min_snps=10, sample_size=downsample_size, pops=['coastal','inland'], show_excluded=False)

# calculate between population metrics
project.calcbpm(table_dir=table_dir, window_size=10000, min_snps=10, pops=['coastal','inland'], show_excluded=False)

# perform genome scans
project.findoutliers(metric_dir=table_dir, percentile=99, metric=None, pops=['coastal','inland'], metrics_prefix=None, to_intersect=True, gff_path=gff_path, to_plot=True, table_dir=table_dir, fai_path=fai_path)