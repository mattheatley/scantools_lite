import argparse, os, sys, subprocess, re, numpy
from ScanToolsLite import CAPTURE, recode

def outliers(bpm_metrics, output_dir, percentile, relevant_metric, intersecting=False, gff_path=None, plotting=False, plot_range=None, pops=None, table=None, fai_path=None):

    if not os.path.isfile(bpm_metrics): print('\nERROR: Problem finding input metrics; set "-i /PATH/TO/INPUT.METRICS" & run again.\n'); sys.exit(0)

    if plotting: 
        intersecting = True # overide to intersect
        r_script = f'{os.path.dirname(__file__)}/plot_afd.r' # specify r script for plotting
        if not os.path.isfile(table): print('\nERROR: Problem finding variant table; set "-table /PATH/TO/INPUT.TABLE" & run again.\n'); sys.exit(0)
        if not os.path.isfile(r_script): print(f'\nERROR: Problem finding {os.path.basename(r_script)} in "{os.path.dirname(r_script)}".\n'); sys.exit(0)

    if intersecting and not os.path.isfile(gff_path): print('\nERROR: Problem finding gff file; set "-gff /PATH/TO/ANNOTATION.GFF" & run again\n'); sys.exit(0)
    
    wrk_path = os.path.dirname(os.path.abspath(bpm_metrics)) # extract working directory path
    output_dir = f'{wrk_path if not os.path.exists(os.path.dirname(output_dir)) else ""}/{output_dir.strip("/")}' # ensure correct path format

    *_, ws, ms, suffix = os.path.basename(bpm_metrics).split('.') # extract bpm info

    coordinates_dir, = output_dirs = [f'{output_dir}/coordinates'] # specify output directories
    
    if intersecting: 
        window_features_dir, window_cds_dir = [ f'{output_dir}/{label}' for label in [f'intersected/{ws}/{ms}', 'genes'] ] # specify output directories
        output_dirs.extend( [window_features_dir, window_cds_dir] ) #  include additional output directories

    if plotting:
        plot_dir = f'{output_dir}/plots/{ws}/{ms}' # specify output directories
        output_dirs.extend([ plot_dir] ) # include additional output directories

    [ os.makedirs(output_dir, exist_ok=True) for output_dir in output_dirs] # create output directories as required

    window_metrics = [] # store window info

    with open(bpm_metrics, 'r') as metrics:
        
        for i, line in enumerate(metrics,0):

            scaffold, window_start, window_end, *_ = line = line.split('\t') # extract window coordinates
            
            if i==0: metric_index = line.index(relevant_metric) # identify relevant metric column from headers
            elif line and scaffold != 'Genome': # process window  

                metric = line[metric_index] # extract relevant metric

                if not metric == '-': # record window info if metric calculated

                    window_info = [
                        scaffold, 
                        int(window_start), 
                        int(window_end)+1, # convert to zero-based coordinates
                        float(metric) ]
                    window_metrics.append(window_info) 
                    
    threshold = numpy.percentile([ value for *_,value in window_metrics ], percentile) # calculate percentile threshold

    outliers = [ window for *window,value in window_metrics if value >= threshold] # extract outliers above threshold

    outlier_window_coordinates_file = f'{coordinates_dir}/{relevant_metric}.pct{percentile}.{ws}.{ms}.bed' # specify outlier window coordinates output

    with open(f'{outlier_window_coordinates_file}', 'w') as coordinates:
        for window in outliers: print(*window, sep='\t', file=coordinates) # write outlier window coordinates
        


    if intersecting:

        relevant_features = ['gene','mRNA']
        genes_file, mrna_file = [ f'{window_cds_dir}/{ws}.{ms}.{feature}'  for feature in relevant_features ] # specify window feature outputs
        
        overlapping_features, overlapping_cds = 0,0 # store summary feature counts
        
        with open(outlier_window_coordinates_file, 'r') as coordinates, open(genes_file, 'w') as genes, open(mrna_file, 'w') as mrna:
            
            for line in coordinates:

                scaffold, window_start, window_end = window = line.strip(' \n').split('\t') # extract window info
                window_features_file = f'{window_features_dir}/{"_".join(map(str,window))}.gff'

                intersect_region_file = f'{outlier_window_coordinates_file}.tmp' # specify temp bed file

                with open(intersect_region_file, 'w') as tmp: print(*window, sep='\t', file=tmp) # write window coordinates
                subprocess.call(f'bedtools intersect -a {gff_path} -b {intersect_region_file} > {window_features_file}', shell=True) # extract window gene features
                os.remove(intersect_region_file) # remove temp bed file

                window_features, *_ = CAPTURE(f'wc -l {window_features_file}').split(' ') # count window features
                
                if int(window_features) > 1: overlapping_features += 1 # check if any window features other than region info from bedtools

                window_cds = CAPTURE(f'grep -c "CDS" {window_features_file}') # count window CDS features
                
                if int(window_cds) > 0: # check if any window cds features
                
                    overlapping_cds += 1
                
                    relevant_features_info = [ CAPTURE(f'''awk -F "\\t" '{{ if ($3 == "{feature}") {{print}} }}' {window_features_file} | cut -f 9''').split('\n') for feature in relevant_features  ] # filter genes & extract attributes

                    id_format = re.compile('ID=(.+?)(?:;|$)') # specify gff ID search term

                    gene_IDs, mrna_IDs = [ [ id_format.search(info).group(1) for info in feature ] for feature in relevant_features_info ] # extract gene & mrna ids

                    for IDs, output in [ (gene_IDs,genes), (mrna_IDs,mrna) ]:

                        for ID in IDs: print(*window, ID, sep='\t', file=output, flush=True) # write window coordinates & feature ids

                else:
                    
                    gene_IDs = "NONCODING"

                if plotting:
                    
                    plot_suffixes = [
                        '.afd',
                        f'.{relevant_metric.lower()}',
                        '.gff',
                        f'{"_".join(gene_IDs)}.png' ]
                    plot_region_afd_file, plot_region_metric_file, plot_region_feature_file, plot_region_file = [ f'{plot_dir}/{scaffold}_{window_start}_{int(window_end)-1}{suffix}' for suffix in plot_suffixes ] # specify plot outputs

                    scaffold_end = int(CAPTURE(f"grep -w '{scaffold}' {fai_path} | cut -f 2")) # extract scaffold end

                    plot_centre = int( int(window_start) + (int(window_end) - int(window_start)) / 2 ) # calculate plot region centre 
                    plot_start, plot_end = int( plot_centre - (plot_range / 2) ), int( plot_centre + (plot_range / 2) ) # calculate plot region start/end

                    offset = 5000 # specify additonal plot parameters
                    offset_start, offset_end = plot_start-offset, plot_end+offset # offset plot regions to extract extra data

                    if plot_start < 0: plot_end += abs(plot_start); plot_start = 0 # adjust plot region to begin at zero as required
                    if offset_start < 0: offset_start = 0 # adjust plot data region as required

                    extract_region_cmd = f'''grep -w "{scaffold}" {{0}} | awk -F "\\t" '{{{{ if (${{1}} == "{scaffold}" && ${{2}} >= {offset_start} && ${{3}} <= {offset_end}) {{{{print}}}} }}}}' '''

                    plot_region_sites = CAPTURE(extract_region_cmd.format( table, 2,3,3 )).split('\n') # extract plot region sites
                    
                    plot_region_metrics = CAPTURE(extract_region_cmd.format( bpm_metrics, 1,2,3) ).split('\n') # extract plot region metrics    
                    
                    with open(intersect_region_file, 'w') as tmp: print(scaffold, offset_start, offset_end, sep='\t', file=tmp) # write plot region coordinates
                    plot_region_features = CAPTURE(f'bedtools intersect -a {gff_path} -b {intersect_region_file}').split('\n') # extract plot region gene features
                    os.remove(intersect_region_file) # remove temp bed file

                    minuend, subtrahend = pops # specify population order in afd calculation
                    afd_order = f'{minuend} - {subtrahend}' # specify y-axis (afd) label

                    plot_region_frequencies = {} # stored window allele frequencies

                    for line in plot_region_sites:

                        individual_table_basename, scaffold, position, ploidy, gts = line = recode(line) # extract site info

                        if not position in plot_region_frequencies: plot_region_frequencies[position] = {} # create position entry as required

                        prefix_matched = re.match(f'^({re.escape(minuend)})\.(.*)\.', individual_table_basename) # check if minuend or subtrahend population site

                        if gts: # process genotypes
                            n = float(len(gts)) # calculate individual number
                            an = float(ploidy) * n # re-calculate allele number
                            ac = sum(gts) # calculate allele count
                            p = round(ac/an, 3) # calculate allele frequency
                            plot_region_frequencies[position].update({ minuend if prefix_matched else subtrahend: p  }) # record allele frequency
                    
                    with open(plot_region_afd_file, 'w') as afds, open(plot_region_metric_file, 'w') as metrics, open(plot_region_feature_file, 'w') as features:
                        
                        print('scaffold', 'position', 'afd', sep='\t', file=afds)
                        for position, p in plot_region_frequencies.items():
                            if len(p) > 1:
                                afd = p[minuend]  -  p[subtrahend] # calculate afd
                                print(scaffold, position, afd, sep='\t', file=afds) # write afd

                        print('scaffold', 'start', 'end', 'middle', relevant_metric.lower(), sep='\t', file=metrics)
                        print(scaffold, 0, 0, 0, 0, sep='\t', file=metrics) # write plot region metric start tether
                        for i, line in enumerate(plot_region_metrics,0):
                            scaffold, start, end, *_ = line = line.split('\t') # extract window info
                            metric = line[metric_index] # extract relevant metric
                            mid_window = ((int(end)-int(start))/2)+int(start) # calculate mid-window coordinate
                            if not metric == '-': print(scaffold, start, end, mid_window, metric, sep='\t', file=metrics) # write plot region metric coordinates
                        print(scaffold, 0, 0, end, 0, sep='\t', file=metrics) # write plot region metric end tether
                        
                        print('scaffold', 'description', 'start', 'end', 'orientation', sep='\t', file=features)
                        for line in plot_region_features:
                            scaffold, source, feature, start, end, score, strand, *_ = line = line.split('\t') # extract feature info
                            print(scaffold, feature, start, end, strand, sep='\t', file=features) # write plot region features

                    subprocess.call(f'Rscript {r_script} {relevant_metric} "{afd_order}" {plot_region_afd_file} {plot_region_metric_file} {plot_region_feature_file} {plot_region_file} {scaffold} {window_start} {int(window_end)-1} {plot_start} {plot_end}', shell=True) # plot window data



    # summary
    print('total', 'features', 'cds')
    print(len(outliers), overlapping_features, overlapping_cds) # write sumnmnary feature counts

    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='</path/to/bpm.metrics>', required=True)
    parser.add_argument('-o', type=str, metavar='<output.directory>', required=True)
    parser.add_argument('-percentile', type=int, metavar='<number>', required=True, help='Outlier percentile threshold')
    parser.add_argument('-metric', type=str, metavar='<name>', required=True, help='Metric to examine')
    parser.add_argument('-intersect', action='store_true', help='Intersect regions with gff')
    parser.add_argument('-gff', type=str, metavar='</path/to/annotation.gff>', default=None, help='Genome annotation file')
    parser.add_argument('-plot', action='store_true', help='Plot regions that include coding')
    parser.add_argument('-range', type=int, metavar='<number>', default=50000, help='Plot range (bp)')
    parser.add_argument('-pops', type=str, nargs=2, metavar='<Minuend Population Name> <Subtrahend Population Name>', default=None, help='Minuend & Subtrahend populations for AFD Calculation')
    parser.add_argument('-table', type=str, metavar='</path/to/input.table>', default=None, help='Variant table')
    parser.add_argument('-fai', type=str, metavar='</path/to/reference.fasta.fai>', default=None, help='Genome index file')

    input_metrics, output_dir, percentile, metric, to_intersect, gff_path, to_plot, plot_range, pops, variant_table, fai_path = vars(parser.parse_args()).values()

    outliers(bpm_metrics=input_metrics, output_dir=output_dir, percentile=percentile, relevant_metric=metric, 
    intersecting=to_intersect, gff_path=gff_path, # intersecting parameters
    plotting=to_plot, plot_range=plot_range, pops=pops, table=variant_table, fai_path=fai_path) # plotting parameters

# testing; -ws 10000 -ms 10 -ni 5
