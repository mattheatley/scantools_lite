import argparse, os, sys, subprocess, numpy
from ScanToolsLite import recode

def PROGRESS(count, total):
    filled = int(round(50 * count / float(total)))
    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled + '-' * (50 - filled)
    sys.stdout.write(f'[{bar}] {percents}% of {int(total):,} sites complete\r')
    sys.stdout.flush()

def calcWPM(variant_table, wpm_metrics, window_size, minimum_snps, ind_num, show_excluded=False, show_progress=False):

    if ind_num <= 3: print(f'ERROR: Unable to calculate within-population metrics on so few individuals; set "-ni" above 3 & run again'); sys.exit(0)

    if not os.path.isfile(variant_table): print('ERROR: Problem finding input file; set "-i /PATH/TO/INPUT.TABLE" & run again'); sys.exit(0)
    
    if not os.path.isdir(os.path.dirname(wpm_metrics)): print('ERROR: Problem finding output directory; set "-o /PATH/TO/OUTPUT.METRICS" & run again.'); sys.exit(0)

    with open(wpm_metrics, 'w') as metrics:
        
        headers = ['ploidy','inds','scaff','start','end','site_count','snp_count','singleton_count','avg_freq','avg_Ehet','Diversity','ThetaW','ThetaW_per_nuc','Pi','ThetaH','ThetaL','D','H','E'] # specify headers
        print(*headers, sep='\t', file=metrics, flush=True) # write metrics header

        sites_to_process, *_ = subprocess.run(f'wc -l {variant_table}', shell=True, capture_output=True).stdout.decode('utf-8').strip(' \n').split(' ') # calculate total sites in table

        windows_processed, windows_excluded = 0, 0 # establish GENOME window counts

        scanning, EOF, summary = True, False, False # switches (see below)

        while scanning:

            with open(variant_table, 'r') as table:

                for i, line in enumerate(table):

                    setup = i == 0 and not EOF # first pass of first loop
                    
                    if EOF: pass # data already processed; skip
                    else: # proceed processing data

                        *_, scaffold, position, ploidy, gts = recode(line) # extract site info
                        
                        gts = [ int(gt) for gt in gts if gt != "-9"] # filter no-call genotypes
                        an = float(ind_num *  float(ploidy)) # downsampled population allele number

                        new_scaffold = setup or scaffold != previous_scaffold # register scaffold change
                        moved_window = setup or int(position) > window_end # register new window


                    # CALCULATE WINDOW/GENOME METRICS AS APPROPRIATE
                    if setup: pass # not relevant yet; skip
                    elif new_scaffold or moved_window or EOF: # calculate metrics if moved scaffold / moved window / final window / summary  

                        region = (previous_scaffold, window_start, window_end) if not summary else ('Genome','-9','-9') # specify relevant region
                    
                        sites, *_ = sites_snps_singletons = (site_count, snp_count, singleton_count) if not summary else (SITE_count, SNP_count, SINGLETON_count) # specify relevant counts

                        region_info = [ploidy, ind_num, *region, *sites_snps_singletons] # specify region info
                        
                        if not summary and snp_count < minimum_snps: pass # not summary & too few snps to process
                        else: # proceed with calculating metrics

                            Afs = afs if not summary else AFS # specify relevant allele frequency spectrum

                            p_i, Ehet = [ numpy.mean(metrics) for metrics in [p_i, Ehet] ] if not summary else ('-9','-9') # specify allele frequency & ehet as appropriate
                            Pi,h,L = 0.0, 0.0, 0.0
                            S = float(sum(Afs[1:-1]))
                            S_per_nuc = S / float(sites) # per nucleotide
                            W = S / aw
                            W_per_nuc = S_per_nuc / aw # per nucleotide
                            W2 = S * (S - 1) / ((aw**2) + a2)

                            for j in range(1, int(an)):
                                Pi += Afs[j] * j * (an - j)
                                h += Afs[j] * (j**2)
                                L += Afs[j] * j

                            Pi = 2 * Pi / (an * (an - 1))
                            h = 2 * h / (an * (an - 1))
                            L = L / (an - 1)
                            div = Pi / float(sites)

                            varPi_W = e1 * S + e2 * S * (S - 1)
                            varPi_L1 = (((an - 2) / (6 * an - 6)) * W)
                            varPi_L2 = (((18 * an**2 * (3 * an + 2) * bw) - (88 * an**3 + 9 * an**2 - 13 * an + 6)) / (9 * an * (an - 1)**2) * W2)
                            varL_W = (((an / (2 * an - 2)) - 1 / aw) * W) + ((a2 / aw**2 + (2 * (an / (an - 1))**2) * a2 - 2 * (an * a2 - an + 1) / ((an - 1) * aw) - (3 * an + 1) / (an - 1)) * W2)
                            varPi_L = varPi_L1 + varPi_L2
                            
                            try:   
                                if not summary: windows_processed += 1 # record processed window

                                D = (Pi - W) / (varPi_W**0.5)
                                H = (Pi - L) / (varPi_L**0.5)  # Normalized H according to Zeng 2006
                                E = (L - W) / (varL_W**0.5)
                                                            
                                region_info.extend([p_i, Ehet, div, W, W_per_nuc, Pi, h, L, D, H, E])

                            except ZeroDivisionError: pass

                        if len(region_info) != len(headers):
                            if not summary: windows_excluded += 1 # record skipped window
                            if show_excluded: [ region_info.append('-9') for column in range(len(region_info),len(headers)) ] # record empty metrics

                        if len(region_info) == len(headers): print(*region_info, sep='\t', file=metrics, flush=True) # write region metrics to output


                    if EOF: break # data alerady processed; exit loop & sequentially trigger remaining switches
                    else: 
                        

                        # ESTABLISH/RESET VARIABLES & CONSTANTS
                        if setup: 
                            SITE_count, SNP_count, SINGLETON_count = 0, 0, 0 # GENOME site counts
                            AFS = [ 0 for category in range(0, int(an)+1) ] # GENOME allele frequency spectrum

                            aw, bw, a2 = 0.0, 0.0, 0.0 # N.B. bw is b sub n+1 as per Zeng
                            for j in range(1, int(an)):
                                aw += 1.0 / float(j)  # a1 in Tajima 1989 and an in Zeng 2006
                                a2 += 1.0 / float(j**2)  # This is bn according to Zeng 2006
                                bw += 1.0 / float(j**2)
                            bw += 1.0 / (an**2)  # This is the n+1 part from Zeng
                            b1 = (an + 1) / (3 * (an - 1))  # b1 through e2 are taken directly rom Tajima 1989
                            b2 = (2 * (an**2 + an + 3)) / (9 * (an * (an - 1)))
                            c1 = b1 - (1 / aw)
                            c2 = b2 - (an + 2) / (aw * an) + a2 / aw**2
                            e1 = c1 / aw
                            e2 = c2 / (aw**2 + a2)

                        if setup or new_scaffold or moved_window:
                            site_count, snp_count, singleton_count = 0, 0, 0 # WINDOW site counts
                            afs = [ 0 for category in range(0, int(an)+1) ] # WINDOW allele frequency spectrum
                            p_i, Ehet = [], [] # WINDOW allele frequencies & ehet 


                        # PROCESS SITE
                        if len(gts) >= ind_num: # proceed if enough individuals

                            site_count += 1; SITE_count += 1
                            sgts = numpy.random.choice(gts, size=ind_num, replace=False) # downsample genotypes
                            sac = sum(sgts) # calculate downsampled alt allele count
                            p = float(sac) / an # calculate downsampled allele frequency

                            if 0 < sac < an: # proceed if polymorphic site
                                                                
                                snp_count += 1; SNP_count += 1 # record snp
                                if sac == 1: singleton_count += 1; SINGLETON_count += 1 # record single polymorphism
                                afs[sac] += 1; AFS[sac] += 1 # record allele frequency in spectrum

                                [ stored.append(info) for stored,info in [(p_i, p), (Ehet, p*(1-p))] ] # store site info


                        # PREPARE FOR NEXT SITE
                        if setup or new_scaffold: window_start, window_end = 0, window_size # establish/reset window coordinates
                        if setup or new_scaffold or moved_window: # specify coordinates of next relevant window as required 
                            while int(position) > window_end: window_end += window_size # shift window end forward
                            window_start = window_end - window_size # calculate window start

                        previous_scaffold = scaffold # record region just examined

                        if show_progress: PROGRESS(i, sites_to_process) # show progress


                # ALL SITES PROCESSED
                if summary: scanning = False # 3rd switch; FINAL pass of FINAL loop completed; genome metrics calculated >>> stop scanning & exit
                if EOF: summary = True # 2nd switch; FIRST pass of FINAL loop completed; final window processed >>> calculate genome metrics
                EOF = True # 1st swtich; FINAL pass of FIRST loop complete; all sites processed >>> process final window

    if show_progress: sys.stdout.write("\033[K"); sys.stdout.flush() # clear progress bar

    return windows_processed, windows_excluded


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='</path/to/input.table>', required=True)
    parser.add_argument('-o', type=str, metavar='<output.metrics>', required=True)
    parser.add_argument('-ws', type=int, metavar='<number>', required=True, default=10000, help='Window size (bp)')
    parser.add_argument('-ms', type=int, metavar='<number>', required=True, default=2, help='Minimum SNPs per window')
    parser.add_argument('-ni', type=int, metavar='<number>', required=True, default=5, help='Number of individuals (for downsampling)')
    parser.add_argument('-excluded', action='store_true', help='Include excluded windows in output')
    parser.add_argument('-progress', action='store_true', help='Show progress bar whilst processing')

    input_table, output_metrics, ws, ms, ni, excluded, progress = vars(parser.parse_args()).values()

    calcWPM(variant_table=input_table, wpm_metrics=output_metrics, window_size=ws, minimum_snps=ms, ind_num=ni, show_excluded=excluded, show_progress=progress)

# testing; -ws 10000 -ms 10 -ni 5