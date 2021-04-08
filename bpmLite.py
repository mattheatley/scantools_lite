import argparse, os, sys, subprocess
from ScanToolsLite import recode

def PROGRESS(count, total):
    filled = int(round(50 * count / float(total)))
    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled + '-' * (50 - filled)
    sys.stdout.write(f'[{bar}] {percents}% of {int(total):,} sites complete\r')
    sys.stdout.flush()

def calcBPM(variant_table, bpm_metrics, window_size, minimum_snps, pop_num, show_excluded=False, show_progress=False):
    
    if not os.path.isfile(variant_table): print('ERROR: Problem finding input file; set "-i /PATH/TO/INPUT.TABLE" & run again.'); sys.exit(0)
        
    if not os.path.isdir(os.path.dirname(bpm_metrics)): print('ERROR: Problem finding output directory; set "-o /PATH/TO/OUTPUT.METRICS" & run again.'); sys.exit(0)

    with open(bpm_metrics, 'w') as metrics:

        headers = ['scaffold','window_start','window_end','site_count','snp_count', 'Fst_WC','Fst_H','Fst_N','Rho','Dxy','AFD','FixedDiff'] # specify headers
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

                        *_, scaffold, position, ploidy, gts = line = recode(line) # extract site info
                        
                        n = float(len(gts)) # calculate individual number
                        an = float(ploidy) * n # re-calculate allele number
                        ac = sum(gts) # re-calculate alt allele count
                        p = ac/an # calculate filtered allele frequency

                        new_scaffold = setup or scaffold != previous_scaffold # register scaffold change
                        new_position = setup or position != previous_position # register site change
                        moved_window = setup or int(position) > window_end # register new window
                        
                        if new_scaffold or new_position or moved_window: locus = [] # establish/reset locus info
                        if gts: locus.append([gts, float(ploidy), n, an, ac, p]) # store site info


                    # PROCESS LOCUS AS APPROPRIATE

                    if EOF: pass # data already processed; skip
                    elif pop_num == len(locus): # process if site available for all populations

                        site_count += 1; SITE_count += 1 # record site
                        x_i, n_i, p_i, p_ij, an_i, ac_i, ac_ij  = [[] for i in range(0,7)] # population (i) ploidies, (ii) sizes, (iii) allele frequencies, (iv) within-individual allele frequencies, (v) allele numbers, (vi) alt allele counts & (vii) within-individual alt allele counts 
                        fst_WC_ss_g, fst_WC_ss_i, fst_WC_ss_w, fst_WC_ss_t = [ 0.0 for i in range(0,4) ] # fst sum of squares (i) between populations (ii) between individuals within populations (iii) within individual & (iv) total
                        rho_ss_g, rho_ss_i, rho_ss_t = [ 0.0 for i in range(0,3) ] # rho sum of squares (i) between populations (ii) between individuals within populations & (iii) total
                        FS_Nij2, FS_SNij2, FS_SNij2_over_SNij, RS_SNij2 = [ 0.0 for i in range(0,4) ] # to calculate sample size coefficients in variance components as per spagedi 

                        r = float(len(locus))  # population number; N.B. should always be 2 for pairwise analysis

                        for pop_site in locus:

                            gts, x, n, an, ac, p = pop_site # unpack site info

                            ac_j, p_j = [], [] # population (i) individual alt allele counts & (ii) individual allele frequencies
                            FSNij2temp = 0.0

                            [ stored.append(info) for stored,info in [(x_i, x), (n_i, n), (an_i, an), (ac_i, ac), (p_i, p)] ] # store population info

                            for gt in gts:
                                
                                [ stored.append(info) for stored,info in [(ac_j, gt), (p_j, gt/x)] ] # store individual info 

                                FS_Nij2 += x**2 # individual allele number squared per individual
                                FSNij2temp += x**2 # individual allele number squared per individual

                            [ stored.append(info) for stored,info in [(p_ij, p_j), (ac_ij, ac_j)] ] # store within-population info

                            FS_SNij2_over_SNij += FSNij2temp / an # squared individual allele numbers divided by population allele number
                            FS_SNij2 += an**2  # population allele number squared per population
                            RS_SNij2 += n**2  # population size squared per population

                        tan, tac = [ sum(info) for info in [an_i, ac_i] ] # total (i) allele counts & (ii) allele numbers

                        p_bar = tac / tan # panmixia allele frequency
                        
                        # degrees of freedom
                        df_i = sum([(n - 1.0) for n in n_i])
                        df_w = sum([(x - 1.0) * n for x,n in zip(x_i,n_i)])
                        df_g = r - 1.0
                        df_t = df_i + df_g + df_w # equates to tan - 1
                        
                        # sample size coefficients
                        fn0bis = (FS_SNij2_over_SNij - (FS_Nij2 / tan)) / df_g
                        fn0 = (tan - FS_SNij2_over_SNij) / df_i
                        fnb0 = (tan - (FS_SNij2 / tan)) / df_g
                        rnb0 = (float(sum(n_i)) - (RS_SNij2 / float(sum(n_i)))) / df_g

                        # fst/rho sum of squares
                        for i, pop in enumerate(ac_ij):
                            for j, ind in enumerate(pop):
                                for k in range(0, int(x_i[i])):
                                    k_u = 1 if k < ind else 0 # 1/0 for alt/ref
                                    fst_WC_ss_g += (p_i[i] - p_bar)**2
                                    fst_WC_ss_i += (p_ij[i][j] - p_i[i])**2
                                    fst_WC_ss_w += (k_u - p_ij[i][j])**2
                                    fst_WC_ss_t += (k_u - p_bar)**2 # not used?!
                                rho_ss_i += (p_ij[i][j] - p_i[i])**2
                                rho_ss_g += (p_i[i] - p_bar)**2
                                rho_ss_t += (p_ij[i][j] - p_bar)**2 # not used?!

                        # fst/rho mean squares
                        fst_WC_ms_w = fst_WC_ss_w / df_w
                        fst_WC_ms_i = fst_WC_ss_i / df_i
                        fst_WC_ms_g = fst_WC_ss_g / df_g
                        rho_ms_i = rho_ss_i / df_i
                        rho_ms_g = rho_ss_g / df_g

                        # fst/rho variance components
                        fst_WC_s2_w = fst_WC_ms_w
                        fst_WC_s2_i = (fst_WC_ms_i - fst_WC_s2_w) / fn0
                        fst_WC_s2_g = (fst_WC_ms_g - fst_WC_s2_w - fn0bis * fst_WC_s2_i) / fnb0
                        rho_s2_i = rho_ms_i
                        rho_s2_g = (rho_ms_g - rho_s2_i) / rnb0

                        # fst/rho numerators & denominators
                        fst_WC_num = fst_WC_s2_g
                        fst_WC_den = fst_WC_s2_w + fst_WC_s2_g + fst_WC_s2_i
                        rho_num = rho_s2_g
                        rho_den = rho_s2_i + rho_s2_g

                        # hudson/nei numerators & denomiantors; pairwise analysis only?
                        (ac1, an1, p1), (ac2, an2, p2) = zip(ac_i, an_i, p_i) # unpack population (i) allele counts, (ii) allele numbers & (iii) allele frequencies
                        h_m1, h_m2 = [ (ac*(an-ac))/(an*(an-1)) for ac,an in zip(ac_i,an_i) ] ###
                        fst_H_num = (p1-p2)**2 - h_m1/an1 - h_m2/an2 ###
                        fst_H_den = fst_H_num + h_m1 + h_m2 ###

                        p_avr = (p1+p2)/2 # average population allele frequency
                        fst_N_num = (((p1)-(p2))**2)/2 # to make maximum 1 instead of 2 (Bhatia, 2013)
                        fst_N_den = 2 * p_avr * (1-p_avr) # to make maximum 1 instead of 2 (Bhatia, 2013)
                                                    
                        if any( 0 < p < 1 for p in p_i): # proceed if polymorphic site

                            snp_count += 1; SNP_count += 1 # record snp

                            if pop_num == 2: # pairwise analysis only?
                                d = (p1 * (1.0 - p2)) + (p2 * (1.0 - p1))
                                dxy += d; DXY += d
                                a = p2 - p1
                                afd += a; AFD += a
                                if a == 1.0: dn += 1; Dn += 1
    
                            current_fractions = [(fst_WC_num, fst_WC_den), (fst_H_num, fst_H_den), (fst_N_num, fst_N_den), (rho_num, rho_den)]
                            window_fractions = [ [sum(values) for values in zip(previous, current)] for previous, current in zip(window_fractions, current_fractions) ] # update running WINDOW totals
                            GENOME_fractions = [ [sum(values) for values in zip(previous, current)] for previous, current in zip(GENOME_fractions, current_fractions) ] # update running GENOME totals


                    # CALCULATE WINDOW/GENOME METRICS AS APPROPRIATE

                    if setup: pass # not relevant yet; skip
                    elif new_scaffold or moved_window or EOF: # calculate metrics if moved scaffold / moved window / final window / summary  

                        region = (previous_scaffold, window_start, window_end) if not summary else ('Genome','-','-') # specify relevant region
                        sites, *_ = sites_snps = (site_count, snp_count) if not summary else (SITE_count, SNP_count) # specify relevant counts
                    
                        region_info = [*region, *sites_snps] # specify region info
                                                
                        if not summary and snp_count < minimum_snps: pass # not summary & too few snps to process
                        else: # proceed with calculating metrics
                            
                            try:
                                if not summary: windows_processed += 1 # record processed window

                                fractions, Dxy, Afd, Dn = (window_fractions, dxy, afd, dn) if not summary else (GENOME_fractions, DXY, AFD, DN) # specify relevant metrics

                                Fst_WC, Fst_H, Fst_N, Fac = [ numerator / denominator for numerator, denominator in  fractions ] # calculate (i) fst WC, (ii) fst H, (iii) fst N & (iv) rho component
                                Rho = Fac / (1 + Fac) # calculate rho
                                if pop_num == 2: Dxy, Afd = [ count / float(sites) for count in [Dxy, Afd] ] # calcualte dxy & afd
                                else: Dxy = Afd = Dn = '-' # dxy & afd irrelevant 

                                region_info.extend([Fst_WC, Fst_H, Fst_N, Rho, Dxy, Afd, Dn]) # record region metrics

                            except ZeroDivisionError: pass # no values to calculate

                        if len(region_info) != len(headers):
                            if not summary: windows_excluded += 1 # record skipped window
                            if show_excluded: [ region_info.append('-') for column in range(len(region_info),len(headers)) ] # record empty metrics

                        if len(region_info) == len(headers): print(*region_info, sep='\t', file=metrics, flush=True) # write region metrics to output


                    if EOF: break # data alerady processed; exit loop & sequentially trigger remaining switches
                    else: # prepare for next data


                        # ESTABLISH/RESET VARIABLES
                        if setup:
                            SNP_count, SITE_count = 0, 0 # GENOME site counts
                            FST_WC, FST_H, FST_N, RHO = GENOME_fractions = [ [0.0, 0.0] for i in range(0,4)] # GENOME fst/rho numerators/denominators
                            DN, DXY, AFD = 0, 0.0, 0.0 # GENOME dxy/afd counts

                        if setup or new_scaffold or moved_window:
                            site_count, snp_count = 0, 0 # WINDOW site counts
                            fst_WC, fst_H, fst_N, rho = window_fractions = [ [0.0, 0.0] for i in range(0,4)] # WINDOW fst/rho numerators/denominators
                            dn, dxy, afd = 0, 0.0, 0.0  # WINDOW dxy/afd counts


                        # PREPARE FOR NEXT SITE
                        if setup or new_scaffold: window_start, window_end = 0, window_size # establish/reset window coordinates
                        if setup or new_scaffold or moved_window: # specify coordinates of next relevant window as required 
                            while int(position) > window_end: window_end += window_size # shift window end forward
                            window_start = window_end - window_size # calculate window start

                        previous_scaffold, previous_position = scaffold, position # record region just examined
                    
                    if show_progress: PROGRESS(i, sites_to_process) # show progress


                # ALL SITES PROCESSED
                if summary: scanning = False # 3rd switch; FINAL pass of FINAL loop completed; genome metrics calculated >>> stop scanning & exit
                if EOF: summary = True # 2nd switch; FIRST pass of FINAL loop completed; final window processed >>> calculate genome metrics
                EOF = True # 1st swtich; FINAL pass of FIRST loop complete; all sites processed >>> process final window

    if show_progress: sys.stdout.write("\033[K"); sys.stdout.flush() # clear progress bar

    return windows_processed, windows_excluded


if __name__ == '__main__':  # Used to run code from command line

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='</path/to/input.table>', required=True, help='Variant table')
    parser.add_argument('-o', type=str, metavar='<output.metrics>', required=True, help='Metrics prefix')
    parser.add_argument('-ws', type=int, metavar='<number>', required=True, default=10000, help='Window size (bp)')
    parser.add_argument('-ms', type=int, metavar='<number>', required=True, default=2, help='Minimum SNPs per window')
    parser.add_argument('-np', type=int, metavar='<number>', required=True, default=2, help='Number of populations')
    parser.add_argument('-excluded', action='store_true', help='Include excluded windows in output')
    parser.add_argument('-progress', action='store_true', help='Show progress bar whilst processing')

    input_table, output_metrics, ws, ms, np, excluded, progress = vars(parser.parse_args()).values()

    calcBPM(variant_table=input_table, bpm_metrics=output_metrics, window_size=ws, minimum_snps=ms, pop_num=np, show_excluded=excluded, show_progress=progress)

# testing; -ws 10000 -ms 10 -np 2