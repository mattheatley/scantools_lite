# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#    _____  ______ ___     _   __  ______ ____   ____   __    _____ #
#   / ___/ / ____//   |   / | / / /_  __// __ \ / __ \ / /   / ___/ #
#   \__ \ / /    / /| |  /  |/ /   / /  / / / // / / // /    \__ \  #
#  ___/ // /___ / ___ | / /|  /   / /  / /_/ // /_/ // /___ ___/ /  #
# /____/ \____//_/  |_|/_/ |_/   /_/   \____/ \____//_____//____/   #
#                                                                   #
#    -. ,-.   ,-. ,-.   ,-. ,-.   .      __    ____ ______ ____     #
#    ||\|G|\ /|T|\|C|\ /|A|\|G|\ /|     / /   /  _//_  __// __/     #
#    |/ \|C|\|A|/ \|G|\|T|/ \|C|\||    / /__ _/ /   / /  / _/       #
#    '   `-' `-'   `-' `-'   `-' `-   /____//___/  /_/  /___/       #
#                                                                   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# February 2021 release
# software; bcftools (v1.9), gatk (v3.8), bedtools (v2.29.2)
import os, sys, subprocess, pandas, math, datetime, time, glob, re


# _-_-_-_-_-_-_-_-_-_-|  G E N E R A L   F U N C T I O N S  |-_-_-_-_-_-_-_-_-_-_ #

def CAPTURE(cmd): return subprocess.run(f'{cmd}', shell=True, capture_output=True).stdout.decode('utf-8').strip(' \n') # capture & format terminal output

def user_input(request, choices=None, integer_expected=False, attempts=3):
    if choices and not isinstance(choices,list): choices = list(choices) # ensure correct format
    proceeding = False; loop = 0    
    while not proceeding:
        if loop == attempts: print('\ni give up...\n'); sys.exit(0) # exit if repeatedly incorrect
        response = input(f'\n{request}:  '); loop += 1 # request user input
        if integer_expected:
            if response.isdigit(): response = int(response) # ensure input appropriate
            else: continue # try again...
        if choices: proceeding = response in choices # ensure input appropriate
        else: break # return any input
    return response

def join_neatly(to_join):
    *leading, trailing = to_join
    if leading: leading, *_ = [(", ").join(leading)] if len(leading) > 1 else leading # join by commas if relevant
    return f'{leading} & {trailing}' if leading else trailing # join by & if relevant

def show_missing(input_type, missing_list):
    print(f'{input_type} missing for:')
    for identifier in missing_list: print(f' - {identifier}')

def find_processed(search_dir, search_prefixes, search_suffix, extract_identifiers=False): # search for existing outputs  
    # define search terms
    raw_seperator = f'\.' if not extract_identifiers else f'\.(.*)\.' # seperator without unintended special characters
    *raw_prefixes, raw_suffix = [ re.escape(affix) for affix in [*search_prefixes, search_suffix] ] # affixes without unintended special characters
    search_term = re.compile( f'^({"|".join(raw_prefixes)}){raw_seperator}{raw_suffix}$') # specify file format
    # search directories
    matches = [ (path, matched.groups()) for path, matched in [ (contents.path, search_term.search(contents.name)) for contents in os.scandir(search_dir) ] if matched ] # extract matching files
    matched_files, *_, = *_, matched_identifiers = list(zip(*matches)) if matches else [[]]
    # extract unique file identifiers
    matched_prefixes, *_ = *_, matched_seperators = [ set(identifiers) for identifiers in zip(*matched_identifiers) ] if any(matched_identifiers) else [set()] 
    missing_prefixes = set(search_prefixes).difference(matched_prefixes)
    if extract_identifiers and len(matched_seperators) > 1: 
        print('\nCombining variant tables accross multiple VCFs:')
        for identifier in sorted(matched_seperators): print(f' - {identifier}')
    return missing_prefixes, matched_files
    
def recode(line):
    file_basename, scaffold, position, reference, an, dp, *raw_gts = line.strip(' \n').split('\t') # extract relevant site info
    raw_gts = [ re.split('\||/', gt) for gt in raw_gts ] # split phased & un-phased genotypes
    ploidy = float(max(len(gt) for gt in raw_gts)) # determine ploidy from genotype
    recoded_gts = [ sum(1 for nucleotide in gt if nucleotide != reference) for gt in raw_gts if not '.' in gt ] # recode genotypes as alt allele counts & filter no-calls
    return file_basename, scaffold, position, ploidy, recoded_gts

def SBATCH(job_id, partition, nodes, ntasks, memory, walltime, out_err, task=None, email=None, conda_env=None):
    labels = ['job-name','partition','nodes','ntasks-per-node','mem','time','output','error', 'parsable'] # specify relevant directives
    inputs = [job_id, partition, nodes, ntasks, memory, walltime, *[ f'{out_err}/%x{f".{task}" if task else ""}.{suffix}' for suffix in ['out','err'] ], None] # organise settings
    if email: labels.extend(['mail-user','mail-type']); inputs.extend([email,'END']) # add optional user email address
    sbatch = ''.join([ f'#SBATCH --{option}{f"={info}" if info else ""} \n' for option,info in zip(labels,inputs) ]) # format directives & settings
    sbatch += '\nsource $HOME/.bash_profile\n'
    if conda_env: sbatch += f'\nconda activate {conda_env}'+'\necho ENVIRONMENT $CONDA_DEFAULT_ENV ACTIVE\n' # add optional conda environment
    return '#!/bin/bash\n'+f'{sbatch}\n'+'echo RUNNING ON `hostname`\n'


# ERRORS
error_bad_values = '\nERROR: Key file contains {}; replace these & run again.\n'
error_bad_name = '\nERROR: {} failed; {}.\n'
error_bad_path = '\nERROR: Problem finding {}; set "{} = /PATH/TO/{}" correctly & run again.\n'
error_no_files = '\nERROR: Problem finding {} in "{}".\n'
error_already_exists = '\nERROR: {} already exist; set "overwrite = True" if reprocessing & run again.\n'
error_no_software = '\nERROR: {} not installed.\n'
error_no_info = '\nERROR: {} not specified; set "{}" & run again.\n'
error_either_or = '\nERROR: {} not specified; set either "{}" or "{}" & run again.\n'
error_bad_conda = '\nERROR: GATK could not be run in "{}"; change environment or set "gatk_path = /PATH/TO/GenomeAnalysisTK.jar" & run again.\n'
error_ind_num = '\nERROR: Unable to calculate within-population metrics with so few individuals; set "sample_size = >3" & run again.\n'
error_pop_num = '\nERROR: At least 2 populations must be provided; set "pops = [pop1, pop2...]" & run again.\n'

# MESSAGES
skipping_existing = '\nExisting outputs identified; only oustanding tasks will be processed.'
task_skipped = '\nSkipping {} for {}; already processed.'
task_processed = '\nProcessing {} for {}'



# _-_-_-_-_-_-_-_-_-_-|  C R E A T E   P R O J E C T  |-_-_-_-_-_-_-_-_-_-_ #

class scantools:

    def __init__(self, wrk_dir, vcf_dir=None, encoding='-9'):

        # setup core elemants of scantools object 
        self.code_dir = os.path.dirname(__file__) # extract & store scantools.py directory
        self.wpm_script, self.bpm_script, self.outlier_script = [ f'{self.code_dir}/{name}Lite.py' for name in ['wpm', 'bpm','outliers'] ] # extract paths of accompanying scantool scripts
        
        wrk_dir = f'/{wrk_dir.strip("/")}' # ensure correct path format
        self.wrk_dir = wrk_dir # store working directory
        
        out_error_dir, popkey, log = [ f'{wrk_dir}/{name}' for name in ['OandE','PopKey.csv', 'LogFile.txt'] ] # specify out & error file directory
        self.oande = out_error_dir # store out & error file directory        
        os.makedirs(out_error_dir, exist_ok=True) # create directories as required
        
        # specify formating of later stages
        self.split_dir_prefix = 'VCF' # specify & store vcf dir prefix
        self.split_dirs = [ contents.path for contents in os.scandir(self.wrk_dir) if contents.name.startswith(self.split_dir_prefix) ] # list vcf directories
        self.vcf_suffix = ('vcf.gz','vcf') # specify & store vcf suffixes (as tuple!)
        
        self.raw_suffix = f'raw.table' # specify & store table suffixes
        
        # create and/or import information from key file
        ind_header, all_header = 'Samples','ALL' # specify key file headers
        
        if not os.path.exists(popkey): # check if key file exists
        
            response = user_input(request='Key file could not be found in the working directory; would you like to generate a template for this now? (y/n)', choices=['y','n']) # offer template creation
            if response != 'y': print(f'\nCreate a key file manually in "{wrk_dir}" & run again (see instructions for further details).\n')
            
            else: # create template
                if vcf_dir: vcf_dir = f'/{vcf_dir.strip("/")}' # ensure correct path format
                if not os.path.isdir(vcf_dir): print(error_bad_path.format('input directory','vcf_dir','VCF_DIRECTORY'))

                else: # input directory found
                    vcf_files = { i:relevant for i,relevant in enumerate([contents.path for contents in os.scandir(vcf_dir) if contents.name.endswith(self.vcf_suffix)], 1) } # find vcfs
                    if not vcf_files: print(error_no_files.format('vcfs',vcf_dir))

                    else: # input files found
                        multiple_vcfs_found = len(vcf_files) > 1 # check if multiple vcfs
                        if multiple_vcfs_found:
                            print('\nMultiple vcfs found in the directory provided.')
                            for key,path in vcf_files.items(): print(f'{key}:  {os.path.basename(path)}') # print vcfs to terminal
                            input_key = user_input(request='Select a file from the above list for extracting sample names (enter an integer)', choices=vcf_files.keys(), integer_expected=True) # request vcf selection                   
                        vcf_file = vcf_files[1 if not multiple_vcfs_found else input_key] # specify relevant vcf

                        bcftool_cmd = CAPTURE('command -v bcftools') # test gatk wrapper shortcuts
                        if not bcftool_cmd: print(error_no_software.format('bcftools'))
                        
                        else: # software found
                            found_inds = CAPTURE(f'bcftools query -l {vcf_file}').split('\n') # extract individual names from vcf
                            input_num = user_input(request='Specify the number of populations required (enter an integer)', integer_expected=True) # request population number

                            with open(popkey, 'w') as template:
                                pop_labels = [ f'Pop{i}' for i in range(1,int(input_num)+1) ] # specify population names
                                print(ind_header, *pop_labels, all_header, sep=',', file=template) # write key header
                                for ind in found_inds:
                                    values = [0]*input_num; values.append(1) # create default values
                                    print(ind, *values, sep=',', file=template) # write individual entry
                            print(f'\nKey file template created in "{wrk_dir}"; manually enter values to indicate the relevant population(s) & run again (see instructions for further details).\n')
                
            sys.exit(0) # EXIT TO FINISH TEMPLATE MANUALLY OR FIX ERRORS

        key_data = pandas.read_csv(popkey, header=0, encoding=encoding if encoding != "-9" else None) # read key file
        pop_names = list(key_data.columns.values)[1:] # extract population names
        ind_names = list(key_data[ind_header]) # extract individual names
        inds, sizes = {}, {} # stored population info

        if not all(pandas.to_numeric(key_data[pop], errors='coerce').notnull().all() for pop in pop_names): print(error_bad_values.format('non-numeric &/or empty values'))
        
        else: # all values are numbers
            if not all((key_data[pop_names] <= 1 ).all()): print(error_bad_values.format('values other than 0/1'))
        
            else: # proceed with extracting info
                for pop in pop_names:
                    pop_inds = [ sample for sample, value in zip(ind_names, list(key_data[pop])) if value == 1 ] # extract relevant individuals
                    inds[pop], sizes[pop] = pop_inds, len(pop_inds) # create entry for population/individuals names & size
                self.pops, self.inds, self.sizes = pop_names, inds, sizes # store population/individuals names & sizes
                self.min_ind = min(self.sizes.values()) # extract & store minimum population size
                self.log_file = open(log, 'a+') # create & open a log file
                self.log_file.write('New Instance at: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.datetime.now())) # log info





# _-_-_-_-_-_-_-_-_-_-|  E D I T   P O P U L A T I O N S  |-_-_-_-_-_-_-_-_-_-_ #


    def removePops(self, to_remove):
        if not isinstance(to_remove, list): to_remove = [to_remove] # ensure correct format if single population
        removed = 0
        if not set(to_remove).issubset(set(self.pops)): 
            absent = sorted(set(to_remove).difference(set(self.pops)))
            print(error_bad_name.format(f'Removing populations', f'{join_neatly(absent)} cannot be found'))
        else: # populations exists
            for pop in to_remove:
                removed += 1
                [ stored.remove(pop) if isinstance(stored, list) else stored.pop(pop, None) for stored in [self.pops, self.inds, self.sizes] ] # remove entries
                print(f'\n{pop} removed.')
                self.log_file.write(f'Removed Population: {pop}\n') # log info
            self.min_ind = min(self.sizes.values()) # update minimum population size
            print(f'\n{removed} population(s) removed in total.\n')


    def removeInds(self, to_remove):
        if not isinstance(to_remove, list): to_remove = [to_remove] # ensure correct format if single population
        all_inds = [name for sublist in self.inds.values() for name in sublist] # extract all individual names
        removed = 0
        if not set(to_remove).issubset(set(all_inds)): 
            absent = sorted(set(to_remove).difference(set(all_inds)))
            print(error_bad_name.format(f'Removing individuals', f'{join_neatly(absent)} cannot be found'))
        else: # individuals exists
            for ind in to_remove:
                removed += 1
                for pop, inds in self.inds.items():
                    if ind in inds:
                        self.inds[pop].remove(ind) # remove individual from entry
                        self.sizes[pop] -= 1 # update population size
                        print(f'\n{ind} removed from {pop}.')
                        self.log_file.write(f'Removed Ind: {ind}\n') # log info
            self.min_ind = min(self.sizes.values()) # update minimum population size
            print(f'\n{removed} individual(s) removed in total.\n')

                
    def combinePops(self, to_combine, new_pop):
        if new_pop in self.pops: print(error_bad_name.format(f'Combining populations', f'{new_pop} already exists'))
        else: # new population name unique
            if len(to_combine) != len(set(to_combine)): 
                duplicates = sorted({ pop for pop in to_combine if to_combine.count(pop) > 1 })
                print(error_bad_name.format('Combining populations', f'{join_neatly(duplicates)} duplicated'))
            else: # population names not duplicated
                if not set(to_combine).issubset(set(self.pops)): 
                    absent = sorted(set(to_combine).difference(set(self.pops)))
                    print(error_bad_name.format(f'Combining populations', f'{join_neatly(absent)} cannot be found'))
                else: # all populations exist
                    new_inds = [ ind for inds in [self.inds[pop] for pop in to_combine] for ind in inds ] # combine relevant individuals
                    self.pops.append(new_pop) # store new population name
                    self.inds[new_pop], self.sizes[new_pop] = new_inds, len(new_inds) # create entry for new population/individuals names & size
                    self.log_file.write(f'Combined Populations: {join_neatly(to_combine)} combined as {new_pop}\n') # log info
                    print(f'\n{join_neatly(to_combine)} combined as {new_pop}.\n')





# _-_-_-_-_-_-_-_-_-_-|  S P L I T   V C F S   B Y   P O P U L A T I O N  |-_-_-_-_-_-_-_-_-_-_ #

    def splitVCFs(self, ref_path, vcf_dir, min_dp, mffg, # required settings
    pops='all', gatk_path=None, conda_env=None, # optional settings
    partition='defq', nodes=1, ntasks=1, memory='4g', walltime='10:00:00', description='split', # slurm settings
    overwrite=False, keep_intermediates=False, use_scratch=False, scratch_path=None, debug=False): # general settings
        
        print('\nPreparing to split vcfs')
        passed_prechecks = False
        if keep_intermediates and use_scratch: use_scratch = False # overide to keep intermediates

        ref_path = f'/{ref_path.strip("/")}' # ensure correct path format
        if not os.path.isfile(ref_path) or not ref_path.endswith(('.fasta','.fa')): print(error_bad_path.format('reference genome','ref_path','reference_genome.[fasta|fa]'))

        else: # reference genome found
            ref_name, *_ = os.path.splitext(os.path.basename(ref_path)) # extract fasta basename
            index_files = [ contents.name for contents in os.scandir(os.path.dirname(ref_path)) if contents.name.endswith(('.dict','.fai')) and contents.name.startswith(ref_name) ] # check reference indexes present            
            if len(index_files) != 2: print(error_no_files.format('fasta index &/or sequence dictionary for reference genome', os.path.dirname(ref_path)))
    
            else: # index & sequence dictionary found
                vcf_dir = f'/{vcf_dir.strip("/")}' # ensure correct path format        
                if not os.path.isdir(vcf_dir):print(error_bad_path.format('input directory','vcf_dir','VCF_DIRECTORY'))

                else: # input directory found
                    vcf_list = [ contents.path for contents in os.scandir(vcf_dir) if contents.name.endswith(self.vcf_suffix) ] # list vcfs found
                    if not vcf_list: print(error_no_files.format('vcfs',vcf_dir))

                    else: # input files found
                        if scratch_path: scratch_path = f'/{scratch_path.strip("/")}' # ensure correct path format        
                        if use_scratch and (not scratch_path or not os.path.isdir(scratch_path)): print(error_bad_path.format('scratch path', 'scratch_path', 'SCRATCH_PATH'))

                        else: # scratch path found or not required
                            self.vcf_dir = vcf_dir # store vcf directory
                            vcf_dir_name = os.path.basename(vcf_dir) # extract vcf directory name
                            
                            split_dir = f'{self.wrk_dir}/{self.split_dir_prefix}_{vcf_dir_name}_DP{min_dp}.M{mffg}' # specify split vcf directory'
                            if split_dir not in self.split_dirs: self.split_dirs.append(split_dir) # store split vcf directory as required
                            os.makedirs(split_dir, exist_ok=True) # create split vcf directory as required
                                                
                            tmp_dir = scratch_path if use_scratch else split_dir # specify appropriate tmp directory

                            if pops == 'all': pops = self.pops # specify relevant pop names
                            if not isinstance(pops, list): pops = [pops] # ensure correct format if single population
                            
                            *_, output_tables = find_processed(search_dir=split_dir, search_prefixes=pops, search_suffix=self.raw_suffix, extract_identifiers=True) # check if available outputs

                            if len(output_tables) == len(vcf_list)*len(pops) and not overwrite: print(error_already_exists.format('Variant tables'))
                            
                            else: # outputs do not exist or being overwritten

                                if output_tables and not overwrite: print(skipping_existing)

                                if not gatk_path and not conda_env: print(error_either_or.format('GATK access','gatk_path = /PATH/TO/GenomeAnalysisTK.jar','conda_env = ENVIRONMENT_NAME'))
                                else: # gatk access specified

                                    if gatk_path: # prioritise gatk path
                                        if conda_env: print(f'\nBoth a GATK path & conda environment were provided; using "{gatk_path}"\n')
                                        gatk_path = f'/{gatk_path.strip("/")}' # ensure correct path format
                                        if not os.path.isfile(gatk_path) or not gatk_path.endswith('.jar'): # check path correct
                                            print(error_bad_path.format('GATK path','gatk_path','GenomeAnalysisTK.jar'))
                                        else: passed_prechecks = True

                                    else: # alternatively use conda environment
                                        print('\nChecking conda environment...')
                                        env_active = conda_env == os.getenv('CONDA_DEFAULT_ENV') # check if active env is specified env
                                        test_shortcut = 'command -v {}' if env_active else f'source ~/.bash_profile; conda activate {conda_env}; if [ "$CONDA_DEFAULT_ENV" == "{conda_env}" ]; then command -v {{}}; fi' # specify test shortcut command
                                        working_shortcuts = [ tested for tested,path in [ (cmd,CAPTURE(test_shortcut.format(cmd))) for cmd in ['gatk3','gatk'] ] if path ] # test gatk wrapper shortcuts

                                        if not working_shortcuts: print(error_bad_conda.format(conda_env))
                                        else: # working shortcut found

                                            print('\nEnvironment OK!\n')
                                            gatk_shortcut, *_ = working_shortcuts # extract first working shortcut
                                            if len(working_shortcuts) > 1: print(f'Multiple GATK shortcuts found; using "{gatk_shortcut}"\n')       
                                            passed_prechecks = True
                                                                        
                                    if passed_prechecks:
                                        
                                        gatk_cmd = f'java -Xmx{memory} -jar {gatk_path}' if gatk_path else f'{gatk_shortcut} -Xmx{memory}' # specify appropriate gatk command

                                        for pop in pops:

                                            ind_flags = ' '.join([ f'-sn {ind}' for ind in self.inds[pop] ]) # specify individuals to be extracted
                                            mfg = int(math.ceil(float(self.sizes[pop]) * float(mffg))) # calculate maximum no calls permitted
                                            
                                            for input_vcf in vcf_list:

                                                basename_vcf, *_ = os.path.basename(input_vcf).rsplit('.vcf',1) # extract vcf basename
                                                basename_pop_vcf = f'{pop}.{basename_vcf}' # specify population vcf name

                                                if any(f'/{basename_pop_vcf}.' in table for table in output_tables) and not overwrite: print(task_skipped.format(os.path.basename(input_vcf), pop)) # only process outstanding vcfs for population
                                                else: # requires processing

                                                    print(task_processed.format(os.path.basename(input_vcf), pop))

                                                    sh_file, raw_table = [ f'{split_dir}/{basename_pop_vcf}.{suffix}' for suffix in ['sh', self.raw_suffix]  ] # specify output formats
                                                    intermediates = [ f'{tmp_dir}/{basename_pop_vcf}.{suffix}' for suffix in ['vcf', f'dp{min_dp}.vcf', f'dp{min_dp}.init.vcf', f'dp{min_dp}.m{mffg}.vcf']  ] # specify intermediates formats
                                                    split_vcf, depth_filtered_vcf, initialised_vcf, nocall_filtered_vcf = intermediates # specify intermediate files

                                                    hpcc_directives = SBATCH(job_id=pop, partition=partition, nodes=nodes, ntasks=ntasks, memory=memory, walltime=walltime, 
                                                    out_err=self.oande, task=description, conda_env=None if gatk_path else conda_env) # specify hpcc directives (slurm)

                                                    with open(sh_file, 'w') as sh:
                                                                                                        
                                                        sh.write(f'{hpcc_directives}\n'+ # hpcc directives (slurm)
 
                                                        f'{gatk_cmd} -T SelectVariants -R {ref_path} -V {input_vcf} {ind_flags} -o {split_vcf}\n'+ # extract population sites

                                                        f'{gatk_cmd} -T VariantFiltration -R {ref_path} -V {split_vcf} --genotypeFilterExpression "DP < {min_dp}" --genotypeFilterName "DP" -o {depth_filtered_vcf}\n') # filter sites below minimum depth
                                                        if not keep_intermediates: sh.write(f'rm {split_vcf}*\n') # remove associated intermediate files  
                                                        
                                                        sh.write(f'{gatk_cmd} -T VariantFiltration -R {ref_path} -V {depth_filtered_vcf} --setFilteredGtToNocall -o {initialised_vcf}\n') # initialise for no-call filtering
                                                        if not keep_intermediates: sh.write(f'rm {depth_filtered_vcf}*\n') # remove associated intermediate files  
                                                        
                                                        sh.write(f'{gatk_cmd} -T SelectVariants -R {ref_path} -V {initialised_vcf} --maxNOCALLnumber {mfg} -o {nocall_filtered_vcf}\n') # filter sites exceeding permitted no-calls
                                                        if not keep_intermediates: sh.write(f'rm {initialised_vcf}*\n') # remove associated intermediate files                     
                                                        
                                                        sh.write(f'{gatk_cmd} -T VariantsToTable -R {ref_path} -V {nocall_filtered_vcf} -F CHROM -F POS -F REF -F AN -F DP -GF GT -o {raw_table}\n') # create variants table
                                                        if not keep_intermediates: sh.write(f'rm {nocall_filtered_vcf}*\n') # remove associated intermediate files                    
                                                    
                                                    if not debug: 
                                                        p1 = subprocess.Popen(f'sbatch -d singleton {sh_file}', shell=True) # submit script to hpcc (slurm)
                                                        sts1 = os.waitpid(p1.pid, 0)[1]
                                                        os.remove(sh_file) # remove associated files
                                                    else: print(open(sh_file, 'r').read()) # print script to terminal

                                        if not debug: # log info
                                            self.log_file.write(f'###  Split VCFs  ###\n"'+
                                            f'Input Directory: {vcf_dir}\n'+
                                            f'Reference Path: {ref_path}\n'+
                                            f'Output Directory: {split_dir}\n'+
                                            f'Min Depth Per Individual: {min_dp}\n'+
                                            f'Max Fraction of Filtered Genotypes: {mffg}\n'+
                                            f'Populations: {", ".join(pops)}\n\n')





# _-_-_-_-_-_-_-_-_-_-|  C A L C U L A T E   W I T H I N   P O P U L A T I O N   M E T R I C S  |-_-_-_-_-_-_-_-_-_-_ #

    def calcwpm(self, table_dir, window_size, min_snps, # required settings    
    sample_size=None, pops='all', show_excluded=False, # optional settings
    partition='defq', nodes=1, ntasks=1, memory='1g', walltime='10:00:00', conda_env=None, description='wpm', # slurm & conda settings
    overwrite=False, keep_intermediates=False, use_scratch=False, scratch_path=None, debug=False): # general settings

        print('\nPreparing to calculate within population metrics')
        if keep_intermediates and use_scratch: use_scratch = False # overide to keep intermediates

        table_dir = f'/{table_dir.strip("/")}' # ensure correct path format

        if not os.path.exists(table_dir): print(error_bad_path.format('input directory','table_dir','TABLE_DIRECTORY'))
        else: # input directory found

            if pops == 'all': pops = self.pops # specify relevant population names
            if not isinstance(pops, list): pops = [pops] # ensure correct format if single population

            if not sample_size: sample_size = self.min_ind if pops == 'all' else min([size for pop,size in self.sizes.items() if pop in pops ]) # extract minimum population size
            print(f'\nPopulations will be downsampled to {sample_size} individuals.')

            if sample_size <= 3: print(error_ind_num)
            else: # downsampling size valid

                if scratch_path: scratch_path = f'/{scratch_path.strip("/")}' # ensure correct path format  
                      
                if use_scratch and (not scratch_path or not os.path.isdir(scratch_path)): print(error_bad_path.format('scratch path', 'scratch_path', 'SCRATCH_PATH'))
                else: # scratch path found or not required

                    tmp_dir = scratch_path if use_scratch else table_dir # specify appropriate tmp directory

                    missing_inputs, table_list = find_processed(search_dir=table_dir, search_prefixes=pops, search_suffix=self.raw_suffix, extract_identifiers=True) # check if available inputs

                    if missing_inputs: 
                        print(error_no_files.format('raw variant tables for all populations', table_dir))
                        show_missing(input_type='Variant tables', missing_list=missing_inputs)
                    else: # input files found                

                        metrics_suffix = f'WS{window_size}.MS{min_snps}.SS{sample_size}.wpm' # specify metrics suffix 

                        missing_outputs, output_tables = find_processed(search_dir=table_dir, search_prefixes=pops, search_suffix=metrics_suffix) # check if already outputs
                        
                        if not missing_outputs and not overwrite: print(error_already_exists.format('Within-population metrics for these populations'))
                        else: # outputs do not exist or being overwritten
                            
                            if output_tables and not overwrite: print(skipping_existing)

                            for pop in pops:
                            
                                metrics_prefix = f'{pop}.{metrics_suffix}' # specify metrics file name

                                if not pop in missing_outputs and not overwrite: print(task_skipped.format('metrics', pop)) # only process outstanding vcfs for population
                                else: # requires processing

                                    print(task_processed.format('metrics', pop))

                                    sh_file, wpm_metrics = [ f'{table_dir}/{metrics_prefix}{suffix}' for suffix in ['.sh',''] ] # specify output formats
                                    concat_table = f'{tmp_dir}/{metrics_prefix}.{self.raw_suffix}.cat' # specify intermediate file

                                    hpcc_directives = SBATCH(job_id=metrics_prefix, partition=partition, nodes=nodes, ntasks=ntasks, memory=memory, walltime=walltime, 
                                    out_err=self.oande, task=description) # specify slurm directives

                                    with open(sh_file, 'w') as sh:
                                        
                                        sh.write(f'{hpcc_directives}\n'+ # hpcc directives (slurm)

                                        f'awk \'FNR > 1 {{OFS="\\t"; print FILENAME,$0}}\' {table_dir}/{pop}.*.{self.raw_suffix} > {concat_table}\n'+ # merge relevant raw tables from multiple vcfs (*) excluding the headers (FNR>1)
                                        f'python3 {self.wpm_script} -i {concat_table} -o {wpm_metrics} -ws {window_size} -ms {min_snps} -ni {sample_size} {"-excluded" if show_excluded else ""}\n')
                                        if not keep_intermediates: sh.write(f'rm {concat_table}\n') # remove associated intermediate files  

                                    if not debug:
                                        p3 = subprocess.Popen(f'sbatch {sh_file}', shell=True) # submit script to hpcc (slurm)
                                        sts3 = os.waitpid(p3.pid, 0)[1]
                                        os.remove(sh_file) # remove associated files  
                                    else: print(open(sh_file, 'r').read()) # print script to terminal

                            if not debug:
                                self.log_file.write( # log info
                                '###  Calculate Within-Population-Metrics  ###\n'+
                                f'Input Directory: {table_dir}\n'+
                                f'Window Size: {window_size}\n'+
                                f'Minimum SNPs: {min_snps}\n'+
                                f'Variant Table: Raw\n'+
                                f'Populations: {", ".join(pops)}\n\n')





# _-_-_-_-_-_-_-_-_-_-|  C A L C U L A T E   B E T W E E N   P O P U L A T I O N   M E T R I C S  |-_-_-_-_-_-_-_-_-_-_ #


    def calcbpm(self, table_dir, window_size, min_snps, pops, # required settings    
    show_excluded=False, metrics_prefix=None, # optional settings
    partition='defq', nodes=1, ntasks=1, memory='1g', walltime='10:00:00', conda_env=None, description='bpm', # slurm & conda settings
    overwrite=False, keep_intermediates=False, use_scratch=False, scratch_path=None, debug=False): # general settings
        
        print('\nPreparing to calculate between population metrics')
        if keep_intermediates and use_scratch: use_scratch = False # overide to keep intermediates

        if not isinstance(pops, list) or len(pops) < 2: print(error_pop_num)
        else: # valid population number

            table_dir = f'/{table_dir.strip("/")}' # ensure correct path format
            
            if not os.path.exists(table_dir): print(error_bad_path.format('input directory','table_dir','TABLE_DIRECTORY'))
            else: # input directory found

                if scratch_path: scratch_path = f'/{scratch_path.strip("/")}' # ensure correct path format        
                
                if use_scratch and (not scratch_path or not os.path.isdir(scratch_path)): print(error_bad_path.format('scratch path', 'scratch_path', 'SCRATCH_PATH'))    
                else: # scratch path found or not required

                    tmp_dir = scratch_path if use_scratch else table_dir # specify appropriate tmp directory

                    missing_inputs, table_list = find_processed(search_dir=table_dir, search_prefixes=pops, search_suffix=self.raw_suffix, extract_identifiers=True) # check if available inputs

                    if missing_inputs: 
                        print(error_no_files.format('raw variant tables for all populations', table_dir))
                        show_missing(input_type='Variant tables', missing_list=missing_inputs)
                    else: # input files found                

                        tables_to_sort = ' '.join(sorted(table_list))
                        metrics_prefix = f'{metrics_prefix if metrics_prefix else "-vs-".join(sorted(pops))}' # specify metrics prefix
                        metrics_suffix = f'WS{window_size}.MS{min_snps}.bpm' # specify metrics suffix

                        missing_outputs, *_ = find_processed(search_dir=table_dir, search_prefixes=[metrics_prefix], search_suffix=metrics_suffix) # check if already outputs
                        
                        if not missing_outputs and not overwrite: print(error_already_exists.format('Between-population metrics for these populations'))
                        else: # outputs do not exist or being overwritten

                            print(task_processed.format('metrics', join_neatly(pops)))

                            sh_file, bpm_metrics = [ f'{table_dir}/{metrics_prefix}.{suffix}' for suffix in ['sh',metrics_suffix] ] # specify output formats
                            sorted_concat_table = f'{tmp_dir}/{metrics_prefix}.{metrics_suffix}.{self.raw_suffix}.cat.sort' # specify intermediate file
                                
                            with open(sh_file, 'w') as sh:

                                hpcc_directives = SBATCH(job_id=metrics_prefix, partition=partition, nodes=nodes, ntasks=ntasks, memory=memory, walltime=walltime, 
                                out_err=self.oande, task=description) # specify hpcc directives (slurm)

                                sh.write(f'{hpcc_directives}\n'+ # hpcc directives
                                f'awk \'FNR > 1 {{OFS="\\t"; print FILENAME,$0}}\' {tables_to_sort} | sort -k 2,2b -k 3,3bn -k 1,1 > {sorted_concat_table}\n'+ # sort by scaffold, site & filename (population)
                                f'python3 {self.bpm_script} -i {sorted_concat_table} -o {bpm_metrics} -ws {window_size} -ms {min_snps} -np {len(pops)} {"-excluded" if show_excluded else ""}\n')
                                if not keep_intermediates: sh.write(f'rm {sorted_concat_table}\n') # remove associated files  

                            if not debug:
                                p4 = subprocess.Popen(f'sbatch {sh_file}', shell=True) # submit script to hpcc (slurm)
                                sts4 = os.waitpid(p4.pid, 0)[1]
                                os.remove(sh_file) # remove associated files  
                            else: print(open(sh_file, 'r').read()) # print script to terminal

                            if not debug:
                                self.log_file.write( # log info
                                '###  Calculate Between-Population-Metrics  ###\n'+
                                f'Input Directory: {table_dir}\n'+
                                f'Window Size: {window_size}\n'+
                                f'Minimum SNPs: {min_snps}\n'+
                                f'Variant Table: Raw\n'+
                                f'Populations: {", ".join(pops)}\n\n')





# _-_-_-_-_-_-_-_-_-_-|  E X T R A C T   O U T L I E R S  |-_-_-_-_-_-_-_-_-_-_ #


    def findoutliers(self, metric_dir, percentile, # required settings    
    metric=None, metrics_prefix=None, # optional settings
    to_intersect=False, gff_path=None, # intersect settings
    to_plot=False, plot_range=50000, pops=None, table_dir=None, fai_path=None, # plot settings
    partition='defq', nodes=1, ntasks=1, memory='4g', walltime='10:00:00', conda_env=None, description='outliers', # slurm & conda settings
    overwrite=False, keep_intermediates=False, use_scratch=False, scratch_path=None, debug=False): # general settings
        
        print('\nPreparing to find outliers')
        if keep_intermediates and use_scratch: use_scratch = False # overide to keep intermediates

        if to_plot: to_intersect = True # overide to intersect
        passed_prechecks = not to_intersect # determine if additional prechecks required

        metric_dir, table_dir, gff_path, fai_path = [ f'/{path.strip("/")}' if path else path for path in [metric_dir, table_dir, gff_path, fai_path] ] # ensure correct path format
        
        if not os.path.exists(metric_dir): print(error_bad_path.format('input directory','metric_dir','METRIC_DIRECTORY'))
        else: # input directory found

            if not pops and not metrics_prefix: print(error_either_or.format('Between-population metric table identifer', 'pops = [Pop1, Pop2]','metrics_prefix = <Metric table prefix>'))
            else: # metric identifier provided

                if pops and not isinstance(pops, list) or (isinstance(pops,list) and len(pops) < 2): print(error_pop_num)
                else: # valid population number

                    metrics_prefix = f'{metrics_prefix if metrics_prefix else "-vs-".join(sorted(pops))}' # specify metrics prefix
            
                    metrics_list = [ contents.path for contents in os.scandir(metric_dir) if contents.name.startswith(metrics_prefix) and contents.name.endswith('.bpm') ] # list relevant tables
        
                    if not metrics_list: print(error_no_files.format(f'Between population metrics for "{metrics_prefix}"',metric_dir))
                    else: # input files found

                        if to_intersect:
                        
                            if not gff_path or not os.path.isfile(gff_path) or not gff_path.endswith(('.gff','.gff3')): print(error_bad_path.format('gff file','gff_path','annotation.[gff|gff3]'))
                            else: # gff file found

                                if not to_plot: passed_prechecks = True
                                else: # plotting prechecks required
                                                    
                                    if not pops: print(error_no_info.format('AFD plot minuend & subtrahend populations', 'pops = [Minuend Pop, Subtrahend Pop]'))
                                    else: # afd plot info provided

                                        if not table_dir or not os.path.isdir(table_dir): print(error_bad_path.format('input directory','table_dir','TABLE_DIRECTORY'))
                                        else: # input directory found

                                            if not fai_path or not os.path.isfile(fai_path) or not fai_path.endswith('.fai'): print(error_bad_path.format('fai file','fai_path','reference.fasta.fai'))
                                            else: # fai file found

                                                if scratch_path: scratch_path = f'/{scratch_path.strip("/")}' # ensure correct path format        
                                                
                                                if use_scratch and (not scratch_path or not os.path.isdir(scratch_path)): print(error_bad_path.format('scratch path', 'scratch_path', 'SCRATCH_PATH'))
                                                else: # scratch path found or not required

                                                    tmp_dir = scratch_path if use_scratch else table_dir # specify appropriate tmp directory

                                                    missing_inputs, table_list = find_processed(search_dir=table_dir, search_prefixes=pops, search_suffix=self.raw_suffix, extract_identifiers=True) # check if available inputs

                                                    if missing_inputs: 
                                                        print(error_no_files.format('raw variant tables for all populations', table_dir))
                                                        show_missing(input_type='Variant tables', missing_list=missing_inputs)
                                                    else: # input files found                

                                                        tables_to_sort = ' '.join(sorted(table_list))
                                                        minuend, subtrahend = pops
                                                        print(f'\nAFD minuend: {minuend}')
                                                        passed_prechecks = True

                        if passed_prechecks:
                    
                            single_metric, *_ = headers = sorted(set().union(*[set(CAPTURE(f'head -n 1 {metrics} | cut -f 6-').split('\t')) for metrics in metrics_list]))

                            if not metric in headers:
                                print('\nMetric not specified or not found within between-population metric table.')
                                multiple_metrics_found = len(headers) > 1 # check if multiple vcfs
                                metrics = { i:relevant for i,relevant in enumerate(headers, 1) } # find vcfs
                                if multiple_metrics_found:
                                    print('\nMultiple metrics found in input files.')
                                    for key,path in metrics.items(): print(f'{key}:  {os.path.basename(path)}') # print vcfs to terminal
                                    input_key = user_input(request='Select a metric from the above list to use (enter an integer)', choices=metrics.keys(), integer_expected=True) # request vcf selection                   
                                else: print(f'\nOnly one metric found; using {single_metric}.')
                                metric = metrics[1 if not multiple_metrics_found else input_key] # specify relevant vcf

                            genome_scan_dir = f'{table_dir}/{metrics_prefix}/{metric}'

                            if os.path.exists(genome_scan_dir) and not overwrite: print(error_already_exists.format('Genome scan outputs'))
                            else: # outputs do not exist or being overwritten

                                description = f'{description}-{metric}'

                                for metrics_file in metrics_list:

                                    bpm_prefix, *_ = os.path.splitext(os.path.basename(metrics_file))
                                    genome_scan_info = f'{bpm_prefix}_{metric}'

                                    print(task_processed.format('genome scan', genome_scan_info))

                                    hpcc_directives = SBATCH(job_id=genome_scan_info, partition=partition, nodes=nodes, ntasks=ntasks, memory=memory, walltime=walltime, 
                                    out_err=self.oande, task=description, conda_env=conda_env) # specify hpcc directives (slurm)
                                    
                                    sh_file = f'{metric_dir}/{genome_scan_info}.sh' # specify output formats

                                    with open(sh_file, 'w') as sh:

                                        sh.write(f'{hpcc_directives}\n') # hpcc directives
                                        if overwrite: sh.write(f'rm -r {genome_scan_dir}\n') # remove any existing genome scan directory
                                        
                                        scan_args = f'-i {metrics_file} -o {genome_scan_dir} -percentile {percentile} -metric {metric}'

                                        if to_intersect: scan_args += f' -intersect -gff {gff_path}'
                                        
                                        if to_plot: 
                                            sorted_concat_table = f'{tmp_dir}/{genome_scan_info}.{self.raw_suffix}.cat.sort' # specify intermediate file
                                            scan_args += f' -plot -range {plot_range} -pops {minuend} {subtrahend} -table {sorted_concat_table} -fai {fai_path}'
                                            sh.write(f'awk \'FNR > 1 {{OFS="\\t"; print FILENAME,$0}}\' {tables_to_sort} | sort -k 2,2b -k 3,3bn -k 1,1 > {sorted_concat_table}\n') # sort by scaffold, site & filename (population)
                                        
                                        sh.write(f'python3 {self.outlier_script} {scan_args}\n')
                                        if to_plot and not keep_intermediates: sh.write(f'rm {sorted_concat_table}\n') # remove associated files  

                                    if not debug:
                                        p4 = subprocess.Popen(f'sbatch {sh_file}', shell=True) # submit script to hpcc (slurm)
                                        sts4 = os.waitpid(p4.pid, 0)[1]
                                        os.remove(sh_file) # remove associated files 
                                    else: print(open(sh_file, 'r').read()) # print script to terminal
