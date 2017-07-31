# -*- coding: utf-8 -*-
#BEGIN_HEADER
from __future__ import print_function
from __future__ import division

import os
import sys
import shutil
import hashlib
import subprocess
import requests
requests.packages.urllib3.disable_warnings()
import re
import traceback
import uuid
from datetime import datetime
from pprint import pprint, pformat

import numpy as np
import math
from Bio import SeqIO

from Workspace.WorkspaceClient import Workspace as workspaceService
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from DataFileUtil.DataFileUtilClient import DataFileUtil as DFUClient
from SetAPI.SetAPIServiceClient import SetAPI
from KBaseReport.KBaseReportClient import KBaseReport

#END_HEADER


class kb_assembly_compare:
    '''
    Module Name:
    kb_assembly_compare

    Module Description:
    ** A KBase module: kb_assembly_compare
**
** This module contains Apps for comparing, combining, and benchmarking assemblies
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/dcchivian/kb_assembly_compare"
    GIT_COMMIT_HASH = "751d420f0c2e542922983f91387e8536bd4bb373"

    #BEGIN_CLASS_HEADER
    workspaceURL     = None
    shockURL         = None
    handleURL        = None
    serviceWizardURL = None
    callbackURL      = None
    scratch          = None

    # wrapped program(s)
    MUMMER_bin = '/usr/local/bin/mummer'
    NUCMER_bin = '/usr/local/bin/nucmer'

    # log
    def log(self, target, message):
        timestamp = str(datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3])
        if target is not None:
            target.append('['+timestamp+'] '+message)
        print('['+timestamp+'] '+message)
        sys.stdout.flush()

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['srv-wiz-url']
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.scratch = os.path.abspath(config['scratch'])

        pprint(config)

        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
        #END_CONSTRUCTOR
        pass


    def run_contig_distribution_compare(self, ctx, params):
        """
        :param params: instance of type "Contig_Distribution_Compare_Params"
           (contig_distribution_compare() ** **  Compare Assembly Contig
           Length Distributions) -> structure: parameter "workspace_name" of
           type "workspace_name" (** The workspace object refs are of form:
           ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_assembly_refs" of type "data_obj_ref"
        :returns: instance of type "Contig_Distribution_Compare_Output" ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN run_contig_distribution_compare

        #### STEP 0: basic init
        ##
        console = []
        invalid_msgs = []
        report_text = ''
        self.log(console, 'Running run_contig_distribution_compare(): ')
        self.log(console, "\n"+pformat(params))

        # Auth
        token = ctx['token']
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # API Clients
        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'
        # wsClient
        try:
            wsClient = workspaceService(self.workspaceURL, token=token)
        except Exception as e:
            raise ValueError('Unable to instantiate wsClient with workspaceURL: '+ self.workspaceURL +' ERROR: ' + str(e))
        # setAPI_Client
        try:
            #setAPI_Client = SetAPI (url=self.callbackURL, token=ctx['token'])  # for SDK local.  local doesn't work for SetAPI
            setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
        except Exception as e:
            raise ValueError('Unable to instantiate setAPI_Client with serviceWizardURL: '+ self.serviceWizardURL +' ERROR: ' + str(e))
        # auClient
        try:
            auClient = AssemblyUtil(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate auClient with callbackURL: '+ self.callbackURL +' ERROR: ' + str(e))
        # dfuClient
        try:
            dfuClient = DFUClient(self.callbackURL)
        except Exception as e:
            raise ValueError('Unable to instantiate dfu_Client with callbackURL: '+ self.callbackURL +' ERROR: ' + str(e))

        # param checks
        required_params = ['workspace_name',
                           'input_assembly_refs'
                          ]
        for arg in required_params:
            if arg not in params or params[arg] == None or params[arg] == '':
                raise ValueError ("Must define required param: '"+arg+"'")

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[]
        for input_ref in params['input_assembly_refs']:
            provenance[0]['input_ws_objects'].append(input_ref)

        # set the output paths
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        html_output_dir = os.path.join(output_dir,'html')
        if not os.path.exists(html_output_dir):
            os.makedirs(html_output_dir)


        #### STEP 1: get assembly refs
        ##
        if len(invalid_msgs) == 0:
            set_obj_type = "KBaseSets.AssemblySet"
            assembly_obj_types = ["KBaseGenomeAnnotations.Assembly", "KBaseGenomes.ContigSet"]
            accepted_input_types = [set_obj_type] + assembly_obj_types
            assembly_refs = []
            assembly_names = []
            assembly_refs_seen = dict()
        
            for i,input_ref in enumerate(params['input_assembly_refs']):
                # assembly obj info
                try:
                    [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
                    input_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_ref}]})[0]
                    input_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_obj_info[TYPE_I])  # remove trailing version
                    input_obj_name = input_obj_info[NAME_I]

                    self.log (console, "GETTING ASSEMBLY: "+str(input_ref)+" "+str(input_obj_name))  # DEBUG

                except Exception as e:
                    raise ValueError('Unable to get object from workspace: (' + input_ref +')' + str(e))
                if input_obj_type not in accepted_input_types:
                    raise ValueError ("Input object of type '"+input_obj_type+"' not accepted.  Must be one of "+", ".join(accepted_input_types))

                # add members to assembly_ref list
                if input_obj_type in assembly_obj_types:
                    try:
                        assembly_seen = assembly_refs_seen[input_ref]
                        continue
                    except:
                        assembly_refs_seen[input_ref] = True
                        assembly_refs.append(input_ref)
                        assembly_names.append(input_obj_name)
                elif input_obj_type != set_obj_type:
                    raise ValueError ("bad obj type for input_ref: "+input_ref)
                else:  # add assembly set members

                    try:
                        assemblySet_obj = setAPI_Client.get_assembly_set_v1 ({'ref':input_ref, 'include_item_info':1})
                    except Exception as e:
                        raise ValueError('Unable to get object from workspace: (' + input_ref +')' + str(e))
                    
                    for assembly_obj in assemblySet_obj['data']['items']:
                        this_assembly_ref = assembly_obj['ref']
                        try:
                            assembly_seen = assembly_refs_seen[this_assembly_ref]
                            continue
                        except:
                            assembly_refs_seen[this_assembly_ref] = True
                            assembly_refs.append(this_assembly_ref)
                            try:
                                [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
                                this_input_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':this_assembly_ref}]})[0]
                                this_input_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_obj_info[TYPE_I])  # remove trailing version
                                this_input_obj_name = input_obj_info[NAME_I]
                            
                                assembly_names.append(this_input_obj_name)
                            except Exception as e:
                                raise ValueError('Unable to get object from workspace: (' + this_assembly_ref +')' + str(e))


        #### STEP 2: Get assemblies to score as fasta files
        ##
        if len(invalid_msgs) == 0:
            self.log (console, "Retrieving Assemblies")  # DEBUG

            #assembly_outdir = os.path.join (output_dir, 'score_assembly')
            #if not os.path.exists(assembly_outdir):
            #    os.makedirs(assembly_outdir)
            score_assembly_file_paths = []

            for ass_i,input_ref in enumerate(assembly_refs):
                self.log (console, "\tAssembly: "+assembly_names[ass_i]+" ("+assembly_refs[ass_i]+")")  # DEBUG
                contig_file = auClient.get_assembly_as_fasta({'ref':assembly_refs[ass_i]}).get('path')
                sys.stdout.flush()
                contig_file_path = dfuClient.unpack_file({'file_path': contig_file})['file_path']
                score_assembly_file_paths.append(contig_file_path)
                #clean_ass_ref = assembly_ref.replace('/','_')
                #assembly_outfile_path = os.join(assembly_outdir, clean_assembly_ref+".fna")
                #shutil.move(contig_file_path, assembly_outfile_path)


        #### STEP 3: Get distributions of contig attributes
        ##
        if len(invalid_msgs) == 0:

            # score fasta lens in contig files
            read_buf_size  = 65536
            #write_buf_size = 65536

            max_len = 0
            lens = []
            for ass_i,assembly_file_path in enumerate(score_assembly_file_paths):
                ass_name = assembly_names[ass_i]
                self.log (console, "Reading contig lengths in assembly: "+ass_name)  # DEBUG

                lens.append([])
                with open (assembly_file_path, 'r', read_buf_size) as ass_handle:
                    seq_buf = ''
                    for fasta_line in ass_handle:
                        if fasta_line.startswith('>'):
                            if seq_buf != '':
                                seq_len = len(seq_buf)
                                if seq_len > max_len:
                                    max_len = seq_len
                                lens[ass_i].append(seq_len)
                                seq_buf = ''
                        else:
                            seq_buf += ''.join(fasta_line.split())
                    if seq_buf != '':
                        seq_len = len(seq_buf)
                        if seq_len > max_len:
                            max_len = seq_len
                        lens[ass_i].append(seq_len)
                        seq_buf = ''

            # count buckets and develop histogram
            hist_window_width = 10000  # make it log scale?
            N_hist_windows = int(max_len % hist_window_width)  
            len_buckets = [ 1000000, 500000, 100000, 50000, 10000, 5000, 1000, 500, 0 ]
            summary_stats = []
            hist = []
            for ass_i,ass_name in enumerate(assembly_names):
                self.log (console, "Building summary and histograms from assembly: "+ass_name)  # DEBUG
                lens[ass_i].sort(key=int, reverse=True)

                # summary stats
                summary_stats.append(dict())
                for bucket in len_buckets:
                    summary_stats[ass_i][bucket] = 0
                    for val in lens[ass_i]:
                        if val >= bucket:
                            summary_stats[ass_i][bucket] += 1
                        else:
                            break

                # histogram
                hist.append([])
                for hist_i in range(N_hist_windows):
                    hist[ass_i].append(0)
                for val in lens[ass_i]:
                    hist_i = int(val / hist_window_width)
                    hist[ass_i][hist_i] += 1

            # determine max height across histograms
            max_hist_height = 0
            for ass_i,ass_name in enumerate(assembly_names):
                for hist_val in hist[ass_i]:
                    if hist_val > max_hist_height:
                        max_hist_height = hist_val


        #### STEP 4: build text report
        ##
        if len(invalid_msgs) == 0:
            for ass_i,ass_name in enumerate(assembly_names):
                report_text += "ASSEMBLY STATS for "+ass_name+"\n"
                for bucket in len_buckets:
                    report_text += "\t"+"Ncontigs >= "+str(bucket)+"bp:\t"+str(summary_stats[ass_i][bucket])+"\n"
                report_text += "\n"
        self.log(console, report_text)  # DEBUG


        #### STEP 5: Build report
        ##
        reportName = 'run_contig_distribution_compare_report_'+str(uuid.uuid4())
        reportObj = {'objects_created': [],
                     #'text_message': '',  # or is it 'message'?
                     'message': '',  # or is it 'text_message'?
                     'direct_html': '',
                     #'direct_html_link_index': 0,
                     'file_links': [],
                     'html_links': [],
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }

        # message
        if len(invalid_msgs) > 0:
            report_text = "\n".join(invalid_msgs)
        reportObj['message'] = report_text

        if len(invalid_msgs) == 0:
            pass
            # html report
            """
            try:
                html_upload_ret = dfuClient.file_to_shock({'file_path': html_output_dir,
                                                     'make_handle': 0,
                                                     'pack': 'zip'})
            except:
                raise ValueError ('error uploading html report to shock')
            reportObj['direct_html_link_index'] = 0
            reportObj['html_links'] = [{'shock_id': html_upload_ret['shock_id'],
                                        'name': html_file,
                                        'label': params['output_name']+' HTML'
                                    }
                                   ]
            """


        # save report object
        #
        SERVICE_VER = 'release'
        reportClient = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        #report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})
        report_info = reportClient.create_extended_report(reportObj)

        returnVal = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END run_contig_distribution_compare

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method run_contig_distribution_compare return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def run_benchmark_assemblies_against_genomes_with_MUMmer4(self, ctx, params):
        """
        :param params: instance of type
           "Benchmark_assemblies_against_genomes_with_MUMmer4_Params"
           (benchmark_assemblies_against_genomes_with_MUMmer4() ** **  Align
           benchmark genomes to assembly contigs) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_genome_refs" of type "data_obj_ref", parameter
           "input_assembly_refs" of type "data_obj_ref", parameter "desc" of
           String
        :returns: instance of type
           "Benchmark_assemblies_against_genomes_with_MUMmer4_Output" ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN run_benchmark_assemblies_against_genomes_with_MUMmer4

        #### STEP 0: basic init
        ##
        console = []
        invalid_msgs = []
        report_text = ''
        self.log(console, 'Running run_benchmark_assemblies_against_genomes_with_MUMmer4(): ')
        self.log(console, "\n"+pformat(params))

        # Auth
        token = ctx['token']
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # API Clients
        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'
        # wsClient
        try:
            wsClient = workspaceService(self.workspaceURL, token=token)
        except Exception as e:
            raise ValueError('Unable to instantiate wsClient with workspaceURL: '+ self.workspaceURL +' ERROR: ' + str(e))
        # setAPI_Client
        try:
            #setAPI_Client = SetAPI (url=self.callbackURL, token=ctx['token'])  # for SDK local.  local doesn't work for SetAPI
            setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
        except Exception as e:
            raise ValueError('Unable to instantiate setAPI_Client with serviceWizardURL: '+ self.serviceWizardURL +' ERROR: ' + str(e))
        # auClient
        try:
            auClient = AssemblyUtil(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate auClient with callbackURL: '+ self.callbackURL +' ERROR: ' + str(e))
        # dfuClient
        try:
            dfuClient = DFUClient(self.callbackURL)
        except Exception as e:
            raise ValueError('Unable to instantiate dfu_Client with callbackURL: '+ self.callbackURL +' ERROR: ' + str(e))

        # param checks
        required_params = ['workspace_name',
                           'input_genome_refs',
                           'input_assembly_refs',
                           'desc'
                          ]
        for arg in required_params:
            if arg not in params or params[arg] == None or params[arg] == '':
                raise ValueError ("Must define required param: '"+arg+"'")

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[]
        for input_ref in params['input_genome_refs']:
            provenance[0]['input_ws_objects'].append(input_ref)
        for input_ref in params['input_assembly_refs']:
            provenance[0]['input_ws_objects'].append(input_ref)

        # set the output paths
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        html_output_dir = os.path.join(output_dir,'html')
        if not os.path.exists(html_output_dir):
            os.makedirs(html_output_dir)


        #### STEP 1: get benchmark genome refs
        ##
        if len(invalid_msgs) == 0:
            set_obj_type = "KBaseSearch.GenomeSet"
            genome_obj_type = "KBaseGenomes.Genome"
            accepted_input_types = [set_obj_type, genome_obj_type]
            genome_refs = []
            genome_refs_seen = dict()
        
            for i,input_ref in enumerate(params['input_genome_refs']):
                # genome obj info
                try:
                    [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
                    input_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_ref}]})[0]
                    input_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_obj_info[TYPE_I])  # remove trailing version

                except Exception as e:
                    raise ValueError('Unable to get object from workspace: (' + input_ref +')' + str(e))
                if input_obj_type not in accepted_input_types:
                    raise ValueError ("Input object of type '"+input_obj_type+"' not accepted.  Must be one of "+", ".join(accepted_input_types))

                # add members to genome_ref list
                if input_obj_type == genome_obj_type:
                    try:
                        genome_seen = genome_refs_seen[input_ref]
                        continue
                    except:
                        genome_refs_seen[input_ref] = True
                    genome_refs.append(input_ref)
                elif input_obj_type != set_obj_type:
                    raise ValueError ("bad obj type for input_ref: "+input_ref)
                else:  # add genome set members
                    try:
                        objects= wsClient.get_objects2({'objects':[{'ref':input_ref}]})['data']
                    except Exception as e:
                        raise ValueError('Unable to get object from workspace: (' + input_ref +')' + str(e))
                    data = objects[0]['data']
                    info = objects[0]['info']
                    set_obj = data
                        
                    # get Genome items from Genome Set
                    for genome_id in set_obj['elements'].keys().sort():
                        genome_ref = set_obj['elements'][genome_id]['ref']
                        try:
                            genome_seen = genome_refs_seen[genome_ref]
                            continue
                        except:
                            genome_refs_seen[genome_ref] = True
                            genome_refs.append(genome_ref)


        #### STEP 2: get benchmark genome assembly seqs and other attributes
        ##
        if len(invalid_msgs) == 0:
            genome_obj_names = []
            genome_sci_names = []
            genome_assembly_refs = []
            
            for i,input_ref in enumerate(genome_refs):
                # genome obj data
                try:
                    [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
                    objects = wsClient.get_objects2({'objects':[{'ref':input_ref}]})['data']
                    genome_obj = objects[0]['data']
                    genome_obj_info = objects[0]['info']
                    genome_obj_names.append(genome_obj_info[NAME_I])
                    genome_sci_names.append(genome_obj['scientific_name'])
                except:
                    raise ValueError ("unable to fetch genome: "+input_ref)

                # Get genome_assembly_refs
                if ('contigset_ref' not in genome_obj or genome_obj['contigset_ref'] == None) \
                   and ('assembly_ref' not in genome_obj or genome_obj['assembly_ref'] == None):
                    msg = "Genome "+genome_obj_names[i]+" (ref:"+input_ref+") "+genome_sci_names[i]+" MISSING BOTH contigset_ref AND assembly_ref.  Cannot process.  Exiting."
                    self.log(console, msg)
                    self.log(invalid_msgs, msg)
                    continue
                elif 'assembly_ref' in genome_obj and genome_obj['assembly_ref'] != None:
                    msg = "Genome "+genome_obj_names[i]+" (ref:"+input_ref+") "+genome_sci_names[i]+" USING assembly_ref: "+str(genome_obj['assembly_ref'])
                    self.log (console, msg)
                    genome_assembly_refs.append(genome_obj['assembly_ref'])
                elif 'contigset_ref' in genome_obj and genome_obj['contigset_ref'] != None:
                    msg = "Genome "+genome_obj_names[i]+" (ref:"+input_ref+") "+genome_sci_names[i]+" USING contigset_ref: "+str(genome_obj['contigset_ref'])
                    self.log (console, msg)
                    genome_assembly_refs.append(genome_obj['contigset_ref'])

        # get fastas for scaffolds
        if len(invalid_msgs) == 0:
            #genomes_outdir = os.path.join (output_dir, 'benchmark_genomes')
            #if not os.path.exists(genomes_outdir):
            #    os.makedirs(genomes_outdir)
            read_buf_size  = 65536
            write_buf_size = 65536
            benchmark_assembly_file_paths = []

            for genome_i,input_ref in enumerate(genome_refs):
                contig_file = auClient.get_assembly_as_fasta({'ref':genome_assembly_refs[genome_i]}).get('path')
                sys.stdout.flush()
                contig_file_path = dfuClient.unpack_file({'file_path': contig_file})['file_path']
                benchmark_assembly_file_paths.append(contig_file_path)
                #clean_genome_ref = genome_ref.replace('/','_')
                #genome_outfile_path = os.join(benchmark_outdir, clean_genome_ref+".fna")
                #shutil.move(contig_file_path, genome_outfile_path)


        #### STEP 3: get assembly refs
        ##
        if len(invalid_msgs) == 0:
            set_obj_type = "KBaseSets.AssemblySet"
            assembly_obj_types = ["KBaseGenomeAnnotations.Assembly", "KBaseGenomes.ContigSet"]
            accepted_input_types = [set_obj_type] + assembly_obj_types
            assembly_refs = []
            assembly_refs_seen = dict()
        
            for i,input_ref in enumerate(params['input_assembly_refs']):
                # assembly obj info
                try:
                    [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
                    input_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_ref}]})[0]
                    input_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_obj_info[TYPE_I])  # remove trailing version

                    # DEBUG
                    input_obj_name = input_obj_info[NAME_I]
                    self.log (console, "GETTING ASSEMBLY: "+str(input_ref)+" "+str(input_obj_name))

                except Exception as e:
                    raise ValueError('Unable to get object from workspace: (' + input_ref +')' + str(e))
                if input_obj_type not in accepted_input_types:
                    raise ValueError ("Input object of type '"+input_obj_type+"' not accepted.  Must be one of "+", ".join(accepted_input_types))

                # add members to assembly_ref list
                if input_obj_type in assembly_obj_types:
                    try:
                        assembly_seen = assembly_refs_seen[input_ref]
                        continue
                    except:
                        assembly_refs_seen[input_ref] = True
                        assembly_refs.append(input_ref)
                elif input_obj_type != set_obj_type:
                    raise ValueError ("bad obj type for input_ref: "+input_ref)
                else:  # add assembly set members

                    try:
                        assemblySet_obj = setAPI_Client.get_assembly_set_v1 ({'ref':input_ref, 'include_item_info':1})
                    except Exception as e:
                        raise ValueError('Unable to get object from workspace: (' + input_ref +')' + str(e))
                    
                    for assembly_obj in assemblySet_obj['data']['items']:
                        this_assembly_ref = assembly_obj['ref']
                        try:
                            assembly_seen = assembly_refs_seen[this_assembly_ref]
                            continue
                        except:
                            assembly_refs_seen[this_assembly_ref] = True
                            assembly_refs.append(this_assembly_ref)


        #### STEP 4: Get assemblies to score as fasta files
        ##
        if len(invalid_msgs) == 0:
            #assembly_outdir = os.path.join (output_dir, 'score_assembly')
            #if not os.path.exists(assembly_outdir):
            #    os.makedirs(assembly_outdir)
            read_buf_size  = 65536
            write_buf_size = 65536
            score_assembly_file_paths = []

            for ass_i,input_ref in enumerate(assembly_refs):
                contig_file = auClient.get_assembly_as_fasta({'ref':assembly_refs[ass_i]}).get('path')
                sys.stdout.flush()
                contig_file_path = dfuClient.unpack_file({'file_path': contig_file})['file_path']
                score_assembly_file_paths.append(contig_file_path)
                #clean_ass_ref = assembly_ref.replace('/','_')
                #assembly_outfile_path = os.join(assembly_outdir, clean_assembly_ref+".fna")
                #shutil.move(contig_file_path, assembly_outfile_path)



        #### STEP 5: Run MUMmer
        ##
        if len(invalid_msgs) == 0:
            cmd = []
            cmd.append (self.NUCMER_bin)
#            # output
#            cmd.append ('-base_name')
#            cmd.append (params['output_name'])
#            cmd.append ('-output_dir')
#            cmd.append (output_dir)
#            # contigs input
#            cmd.append ('-reference_file')
#            cmd.append (genomes_src_db_file_path)


            # RUN
            cmd_str = " ".join(cmd)
            self.log (console, "===========================================")
            self.log (console, "RUNNING: "+cmd_str)
            self.log (console, "===========================================")

            """
            cmdProcess = subprocess.Popen(cmd_str, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            outputlines = []
            while True:
                line = cmdProcess.stdout.readline()
                outputlines.append(line)
                if not line: break
                self.log(console, line.replace('\n', ''))
                
            cmdProcess.stdout.close()
            cmdProcess.wait()
            self.log(console, 'return code: ' + str(cmdProcess.returncode) + '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running run_benchmark_assemblies_against_genomes_with_MUMmer4, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

            #report_text += "\n".join(outputlines)
            #report_text += "cmdstring: " + cmdstring + " stdout: " + stdout + " stderr " + stderr
            """


        #### STEP 5: Build report
        ##
        reportName = 'run_benchmark_assemblies_against_genomes_with_MUMmer4_report_'+str(uuid.uuid4())
        reportObj = {'objects_created': [],
                     #'text_message': '',  # or is it 'message'?
                     'message': '',  # or is it 'text_message'?
                     'direct_html': '',
                     #'direct_html_link_index': 0,
                     'file_links': [],
                     'html_links': [],
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }

        # message
        if len(invalid_msgs) > 0:
            report_text = "\n".join(invalid_msgs)
        reportObj['message'] = report_text

        if len(invalid_msgs) == 0:

            # html report
            """
            try:
                html_upload_ret = dfuClient.file_to_shock({'file_path': html_output_dir,
                                                     'make_handle': 0,
                                                     'pack': 'zip'})
            except:
                raise ValueError ('error uploading html report to shock')
            reportObj['direct_html_link_index'] = 0
            reportObj['html_links'] = [{'shock_id': html_upload_ret['shock_id'],
                                        'name': html_file,
                                        'label': params['output_name']+' HTML'
                                    }
                                   ]
            """


        # save report object
        #
        SERVICE_VER = 'release'
        reportClient = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        #report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})
        report_info = reportClient.create_extended_report(reportObj)

        returnVal = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END run_benchmark_assemblies_against_genomes_with_MUMmer4

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method run_benchmark_assemblies_against_genomes_with_MUMmer4 return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
