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
#from Bio import SeqIO
#import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
#from matplotlib.patches import Arc
#from matplotlib.patches import Rectangle

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
                                this_input_obj_name = this_input_obj_info[NAME_I]
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
                                lens[ass_i].append(seq_len)
                                seq_buf = ''
                        else:
                            seq_buf += ''.join(fasta_line.split())
                    if seq_buf != '':
                        seq_len = len(seq_buf)
                        lens[ass_i].append(seq_len)
                        seq_buf = ''

            # sort lens (absolutely critical to subsequent steps)
            for ass_i,ass_name in enumerate(assembly_names):
                self.log (console, "Sorting contig lens for "+ass_name)  # DEBUG
                lens[ass_i].sort(key=int, reverse=True)  # sorting is critical

            # get min_max ranges
            huge_val = 100000000000000000
            len_buckets = [ 1000000, 100000, 10000, 1000, 500, 1 ]
            percs = [50, 75, 90]
            max_len = 0
            max_lens = []
            best_val = { 'len': 0,
                         'N': {},
                         'L': {},
                         'summary_stats': {},
                         'cumulative_len_stats': {}
                         }
            worst_val = { 'len': huge_val,
                          'N': {},
                          'L': {},
                          'summary_stats': {},
                          'cumulative_len_stats': {}
                          }
            for perc in percs:
                best_val['N'][perc] = 0
                best_val['L'][perc] = huge_val
                worst_val['N'][perc] = huge_val
                worst_val['L'][perc] = 0
            for bucket in len_buckets:
                best_val['summary_stats'][bucket] = 0
                best_val['cumulative_len_stats'][bucket] = 0
                worst_val['summary_stats'][bucket] = huge_val
                worst_val['cumulative_len_stats'][bucket] = huge_val

            for ass_i,ass_name in enumerate(assembly_names):
                self.log (console, "Getting max lens "+ass_name)  # DEBUG
                this_max_len = lens[ass_i][0]
                max_lens.append(this_max_len)
                if this_max_len > max_len:
                    max_len = this_max_len

            # Sum cumulative plot data
            max_total = 0
            total_lens = []
            cumulative_lens = []
            for ass_i,ass_name in enumerate(assembly_names):
                self.log (console, "Summing cumulative lens for "+ass_name)  # DEBUG
                total_lens.append(0)
                cumulative_lens.append([])
                for val in lens[ass_i]:
                    cumulative_lens[ass_i].append(total_lens[ass_i]+val)
                    total_lens[ass_i] += val
            for ass_i,ass_name in enumerate(assembly_names):
                if total_lens[ass_i] > max_total:
                    max_total = total_lens[ass_i]
                
            """
            # DEBUG
            self.log (console, "SORTED VALS\n================================")
            for ass_i,ass_name in enumerate(assembly_names):
                self.log (console, ass_name)
                self.log (console, "\t\t"+"TOTAL_LEN: "+str(total_lens[ass_i]))
                for val_i,val in enumerate(lens[ass_i]):
                    self.log (console, "\t\t\t"+"VAL: "+str(val))
                for len_i,length in enumerate(cumulative_lens[ass_i]):
                    self.log (console, "\t\t\t"+"LEN: "+str(length))
            # END DEBUG
            """

            # get N50 and L50 (and 75s, and 90s)
            N = { 50: [],
                  75: [],
                  90: []
              }
            L = { 50: [],
                  75: [],
                  90: []
              }
            for perc in N.keys():
                frac = perc/100.0
                for ass_i,ass_name in enumerate(assembly_names):
                    self.log (console, "Getting N/L"+str(perc)+" for "+ass_name)  # DEBUG
                    for val_i,val in enumerate(lens[ass_i]):
                        if cumulative_lens[ass_i][val_i] >= frac * total_lens[ass_i]:
                            N[perc].append(val)
                            L[perc].append(val_i+1)
                            break

            """
            # DEBUG
            self.log (console, "N50, etc.\n================================")
            for ass_i,ass_name in enumerate(assembly_names):
                self.log (console, ass_name)
                for perc in N.keys():
                    self.log (console, "\t"+"N"+str(perc)+": "+str(N[perc][ass_i]))
                    self.log (console, "\t"+"L"+str(perc)+": "+str(L[perc][ass_i]))
            # END DEBUG
            """

            # count buckets and develop histogram
            hist_window_width = 10000  # make it log scale?
            N_hist_windows = int(max_len % hist_window_width)  
            #len_buckets = [ 1000000, 500000, 100000, 50000, 10000, 5000, 1000, 500, 0 ]
            summary_stats = []
            cumulative_len_stats = []
            hist = []
            for ass_i,ass_name in enumerate(assembly_names):
                self.log (console, "Building summary and histograms from assembly: "+ass_name)  # DEBUG
                lens[ass_i].sort(key=int, reverse=True)  # sorting is critical

                # summary stats
                summary_stats.append(dict())
                cumulative_len_stats.append(dict())
                for bucket in len_buckets:
                    summary_stats[ass_i][bucket] = 0
                    cumulative_len_stats[ass_i][bucket] = 0
                curr_bucket_i = 0
                for val in lens[ass_i]:
                    for bucket_i in range(curr_bucket_i,len(len_buckets)):
                        bucket = len_buckets[bucket_i]
                        if val >= bucket:
                            summary_stats[ass_i][bucket] += 1
                            cumulative_len_stats[ass_i][bucket] += val
                            #curr_bucket_i = bucket_i
                            #break
                        else:
                            curr_bucket_i = bucket_i + 1

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

            # adjust best and worst values
            for ass_i,ass_name in enumerate(assembly_names):
                if max_lens[ass_i] > best_val['len']:
                    best_val['len'] = max_lens[ass_i]
                if max_lens[ass_i] < worst_val['len']:
                    worst_val['len'] = max_lens[ass_i]

                for perc in percs:
                    if N[perc][ass_i] > best_val['N'][perc]:
                        best_val['N'][perc] = N[perc][ass_i]
                    if N[perc][ass_i] < worst_val['N'][perc]:
                        worst_val['N'][perc] = N[perc][ass_i]
                    if L[perc][ass_i] < best_val['L'][perc]:
                        best_val['L'][perc] = L[perc][ass_i]
                    if L[perc][ass_i] > worst_val['L'][perc]:
                        worst_val['L'][perc] = L[perc][ass_i]

                for bucket in len_buckets:
                    if summary_stats[ass_i][bucket] > best_val['summary_stats'][bucket]:
                        best_val['summary_stats'][bucket] = summary_stats[ass_i][bucket]
                    if summary_stats[ass_i][bucket] < worst_val['summary_stats'][bucket]:
                        worst_val['summary_stats'][bucket] = summary_stats[ass_i][bucket]
                    if cumulative_len_stats[ass_i][bucket] > best_val['cumulative_len_stats'][bucket]:
                        best_val['cumulative_len_stats'][bucket] = cumulative_len_stats[ass_i][bucket]
                    if cumulative_len_stats[ass_i][bucket] < worst_val['cumulative_len_stats'][bucket]:
                        worst_val['cumulative_len_stats'][bucket] = cumulative_len_stats[ass_i][bucket]


        #### STEP 4: build text report
        ##
        if len(invalid_msgs) == 0:
            for ass_i,ass_name in enumerate(assembly_names):
                report_text += "ASSEMBLY STATS for "+ass_name+"\n"

                report_text += "\t"+"Len longest contig: "+str(max_lens[ass_i])+" bp"+"\n"
                for perc in N.keys():
                    report_text += "\t"+"N"+str(perc)+" (L"+str(perc)+"):\t"+str(N[perc][ass_i])+" ("+str(L[perc][ass_i])+")"+"\n"
                for bucket in len_buckets:
                    report_text += "\t"+"Num contigs >= "+str(bucket)+" bp:\t"+str(summary_stats[ass_i][bucket])+"\n"
                report_text += "\n"

                for bucket in len_buckets:
                    report_text += "\t"+"Len contigs >= "+str(bucket)+" bp:\t"+str(cumulative_len_stats[ass_i][bucket])+" bp"+"\n"
                report_text += "\n"

        self.log(console, report_text)  # DEBUG


        #### STEP 5: Make figures with matplotlib
        ##
        file_links = [] 

        # Cumulative len plot
        plot_name = "cumulative_len_plot"
        plot_name_desc = "Cumulative Length (in bp)"
        self.log (console, "GENERATING PLOT "+plot_name_desc)
        img_dpi = 200
        img_units = "in"
        img_in_width  = 6.0
        img_in_height = 3.0
        x_margin = 0.01
        y_margin = 0.01
        title_fontsize = 12
        text_color = "#606060"
        fig = plt.figure()
        fig.set_size_inches(img_in_width, img_in_height)
        ax = plt.subplot2grid ( (1,1), (0,0), rowspan=1, colspan=1)
        """
        # Let's turn off visibility of all tic labels and boxes here
        for ax in fig.axes:
            ax.xaxis.set_visible(False)  # remove axis labels and tics
            ax.yaxis.set_visible(False)
            for t in ax.get_xticklabels()+ax.get_yticklabels():  # remove tics
                t.set_visible(False)
            ax.spines['top'].set_visible(False)     # Get rid of top axis line
            ax.spines['bottom'].set_visible(False)  # bottom axis line
            ax.spines['left'].set_visible(False)    # left axis line
            ax.spines['right'].set_visible(False)   # right axis line
        """
        ax = fig.axes[0]
        ax.set_title (plot_name_desc)
        #ax.text (x_margin, 1.0-(y_margin), plot_name, verticalalignment="bottom", horizontalalignment="left", color=text_color, fontsize=title_fontsize, zorder=1)

        # build x and y coord lists
        for ass_i,ass_name in enumerate(assembly_names):
            x_coords = []
            y_coords = []
            for val_i,val in enumerate(cumulative_lens[ass_i]):
                x_coords.append(val_i+1)
                y_coords.append(val)
            plt.plot(x_coords, y_coords, lw=2)

        """
        # capture data into pandas
        cumulative_lens_by_ass_name = dict()
        for ass_i,ass_name in enumerate(assembly_names):
            cumulative_lens_by_ass_name[ass_name] = pd.Series(cumulative_lens[ass_i])
        cumulative_lens_df = pd.DataFrame(cumulative_lens_by_ass_name)
        cumulative_lens_plot = cumulative_lens_df \
                                .plot(kind="line", figsize=(15,5), ylim=(0,max_total), fontsize=15, lw=5)
        cumulative_lens_plot.xaxis.grid(True)
        cumulative_lens_plot.yaxis.grid(True) 
        fig = cumulative_lens_plot.get_figure()
        """

        # save plot
        self.log (console, "SAVING PLOT "+plot_name_desc)
        cumulative_lens_png_file = png_file = plot_name+".png"
        cumulative_lens_pdf_file = pdf_file = pdf_file = plot_name+".pdf"
        output_png_file_path = os.path.join (html_output_dir, png_file)
        output_pdf_file_path = os.path.join (html_output_dir, pdf_file)
        fig.savefig (output_png_file_path, dpi=img_dpi)
        fig.savefig (output_pdf_file_path, format='pdf')

        # upload PNG
        try:
            upload_ret = dfuClient.file_to_shock({'file_path': output_png_file_path,
                                                  'make_handle': 0,
                                                  'pack': 'zip'})
            file_links.append({'shock_id': upload_ret['shock_id'],
                               'name': png_file,
                               'label': plot_name_desc+' PNG'
                               }
                              )
        except:
            raise ValueError ('Logging exception loading png_file '+png_file+' to shock')
        # upload PDF
        try:
            upload_ret = dfuClient.file_to_shock({'file_path': output_pdf_file_path,
                                                  'make_handle': 0,
                                                  'pack': 'zip'})
            file_links.append({'shock_id': upload_ret['shock_id'],
                               'name': pdf_file,
                               'label': plot_name_desc+' PDF'
                               }
                              )
        except:
            raise ValueError ('Logging exception loading pdf_file '+pdf_file+' to shock')


        # Sorted Contig len plot
        plot_name = "sorted_contig_lengths"
        plot_name_desc = "Sorted Contig Lengths (in bp)"
        self.log (console, "GENERATING PLOT "+plot_name_desc)
        img_dpi = 200
        img_units = "in"
        img_in_width  = 6.0
        img_in_height = 3.0
        x_margin = 0.01
        y_margin = 0.01
        title_fontsize = 12
        text_color = "#606060"
        fig = plt.figure()
        fig.set_size_inches(img_in_width, img_in_height)
        ax = plt.subplot2grid ( (1,1), (0,0), rowspan=1, colspan=1)
        """
        # Let's turn off visibility of all tic labels and boxes here
        for ax in fig.axes:
            ax.xaxis.set_visible(False)  # remove axis labels and tics
            ax.yaxis.set_visible(False)
            for t in ax.get_xticklabels()+ax.get_yticklabels():  # remove tics
                t.set_visible(False)
            ax.spines['top'].set_visible(False)     # Get rid of top axis line
            ax.spines['bottom'].set_visible(False)  # bottom axis line
            ax.spines['left'].set_visible(False)    # left axis line
            ax.spines['right'].set_visible(False)   # right axis line
        """
        ax = fig.axes[0]
        ax.set_title (plot_name_desc)
        #ax.text (x_margin, 1.0-(y_margin), plot_name, verticalalignment="bottom", horizontalalignment="left", color=text_color, fontsize=title_fontsize, zorder=1)

        # build x and y coord lists
        mini_delta = .000001
        for ass_i,ass_name in enumerate(assembly_names):
            x_coords = []
            y_coords = []
            running_sum = 0
            for val_i,val in enumerate(lens[ass_i]):
                x_coords.append(running_sum + mini_delta)
                y_coords.append(val)
                running_sum += val
                x_coords.append(running_sum)
                y_coords.append(val)
            plt.plot(x_coords, y_coords, lw=2)

        # save plot
        self.log (console, "SAVING PLOT "+plot_name_desc)
        sorted_lens_png_file = png_file = plot_name+".png"
        sorted_pens_pdf_file = pdf_file = plot_name+".pdf"
        output_png_file_path = os.path.join (html_output_dir, png_file)
        output_pdf_file_path = os.path.join (html_output_dir, pdf_file)
        fig.savefig (output_png_file_path, dpi=img_dpi)
        fig.savefig (output_pdf_file_path, format='pdf')

        # upload PNG
        try:
            upload_ret = dfuClient.file_to_shock({'file_path': output_png_file_path,
                                                  'make_handle': 0,
                                                  'pack': 'zip'})
            file_links.append({'shock_id': upload_ret['shock_id'],
                               'name': png_file,
                               'label': plot_name_desc+' PNG'
                               }
                              )
        except:
            raise ValueError ('Logging exception loading png_file '+png_file+' to shock')
        # upload PDF
        try:
            upload_ret = dfuClient.file_to_shock({'file_path': output_pdf_file_path,
                                                  'make_handle': 0,
                                                  'pack': 'zip'})
            file_links.append({'shock_id': upload_ret['shock_id'],
                               'name': pdf_file,
                               'label': plot_name_desc+' PDF'
                               }
                              )
        except:
            raise ValueError ('Logging exception loading pdf_file '+pdf_file+' to shock')


        #### STEP 6: Create and Upload HTML Report
        ##
        self.log (console, "CREATING HTML REPORT")
        def get_cell_color (val, best, worst, low_good=False):
            self.log (console, "VAL: "+str(val)+" BEST: "+str(best)+" WORST: "+str(worst))  # DEBUG

            if best == worst:
                return '#ffffff'
            val_color_map = { 0: '00',
                              1: '11',
                              2: '22',
                              3: '33',
                              4: '44',
                              5: '55',
                              6: '66',
                              7: '77',
                              8: '88',
                              9: '99',
                              10: 'aa',
                              11: 'bb',
                              12: 'cc',
                              13: 'dd',
                              14: 'ee',
                              15: 'ff'
                             }
            base_intensity = 0  # 5
            top = 15 - base_intensity
            mid = 0.5 * (best + worst)
            if val == mid:
                return '#ffffff'
            if low_good:
                if val < mid:
                    rescaled_val = int(0.5 + top * (val-best) / (mid-best))
                    self.log (console, "A, MID: "+str(mid)+" RESCALED_VAL: "+str(rescaled_val))  # DEBUG
                    r = val_color_map[base_intensity + rescaled_val]
                    g = val_color_map[base_intensity + rescaled_val]
                    b = 'ff'
                else:
                    rescaled_val = int(0.5 + top * (val-worst) / (mid-worst))
                    self.log (console, "B, MID: "+str(mid)+" RESCALED_VAL: "+str(rescaled_val))  # DEBUG
                    r = 'ff'
                    g = val_color_map[base_intensity + rescaled_val]
                    b = val_color_map[base_intensity + rescaled_val]
            else:
                if val > mid:
                    rescaled_val = int(0.5 + top * (val-best) / (mid-best))
                    self.log (console, "C, MID: "+str(mid)+" RESCALED_VAL: "+str(rescaled_val))  # DEBUG
                    r = val_color_map[base_intensity + rescaled_val]
                    g = val_color_map[base_intensity + rescaled_val]
                    b = 'ff'
                else:
                    rescaled_val = int(0.5 + top * (val-worst) / (mid-worst))
                    self.log (console, "D, MID: "+str(mid)+" RESCALED_VAL: "+str(rescaled_val))  # DEBUG
                    r = 'ff'
                    g = val_color_map[base_intensity + rescaled_val]
                    b = val_color_map[base_intensity + rescaled_val]
            self.log (console, "RGB: "+r+g+b)  # DEBUG
            return '#'+r+g+b

        hist_colspan = 5
        col_width   = 5 + hist_colspan  # in cells
        half_col_width = col_width // 2
        img_height = 300  # in pixels
        head_color = "#eeeeff"
        border_head_color = "#ffccff"
        text_fontsize = "2"
        text_color = '#606060'
        border_body_color = "#cccccc"
        base_cell_color = "#eeeeee"
        cellpadding = "3"
        cellspacing = "2"
        subtab_cellpadding = "1"
        subtab_cellspacing = "2"
        border = "0"
        sp = '&nbsp;'

        html_report_lines = []
        html_report_lines += ['<html>']
        html_report_lines += ['<head>']
        html_report_lines += ['<title>KBase Assembled Contig Distributions</title>']
#        html_report_lines += ['<style>']
#        html_report_lines += [".vertical-text {\ndisplay: inline-block;\noverflow: hidden;\nwidth: 0.65em;\n}\n.vertical-text__inner {\ndisplay: inline-block;\nwhite-space: nowrap;\nline-height: 1.1;\ntransform: translate(0,100%) rotate(-90deg);\ntransform-origin: 0 0;\n}\n.vertical-text__inner:after {\ncontent: \"\";\ndisplay: block;\nmargin: 0.0em 0 100%;\n}"]
#        html_report_lines += [".vertical-text_title {\ndisplay: inline-block;\noverflow: hidden;\nwidth: 1.0em;\n}\n.vertical-text__inner_title {\ndisplay: inline-block;\nwhite-space: nowrap;\nline-height: 1.0;\ntransform: translate(0,100%) rotate(-90deg);\ntransform-origin: 0 0;\n}\n.vertical-text__inner_title:after {\ncontent: \"\";\ndisplay: block;\nmargin: 0.0em 0 100%;\n}"]
#        html_report_lines += ['</style>']
        html_report_lines += ['</head>']
        html_report_lines += ['<body bgcolor="white">']

        #html_report_lines += ['<tr><td valign=top align=left rowspan=1><div class="vertical-text_title"><div class="vertical-text__inner_title"><font color="'+text_color+'">'+label+'</font></div></div></td>']

        html_report_lines += ['<table cellpadding='+str(cellpadding)+' cellspacing='+str(cellspacing)+' border='+str(border)+'>']
        html_report_lines += ['<tr><td valign=top align=left rowspan=1 colspan='+str(half_col_width)+'><img src="'+cumulative_lens_png_file+'" height='+str(img_height)+'></td>']
        html_report_lines += ['<td valign=top align=left rowspan=1 colspan='+str(half_col_width)+'><img src="'+sorted_lens_png_file+'" height='+str(img_height)+'></td></tr>']

        # header
        html_report_lines += ['<tr><td>'+sp+'</td></tr>']
        html_report_lines += ['<tr><td>'+sp+'</td></tr>']
        html_report_lines += ['<tr bgcolor="'+head_color+'">']
        # name
        html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+' align="left">'+'ASSEMBLY'+'</font></td>']
        # Longest Len
        html_report_lines += ['<td align="center" style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'Longest<br>Contig<br>(bp)'+'</font></td>']
        # N50,L50 etc.
        html_report_lines += ['<td align="center" style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'" colspan=2><font color="'+text_color+'" size='+text_fontsize+'>'+'Nx (Lx)'+'</font></td>']
        # Summary Stats
        html_report_lines += ['<td align="center" style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'LENGTH<br>(bp)'+'</font></td>']
        html_report_lines += ['<td align="center" style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'NUM CONTIGS'+'</font></td>']
        html_report_lines += ['<td align="center" style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'<nobr>SUM LENGTH</nobr><br>(bp)'+'</font></td>']
        # hist
        #html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'" colspan='+str(hist_colspan)+'><font color="'+text_color+'" size='+text_fontsize+'>'+'Contig Length Histogram'+'</font></td>']
        html_report_lines += ['</tr>']

        # report stats
        for ass_i,ass_name in enumerate(assembly_names):
            html_report_lines += ['<tr>']
            # name
            html_report_lines += ['<td align="left" bgcolor="'+base_cell_color+'" rowspan=6><font color="'+text_color+'" size='+text_fontsize+'>'+ass_name+'</font></td>']

            # longest contig
            cell_color = get_cell_color (max_lens[ass_i], best_val['len'], worst_val['len'])
            html_report_lines += ['<td bgcolor="'+cell_color+'" align="center" valign="middle" rowspan=6><font color="'+text_color+'" size='+text_fontsize+'>'+str(max_lens[ass_i])+'</font></td>']

            # subtable
            N_rows = 6
            edges = ' style="border-right:solid 2px '+border_body_color+'"'
            bottom_edge = ''
            for sub_i in range(N_rows):
                perc = sorted(N.keys(), key=int)[sub_i // 2]
                bucket = len_buckets[sub_i]
                if sub_i == N_rows-1:
                    edges = ' style="border-right:solid 2px '+border_body_color+'; border-bottom:solid 2px '+border_body_color+'"'
                    bottom_edge = ' style="border-bottom:solid 2px '+border_body_color+'"'

                # N50, L50, etc.
                if sub_i > 0:
                    html_report_lines += ['<tr>']
                
                if (sub_i % 2) == 0:
                    cell_color = get_cell_color (N[perc][ass_i], best_val['N'][perc], worst_val['N'][perc])
                    html_report_lines += ['<td align="right"'+bottom_edge+'>'+'<font color="'+text_color+'" size='+text_fontsize+'>'+'N'+str(perc)+':</font></td><td bgcolor="'+cell_color+'" align="right"'+edges+'>'+'<font color="'+text_color+'" size='+text_fontsize+'>'+sp+str(N[perc][ass_i])+'</font></td>']
                else:
                    cell_color = get_cell_color (L[perc][ass_i], best_val['L'][perc], worst_val['L'][perc], low_good=True)
                    html_report_lines += ['<td align="right"'+bottom_edge+'>'+'<font color="'+text_color+'" size='+text_fontsize+'>'+'L'+str(perc)+':</font></td><td bgcolor="'+cell_color+'" align="right"'+edges+'>'+'<font color="'+text_color+'" size='+text_fontsize+'>'+sp+'('+str(L[perc][ass_i])+')'+'</font></td></tr>']

                # Summary Stats
                html_report_lines += ['<td align="right"'+bottom_edge+'>'+'<font color="'+text_color+'" size='+text_fontsize+'>']
                if bucket >= 1000:
                    html_report_lines += ['<nobr>'+'&gt;= '+'10'+'<sup>'+str(int(math.log(bucket,10)+0.1))+'</sup>'+'</nobr>']
                else:
                    html_report_lines += ['<nobr>'+'&gt;= '+str(bucket)+'</nobr>']
                html_report_lines += ['</font></td>']
                
                cell_color = get_cell_color (summary_stats[ass_i][bucket], best_val['summary_stats'][bucket], worst_val['summary_stats'][bucket])
                html_report_lines += ['<td bgcolor="'+cell_color+'" align="right"'+bottom_edge+'>'+'<font color="'+text_color+'" size='+text_fontsize+'>'+str(summary_stats[ass_i][bucket])+'</font></td>']

                cell_color = get_cell_color (cumulative_len_stats[ass_i][bucket], best_val['cumulative_len_stats'][bucket], worst_val['cumulative_len_stats'][bucket])
                html_report_lines += ['<td bgcolor="'+cell_color+'" align="right"'+edges+'>'+'<font color="'+text_color+'" size='+text_fontsize+'>'+str(cumulative_len_stats[ass_i][bucket])+'</font></td>']
                html_report_lines += ['</tr>']

            # Hist

        html_report_lines += ['</table>']
        html_report_lines += ['</body>']
        html_report_lines += ['</html>']

        # write html to file and upload
        self.log (console, "SAVING AND UPLOADING HTML REPORT")
        html_report_str = "\n".join(html_report_lines)
        html_file = 'contig_distribution_report.html'
        html_file_path = os.path.join (html_output_dir, html_file)
        with open (html_file_path, 'w', 0) as html_handle:
            html_handle.write(html_report_str)
        try:
            html_upload_ret = dfuClient.file_to_shock({'file_path': html_output_dir,
                                                       'make_handle': 0,
                                                       'pack': 'zip'})
        except:
            raise ValueError ('Logging exception loading html_report to shock')


        #### STEP N: Build report
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

            # html report
            reportObj['direct_html_link_index'] = 0
            reportObj['html_links'] = [{'shock_id': html_upload_ret['shock_id'],
                                        'name': html_file,
                                        'label': 'Contig Distribution Report'+' HTML'
                                    }
                                   ]
            reportObj['file_links'] = file_links


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
