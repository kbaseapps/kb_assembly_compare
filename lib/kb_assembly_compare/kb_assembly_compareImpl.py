# -*- coding: utf-8 -*-
#BEGIN_HEADER
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
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        #END_CONSTRUCTOR
        pass


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
