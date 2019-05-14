# -*- coding: utf-8 -*-
import os  # noqa: F401
import shutil
import time
import unittest
from configparser import ConfigParser  # py3
from os import environ
from pprint import pprint  # noqa: F401

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService
from kb_assembly_compare.authclient import KBaseAuth as _KBaseAuth
from kb_assembly_compare.kb_assembly_compareImpl import kb_assembly_compare
from kb_assembly_compare.kb_assembly_compareServer import MethodContext


class kb_assembly_compareTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_assembly_compare'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_assembly_compare',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_assembly_compare(cls.cfg)
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.scratch = cls.cfg['scratch']
        if not os.path.exists(cls.scratch):
            os.makedirs(cls.scratch)        

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_assembly_compare_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx


    ##############
    # UNIT TESTS #
    ##############

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa

    #### test_filter_contigs_by_length_01()
    ##
    def test_filter_contigs_by_length_01 (self):
        method = 'filter_contigs_by_length_01'
        
        print ("\n\nRUNNING: test_filter_contigs_by_length_01()")
        print ("===========================================\n\n")

        # upload test data
        try:
            auClient = AssemblyUtil(self.callback_url, token=self.getContext()['token'])
        except Exception as e:
            raise ValueError('Unable to instantiate auClient with callbackURL: '+ self.callback_url +' ERROR: ' + str(e))
        ass_file_1 = 'assembly_1.fa'
        ass_file_2 = 'assembly_2.fa'
        ass_path_1 = os.path.join(self.scratch, ass_file_1)
        ass_path_2 = os.path.join(self.scratch, ass_file_2)
        shutil.copy(os.path.join("data", ass_file_1), ass_path_1)
        shutil.copy(os.path.join("data", ass_file_2), ass_path_2)
        ass_ref_1 = auClient.save_assembly_from_fasta({
            'file': {'path': ass_path_1},
            'workspace_name': self.getWsName(),
            'assembly_name': 'assembly_1'
        })
        ass_ref_2 = auClient.save_assembly_from_fasta({
            'file': {'path': ass_path_2},
            'workspace_name': self.getWsName(),
            'assembly_name': 'assembly_2'
        })

        # run method
        input_refs = [ ass_ref_1, ass_ref_2 ]
        base_output_name = method+'_output'
        params = {
            'workspace_name': self.getWsName(),
            'input_assembly_refs': input_refs,
            'min_contig_length': 1000,
            'output_name': 'test_filtered'
        }
        result = self.getImpl().run_filter_contigs_by_length(self.getContext(),params)
        print('RESULT:')
        pprint(result)
        pass


    #### test_contig_distribution_compare_01()
    ##
    def test_contig_distribution_compare_01 (self):
        method = 'contig_distribution_compare_01'
        
        print ("\n\nRUNNING: test_contig_distribution_compare_01()")
        print ("==============================================\n\n")

        # upload test data
        try:
            auClient = AssemblyUtil(self.callback_url, token=self.getContext()['token'])
        except Exception as e:
            raise ValueError('Unable to instantiate auClient with callbackURL: '+ self.callback_url +' ERROR: ' + str(e))
        ass_file_1 = 'assembly_1.fa'
        ass_file_2 = 'assembly_2.fa'
        ass_path_1 = os.path.join(self.scratch, ass_file_1)
        ass_path_2 = os.path.join(self.scratch, ass_file_2)
        shutil.copy(os.path.join("data", ass_file_1), ass_path_1)
        shutil.copy(os.path.join("data", ass_file_2), ass_path_2)
        ass_ref_1 = auClient.save_assembly_from_fasta({
            'file': {'path': ass_path_1},
            'workspace_name': self.getWsName(),
            'assembly_name': 'assembly_1'
        })
        ass_ref_2 = auClient.save_assembly_from_fasta({
            'file': {'path': ass_path_2},
            'workspace_name': self.getWsName(),
            'assembly_name': 'assembly_2'
        })

        # run method
        input_refs = [ ass_ref_1, ass_ref_2 ]
        base_output_name = method+'_output'
        params = {
            'workspace_name': self.getWsName(),
            'input_assembly_refs': input_refs
        }
        result = self.getImpl().run_contig_distribution_compare(self.getContext(),params)
        print('RESULT:')
        pprint(result)
        pass


    def HIDE_run_benchmark_assemblies_against_genomes_with_MUMmer4_01 (self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods

        # input_data
        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #reference_prok_genomes_WS = '19217'  # PROD
        #reference_prok_genomes_WS = '15792'  # CI

        genome_ref_1 = reference_prok_genomes_WS+'/GCF_001566335.1/1'  # E. coli K-12 MG1655
        genome_ref_2 = reference_prok_genomes_WS+'/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'
        genome_ref_3 = reference_prok_genomes_WS+'/GCF_000006805.1/1'  # Halobacterium salinarum NRC-1 (3 contigs)

        sister_genome_ref_1 = reference_prok_genomes_WS+'/GCF_000357445.2/1'  # Escherichia coli P0298942.3 (207 contigs)
        sister_genome_ref_2 = reference_prok_genomes_WS+'/GCF_000195755.1/1'  # DvH (2 contigs)
        sister_genome_ref_3 = reference_prok_genomes_WS+'/GCF_000069025.1/1'  # Halobacterium salinarum R1 (5 contigs)

        # PROD refs
        #sister_assembly_ref_1 = reference_prok_genomes_WS+'/442/1'
        #sister_assembly_ref_2 = reference_prok_genomes_WS+'/48316/1'
        #sister_assembly_ref_3 = reference_prok_genomes_WS+'/194909/1'

        # CI refs
        sister_assembly_ref_1 = reference_prok_genomes_WS+'/29930/1'
        sister_assembly_ref_2 = reference_prok_genomes_WS+'/62326/1'
        sister_assembly_ref_3 = reference_prok_genomes_WS+'/114529/1'

        parameters = { 'workspace_name': self.getWsName(),
                       'desc': 'test assembly benchmark',
                       'input_genome_refs': [genome_ref_1, genome_ref_2, genome_ref_3],
                       'input_assembly_refs': [sister_assembly_ref_1, sister_assembly_ref_2, sister_assembly_ref_3]
                     }

        ret = self.getImpl().run_benchmark_assemblies_against_genomes_with_MUMmer4(self.getContext(), parameters)

        pass


