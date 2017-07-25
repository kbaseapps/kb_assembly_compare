# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from kb_assembly_compare.kb_assembly_compareImpl import kb_assembly_compare
from kb_assembly_compare.kb_assembly_compareServer import MethodContext
from kb_assembly_compare.authclient import KBaseAuth as _KBaseAuth


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
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

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
    def test_run_benchmark_assemblies_against_genomes_with_MUMmer4_test_01 (self):
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
